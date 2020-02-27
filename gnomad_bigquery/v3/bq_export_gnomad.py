import argparse
import sys

import hail as hl

from gnomad_bigquery.v3.bq_utils import export_ht_for_bq, get_gnomad_data_for_bq, logger
from gnomad_hail.utils.constants import CSQ_ORDER
from gnomad_hail.utils.slack import try_slack
from gnomad_hail.resources.grch38.reference_data import lcr_intervals
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources import annotations

def export_genotypes(
        data_type: str,
        export_missing_genotypes: bool,
        output_dir: str,
        least_consequence: str = None,
        variant_index: bool = True,
) -> None:
    """
    Export gnomAD genotypes, present or missing, within variants where
    there is a transcript consequence that is at least the least_consequence passed
    :param data_type: exomes or genomes (used to write out a path)
    :param export_missing_genotypes: whether to export only missing GTs
    :param output_dir: Directory where to write parquet file
    :param least_consequence: Least consequence admitted, default is `3_prime_UTR_variant`
    :param variant_index: Index variant to save on storage cost
    :return:
    """
    mt = get_gnomad_data_for_bq(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=False,
        adjust_sex_and_annotate_adj=True,
    )
    mt = mt.select_cols().select_rows().add_row_index()

    vep = annotations.vep.ht()
    if least_consequence is not None:
        vep_consequences = hl.literal(
            set(CSQ_ORDER[0 : CSQ_ORDER.index(least_consequence) + 1])
        )
        vep = vep.filter(
            vep.vep.transcript_consequences.any(
                lambda x: (x.biotype == "protein_coding")
                & hl.any(lambda csq: vep_consequences.contains(csq), x.consequence_terms)
            )
        )
    else:
        vep = vep.filter(
            vep.vep.transcript_consequences.any(lambda x: x.biotype == "protein_coding")
        )
    vep = vep.persist()
    logger.info(
        f"Found {vep.count()} variants with a VEP coding transcript and a consequence worst or equal to {least_consequence}."
    )

    select_expr = hl.is_defined(vep[mt.row_key])

    freq = annotations.freq.ht()
    select_expr = select_expr & hl.is_defined(freq[mt.row_key])

    mt = mt.filter_rows(select_expr)
    ht = mt.entries()
    if export_missing_genotypes:
        ht = ht.filter(ht.is_missing)
    else:
        ht = ht.filter(hl.is_defined(ht.GT))

    ht = ht.key_by()
    if variant_index:
        select_expr = {"v": ht.row_idx}
    else:
        select_expr = {
            "chrom": ht.locus.contig,
            "pos": ht.locus.position,
            "ref": ht.alleles[0],
            "alt": ht.alleles[1],
        }
    select_expr.update(
        {
            "s": ht.s,
            "is_het": ht.GT.is_het(),
            "is_adj": ht.adj,
            "dp": ht.DP,
            "gq": ht.GQ,
            "ad0": ht.AD[0],
            "ad1": ht.AD[1],
            "pl0": ht.PL[0],
            "pl1": ht.PL[1],
            "pl2": ht.PL[2],
            "pid": ht.PID,
            "pgt": ht.PGT,
        }
    )
    ht = ht.select(**select_expr)

    ht.to_spark().write.parquet(
        "{}/gnomad_{}{}_genotypes.parquet".format(
            output_dir, data_type, "missing_" if export_missing_genotypes else ""
        ),
        mode="overwrite",
    )


def export_variants(
        data_type: str,
        export_subsets_freq: bool,
        export_all_sex_freq: bool,
        output_dir: str,
        add_vep: bool
) -> None:
    """
    Exports all variants in gnomAD with their frequencies and optional vep transcript consequences
    :param data_type: exomes or genomes (used to write out a path)
    :param export_subsets_freq: Whether to export subset frequencies alongside overall frequencies
    :param export_all_sex_freq: Whether to export sex frequencies alongside overall frequencies
    :param output_dir: Directory where to write parquet file
    :param add_vep: Whether to annotate with vep's transcript_consequences and most_severe_consequence
    :return:
    """
    def format_freq(ht: hl.Table, subset: str):
        def filter_freqs(x: hl.expr.ArrayExpression):
            if export_all_sex_freq:
                return ~x[0].contains("platform") & ~x[0].contains("downsampling")
            else:
                return (
                    ~x[0].contains("platform")
                    & ~x[0].contains("downsampling")
                    & ~(x[0].contains("sex") & x[0].contains("pop"))
                )

        def get_freq_struct(x: hl.expr.ArrayExpression):
            freq_struct = hl.struct(
                adj=x[0]["group"] == "adj",
                pop=x[0].get("subpop", x[0].get("pop", "all")),
                sex=x[0].get("sex", "all"),
                ac=x[1].AC,
                an=x[1].AN,
                af=x[1].AF,
                hom=x[1].homozygote_count,
            )

            if export_subsets_freq:
                freq_struct = freq_struct.annotate(subset=subset)
            return freq_struct

        return ht.annotate(
            freq=hl.zip(freq_meta, freq[ht.key].freq)
            .filter(filter_freqs)
            .map(get_freq_struct)
        )

    mt = get_gnomad_data_for_bq(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=False,
        adjust_sex_and_annotate_adj=False,
    )
    ht = mt.rows()
    ht = ht.select().add_index()  # TODO: Add interesting annotations
    freq = annotations.freq.ht()
    freq_meta = hl.literal(freq.freq_meta.collect()[0])
    ht = format_freq(ht, "all")

    if add_vep:
        vep = annotations.vep.ht()
        vep = vep.select(
            transcript_consequences=vep.vep.transcript_consequences,
            most_severe_consequence=vep.vep.most_severe_consequence,
        )
        ht = ht.annotate(**vep[ht.key])

    ht = ht.annotate(lcr=hl.is_defined(lcr_intervals.ht()[ht.locus]))
    ht = ht.annotate(nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar()))

    export_ht_for_bq(ht, f"{output_dir}/gnomad_{data_type}_variants.parquet")


def main(args):

    hl.init(log="/bq.log", default_reference='GRCh38')
    data_types = []
    if args.exomes:
        data_types.append("exomes")
    if args.genomes:
        data_types.append("genomes")

    for data_type in data_types:

        if args.export_metadata:
            metadata = meta.ht()
            export_ht_for_bq(
                metadata, f"{args.output_dir}/gnomad_{data_type}_meta.parquet"
            )

        if args.export_genotypes:
            export_genotypes(
                data_type,
                False,
                args.output_dir,
                args.least_consequence,
                True,
            )

        if args.export_missing_genotypes:
            export_genotypes(
                data_type,
                True,
                args.output_dir,
                args.least_consequence,
                True,
            )

        if args.export_variants:
            export_variants(
                data_type,
                args.export_subset_freq,
                args.export_all_sex_freq,
                args.output_dir,
                True,
            )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--exomes",
        help="Run on exomes. At least one of --exomes or --genomes is required.",
        action="store_true",
    )
    parser.add_argument(
        "--genomes",
        help="Run on genomes. At least one of --exomes or --genomes is required.",
        action="store_true",
    )
    parser.add_argument(
        "--export_metadata", help="Export samples metadata", action="store_true"
    )
    parser.add_argument(
        "--export_genotypes", help="Export non-ref genotypes", action="store_true"
    )
    parser.add_argument(
        "--export_variants", help="Export variants", action="store_true"
    )
    parser.add_argument(
        "--output_dir",
        help="Output root. default: gs://gnomad-tmp/bq",
        default="gs://gnomad-tmp/bq",
    )

    var_exp = parser.add_argument_group(
        "Export variants", description="Options related to exporting gnomAD variants"
    )
    var_exp.add_argument(
        "--export_subset_freq",
        help="If set, exports subset frequencies",
        action="store_true",
    )
    var_exp.add_argument(
        "--export_all_sex_freq",
        help="If set, exports sex-specific frequencies for each pop",
        action="store_true",
    )

    gt_exp = parser.add_argument_group(
        "Export genotypes", description="Options related to exporting gnomAD genotypes"
    )
    gt_exp.add_argument(
        "--least_consequence",
        help="When exporting genotypes, includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: 3_prime_UTR_variant)",
        default="3_prime_UTR_variant",
    )
    gt_exp.add_argument(
        "--export_missing_genotypes",
        help="Export missing genotypes (missing genotypes export to a different (_missing) file.",
        action="store_true",
    )

    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit("Error: At least one of --exomes or --genomes must be specified.")

    if (
        not args.export_genotypes
        and not args.export_variants
        and not args.export_metadata
    ):
        sys.exit(
            "Error: At least one of --export_metadata, --export_variants or --export_genotypes must be specified"
        )

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
