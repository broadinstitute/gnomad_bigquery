import argparse
import sys
from typing import Optional

import hail as hl
from gnomad.resources.grch37.reference_data import decoy_intervals as decoy_intervals_37
from gnomad.resources.grch37.reference_data import lcr_intervals as lcr_intervals_37
from gnomad.resources.grch37.reference_data import seg_dup_intervals as seg_dup_intervals_37
from gnomad.resources.grch38.reference_data import lcr_intervals as lcr_intervals_38
from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import try_slack
from gnomad.utils.vep import CSQ_ORDER

from gnomad_bigquery.bq_utils import (GNOMAD_VERSIONS, export_ht_for_bq,
                                      get_gnomad_v3_data_for_bq, logger)
from gnomad_qc.v2.resources.annotations import frequencies as v2_freq
from gnomad_qc.v2.resources.annotations import rf as v2_rf
from gnomad_qc.v2.resources.annotations import vep as v2_vep
from gnomad_qc.v2.resources.basics import get_gnomad_data, get_gnomad_meta
from gnomad_qc.v3.resources.annotations import freq as v3_freq
from gnomad_qc.v3.resources.annotations import vep as v3_vep
from gnomad_qc.v3.resources.meta import meta as v3_meta


def export_genotypes(
    data_type: str,
    version: int,
    output_dir: str,
    max_freq: Optional[float] = None,
    least_consequence: str = None,
    variant_index: bool = True,
    ) -> None:
    """
    Export gnomAD genotypes, present or missing, within variants where
    there is a transcript consequence that is at least the least_consequence passed
    :param data_type: exomes or genomes
    :param version: Version of gnomAD
    :param output_dir: Directory where to write parquet file
    :param max_freq: Maximum AF to allow in export
    :param least_consequence: Least consequence admitted, default is `3_prime_UTR_variant`
    :param variant_index: Index variant to save on storage cost
    :return:
    """
    if version == 2:
        mt = get_gnomad_data(data_type, non_refs_only=True)
        vep = v2_vep(data_type).ht()
        freq = v2_freq(data_type).ht()
    elif version == 3 and data_type == "genomes":
        mt = get_gnomad_v3_data_for_bq(
            key_by_locus_and_alleles=True,
            remove_hard_filtered_samples=False,
            adjust_sex_and_annotate_adj=True,
        )
        vep = v3_vep.ht()
        freq = v3_freq.ht()
    else:
        raise DataException(f"There is no version {version} of gnomAD {data_type}. Please choose one from {GNOMAD_VERSIONS}")

    mt = mt.select_cols().select_rows().add_row_index()

    if least_consequence is not None:
        vep_consequences = hl.literal(set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1]))
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
    logger.debug(
        f"Found {vep.count()} variants with a VEP coding transcript and a consequence worst or equal to {least_consequence}."
        )

    select_expr = hl.is_defined(vep[mt.row_key])

    if max_freq is not None:
        freq = freq.filter(freq.freq[0].AF < max_freq)
        freq = freq.persist()
        logger.debug(f"Found {freq.count()} variants with AF < {max_freq}")
        select_expr = select_expr & hl.is_defined(freq[mt.row_key])  #TODO: number of variants in gnomAD without vep, is there a pattern

    mt = mt.filter_rows(select_expr)
    ht = mt.entries()
    logger.debug(
        f"Found {ht.count()} variants with a VEP coding transcript and a consequence worst or equal to {least_consequence}."
        )
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
            "pgt1_non_ref": ht.PGT[0] != 0,
            "pgt2_non_ref": ht.PGT[1] != 0,
        }
    )
    ht = ht.select(**select_expr)
    ht.to_spark().write.parquet(f'{output_dir}/gnomad_v{version}_{data_type}_genotypes.parquet', mode='overwrite')


def export_variants(data_type: str, version: int, export_all_sex_freq: bool, export_filters: bool, output_dir: str, add_vep: bool) -> None:
    """
    Exports all variants in gnomAD with their frequencies and optional vep transcript consequences
    :param data_type: exomes or genomes
    :param version: Version of gnomAD
    :param export_all_sex_freq: Whether to export sex frequencies alongside overall frequencies
    :param export_filters: Whether to add variant filters such as random forest annotations and vqsr 
    :param output_dir: Directory where to write parquet file
    :param add_vep: Whether to annotate with vep's transcript_consequences and most_severe_consequence
    :return:
    """
    def format_freq(ht: hl.Table):
        def filter_freqs(x: hl.expr.ArrayExpression):
            if export_all_sex_freq:
                return ~x[0].contains('platform') & ~x[0].contains('downsampling')
            else:
                return ~x[0].contains('platform') & ~x[0].contains('downsampling') & ~(x[0].contains('sex') & x[0].contains('pop'))

        def get_freq_struct(x: hl.expr.ArrayExpression):
            freq_struct = hl.struct(
                adj=x[0]['group'] == 'adj',
                pop=x[0].get('subpop', x[0].get('pop', 'all')),
                sex=x[0].get('sex', 'all'),
                ac=x[1].AC,
                an=x[1].AN,
                af=x[1].AF,
                hom=x[1].homozygote_count
            )
            return freq_struct

        return  ht.annotate(freq=hl.zip(freq_meta, freq[ht.key].freq)
                            .filter(filter_freqs)
                            .map(get_freq_struct))

    if version == 2:
        mt = get_gnomad_data(data_type, non_refs_only=True)
        vep = v2_vep(data_type).ht()
        freq = v2_freq(data_type).ht()
        rf = v2_rf(data_type).ht()
        lcr = lcr_intervals_37.ht()
        decoy = decoy_intervals_37.ht()
        segdup = seg_dup_intervals_37.ht()
    elif version == 3 and data_type == "genomes":
        mt = get_gnomad_v3_data_for_bq(
            key_by_locus_and_alleles=True,
            remove_hard_filtered_samples=False,
            adjust_sex_and_annotate_adj=True,
        )
        vep = v3_vep.ht()
        freq = v3_freq.ht()
        lcr = lcr_intervals_38.ht()
        vqsr = hl.read_table("gs://gnomad/annotations/hail-0.2/ht/genomes_v3/filtering/vqsr_alleleSpecificTrans.finalized.split.ht") #TODO: Update once variant QC is in gnomad_qc repo
    else:
        raise DataException(f"There is no version {version} of gnomAD {data_type}. Please choose one from {GNOMAD_VERSIONS}")

    ht = mt.rows()
    ht = ht.select().add_index()
    freq_meta = hl.literal(freq.freq_meta.collect()[0])
    ht = format_freq(ht)

    if add_vep:
        vep = vep.select(
            transcript_consequences=vep.vep.transcript_consequences,
            most_severe_consequence=vep.vep.most_severe_consequence
        )
        ht = ht.annotate(**vep[ht.key])

    ht = ht.annotate(nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar()))
    ht = ht.annotate(lcr=hl.is_defined(lcr[ht.locus]))

    if version == 2:
        ht = ht.annotate(decoy=hl.is_defined(decoy[ht.locus]), segdup=hl.is_defined(segdup[ht.locus]))
        if export_filters:
            rf_expr = {f: rf[ht.key][f] for f in rf.row_value if not f.endswith('rank')}
            ht = ht.annotate(**rf_expr)

    if version == 3 and export_filters:
        ht = ht.annotate(vqsr_filters=vqsr[ht.key].filters)

    export_ht_for_bq(ht, f'{output_dir}/gnomad_v{version}_{data_type}_variants.parquet')


def main(args):

    hl.init(log='/bq.log')
    data_types = []
    version = args.version

    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:

        if args.export_metadata:
            if version == 2:
                meta = get_gnomad_meta(data_type=data_type, full_meta=True)
            elif version == 3 and data_type == "genomes":
                meta = v3_meta.ht()
            else:
                raise DataException(f"There is no version {version} of gnomAD {data_type}. Please choose one from {GNOMAD_VERSIONS}")
            export_ht_for_bq(meta, f'{args.output_dir}/gnomad_v{version}_{data_type}_meta.parquet')

        if args.export_genotypes:
            export_genotypes(data_type,
                            version,
                            args.output_dir,
                            args.max_freq,
                            args.least_consequence,
                            )

        if args.export_variants:
            export_variants(data_type, version, args.export_all_sex_freq, True, args.output_dir, True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--version', help="gnomAD version", choices=GNOMAD_VERSIONS, type=int, required=True)
    parser.add_argument('--exomes', help='Run on exomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--export_metadata', help='Export samples metadata', action='store_true')
    parser.add_argument('--export_genotypes', help='Export non-ref genotypes', action='store_true')
    parser.add_argument('--export_variants', help='Export variants', action='store_true')
    parser.add_argument('--output_dir', help='Output root. default: gs://gnomad-tmp/bq', default='gs://gnomad-tmp/bq')

    var_exp = parser.add_argument_group('Export variants', description='Options related to exporting gnomAD variants')
    var_exp.add_argument('--export_all_sex_freq', help='If set, exports sex-specific frequencies for each pop', action='store_true')

    gt_exp = parser.add_argument_group('Export genotypes', description='Options related to exporting gnomAD genotypes')
    gt_exp.add_argument('--max_freq', help='If specified, maximum global adj AF for genotypes table to emit. (default: 0.02)', default=0.02, type=float)
    gt_exp.add_argument('--least_consequence', help='When exporting genotypes, includes all variants for which the worst_consequence is at least as bad as the specified consequence. The order is taken from gnomad_hail.constants. (default: 3_prime_UTR_variant)',
                        default='3_prime_UTR_variant')

    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if args.exomes and args.version == 3:
        sys.exit('gnomAD v3 does not contain exomes')

    if not args.export_genotypes and not args.export_variants and not args.export_metadata and not args.export_transcripts:
        sys.exit('Error: At least one of --export_metadata, --export_variants or --export_genotypes or --export_transcripts must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
