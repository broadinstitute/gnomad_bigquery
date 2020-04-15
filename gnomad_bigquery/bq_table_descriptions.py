from gnomad.utils.vep import CSQ_ORDER

LEAST_CSQ = '3_prime_UTR_variant'

def get_variants_table_desc(version: int, data_type: str, split=False, in_gt=False):
    dtype = '' if data_type is None else f' {data_type}'
    if split:
        if in_gt:
            return f"""This table contains gnomad v{version} {dtype} variants that are filtered similarly to the
             genotypes table, optionally filtered by a least consequence term and max AF. 
            Notes:
            * All variants were split, so multi-allelic sites are represented as multiple rows (one per non-ref allele)
            * Only VEP transcript annotations are available (non-coding annotations aren't)
            """
        else:
            return f"""This table contains gnomad v{version} {dtype} variants that are filtered out of the
             genotypes table by not meeting the threshold for the least consequence term and max AF. 
            Notes:
            * All variants were split, so multi-allelic sites are represented as multiple rows (one per non-ref allele)
            * Only VEP transcript annotations are available (non-coding annotations aren't)
            """
    else:
        return f"""This table contains all gnomad v{version} {dtype} variants.
        Notes:
        * All variants were split, so multi-allelic sites are represented as multiple rows (one per non-ref allele)
        * Only VEP transcript annotations are available (non-coding annotations aren't)
        """


def get_meta_table_desc(version: int, data_type: str = None):
    dtype = '' if data_type is None else f' {data_type}'

    return f"""This table contains all gnomad v{version} {dtype} samples metadata.
    """


def get_genotypes_table_desc(version: int, data_type: str = None):
    dtype = '' if data_type is None else f' {data_type}'
    all_str = 'all' if data_type is None else data_type
    return f"""This table contains all gnomad v{version} {dtype} non-reference genotypes for variants meeting the following criteria:
    1. Variant VEP annotation is at least `{LEAST_CSQ}` (i.e. one of: {','.join(CSQ_ORDER[0:CSQ_ORDER.index(LEAST_CSQ) + 1])})
    
    A few notes about data representation:
    * All variants were split beforehand, so all genotypes are bi-allelic.
    * To save storage space, variants are referred using a integer index. To get variants information, including chrom, pos, ref, alt, this table needs to be joined with the{dtype}_variants table. (The {all_str} view provides an easy representation of this join)
    * If `is_het` if `True` then the sample is het at the site, if it is `False` the sample is homozygous non-reference (hom-ref samples are  excluded from this table)
    * `ad0` and `ad1` store the allelic-depth for the reference and non-reference alleles resp.
    * `pl0`, `pl1` and `pl2` are the genotype likelihoods for the hom-ref, het and hom-var genotypes resp.
    """


def get_data_view_desc(version: int, data_type: str = None, split_variant_table=False):
    if data_type is not None:
        view_tables_str = f'{data_type}_variants{"_in_gt" if split_variant_table else ""}, {data_type}_genotypes and {data_type}_meta tables'
        n_tables_str = '3 tables.'
    else:
        view_tables_str = 'exomes and genomes views.'
        n_tables_str = '2 views'

    return f"""This view contains a flat version of joining the gnomAD v{version} {view_tables_str}.
    To get more information about the content of this view, take a look a the description of those {n_tables_str}."""


def get_variants_view_desc(version: int, data_type: str):
    dtype = '' if data_type is None else f' {data_type}'
    return f"""This table is a VARIANTS ONLY view of all gnomad v{version} {dtype} variants, 
    a union of variants in the genotypes table and variants not in the genotypes table.
        Notes:
        * All variants were split, so multi-allelic sites are represented as multiple rows (one per non-ref allele)
        * Only VEP transcript annotations are available (non-coding annotations aren't)
        """