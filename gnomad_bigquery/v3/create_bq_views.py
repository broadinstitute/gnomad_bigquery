import argparse

POPMAX_SQL = """
            (select struct(element.pop, element.ac, element.an, element.af, element.hom) from
            unnest(freq.list) as x
            where x.element.pop in ('nfe', 'eas', 'amr', 'afr', 'sas')
            and x.element.sex = 'all'
            and x.element.af  > 0
            order by element.af desc limit 1) as popmax
            """

def get_data_view(client: bigquery.Client, data_type: str, dataset: bigquery.DatasetReference) -> str:
    """
    This creates a view over one of the 'exomes' or 'genomes' data regrouping the
    variants, genotypes and sample metadata into a single view.

    :param CLient client: BQ client
    :param str data_type: One of 'exomes' or 'genomes'
    :param DatasetReference dataset: BQ Dataset
    :return: SQL for the view
    :rtype: str
    """

    # Exclude data_type from meta, so we can use * in main  query
    meta_table = client.get_table(dataset.table(f'{data_type}_meta'))
    variants_table = client.get_table(dataset.table(f'{data_type}_variants'))
    genotypes_table = client.get_table(dataset.table(f'{data_type}_genotypes'))
    meta_cols = [f.name for f in meta_table.schema if f.name != 'data_type']
    genotypes_cols = [f.name for f in genotypes_table.schema]
    variants_cols = [f.name for f in variants_table.schema]
    first_cols = [f"'{data_type}' as data_type", "chrom", "pos", "ref", "alt"]

    return f"""

    SELECT {",".join(first_cols)},
           {",".join([f'gt.{f}' for f in genotypes_cols if f != 'v'])}, 
           {",".join([f'v.{f}' for f in variants_cols if f not in first_cols])},
           {",".join([f'meta.{f}' for f in meta_cols if f not in genotypes_cols])},
            {POPMAX_SQL}
            
           FROM `{dataset.project}.{dataset.dataset_id}.{data_type}_variants` as v
    LEFT JOIN `{dataset.project}.{dataset.dataset_id}.{data_type}_genotypes` as gt ON v.idx = gt.v
    LEFT JOIN `{dataset.project}.{dataset.dataset_id}.{data_type}_meta` as meta ON gt.s = meta.s    
    """


def main(args):

    client = bigquery.Client()
    dataset = client.dataset(args.dataset)


    logger.info("Creating all_meta view")
    create_table(client,
                 dataset.table('all_meta'),
                 sql=create_union_query(client,
                                                  [client.get_table(dataset.table('genomes_meta')),
                                                   client.get_table(dataset.table('exomes_meta'))],
                                                  True
                                                  ),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_meta_table_desc()
                 )

    logger.info("Creating exomes view")
    create_table(client,
                 dataset.table('exomes'),
                 sql=get_data_view(client, 'exomes', dataset),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc('exomes')
                 )

    logger.info("Creating genomes view")
    create_table(client,
                 dataset.table('genomes'),
                 sql=get_data_view(client, 'genomes', dataset),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc('genomes')
                 )

    logger.info("Creating all view")
    create_table(client,
                 dataset.table('all'),
                 sql=create_union_query(client, [dataset.table('exomes'), dataset.table('genomes')]),
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc()
                 )

    logger.info("Creating all variants view")
    all_variants_sql = create_union_query(client, [dataset.table('exomes_variants'), dataset.table('genomes_variants')])
    all_variants_sql = all_variants_sql.replace('FROM', f', {POPMAX_SQL}, "exomes" as data_type from', 1)
    all_variants_sql = all_variants_sql.replace('FROM', f', {POPMAX_SQL}, "genomes" as data_type from', 1)
    create_table(client,
                 dataset.table('all_variants'),
                 sql=all_variants_sql,
                 view=True,
                 overwrite=args.overwrite,
                 description=get_data_view_desc()
                 )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--overwrite', help='If set, will overwrite all existing tables.', action='store_true')
    args = parser.parse_args()

    main(args)