import logging
import sys
from typing import Dict, List

import hail as hl
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import get_adj_expr
from google.cloud import bigquery

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("bq")
logger.setLevel(logging.INFO)

GNOMAD_VERSIONS = [2, 3]


def drop_columns(client: bigquery.Client, table: bigquery.TableReference, columns_to_exclude: List[str]):
    t = client.get_table(table)
    # Doesn't work with views
    if t.view_query:
        sys.exit("Cannot drop column(s) from view.")

    cols_to_keep = [f.name for f in t.schema if f.name not in columns_to_exclude]
    create_table(client,
                 table,
                 sql=f"select {','.join(cols_to_keep)} from `{table.project}.{table.dataset_id}.{table.table_id}`",
                 overwrite=True,
                 description=t.description)


def create_table(client: bigquery.Client,
                 destination_table: bigquery.TableReference,
                 sql: str,
                 overwrite: bool = False,
                 use_cache: bool = True,
                 view: bool = False,
                 description: str = None) -> None:
    dataset_ref = client.dataset(destination_table.dataset_id)
    if destination_table in [t.reference for t in list(client.list_tables(dataset_ref))]:
        if not overwrite:
            logger.info(f"Table {destination_table.path} exists; using existing version. Use --overwrite to re-generate.")
            return
        elif view:
            client.delete_table(destination_table)

    if view:
        table = bigquery.Table(destination_table)
        table.view_query = sql
        if description is not None:
            table.description = description
        client.create_table(table)
        logger.info(f"View {destination_table.path} created.")
    else:
        job_config = bigquery.QueryJobConfig()
        job_config.destination = destination_table
        job_config.use_query_cache = use_cache
        job_config.allow_large_results = True
        job_config.write_disposition = 'WRITE_TRUNCATE'
        job_config.use_legacy_sql = False

        query_job = client.query(
            sql,
            job_config=job_config
        )
        query_job.result()
        logger.info(f'{client.get_table(destination_table).num_rows} query results loaded to table {destination_table.path}')
        if description is not None:
            table = bigquery.Table(destination_table)
            table.description = description
            client.update_table(table, ['description'])


TYPE_CONVERSION = {
   ('INTEGER', 'FLOAT'): 'INTEGER',
   ('FLOAT', 'INTEGER'): 'INTEGER',
   ('STRING', 'INTEGER'): 'STRING',
   ('INTEGER', 'STRING'): 'STRING',
   ('STRING', 'FLOAT'): 'STRING',
   ('FLOAT', 'STRING'): 'STRING'
}


def create_union_query(client: bigquery.Client, tables: List[bigquery.TableReference], convert_types: bool = True):

    def get_field_sql(schema: Dict[str, str], field: str, field_type: str):
        if field not in schema:
            return f'null as {field}'
        elif field_type != schema[field]:
            return f'cast({field} as {TYPE_CONVERSION[field_type, schema[field]]}) as {field}'
        return f'{field} as {field}'

    schemas = [(f'`{t.project}.{t.dataset_id}.{t.table_id}`', {f.name: f.field_type for f in client.get_table(t).schema}) for t in tables]

    unified_schema = {}

    for schema in schemas:
        for name, field_type in schema[1].items():
            if name in unified_schema:
                if unified_schema[name] != field_type:
                    if convert_types and (field_type, unified_schema[name]) in TYPE_CONVERSION:
                        unified_schema[name] = TYPE_CONVERSION[(field_type, unified_schema[name])]
                        print(f"Found different types for field {name}. Resulting column type will be {unified_schema[name]}.")
                    else:
                        sys.exit(f"Cannot merge schemas: field {name} has different type ({unified_schema[name]} and {field_type})")
            else:
                unified_schema[name] = field_type

    return " UNION ALL \n".join(['SELECT ' + ", \n".join([get_field_sql(s, f, ft) for f, ft in unified_schema.items()]) + f' FROM {table}\n' for table, s in schemas])


def export_ht_for_bq(ht: hl.Table, destination: str, lower: bool = True, mode: str = 'overwrite'):

    # Special case locus and alleles
    ht = ht.key_by()
    if 'locus' in ht.row:
        ht = ht.transmute(
            chrom=ht.locus.contig,
            pos=ht.locus.position
        )

    if 'alleles' in ht.row:
        ht = ht.transmute(
            ref=ht.alleles[0],
            alt=hl.delimit(ht.alleles[1:], ",")
        )

    # Flatten all structs
    while [c for c in ht.row.items() if isinstance(c[1].dtype, hl.tstruct)]:
        ht = ht.flatten()

    # Rename any field with '.' (and lower case if set)
    names = {name: name.replace('.', '_') for name in ht.row}
    if lower:
        names = {old_name: new_name.lower() for old_name, new_name in names.items()}

    ht = ht.rename(names)

    ht.to_spark().write.parquet(destination, mode=mode)
    logger.info(f"Successfully exported {destination}")


def get_gnomad_v3_data_for_bq(key_by_locus_and_alleles=True, remove_hard_filtered_samples=False, adjust_sex_and_annotate_adj=False) -> hl.MatrixTable:
    """
    This script preps gnomAD v3 data for bigquery export by splitting a sparse MT,
    removing all non_ref entries, and filtering out monomorphic sites. Optionally,
    it will adjust sex ploidy and annotate adj.
    :param key_by_locus_and_alleles: Key the rows of the sparse MT by locus and alleles
    :param remove_hard_filtered_samples: Removes hard filtered samples
    :param adjust_sex_and_annotate_adj: Adjusts sex ploidy and annotates GT with adj
    :return: MatrixTable
    """
    from gnomad_qc.v3.resources import get_gnomad_v3_mt, meta

    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=key_by_locus_and_alleles, remove_hard_filtered_samples=remove_hard_filtered_samples)
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_entries(mt.GT.is_non_ref())

    if adjust_sex_and_annotate_adj:
        mt = mt.annotate_cols(**meta.ht()[mt.col_key])
        mt = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex),
            adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD)
        )
    return mt
