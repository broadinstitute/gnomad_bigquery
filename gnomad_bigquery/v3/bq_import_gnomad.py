import argparse
import sys

from gnomad_bigquery.v3 import bq_import
from gnomad_bigquery.v3.bq_utils import logger
from gnomad_bigquery.v3.bq_table_descriptions import get_genotypes_table_desc, get_meta_table_desc, get_variants_table_desc


def main(args):

    def import_data(input_dir: str, data: str, description: str):
        parquet_files = f"{input_dir}_{data}.parquet/*.parquet"
        logger.info(f"Importing {data} from {parquet_files}")
        parser = bq_import.get_parser()
        imp_args = parser.parse_args(['--dataset', args.dataset,
                                      '--parquet_files', parquet_files,
                                      '--table', f'{data_type}_{data}',
                                      '--write_disposition',
                                      'WRITE_TRUNCATE' if args.overwrite else 'WRITE_EMPTY',
                                      '--description', description])
        print(imp_args)
        bq_import.main(imp_args)

    data_types = []
    if args.exomes:
        data_types.append("exomes")
    if args.genomes:
        data_types.append("genomes")

    version = args.version

    for data_type in data_types:

        input_dir = f"{args.input_dir}/gnomad_{data_type}_v{version}"

        if args.import_meta:
            import_data(input_dir, 'meta', get_meta_table_desc(data_type, version))

        if args.import_variants:
            import_data(input_dir, 'variants', get_variants_table_desc(data_type, version))

        if args.import_genotypes:
            import_data(input_dir, 'genotypes', get_genotypes_table_desc(data_type, version))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Will import exomes data. At least one of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Will import genomes data. At least one of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--version', help='Version of gnomAD to import.', type=int)
    parser.add_argument('--dataset', help='Dataset to create the table in. (default: gnomad)', default='gnomad')
    parser.add_argument('--import_meta', help='Imports samples metadata.', action='store_true')
    parser.add_argument('--import_variants', help='Imports variants.', action='store_true')
    parser.add_argument('--import_genotypes', help='Imports genotypes.', action='store_true')
    parser.add_argument('--input_dir', help='Input root directory (assumes data was produced with `bq_export.py`). default: gs://gnomad-tmp/bq', default='gs://gnomad-tmp/bq')
    parser.add_argument('--overwrite', help='If set, overwrites existing table.', action='store_true')
    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit("At least one of --exomes or --genomes needs to be specified.")

    main(args)
