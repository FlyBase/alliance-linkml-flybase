# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase gene expression for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_expression.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE]

Example:
    python AGR_data_retrieval_curation_expression.py -v -t -c /path/to/config.cfg
    -l v1.1.2

Notes:
    This script exports FlyBase gene expression data as a JSON file conforming to the
    Expression LinkML specs (GeneExpressionExperiment) for the Alliance persistent
    curation database. FlyBase feature_expression annotations are grouped into
    GeneExpressionExperiment objects (by assayed gene, reference and assay), each
    holding one to many GeneExpressionAnnotationDTOs. A chado database is required.

"""

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from expression_handler import ExpressionHandler
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
REPORT_LABEL = 'expression_curation'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
output_filename = set_up_dict['output_filename'].replace('tsv', 'json')
log = set_up_dict['log']
testing = set_up_dict['testing']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(
    description='Export FlyBase gene expression data to Alliance LinkML JSON.',
    epilog="""
Environment variables:
  SERVER              Database server (e.g. flysql25)
  DATABASE            Database name (e.g. production_chado)
""",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('-l', '--linkml_release', help='The "agr_curation_schema" LinkML release number.', required=True)

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
# insp = inspect(engine)    # I always have this line, but I do not know what it does.
Session = sessionmaker(bind=engine)
session = Session()


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase gene expression data."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    expression_handler = ExpressionHandler(log, testing)
    export_chado_data(session, log, expression_handler)

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    export_dict[expression_handler.primary_export_set] = expression_handler.export_data[expression_handler.primary_export_set]
    if len(export_dict[expression_handler.primary_export_set]) == 0:
        log.error('The "gene_expression_experiment_ingest_set" is unexpectedly empty.')
        raise ValueError('The "gene_expression_experiment_ingest_set" is unexpectedly empty.')
    else:
        generate_export_file(export_dict, log, output_filename)

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
