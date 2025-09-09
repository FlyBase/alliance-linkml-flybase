# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Report curated expression as a bulk TSV file.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_xprn.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python report_xprn.py -v -t -c /path/to/config.cfg

Notes:
    This script generates a bulk file for feature_expression chado data.

"""

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from expression_handler import ExpressionHandler
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)

# Data types handled by this script.
REPORT_LABEL = 'curated_expression'
REPORT_TITLE = 'FlyBase Expression Report'
HEADER_LIST = [
    'feature_id',
    'feature_symbol',
    'reference_id',
    'expression_type',
    'expression_id',    # BOB: Internal expression.expression_id, only for debugging.
    'assay_term',
    'stage_start',
    'stage_end',
    'stage_qualifiers',
    'stage_slim_terms',
    'anatomical_structure_term',
    'anatomical_structure_qualifiers',
    'anatomical_structure_slim_terms',
    'anatomical_substructure_term',
    'anatomical_substructure_qualifiers',
    'anatomical_substructure_slim_terms',
    'cellular_component_term',
    'cellular_component_qualifiers',
    'notes',
]

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
SERVER = set_up_dict['server']
DATABASE = set_up_dict['database']
USERNAME = set_up_dict['username']
PASSWORD = set_up_dict['password']
DATABASE_RELEASE = set_up_dict['database_release']
OUTPUT_FILENAME = set_up_dict['output_filename']
TESTING = set_up_dict['testing']
log = set_up_dict['log']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + USERNAME + ":" + PASSWORD + '@' + SERVER + '/' + DATABASE
engine = create_engine(engine_var_rep)
# insp = inspect(engine)    # I always have this line, but I do not know what it does.
Session = sessionmaker(bind=engine)
session = Session()


# The main process.
def main():
    """Run the steps for exporting feature_expression data to a TSV file."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {DATABASE_RELEASE}')

    # Get the data and process it.
    expression_handler = ExpressionHandler(log, TESTING)
    expression_handler.get_general_data(session)
    expression_handler.get_datatype_data(session)
    expression_handler.synthesize_info(session)
    expression_handler.process_for_tsv_export()

    # Export the data.
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['data'] = expression_handler.export_data_for_tsv
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
