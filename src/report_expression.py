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
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
REPORT_LABEL = 'curated_expression'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
input_dir = set_up_dict['input_dir']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']
testing = set_up_dict['testing']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
# insp = inspect(engine)    # I always have this line, but I do not know what it does.
Session = sessionmaker(bind=engine)
session = Session()


# The main process.
def main():
    """Run the steps for exporting feature_expression data to a TSV file."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')

    # Get the data and process it.
    expression_handler = ExpressionHandler(log, testing)
    expression_handler.get_general_data(session)
    expression_handler.get_datatype_data(session)
    expression_handler.synthesize_info()

    # Export the data.
    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
