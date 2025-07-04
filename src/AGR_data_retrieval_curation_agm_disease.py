# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Conversion of FlyBase allele disease annotations into AGM Alliance annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_agm_disease.py [-h]
    [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG] [-t TESTING]

Example:
    python AGR_data_retrieval_curation_agm_disease.py -v -l v2.9.0
    -c /path/to/config.cfg

Notes:
    This script exports FlyBase allele-based disease annotations as a JSON file
    conforming to the AGM LinkML specs for the Alliance persistent curation
    database. A chado database with a full "audit_chado" table is required.

"""

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from disease_handlers import AGMDiseaseHandler
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
REPORT_LABEL = 'agm_disease_curation'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
reference_assembly = set_up_dict['assembly']
input_dir = set_up_dict['input_dir']
output_filename = set_up_dict['output_filename'].replace('tsv', 'json')
log = set_up_dict['log']
testing = set_up_dict['testing']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
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
    """Run the steps for exporting LinkML-compliant FlyBase AGM disease annotations."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    agm_disease_handler = AGMDiseaseHandler(log, testing)
    export_chado_data(session, log, agm_disease_handler, testing=testing)

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    export_dict[agm_disease_handler.primary_export_set] = agm_disease_handler.export_data[agm_disease_handler.primary_export_set]
    generate_export_file(export_dict, log, output_filename)

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
