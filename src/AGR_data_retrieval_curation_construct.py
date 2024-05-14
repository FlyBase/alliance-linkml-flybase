# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase construct for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_construct.py [-h]
    [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_curation_construct.py -v -l v1.1.2
    -c /path/to/config.cfg

Notes:
    This script exports FlyBase construct data as a JSON file conforming to the
    construct LinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.psycopg_functions import set_up_db_reading
from utils import get_handler, db_query_transaction, generate_export_file

# Data types handled by this script.
FB_DATA_TYPE = 'construct'
REPORT_LABEL = 'construct_curation'

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
    """Run the steps for exporting LinkML-compliant FlyBase AGM."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    cons_handler = get_handler(log, FB_DATA_TYPE, testing)
    db_query_transaction(session, log, cons_handler)

    # Export the construct data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    export_dict['construct_ingest_set'] = cons_handler.export_data['construct_ingest_set']
    generate_export_file(export_dict, log, output_filename)

    # Export the construct associations to a separate file.
    association_output_filename = output_filename.replace('construct', 'construct_association')
    association_export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    association_export_dict['construct_genomic_entity_association_ingest_set'] = cons_handler.export_data['construct_genomic_entity_association_ingest_set']
    generate_export_file(association_export_dict, log, association_output_filename)

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
