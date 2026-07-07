# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase antibody data for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_antibody.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE]

Example:
    python AGR_data_retrieval_curation_antibody.py -v -t -c /path/to/config.cfg
    -l v1.1.2

Notes:
    This script exports FlyBase antibody data as a JSON file conforming to the
    Antibody LinkML specs for the Alliance persistent curation database.
    FlyBase has no first class antibody entities; antibodies are synthesized from
    gene-level antibody data: "reported_antibod_gen" featureprops (lab-generated),
    and DSHB/Cell Signaling Technology feature_dbxrefs (commercial).

"""

import argparse
from os import environ

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from antibody_handler import AntibodyHandler
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
REPORT_LABEL = 'antibody_curation'

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

output_filename = environ.get('ALT_OUTPUT', output_filename)

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(
    description='Export FlyBase antibody data to Alliance LinkML JSON.',
    epilog="""
Environment variables:
  SERVER              Database server (e.g. flysql25)
  DATABASE            Database name (e.g. production_chado)
  SQL_PORT            Database port (default: 5432)
  ALT_OUTPUT          Override default output file path
""",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('-l', '--linkml_release',
                    help='The "agr_curation_schema" LinkML release number.', required=True)

# Use parse_known_args(), not parse_args(),
# to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info(f'Parsing args specific to this script; ignoring these: {extra_args}')
linkml_release = args.linkml_release

port = environ.get('SQL_PORT', '5432')

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + ':' + port + '/' + database

print(f"Connecting to server:{server} port:{port} database:{database} username:{username}")

engine = create_engine(engine_var_rep)
Session = sessionmaker(bind=engine)
session = Session()


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase antibody data."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    antibody_handler = AntibodyHandler(log, testing)
    export_chado_data(session, log, antibody_handler)

    # Export the antibody data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    set_name = antibody_handler.primary_export_set
    export_dict[set_name] = antibody_handler.export_data[set_name]
    if len(export_dict[set_name]) == 0:
        log.error(f'The "{set_name}" is unexpectedly empty.')
        raise ValueError(f'The "{set_name}" is unexpectedly empty.')
    generate_export_file(export_dict, log, output_filename)

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
