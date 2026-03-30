# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase gene group data for Alliance curation database.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_gene_group.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python AGR_data_retrieval_curation_gene_group.py -v -t -c /path/to/config.cfg
    -l v1.1.2
    -r fb_2024_06_reporting

Notes:
    This script exports FlyBase gene group (FBgg) data as a JSON file conforming to the
    FunctionalGeneSet LinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
from os import environ

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from gene_group_handler import GeneGroupHandler
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
REPORT_LABEL = 'gene_group_curation'
REPORT_TITLE = 'FlyBase Gene Group Report'
TSV_HEADERS = [
    'gene_group_id',
    'symbol',
    'full_name',
    'synonyms',
    'description',
    'go_molecular_function',
    'go_biological_process',
    'go_cellular_component',
    'parent_groups',
    'related_groups',
    'gene_members',
]

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
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-l', '--linkml_release',
                    help='The "agr_curation_schema" LinkML release number.', required=True)
parser.add_argument('-r', '--reference_db',
                    help='The name of a previous reference db for incremental exports.',
                    required=False)

# Use parse_known_args(), not parse_args(),
# to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()

log.info(f'Parsing args specific to this script; ignoring these: {extra_args}')
linkml_release = args.linkml_release
reference_db = args.reference_db

port = environ.get('SQL_PORT', '5432')

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + ':' + port + '/' + database

print(f"Connecting to server:{server} port:{port} database:{database} username:{username}")

engine = create_engine(engine_var_rep)
Session = sessionmaker(bind=engine)
session = Session()

# Create a session to the reference db.
if reference_db:
    engine_var_ref = 'postgresql://' + username + ":" + password + '@' + 'flysql23' + '/' + reference_db
    ref_engine = create_engine(engine_var_ref)
    RefSession = sessionmaker(bind=ref_engine)
    reference_session = RefSession()
else:
    reference_session = None


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase gene group data."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    gg_handler = GeneGroupHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, gg_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, gg_handler)

    # Export the gene group data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    set_name = gg_handler.primary_export_set
    export_dict[set_name] = gg_handler.export_data[set_name]
    if len(export_dict[set_name]) == 0:
        if reference_session:
            log.info('No updates to report.')
        else:
            log.error(f'The "{set_name}" is unexpectedly empty.')
            raise ValueError(f'The "{set_name}" is unexpectedly empty.')
    else:
        generate_export_file(export_dict, log, output_filename)

    # Export the gene member associations to a separate file.
    if not reference_session:
        assoc_set_name = 'gene_functional_gene_set_association_ingest_set'
        if assoc_set_name in gg_handler.export_data:
            association_export_dict = {
                'linkml_version': linkml_release,
                'alliance_member_release_version': database_release,
            }
            association_export_dict[assoc_set_name] = gg_handler.export_data[assoc_set_name]
            if len(association_export_dict[assoc_set_name]) > 0:
                association_output_filename = output_filename.replace(
                    'gene_group', 'gene_functional_gene_set_association')
                generate_export_file(association_export_dict, log, association_output_filename)
            else:
                log.warning(f'The "{assoc_set_name}" is empty.')

    # Export the gene group TSV report.
    if not reference_session:
        gg_handler.process_for_tsv_export()
        tsv_data = generic_FB_tsv_dict(REPORT_TITLE, database)
        tsv_data['data'] = gg_handler.export_data_for_tsv
        tsv_output_filename = output_filename.replace('.json', '.tsv')
        if tsv_data['data']:
            tsv_report_dump(tsv_data, tsv_output_filename, headers=TSV_HEADERS)
            log.info(f'Wrote TSV to {tsv_output_filename}')
        else:
            log.warning('No gene group data for TSV export.')

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
