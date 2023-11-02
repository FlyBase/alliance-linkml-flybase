# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase PLACEHOLDER_DATA_CLASS for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_template.py [-h] [-r FLYBASE_RELEASE] [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_curation_template.py -v -r 2023_05 -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script exports FlyBase PLACEHOLDER_DATA_CLASS data as a JSON file conforming to the
    PLACEHOLDER_DATA_CLASS LinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
import datetime
import strict_rfc3339
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, OrganismDbxref, Pub, PubDbxref, Strain, StrainDbxref,
    StrainPub, StrainSynonym, Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading
from utils import DataHandler, db_query_transaction, generate_export_file

# Now proceed with generic setup.
report_label = 'PLACEHOLDER_DATA_CLASS'
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
reference_assembly = set_up_dict['assembly']
input_dir = set_up_dict['input_dir']
output_filename = set_up_dict['output_filename'].replace('tsv', 'json')
log = set_up_dict['log']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-l', '--linkml_release', help='The "agr_curation_schema" LinkML release number.', required=True)
parser.add_argument('-r', '--fb_release', help='The FlyBase data release from which data was obtained.', required=True)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release
fb_release = args.fb_release

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
# insp = inspect(engine)    # I always have this line, but I do not know what it does.
Session = sessionmaker(bind=engine)
session = Session()


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase PLACEHOLDER_DATA_CLASS."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {fb_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    gene_handler = GeneHandler(log, 'gene', 'GENE')
    db_query_transaction(session, log, gene_handler)
    gene_handler.export_data = [{'bob': 'cool'}]

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': fb_release,
    }
    export_dict[gene_handler.agr_data_type] = gene_handler.export_data
    generate_export_file(export_dict, log, output_filename)

    log.info('Ended main function.\n')


class GeneHandler(DataHandler):
    """This object gets strains, synthesizes/filters the data, then exports it as LinkML JSON."""


if __name__ == "__main__":
    main()
