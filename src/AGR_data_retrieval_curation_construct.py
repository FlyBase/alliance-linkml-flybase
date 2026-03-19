# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase construct for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_construct.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python AGR_data_retrieval_curation_construct.py -v -t -c /path/to/config.cfg
    -l v1.1.2
    -r fb_2024_06_reporting

Notes:
    This script exports FlyBase construct data as a JSON file conforming to the
    construct LinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
import os
from os import getenv

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from construct_handler import ConstructHandler
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
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
parser.add_argument('-r', '--reference_db', help='The name of a previous reference db for incremental exports.', required=False)

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release
reference_db = args.reference_db

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
# insp = inspect(engine)    # I always have this line, but I do not know what it does.
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


def generate_association_tsv_file(export_dict, ingest_name, filename):
    """Generate a TSV file for an association ingest set."""
    first_entity = 'construct_identifier'
    if ingest_name == 'construct_cassette_association_ingest_set':
        second_entity = 'cassette_identifier'
    else:
        second_entity = 'genomic_entity_identifier'
    with open(filename, 'w') as outfile:
        outfile.write(f"#{first_entity}\tRelationship\t{second_entity}\tEvidence\n")
        for entity_dict in export_dict[ingest_name]:
            sub = entity_dict[first_entity]
            obj = entity_dict[second_entity]
            rel_type = entity_dict['relation_name']
            if 'evidence_curies' in entity_dict:
                pubs = "|".join(entity_dict['evidence_curies'])
            else:
                pubs = ""
            outfile.write(f"{sub}\t{rel_type}\t{obj}\t{pubs}\n")


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase AGM."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    cons_handler = ConstructHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, cons_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, cons_handler)

    # Export the construct data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    export_dict['construct_ingest_set'] = cons_handler.export_data['construct_ingest_set']
    if len(export_dict['construct_ingest_set']) == 0:
        if reference_session:
            log.info('No updates to report.')
        else:
            log.error('The "construct_ingest_set" is unexpectedly empty.')
            raise ValueError('The "construct_ingest_set" is unexpectedly empty.')
    else:
        generate_export_file(export_dict, log, output_filename)

    if not reference_session:
        # Export the construct associations to a separate file.
        association_output_filename = output_filename.replace('construct', 'construct_association')
        association_export_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
        }
        association_export_dict['construct_genomic_entity_association_ingest_set'] = cons_handler.export_data['construct_genomic_entity_association_ingest_set']
        if len(association_export_dict['construct_genomic_entity_association_ingest_set']) == 0:
            log.error('The "construct_genomic_entity_association_ingest_set" is unexpectedly empty.')
            raise ValueError('The "construct_genomic_entity_association_ingest_set" is unexpectedly empty.')

        # Because the Alliance is not yet abe to handle cassettes we do not want to add these
        # associations. For testing set the env ADD_CASS_TO_CONSTRUCT which will then do this
        dump_cass_assoc = getenv('ADD_CASS_TO_CONSTRUCT', None)
        if dump_cass_assoc and dump_cass_assoc == 'YES':
            association_export_dict['construct_cassette_association_ingest_set'] = \
                cons_handler.export_data['construct_cassette_association_ingest_set']
        else:
            log.warning('The ADD_CASS_TO_CONSTRUCT environment variable is not set to "YES". '
                        'So no assoc to cassettes added.')

        generate_export_file(association_export_dict, log, association_output_filename)

        # Generate TSV files for each association type.
        tsv_dir = os.path.dirname(set_up_dict['output_filename'])
        for ingest_name in association_export_dict:
            if not ingest_name.endswith('_ingest_set'):
                continue
            if len(association_export_dict[ingest_name]) == 0:
                continue
            tsv_filename = os.path.join(tsv_dir, f"{ingest_name.replace('_ingest_set', '')}.tsv")
            try:
                generate_association_tsv_file(association_export_dict, ingest_name, tsv_filename)
                log.info(f'Generated TSV: {tsv_filename}')
            except KeyError as e:
                log.error(f'TSV generation for "{ingest_name}" failed: KeyError {e}')

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
