# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase sequence targeting reagents (STR) for Alliance curation database.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_str.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python AGR_data_retrieval_curation_str.py -v -t -c /path/to/config.cfg
    -l v2.17.0
    -r fb_2026_01_reporting
Notes:
    This script exports the FlyBase sequence targeting reagents - the subset of FBsf
    sequence features whose feature.type is "RNAi_reagent" or "sgRNA" (FTA-224) - as a
    JSON file conforming to the SequenceTargetingReagentDTO LinkML specs for the Alliance
    persistent curation database. A chado database with a full "audit_chado" table is
    required.

"""

import argparse
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from str_handler import SequenceTargetingReagentHandler
from utils import export_chado_data, generate_export_file
import curation_tsv

# Data types handled by this script.
REPORT_LABEL = 'str_curation'

# SequenceTargetingReagentDTO holds its symbol in a plain "name" string and its synonyms in a plain
# list of strings, not in the "{datatype}_symbol_dto"/"{datatype}_synonym_dtos" slots that
# write_primary_tsv() reads by convention, so point it at the right keys. There is no full name
# slot on this class at all, so that column is left empty.


def _str_symbol(entity_dict):
    """Return the STR symbol for the curator TSV."""
    return entity_dict.get('name')


def _str_synonyms(entity_dict):
    """Return the STR synonym strings for the curator TSV."""
    return entity_dict.get('synonyms', [])


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

output_filename = os.environ.get('ALT_OUTPUT', output_filename)

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(
    description='Export FlyBase sequence targeting reagent data to Alliance LinkML JSON.',
    epilog="""
Environment variables:
  SERVER              Database server (e.g. flysql25)
  DATABASE            Database name (e.g. production_chado)
  SQL_PORT            Database port (default: 5432)
  ALT_OUTPUT          Override default output file path
  ADD_OBSOLETE        Set to 'NO' to exclude obsolete/internal rows from the TSVs only; JSON output is unaffected
""",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
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

port = os.environ.get('SQL_PORT', '5432')

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + ':' + port + '/' + database

print(f"Connecting to server:{server} port:{port} database:{database} username:{username}")

engine = create_engine(engine_var_rep)
Session = sessionmaker(bind=engine)
session = Session()

# Create a session to the reference db.
if reference_db:
    # Reference DB host defaults to the primary SERVER unless REFERENCE_SERVER overrides it.
    reference_server = os.environ.get('REFERENCE_SERVER', server)
    engine_var_ref = 'postgresql://' + username + ":" + password + '@' + reference_server + '/' + reference_db
    ref_engine = create_engine(engine_var_ref)
    RefSession = sessionmaker(bind=ref_engine)
    reference_session = RefSession()
else:
    reference_session = None


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase STR data."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    str_handler = SequenceTargetingReagentHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, str_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, str_handler)

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    set_name = str_handler.primary_export_set
    export_dict[set_name] = str_handler.export_data[set_name]
    if len(export_dict[set_name]) == 0:
        if reference_session:
            log.info('No updates to report.')
            return
        log.error(f'The "{set_name}" is unexpectedly empty.')
        raise ValueError(f'The "{set_name}" is unexpectedly empty.')
    generate_export_file(export_dict, log, output_filename)
    tsv_filename = set_up_dict['output_filename']
    curation_tsv.write_primary_tsv(
        log=log, filename=tsv_filename, entities=export_dict[set_name],
        datatype='str', symbol_source=_str_symbol, synonym_source=_str_synonyms,
    )
    # FTA-240: the notes themselves, plus the diagnostic report of note text that
    # clean_free_text() could not handle, so curators can find the rows to fix.
    curation_tsv.write_notes_tsv(
        filename=tsv_filename.replace('.tsv', '_notes.tsv'), entities=export_dict[set_name],
    )
    curation_tsv.write_note_clean_failures_tsv(
        filename=tsv_filename.replace('.tsv', '_internal_note_clean_failures.tsv'),
        failures=str_handler.internal_note_clean_failures,
    )

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
