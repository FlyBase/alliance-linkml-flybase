# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase cassette for Alliance curation database.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_cassette.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python AGR_data_retrieval_curation_cassette.py -v -t -c /path/to/config.cfg
    -l v1.1.2
    -r fb_2024_06_reporting
Notes:
    This script exports FlyBase cassette data as a JSON file conforming to the
    cassetteLinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from cassette_handler import CassetteHandler
from utils import export_chado_data, generate_export_file
import curation_tsv

# Data types handled by this script.
REPORT_LABEL = 'cassette_curation'

# Map each cassette association ingest set to its 'object' identifier field.
_CASSETTE_ASSOC_SECOND_FIELDS = {
    'cassette_transgenic_tool_association_ingest_set': 'transgenic_tool_identifier',
    'cassette_genomic_entity_association_ingest_set': 'genomic_entity_identifier',
}
# Extra TSV columns specific to cassette associations.
_CASSETTE_ASSOC_EXTRAS = [('Comp type curie', 'component_type_curies', '')]

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
    description='Export FlyBase cassette data to Alliance LinkML JSON.',
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


def _write_cassette_tsvs(entities, base_filename):
    """Write the four cassette curator TSVs (primary, notes, component slots, tool uses)."""
    curation_tsv.write_primary_tsv(
        log=log, filename=base_filename, entities=entities, datatype='cassette',
    )
    curation_tsv.write_notes_tsv(
        filename=base_filename.replace('.tsv', '_notes.tsv'), entities=entities,
    )
    curation_tsv.write_components_tsv(
        filename=base_filename.replace('.tsv', '_component_slots.tsv'),
        entities=entities,
        datatype='cassette',
    )
    curation_tsv.write_tool_uses_tsv(
        filename=base_filename.replace('.tsv', '_tool_uses.tsv'),
        entities=entities,
        datatype='cassette',
    )


def _export_primary(cassette_handler):
    """Write the primary cassette JSON and curator TSVs.

    Raises ValueError on unexpected empty exports for non-reference-DB runs;
    reference-DB runs treat an empty export as 'no updates' and skip writing.
    """
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    set_name = cassette_handler.primary_export_set
    export_dict[set_name] = cassette_handler.export_data[set_name]
    if len(export_dict[set_name]) == 0:
        if reference_session:
            log.info('No updates to report.')
            return
        log.error(f'The "{set_name}" is unexpectedly empty.')
        raise ValueError(f'The "{set_name}" is unexpectedly empty.')
    generate_export_file(export_dict, log, output_filename)
    _write_cassette_tsvs(export_dict[set_name], set_up_dict['output_filename'])
    # FTA-211: diagnostic report of internal_notes whose text failed clean_free_text.
    curation_tsv.write_note_clean_failures_tsv(
        filename=set_up_dict['output_filename'].replace('.tsv', '_internal_note_clean_failures.tsv'),
        failures=cassette_handler.internal_note_clean_failures,
    )


def _export_associations(cassette_handler):
    """Write per-sub-type association TSVs and a combined association JSON."""
    association_export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    for sub_type in ('transgenic_tool', 'genomic_entity'):
        set_name = f"cassette_{sub_type}_association"
        ingest_name = f"{set_name}_ingest_set"
        association_export_dict[ingest_name] = []
        association_export_dict[ingest_name].extend(cassette_handler.export_data[ingest_name])
        if len(association_export_dict[ingest_name]) == 0:
            if os.getenv('ADD_CASS_TO_CONSTRUCT') == 'YES':
                raise ValueError(f'The "{set_name}" is unexpectedly empty.')
            log.info(
                f'The "{set_name}" is empty (ADD_CASS_TO_CONSTRUCT not set to YES); skipping.'
            )
            continue
        association_output_filename = output_filename.replace('cassette', f'{set_name}')
        tsv_filename = association_output_filename.replace('.json', '_associations.tsv')
        try:
            curation_tsv.write_association_tsv(
                filename=tsv_filename,
                rows=association_export_dict[ingest_name],
                first_field='cassette_identifier',
                second_field=_CASSETTE_ASSOC_SECOND_FIELDS[ingest_name],
                extra_fields=_CASSETTE_ASSOC_EXTRAS,
            )
        except KeyError as e:
            log.error(f'The "{sub_type}" blew up on tsv generation. keyError {e}')
            raise

    combined_filename = output_filename.replace('cassette', 'cassette_association')
    generate_export_file(association_export_dict, log, combined_filename)


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase cassette data."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    cassette_handler = CassetteHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, cassette_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, cassette_handler)

    if os.getenv('ADD_CASS_TO_CONSTRUCT') == 'YES':
        cassette_handler.populate_anon_cassettes_from_constructs(session)
    else:
        log.warning('ADD_CASS_TO_CONSTRUCT not set to "YES". '
                    'Skipping anonymous cassette creation.')

    _export_primary(cassette_handler)
    if not reference_session:
        _export_associations(cassette_handler)

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
