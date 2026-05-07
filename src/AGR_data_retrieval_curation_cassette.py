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

# Data types handled by this script.
REPORT_LABEL = 'cassette_curation'

# TSV format constants.
CASSETTE_PRIMARY_HEADER = (
    "# Primary FBid\tValid symbol\tValid full name\tsecondary FBid(s)\tsynonyms\tinternal\n"
)
CASSETTE_NOTES_HEADER = "# Primary FBid\ttype\tcomment\n"
CASSETTE_COMPONENTS_HEADER = "# Primary FBid\tsymbol\trelation\ttaxon\tevidence\n"
CASSETTE_TOOL_USES_HEADER = "# Primary FBid\tevidence\ttool_uses\n"
NO_PUBS_SENTINEL = "NO PUBS"
EVIDENCE_DELIMITER = "|"

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


def generate_tsv_file(export_dict, filename):
    """Generate tsv files for curators to read more easily. This can be commented out later.

    ADD_OBSOLETE=NO suppresses obsolete/internal rows in the TSVs only; the
    JSON output is unaffected.
    """
    skip_obsolete = os.environ.get('ADD_OBSOLETE') == 'NO'
    if skip_obsolete:
        log.info('ADD_OBSOLETE=NO: excluding obsolete/internal cassettes from TSV.')

    def _skip(d):
        return skip_obsolete and (d.get('internal') or d.get('obsolete'))

    with open(filename, 'w') as outfile:
        outfile.write(CASSETTE_PRIMARY_HEADER)
        for entity_dict in export_dict["cassette_ingest_set"]:
            if _skip(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            symbol = ''
            name = ''
            secondary = []
            syns = []
            if "cassette_full_name_dto" in entity_dict:
                name = entity_dict["cassette_full_name_dto"]["format_text"]
            if "cassette_symbol_dto" in entity_dict:
                symbol = entity_dict["cassette_symbol_dto"]["format_text"]
            if "cassette_synonym_dtos" in entity_dict:
                for synonym in entity_dict["cassette_synonym_dtos"]:
                    syns.append(synonym["format_text"])
            if "secondary_identifiers" in entity_dict:
                secondary = entity_dict["secondary_identifiers"]
            internal = entity_dict.get("internal", False)
            secondary_str = EVIDENCE_DELIMITER.join(secondary)
            syns_str = EVIDENCE_DELIMITER.join(syns)
            try:
                outfile.write(
                    f"{primary}\t{symbol}\t{name}\t{secondary_str}\t{syns_str}\t{internal}\n")
            except TypeError:
                log.error(f"entity_dict: {entity_dict}")
                log.error(f"primary: {primary}")
                log.error(f"secondary {secondary}")
                log.error(f"symbol: {symbol}")
                log.error(f"name: {name}")
                log.error(f"syns: {syns}")
                log.error(f"internal: {internal}")
                raise

    filename = filename.replace('.tsv', '_notes.tsv')
    with open(filename, 'w') as outfile:
        outfile.write(CASSETTE_NOTES_HEADER)
        for entity_dict in export_dict["cassette_ingest_set"]:
            if _skip(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            if "note_dtos" in entity_dict:
                for note in entity_dict["note_dtos"]:
                    ntype = note["note_type_name"]
                    txt = note['free_text']
                    outfile.write(f"{primary}\t{ntype}\t{txt}\n")

    filename = filename.replace('_notes.tsv', '_component_slots.tsv')
    with open(filename, 'w') as outfile:
        outfile.write(CASSETTE_COMPONENTS_HEADER)
        for entity_dict in export_dict["cassette_ingest_set"]:
            if _skip(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            if "cassette_component_dtos" in entity_dict:
                for comp in entity_dict["cassette_component_dtos"]:
                    symbol = comp["component_symbol"]
                    relation = comp['relation_name']
                    taxon = comp['taxon_curie']
                    if 'evidence_curies' in comp:
                        evidence = EVIDENCE_DELIMITER.join(comp['evidence_curies'])
                    else:
                        evidence = ""
                    outfile.write(f"{primary}\t{symbol}\t{relation}\t{taxon}\t{evidence}\n")

    filename = filename.replace('_component_slots.tsv', '_tool_uses.tsv')
    with open(filename, 'w') as outfile:
        outfile.write(CASSETTE_TOOL_USES_HEADER)
        for entity_dict in export_dict["cassette_ingest_set"]:
            if _skip(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            if "cassette_use_dtos" in entity_dict:
                for comp in entity_dict["cassette_use_dtos"]:
                    if 'evidence_curies' in comp:
                        evidence = EVIDENCE_DELIMITER.join(comp['evidence_curies'])
                    else:
                        evidence = NO_PUBS_SENTINEL
                    tools = EVIDENCE_DELIMITER.join(comp["use_curies"])
                    outfile.write(f"{primary}\t{tools}\t{evidence}\n")


def generate_association_tsv_file(export_dict, ingest_name, filename):
    """ADD_OBSOLETE=NO suppresses obsolete/internal rows from the TSV only."""
    skip_obsolete = os.environ.get('ADD_OBSOLETE') == 'NO'
    filename = filename.replace('.tsv', '_associations.tsv')
    # To help in debugging, the 'first_entity' and 'second_entity' variables are used:
    # - to get the entities involved in the association out of the relevant 'export_dict[ingest_name]'
    # - AND as the column headers in the tsv output file
    first_entity = 'cassette_identifier'
    if ingest_name == 'cassette_transgenic_tool_association_ingest_set':
        second_entity = 'transgenic_tool_identifier'
    elif ingest_name == 'cassette_genomic_entity_association_ingest_set':
        second_entity = 'genomic_entity_identifier'
    else:
        second_entity = 'sequence_targeting_reagent_identifier'
    with open(filename, 'w') as outfile:
        outfile.write(f"#{first_entity}\tRelationship\t{second_entity}\tEvidence\tComp type curie\n")
        for entity_dict in export_dict[ingest_name]:
            if skip_obsolete and (entity_dict.get('internal') or entity_dict.get('obsolete')):
                continue
            sub = entity_dict[first_entity]
            obj = entity_dict[second_entity]
            rel_type = entity_dict['relation_name']
            if 'evidence_curies' in entity_dict:
                pubs = EVIDENCE_DELIMITER.join(entity_dict['evidence_curies'])
            else:
                pubs = NO_PUBS_SENTINEL
            if 'component_type_curies' in entity_dict:
                comp = EVIDENCE_DELIMITER.join(entity_dict['component_type_curies'])
            else:
                comp = ""
            outfile.write(f"{sub}\t{rel_type}\t{obj}\t{pubs}\t{comp}\n")


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase AGM."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    cassette_handler = CassetteHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, cassette_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, cassette_handler)

    # Optionally pull anonymous cassette data from the construct pipeline.
    if os.getenv('ADD_CASS_TO_CONSTRUCT') == 'YES':
        cassette_handler.populate_anon_cassettes_from_constructs(session)
    else:
        log.warning('ADD_CASS_TO_CONSTRUCT not set to "YES". '
                    'Skipping anonymous cassette creation.')

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    set_name = cassette_handler.primary_export_set
    export_dict[set_name] = cassette_handler.export_data[set_name]
    if len(export_dict[set_name]) == 0:
        if reference_session:
            log.info('No updates to report.')
        else:
            log.error(f'The "{set_name}" is unexpectedly empty.')
            raise ValueError(f'The "{set_name}" is unexpectedly empty.')
    else:
        generate_export_file(export_dict, log, output_filename)
        generate_tsv_file(export_dict, set_up_dict['output_filename'])

    if not reference_session:
        # Export cassette_associations to a separate files.
        association_export_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
        }
        # add each set to association export dict
        # and output tsv's to separate files.
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
            # Print the output tsv file.
            association_output_filename = output_filename.replace('cassette', f'{set_name}')
            tsv_filename = association_output_filename.replace('.json', '.tsv')
            try:
                generate_association_tsv_file(association_export_dict, ingest_name, tsv_filename)
            except KeyError as e:
                log.error(f'The "{sub_type}" blew up on tsv generation. keyError {e}')
                raise

        # output all association in one file.
        association_output_filename = output_filename.replace('cassette', 'cassette_association')
        generate_export_file(association_export_dict, log, association_output_filename)
    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
