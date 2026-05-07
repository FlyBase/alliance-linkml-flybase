# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase experimental tool for Alliance curation database.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_transgenic_tool.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python AGR_data_retrieval_curation_transgenic_tool.py -v -t -c /path/to/config.cfg
    -l v1.1.2
    -r fb_2024_06_reporting
Notes:
    This script exports FlyBase experimental tool data as a JSON file conforming to the
    transgenic_tool LinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from transgenic_tool_handler import ExperimentalToolHandler
from utils import export_chado_data, generate_export_file
import curation_tsv

# Data types handled by this script.
REPORT_LABEL = 'transgenic_tool_curation'

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
parser = argparse.ArgumentParser(
    description='Export FlyBase experimental tool data to Alliance LinkML JSON.',
    epilog="""
Environment variables:
  SERVER              Database server (e.g. flysql25)
  DATABASE            Database name (e.g. production_chado)
  ADD_OBSOLETE        Set to 'NO' to exclude obsolete/internal rows from the TSVs only; JSON output is unaffected
""",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument('-l', '--linkml_release',
                    help='The "agr_curation_schema" LinkML release number.', required=True)
parser.add_argument('-r', '--reference_db',
                    help='The name of a previous reference db for incremental exports.', required=False)

# Use parse_known_args(), not parse_args(),
# to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release
reference_db = args.reference_db

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
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
    """Run the steps for exporting LinkML-compliant FlyBase transgenic tool data."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    tool_handler = ExperimentalToolHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, tool_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, tool_handler)

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    export_dict[tool_handler.primary_export_set] = tool_handler.export_data[tool_handler.primary_export_set]
    if len(export_dict[tool_handler.primary_export_set]) == 0:
        if reference_session:
            log.info('No updates to report.')
        else:
            log.error(f'The "{tool_handler.primary_export_set}" is unexpectedly empty.')
            raise ValueError(f'The "{tool_handler.primary_export_set}" is unexpectedly empty.')
    else:
        generate_export_file(export_dict, log, output_filename)
        tsv_filename = set_up_dict['output_filename']
        entities = export_dict[tool_handler.primary_export_set]
        curation_tsv.write_primary_tsv(
            log=log, filename=tsv_filename, entities=entities, datatype='transgenic_tool',
        )
        curation_tsv.write_notes_tsv(
            filename=tsv_filename.replace('.tsv', '_notes.tsv'), entities=entities,
        )

    if not reference_session:
        # Export tool associations to a separate file.
        association_output_filename = output_filename.replace('tool', 'tool_association')
        association_export_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
        }
        # tool_tool associations.
        assoc = 'transgenic_tool_transgenic_tool_association_ingest_set'
        association_export_dict[assoc] = []
        association_export_dict[assoc].extend(tool_handler.export_data['tool_association_ingest_set'])
        if len(association_export_dict[assoc]) == 0:
            log.error(f'The "{assoc}" is unexpectedly empty.')
            raise ValueError(f'The "{assoc}" is unexpectedly empty.')
        # Print the output file.
        generate_export_file(association_export_dict, log, association_output_filename)
        # Migrate to standard subject/relation/object/evidence layout (was previously
        # a 3-column file with relation_name written under a misleading 'Pub' header
        # and evidence_curies suppressed).
        assoc_tsv_filename = set_up_dict['output_filename'].replace('.tsv', '_associations.tsv')
        curation_tsv.write_association_tsv(
            filename=assoc_tsv_filename,
            rows=association_export_dict[assoc],
            first_field='transgenic_tool_subject_identifier',
            second_field='transgenic_tool_object_identifier',
        )
    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
