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
parser = argparse.ArgumentParser(description='inputs')
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


def generate_tsv_file(export_dict, filename):

    with open(filename, 'w') as outfile:
        outfile.write("# Primary FBid\tValid symbol\tValid full name\tsecondary FBid(s)\tsynonyms\n")
        for entity_dict in export_dict["transgenic_tool_ingest_set"]:
            primary = entity_dict["primary_external_id"]
            symbol = ''
            name = ''
            secondary = []
            syns = []
            if "transgenic_tool_full_name_dto" in entity_dict:
                name = entity_dict["transgenic_tool_full_name_dto"]["format_text"]
            if "transgenic_tool_symbol_dto" in entity_dict:
                symbol = entity_dict["transgenic_tool_symbol_dto"]["format_text"]
            if "transgenic_tool_synonym_dtos" in entity_dict:
                for synonym in entity_dict["transgenic_tool_synonym_dtos"]:
                    syns.append(synonym["format_text"])
            if "secondary_identifiers" in entity_dict:
                secondary = entity_dict["secondary_identifiers"]
            outfile.write(f"{primary}\t{symbol}\t{name}\t{'|'.join(secondary)}\t{'|'.join(syns)}\n")


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase AGM."""
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
        generate_tsv_file(export_dict, set_up_dict['output_filename'])

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
