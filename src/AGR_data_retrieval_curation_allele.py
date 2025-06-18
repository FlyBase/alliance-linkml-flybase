# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase allele for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_allele.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]
    [-l LINKML_RELEASE] [-r REFERENCE_DB] (OPTIONAL)

Example:
    python AGR_data_retrieval_curation_allele.py -v -t -c /path/to/config.cfg
    -l v1.1.2
    -r fb_2024_06_reporting

Notes:
    This script exports FlyBase allele data as a JSON file conforming to the
    Allele LinkML specs for the Alliance persistent curation database.
    A chado database with a full "audit_chado" table is required.

"""

import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from allele_handlers import AlleleHandler, AberrationHandler    #, BalancerHandler
from utils import export_chado_data, generate_export_file

# Data types handled by this script.
REPORT_LABEL = 'allele_curation'

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


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase AGM."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    log.info(f'Exporting data from FlyBase release: {database_release}')
    log.info(f'Output JSON file corresponds to "agr_curation_schema" release: {linkml_release}')

    # Get the data and process it.
    allele_handler = AlleleHandler(log, testing)
    aberration_handler = AberrationHandler(log, testing)
    # balancer_handler = BalancerHandler(log, testing)
    if reference_session:
        export_chado_data(session, log, allele_handler, reference_session=reference_session)
        export_chado_data(session, log, aberration_handler, reference_session=reference_session)
        # export_chado_data(session, log, balancer_handler, reference_session=reference_session)
    else:
        export_chado_data(session, log, allele_handler)
        export_chado_data(session, log, aberration_handler)
        # export_chado_data(session, log, balancer_handler)

    # Export the data.
    export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
    }
    export_dict['allele_ingest_set'] = []
    export_dict['allele_ingest_set'].extend(allele_handler.export_data[allele_handler.primary_export_set])
    export_dict['allele_ingest_set'].extend(aberration_handler.export_data[aberration_handler.primary_export_set])
    # export_dict['allele_ingest_set'].extend(balancer_handler.export_data[balancer_handler.primary_export_set])
    generate_export_file(export_dict, log, output_filename)

    if not reference_session:
        # Export the gene-allele associations to a separate file.
        association_output_filename = output_filename.replace('allele', 'allele_gene_association')
        association_export_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
        }
        association_export_dict['allele_gene_association_ingest_set'] = []
        association_export_dict['allele_gene_association_ingest_set'].extend(allele_handler.export_data['allele_gene_association_ingest_set'])
        association_export_dict['allele_gene_association_ingest_set'].extend(aberration_handler.export_data['allele_gene_association_ingest_set'])
        generate_export_file(association_export_dict, log, association_output_filename)

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
