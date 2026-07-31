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
from os import getenv
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from harvdev_utils.psycopg_functions import set_up_db_reading
from allele_handlers import AlleleHandler, AberrationHandler    # BalancerHandler
from utils import export_chado_data, generate_export_file
import curation_tsv


# Map each allele association ingest set to its 'object' identifier field.
_ALLELE_ASSOC_SECOND_FIELDS = {
    'allele_gene_association_ingest_set': 'gene_identifier',
    'allele_construct_association_ingest_set': 'construct_identifier',
    'allele_allele_association_ingest_set': 'object_allele_identifier',
}
# Allele association TSVs include an extra 'internal' boolean column.
_ALLELE_ASSOC_EXTRAS = [('internal', 'internal', False)]


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
parser = argparse.ArgumentParser(
    description='Export FlyBase allele data to Alliance LinkML JSON.',
    epilog="""
Environment variables:
  SERVER                    Database server (e.g. flysql25)
  DATABASE                  Database name (e.g. production_chado)
  ADD_OBSOLETE              Set to 'NO' to exclude obsolete/internal rows from the TSVs only; JSON output is unaffected
  ADD_ALLELE_ALLELE_ASSOC   Set to 'YES' to emit the FTA-218 'allele_allele_association_ingest_set' (aberration
                            'carries'/'breakpoint_allele' relations). Off by default: the Alliance schema has no
                            such ingest set and lacks the two CV terms, so the data cannot yet be loaded.
""",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
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
    """Run the steps for exporting LinkML-compliant FlyBase allele data."""
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
    if len(export_dict['allele_ingest_set']) == 0:
        if reference_session:
            log.info('No updates to report.')
        else:
            log.error('The "allele_ingest_set" is unexpectedly empty.')
            raise ValueError('The "allele_ingest_set" is unexpectedly empty.')
    else:
        generate_export_file(export_dict, log, output_filename)
        tsv_filename = set_up_dict['output_filename']
        entities = export_dict['allele_ingest_set']
        curation_tsv.write_primary_tsv(
            log=log, filename=tsv_filename, entities=entities, datatype='allele',
        )
        curation_tsv.write_notes_tsv(
            filename=tsv_filename.replace('.tsv', '_notes.tsv'), entities=entities,
        )

    if not reference_session:
        # Export the gene-allele associations to a separate file.
        association_output_filename = output_filename.replace('allele', 'allele_association')
        association_export_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
        }
        # Allele-gene associations.
        association_export_dict['allele_gene_association_ingest_set'] = []
        association_export_dict['allele_gene_association_ingest_set'].extend(allele_handler.export_data['allele_gene_association_ingest_set'])
        association_export_dict['allele_gene_association_ingest_set'].extend(aberration_handler.export_data['allele_gene_association_ingest_set'])
        if len(association_export_dict['allele_gene_association_ingest_set']) == 0:
            log.error('The "allele_gene_association_ingest_set" is unexpectedly empty.')
            raise ValueError('The "allele_gene_association_ingest_set" is unexpectedly empty.')
        # Allele-construct associations.
        association_export_dict['allele_construct_association_ingest_set'] = []
        association_export_dict['allele_construct_association_ingest_set'].extend(allele_handler.export_data['allele_construct_association_ingest_set'])
        if len(association_export_dict['allele_construct_association_ingest_set']) == 0:
            log.error('The "allele_construct_association_ingest_set" is unexpectedly empty.')
            raise ValueError('The "allele_construct_association_ingest_set" is unexpectedly empty.')
        # Allele-allele associations (FTA-218). Gated by ADD_ALLELE_ALLELE_ASSOC=YES because the Alliance schema has
        # no "allele_allele_association_ingest_set" and lacks the "carries"/"breakpoint_allele" CV terms, so including
        # this set would make the file unloadable.
        association_ingest_names = ['allele_gene_association_ingest_set',
                                    'allele_construct_association_ingest_set']
        if getenv('ADD_ALLELE_ALLELE_ASSOC', None) == 'YES':
            association_export_dict['allele_allele_association_ingest_set'] = []
            association_export_dict['allele_allele_association_ingest_set'].extend(
                aberration_handler.export_data['allele_allele_association_ingest_set'])
            if len(association_export_dict['allele_allele_association_ingest_set']) == 0:
                log.error('The "allele_allele_association_ingest_set" is unexpectedly empty.')
                raise ValueError('The "allele_allele_association_ingest_set" is unexpectedly empty.')
            association_ingest_names.append('allele_allele_association_ingest_set')
        else:
            log.info('ADD_ALLELE_ALLELE_ASSOC not set to "YES"; omitting the "allele_allele_association_ingest_set".')
        # Print the output file.
        generate_export_file(association_export_dict, log, association_output_filename)
        # Per-association-set TSV files for easy curator review.
        for ingest_name in association_ingest_names:
            set_name = ingest_name.replace('_ingest_set', '')
            assoc_tsv_filename = set_up_dict['output_filename'].replace('allele', set_name)
            assoc_tsv_filename = assoc_tsv_filename.replace('.tsv', '_associations.tsv')
            try:
                curation_tsv.write_association_tsv(
                    filename=assoc_tsv_filename,
                    rows=association_export_dict[ingest_name],
                    first_field='allele_identifier',
                    second_field=_ALLELE_ASSOC_SECOND_FIELDS[ingest_name],
                    extra_fields=_ALLELE_ASSOC_EXTRAS,
                )
            except KeyError as e:
                log.error(f'The "{ingest_name}" blew up on tsv generation. KeyError {e}')
                raise

    log.info('Ended main function.\n')


if __name__ == "__main__":
    main()
