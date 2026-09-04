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


def _is_aberration_cell(entity_dict):
    """Return 'true' for FBab aberrations, else '' (FTA-217).

    Derived from the primary ID rather than read from the exported 'is_aberration' key so that the
    curator TSV is populated even when ADD_IS_ABERRATION is unset. That env var gates the JSON only,
    to protect Alliance schema validation; it says nothing about which entities are aberrations. The
    two can never disagree: FTA-217 defines an aberration as an entry whose primary ID starts 'FBab'.
    """
    return 'true' if entity_dict.get('primary_external_id', '').startswith('FB:FBab') else ''


def _make_is_balancer_cell(balancer_ids):
    """Return a TSV cell function reporting 'true' for balancer-flagged FBab aberrations (FTA-235).

    Reads the AberrationHandler's detected ID set rather than the exported 'is_balancer' key, so the
    curator TSV is populated even when ADD_IS_ABERRATION is unset. That env var gates the JSON only, to
    protect Alliance schema validation; it says nothing about which aberrations are balancers.
    """
    def _is_balancer_cell(entity_dict):
        return 'true' if entity_dict.get('primary_external_id', '') in balancer_ids else ''
    return _is_balancer_cell


# Allele primary TSV carries the FTA-217 aberration flag for curator review; the FTA-235 balancer flag is
# appended in main(), where the AberrationHandler that detected the balancers is available.
_ALLELE_PRIMARY_EXTRAS = [
    ('is_aberration', _is_aberration_cell, ''),
]


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
                            'carries'/'breakpoint_allele' relations, plus the FTA-236 alleles and insertions
                            carried by a balancer and moved to its parent aberration). Off by default: the
                            Alliance schema has no such ingest set and lacks the two CV terms, so the data
                            cannot yet be loaded.
  ADD_IS_ABERRATION         Set to 'YES' to emit the 'is_aberration' boolean for FBab entities, and the
                            'is_balancer' boolean for those FBab entities flagged as balancers (FTA-235).
                            Both slots come from the same schema PR (#327), so one gate covers them.
                            Requires a LinkML release containing the slots (absent from v2.17.0).

Notes:
  FTA-236: the 38 FBba balancers carrying a "FTA: Balancer - merge with parent ..." internal note have
  their synonyms, secondary IDs, comments, internal notes and references folded into the parent FBab
  aberration. Those parts are NOT gated (the slots are long-released); only the carried alleles and
  insertions wait on ADD_ALLELE_ALLELE_ASSOC. Three balancers are excluded from the 'carries' move by
  curator decision: FM1 (FBba0000011), SM6a (FBba0000039) and SM6b (FBba0000040). See the
  *_balancer_merges.tsv output for what moved.

  FTA-237: 24 aberrations are exported under their balancer's symbol and full name, driven by curated
  "FTA: Balancer - use balancer symbol and fullname for parent ..." internal notes, plus In(2LR)SM6
  (FBab0004818), which is renamed to "SM6" / "Second Multiple 6" from a hard-coded table because no
  flag exists for it. The aberration's own names are kept as synonyms. Not gated: the slots involved
  are long-released. See the *_balancer_renames.tsv output for the full list.
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
        primary_extras = _ALLELE_PRIMARY_EXTRAS + [
            ('is_balancer', _make_is_balancer_cell(aberration_handler.balancer_ids), ''),
        ]
        curation_tsv.write_primary_tsv(
            log=log, filename=tsv_filename, entities=entities, datatype='allele',
            extra_fields=primary_extras,
        )
        curation_tsv.write_notes_tsv(
            filename=tsv_filename.replace('.tsv', '_notes.tsv'), entities=entities,
        )
        # FTA-236: one row per FBba balancer merged into a parent FBab aberration, so curators can check
        # what moved. Written unconditionally; the 'carries' counts stay 0 unless ADD_ALLELE_ALLELE_ASSOC
        # is 'YES', since that gate controls the association export. write_association_tsv() needs a
        # 'relation_name' cell, hence the constant, and the evidence sentinel is blanked because a merge
        # carries no evidence of its own.
        merges_filename = tsv_filename.replace('.tsv', '_balancer_merges.tsv')
        curation_tsv.write_association_tsv(
            filename=merges_filename,
            rows=[dict(row, relation_name='merged_from_balancer') for row in aberration_handler.balancer_merge_report],
            first_field='fbab_id',
            second_field='fbba_id',
            extra_fields=[
                ('fbba_symbol', 'fbba_symbol', ''),
                ('synonyms', 'synonyms', 0),
                ('secondary_ids', 'secondary_ids', 0),
                ('comments', 'comments', 0),
                ('internal_notes', 'internal_notes', 0),
                ('references', 'references', 0),
                ('carries_alleles', 'carries_alleles', 0),
                ('carries_excluded', 'carries_excluded', 0),
            ],
            no_pubs_sentinel='',
        )
        log.info(f'Generated TSV: {merges_filename} ({len(aberration_handler.balancer_merge_report)} balancer merges)')
        # FTA-237: one row per renamed aberration, so curators can check the 24 flag-driven renames and
        # the hard-coded In(2LR)SM6 case against their spreadsheet. write_association_tsv() is reused
        # rather than adding another writer; it needs a 'relation_name' cell, hence the constant below,
        # and no_pubs_sentinel is blanked because a rename carries no evidence of its own.
        renames_filename = tsv_filename.replace('.tsv', '_balancer_renames.tsv')
        curation_tsv.write_association_tsv(
            filename=renames_filename,
            rows=[dict(row, relation_name='renamed_after_balancer') for row in aberration_handler.balancer_rename_report],
            first_field='fbab_id',
            second_field='fbba_id',
            extra_fields=[
                ('source', 'source', ''),
                ('new_symbol', 'new_symbol', ''),
                ('old_symbol', 'old_symbol', ''),
                ('new_full_name', 'new_full_name', ''),
                ('old_full_name', 'old_full_name', ''),
            ],
            no_pubs_sentinel='',
        )
        log.info(f'Generated TSV: {renames_filename} ({len(aberration_handler.balancer_rename_report)} renames)')
        # FTA-221: diagnostic TSV of note props whose text could not be cleaned (NULL, blank, or bad
        # SGML), so curators can find and fix the offending rows. Mirrors the construct/cassette
        # scripts. Combined across both handlers, since each collects its own failures.
        note_failures = allele_handler.internal_note_clean_failures + aberration_handler.internal_note_clean_failures
        failures_filename = tsv_filename.replace('.tsv', '_internal_note_clean_failures.tsv')
        curation_tsv.write_note_clean_failures_tsv(filename=failures_filename, failures=note_failures)
        log.info(f'Generated TSV: {failures_filename} ({len(note_failures)} uncleanable note props)')

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
