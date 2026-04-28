# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cross-reference chado cassette-eligible features against curation TSV files.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    audit_cassettes_in_curation_tsvs.py [-h] [-v VERBOSE] [-c CONFIG]
    [--tsv-dir TSV_DIR]

Example:
    python audit_cassettes_in_curation_tsvs.py -v -c /path/to/config.cfg

Notes:
    Cassettes are not their own chado feature class - they are derived
    from FBal alleles plus the FBti-suffixed anon cassettes built by the
    FTA-136 generic-TI pipeline. So the chado universe for this audit
    is the union of:

        - non-obsolete FBal alleles (cvterm.name = 'allele'), and
        - non-obsolete FBti insertions (cvterm.name in the insertion
          feature_subtypes from handler.py:252).

    The script walks every *_curation*.tsv file under the TSV directory
    (default ./tsvs) and counts FBal/FBti occurrences in column 0 only.
    Anon cassettes are stored in the cassette TSV as 'FB:FBti0xxxxxx_cas';
    the regex \\bFBti\\d{7}\\b matches the FBti portion, so each anon
    cassette counts under its FBti id (no double-counting vs the suffixed
    form).

    Output format mirrors src/audit_alleles_in_curation_tsvs.py: long
    rows of (uniquename, file, count) plus sentinel rows for "in chado,
    not in any TSV" and "in TSVs, not in chado".
"""

import argparse
import re
from collections import Counter
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, aliased
from harvdev_utils.psycopg_functions import set_up_db_reading
from harvdev_utils.reporting import Cvterm, Feature

# Data types handled by this script.
REPORT_LABEL = 'audit_cassettes_in_curation_tsvs'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
output_filename = set_up_dict['output_filename']
log = set_up_dict['log']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(
    description='Audit chado cassette-eligible features (FBal + FBti) against the curation TSV files.')
parser.add_argument('--tsv-dir', default='tsvs',
                    help='Directory containing the *_curation*.tsv files (default: tsvs).')
args, extra_args = parser.parse_known_args()
log.info(f'These args are handled by this specific script: {args}')
log.info(f'These args are handled by modules: {extra_args}')

# Create SQL Alchemy engine from environmental variables.
engine_var_rep = 'postgresql://' + username + ':' + password + '@' + server + '/' + database
ENGINE = create_engine(engine_var_rep)
Session = sessionmaker(bind=ENGINE)

ALLELE_REGEX = r'^FBal[0-9]{7}$'
INSERTION_REGEX = r'^FBti[0-9]{7}$'
# Match FBal*** or FBti*** anywhere in the cell; works for 'FB:FBal0000010'
# and 'FB:FBti0167947_cas' alike.
CASSETTE_RE = re.compile(r'\b(?:FBal|FBti)\d{7}\b')
# Insertion feature_subtypes per src/handler.py:252.
INSERTION_TYPES = (
    'insertion_site',
    'transposable_element',
    'transposable_element_insertion_site',
)
TSV_GLOB = '*_curation*.tsv'
EXCLUDE_SUBSTRINGS = ('_notes', 'association', 'slots', 'tool_uses', 'curation_tsvs')
SENTINEL_NO_TSV = '<not in any tsv>'
SENTINEL_NOT_IN_CHADO = '<not_in_chado>'


def query_chado_cassettes_universe(session):
    """Return the union of non-obsolete FBal alleles and FBti insertions from chado."""
    log.info('Query chado for non-obsolete FBal alleles.')
    allele_type = aliased(Cvterm, name='allele_type')
    fbal_results = session.query(Feature.uniquename).\
        join(allele_type, allele_type.cvterm_id == Feature.type_id).\
        filter(
            Feature.uniquename.op('~')(ALLELE_REGEX),
            Feature.is_obsolete.is_(False),
            allele_type.name == 'allele').\
        distinct()
    fbals = {row.uniquename for row in fbal_results}
    log.info(f'  loaded {len(fbals):,} non-obsolete FBal alleles.')

    log.info('Query chado for non-obsolete FBti insertions.')
    insertion_type = aliased(Cvterm, name='insertion_type')
    fbti_results = session.query(Feature.uniquename).\
        join(insertion_type, insertion_type.cvterm_id == Feature.type_id).\
        filter(
            Feature.uniquename.op('~')(INSERTION_REGEX),
            Feature.is_obsolete.is_(False),
            insertion_type.name.in_(INSERTION_TYPES)).\
        distinct()
    fbtis = {row.uniquename for row in fbti_results}
    log.info(f'  loaded {len(fbtis):,} non-obsolete FBti insertions.')

    universe = fbals | fbtis
    log.info(f'Cassette universe = {len(universe):,} ids ({len(fbals):,} FBal + {len(fbtis):,} FBti).')
    return universe


def collect_tsv_counts(tsv_dir):
    """Walk *_curation*.tsv files and count FBal/FBti occurrences per file (column 0 only)."""
    log.info(f'Walk TSV files matching "{TSV_GLOB}" under {tsv_dir}.')
    counts_by_file = {}
    all_paths = sorted(Path(tsv_dir).glob(TSV_GLOB))
    paths = [p for p in all_paths
             if not any(sub in p.name for sub in EXCLUDE_SUBSTRINGS)]
    skipped = len(all_paths) - len(paths)
    if skipped:
        log.info(
            f'Skipping {skipped} files whose names contain any of '
            f'{EXCLUDE_SUBSTRINGS}.')
    if not paths:
        log.warning(f'No files matched "{TSV_GLOB}" under {tsv_dir}.')
    for path in paths:
        counter = Counter()
        with path.open() as fh:
            for line in fh:
                if not line or line.startswith('#'):
                    continue
                first_col = line.split('\t', 1)[0]
                for m in CASSETTE_RE.findall(first_col):
                    counter[m] += 1
        counts_by_file[path.name] = counter
        log.info(
            f'  {path.name}: {len(counter):,} unique FBal/FBti ids, '
            f'{sum(counter.values()):,} occurrences.')
    return counts_by_file


def write_audit(chado_universe, counts_by_file, output_path):
    """Write the long-format audit TSV plus a summary footer."""
    log.info(f'Write audit output to {output_path}.')
    tsv_ids = set()
    for counter in counts_by_file.values():
        tsv_ids.update(counter)
    universe = sorted(chado_universe | tsv_ids)
    file_names = sorted(counts_by_file.keys())

    in_chado_and_any_tsv = 0
    in_chado_no_tsv = 0
    in_tsv_not_in_chado = 0
    per_file_unique = {fn: 0 for fn in file_names}

    with open(output_path, 'w') as out:
        out.write('# uniquename\tfile\tcount\n')
        for uname in universe:
            in_chado = uname in chado_universe
            found_anywhere = False
            for fn in file_names:
                count = counts_by_file[fn].get(uname, 0)
                if count == 0:
                    continue
                found_anywhere = True
                tag = fn if in_chado else f'{fn} {SENTINEL_NOT_IN_CHADO}'
                out.write(f'{uname}\t{tag}\t{count}\n')
                per_file_unique[fn] += 1
            if in_chado and not found_anywhere:
                out.write(f'{uname}\t{SENTINEL_NO_TSV}\t0\n')
                in_chado_no_tsv += 1
            if in_chado and found_anywhere:
                in_chado_and_any_tsv += 1
            if not in_chado and found_anywhere:
                in_tsv_not_in_chado += 1
        # Summary footer.
        out.write('# --- summary ---\n')
        out.write(f'# total cassette-eligible ids in chado: {len(chado_universe):,}\n')
        out.write(f'# total cassette-eligible ids in any TSV: {len(tsv_ids):,}\n')
        out.write(f'# ids in chado AND in >=1 TSV: {in_chado_and_any_tsv:,}\n')
        out.write(
            f'# ids in chado but absent from every TSV: {in_chado_no_tsv:,}\n')
        out.write(
            f'# ids in TSVs but absent from chado: {in_tsv_not_in_chado:,} '
            f'(orphaned references)\n')
        out.write('# per-file unique-id totals:\n')
        for fn in file_names:
            out.write(f'#   {fn}: {per_file_unique[fn]:,}\n')

    log.info(f'Wrote audit with {len(universe):,} ids to {output_path}.')
    log.info(
        f'  in-chado-and-any-tsv={in_chado_and_any_tsv:,}, '
        f'in-chado-no-tsv={in_chado_no_tsv:,}, '
        f'in-tsv-not-in-chado={in_tsv_not_in_chado:,}.')


def main():
    """Run the audit."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    session = Session()
    try:
        chado_universe = query_chado_cassettes_universe(session)
    finally:
        session.close()
    counts_by_file = collect_tsv_counts(args.tsv_dir)
    write_audit(chado_universe, counts_by_file, output_filename)
    log.info('Ended main function.\n')


if __name__ == '__main__':
    main()
