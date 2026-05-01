# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cross-reference chado constructs against curation TSV files.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    audit_constructs_in_curation_tsvs.py [-h] [-v VERBOSE] [-c CONFIG]
    [--tsv-dir TSV_DIR]

Example:
    python audit_constructs_in_curation_tsvs.py -v -c /path/to/config.cfg

Notes:
    Queries chado for the set of all non-obsolete FBtp constructs, walks
    every *_curation*.tsv file under the TSV directory (default: ./tsvs),
    and writes a long-format audit TSV listing each (construct, file, count)
    tuple. Constructs in chado but absent from every TSV are emitted with a
    sentinel file value '<not in any tsv>' and count=0; FBtps appearing in
    TSVs but absent from chado are emitted with sentinel file
    '<not_in_chado>' so typos and stale references are surfaced.

    Anonymous <FBti>_con constructs use FBti-prefixed primary_external_id
    values and so are intentionally not matched here -- this audit covers
    the canonical FBtp set in chado.
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
REPORT_LABEL = 'audit_constructs_in_curation_tsvs'

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
    description='Audit chado FBtp constructs against the curation TSV files.')
parser.add_argument('--tsv-dir', default='tsvs',
                    help='Directory containing the *_curation*.tsv files (default: tsvs).')
args, extra_args = parser.parse_known_args()
log.info(f'These args are handled by this specific script: {args}')
log.info(f'These args are handled by modules: {extra_args}')

# Create SQL Alchemy engine from environmental variables.
engine_var_rep = 'postgresql://' + username + ':' + password + '@' + server + '/' + database
ENGINE = create_engine(engine_var_rep)
Session = sessionmaker(bind=ENGINE)

CONSTRUCT_REGEX = r'^FBtp[0-9]{7}$'
FBTP_RE = re.compile(r'\bFBtp\d{7}\b')
# Construct feature subtypes per src/handler.py:250.
CONSTRUCT_TYPES = (
    'engineered_transposable_element',
    'engineered_region',
    'transgenic_transposable_element',
)
TSV_GLOB = '*_curation*.tsv'
# Filter out filenames whose names contain any of these substrings; these
# are auxiliary curation TSVs (notes, association, component slots, etc.)
# we don't want in the cross-reference audit. 'curation_tsvs' self-excludes
# both this script's own output AND the sibling allele-audit output.
EXCLUDE_SUBSTRINGS = ('_notes', 'association', 'slots', 'tool_uses', 'curation_tsvs')
SENTINEL_NO_TSV = '<not in any tsv>'
SENTINEL_NOT_IN_CHADO = '<not_in_chado>'
# Used in place of SENTINEL_NOT_IN_CHADO when the TSV row carried
# internal=True (likely an obsolete-in-chado entity that the chado query
# filtered out via is_obsolete=False).
SENTINEL_OBSOLETE_IN_CHADO = '<obsolete_in_chado>'


def query_chado_constructs(session):
    """Return the set of non-obsolete FBtp construct uniquenames from chado."""
    log.info('Query chado for non-obsolete FBtp constructs.')
    construct_type = aliased(Cvterm, name='construct_type')
    results = session.query(Feature.uniquename).\
        join(construct_type, construct_type.cvterm_id == Feature.type_id).\
        filter(
            Feature.uniquename.op('~')(CONSTRUCT_REGEX),
            Feature.is_obsolete.is_(False),
            construct_type.name.in_(CONSTRUCT_TYPES)).\
        distinct()
    chado_constructs = {row.uniquename for row in results}
    log.info(f'Loaded {len(chado_constructs):,} non-obsolete FBtp constructs from chado.')
    return chado_constructs


def collect_tsv_counts(tsv_dir):
    """Walk *_curation*.tsv files and count FBtp occurrences per file.

    Skips auxiliary files whose names contain any of EXCLUDE_SUBSTRINGS
    (notes, association, component slots, tool_uses, prior audit output).
    Only the first tab-separated column is scanned -- so a construct is
    only counted when the row's primary id is that construct.
    """
    log.info(f'Walk TSV files matching "{TSV_GLOB}" under {tsv_dir}.')
    counts_by_file = {}
    internal_ids = set()
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
                parts = line.rstrip('\n').split('\t')
                first_col = parts[0]
                # The 'internal' flag is the last column on rows produced by
                # the curation retrieval scripts; older TSVs without that
                # column will fall through as 'not internal'.
                row_is_internal = (parts[-1].strip() == 'True' and len(parts) > 1)
                for m in FBTP_RE.findall(first_col):
                    counter[m] += 1
                    if row_is_internal:
                        internal_ids.add(m)
        counts_by_file[path.name] = counter
        log.info(
            f'  {path.name}: {len(counter):,} unique FBtps, '
            f'{sum(counter.values()):,} occurrences.')
    return counts_by_file, internal_ids


def write_audit(chado_constructs, counts_by_file, internal_ids, output_path):
    """Write the long-format audit TSV plus a summary footer."""
    log.info(f'Write audit output to {output_path}.')
    # Build the union universe so we cover both chado-only and TSV-only constructs.
    tsv_constructs = set()
    for counter in counts_by_file.values():
        tsv_constructs.update(counter)
    universe = sorted(chado_constructs | tsv_constructs)
    file_names = sorted(counts_by_file.keys())

    in_chado_and_any_tsv = 0
    in_chado_no_tsv = 0
    in_tsv_not_in_chado = 0
    in_tsv_obsolete_in_chado = 0
    per_file_unique = {fn: 0 for fn in file_names}

    with open(output_path, 'w') as out:
        out.write('# uniquename\tfile\tcount\n')
        for construct in universe:
            in_chado = construct in chado_constructs
            orphan_sentinel = (SENTINEL_OBSOLETE_IN_CHADO
                               if construct in internal_ids
                               else SENTINEL_NOT_IN_CHADO)
            found_anywhere = False
            for fn in file_names:
                count = counts_by_file[fn].get(construct, 0)
                if count == 0:
                    continue
                found_anywhere = True
                tag = fn if in_chado else f'{fn} {orphan_sentinel}'
                out.write(f'{construct}\t{tag}\t{count}\n')
                per_file_unique[fn] += 1
            if in_chado and not found_anywhere:
                out.write(f'{construct}\t{SENTINEL_NO_TSV}\t0\n')
                in_chado_no_tsv += 1
            if in_chado and found_anywhere:
                in_chado_and_any_tsv += 1
            if not in_chado and found_anywhere:
                if construct in internal_ids:
                    in_tsv_obsolete_in_chado += 1
                else:
                    in_tsv_not_in_chado += 1
        # Summary footer.
        out.write('# --- summary ---\n')
        out.write(f'# total constructs in chado: {len(chado_constructs):,}\n')
        out.write(f'# total constructs found in any TSV: {len(tsv_constructs):,}\n')
        out.write(f'# constructs in chado AND in >=1 TSV: {in_chado_and_any_tsv:,}\n')
        out.write(
            f'# constructs in chado but absent from every TSV: {in_chado_no_tsv:,}\n')
        out.write(
            f'# constructs in TSVs but absent from chado (orphans): {in_tsv_not_in_chado:,}\n')
        out.write(
            f'# constructs in TSVs marked internal (obsolete-in-chado): '
            f'{in_tsv_obsolete_in_chado:,}\n')
        out.write('# per-file unique-construct totals:\n')
        for fn in file_names:
            out.write(f'#   {fn}: {per_file_unique[fn]:,}\n')

    log.info(f'Wrote audit with {len(universe):,} constructs to {output_path}.')
    log.info(
        f'  in-chado-and-any-tsv={in_chado_and_any_tsv:,}, '
        f'in-chado-no-tsv={in_chado_no_tsv:,}, '
        f'in-tsv-not-in-chado={in_tsv_not_in_chado:,}, '
        f'in-tsv-obsolete-in-chado={in_tsv_obsolete_in_chado:,}.')


def main():
    """Run the audit."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    session = Session()
    try:
        chado_constructs = query_chado_constructs(session)
    finally:
        session.close()
    counts_by_file, internal_ids = collect_tsv_counts(args.tsv_dir)
    write_audit(chado_constructs, counts_by_file, internal_ids, output_filename)
    log.info('Ended main function.\n')


if __name__ == '__main__':
    main()
