# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge FlyBase construct, cassette and transgenic tool exports into two combined files.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

Usage:
    merge_construct_family.py [-h] -l LINKML_RELEASE -i INPUT_DIR [-o OUTPUT_DIR]

Example:
    python merge_construct_family.py -l v2.19.0 -i /data/alliance/PERSISTENT_main_construct_family

Notes:
    The Alliance takes construct-related data as two uploads: one file of all entities
    (Construct, Cassette, TransgenicTool) and one file of their associations, both among
    themselves and with other entities such as genes and STRs.
    This script reads the six JSON files written by the construct, cassette and
    transgenic_tool retrieval scripts and combines them, stamping the given LinkML release.

    Those retrieval scripts must have been run with ADD_CASS_TO_CONSTRUCT=YES (from v2.19.0 the
    schema has no construct components or construct-genomic entity associations) and
    ADD_TOOL_USES=YES (transgenic_tool_use_dtos, from v2.18.0).

"""

import argparse
import glob
import json
import logging
import os

from utils import generate_export_file

# Primary entity ingest sets, keyed by the stem of the retrieval script's output file.
ENTITY_SETS = {
    'construct_curation_': 'construct_ingest_set',
    'cassette_curation_': 'cassette_ingest_set',
    'transgenic_tool_curation_': 'transgenic_tool_ingest_set',
}

# Association ingest sets, keyed by the stem of the retrieval script's output file: associations among
# the three entity types, plus cassette associations to genomic entities and STRs.
ASSOC_SETS = {
    'construct_association_curation_': ['construct_cassette_association_ingest_set'],
    'cassette_association_curation_': [
        'cassette_transgenic_tool_association_ingest_set',
        'cassette_genomic_entity_association_ingest_set',
        'cassette_str_association_ingest_set',
    ],
    'transgenic_tool_association_curation_': ['transgenic_tool_transgenic_tool_association_ingest_set'],
}

log = logging.getLogger(__name__)


def find_input_file(input_dir, stem):
    """Return the single JSON file in input_dir whose name starts with stem."""
    matches = glob.glob(os.path.join(input_dir, f'{stem}*.json'))
    if len(matches) != 1:
        raise ValueError(f'Expected one "{stem}*.json" file in {input_dir}, found {len(matches)}: {matches}')
    return matches[0]


def load_sets(filename, wanted_sets, member_releases):
    """Load the wanted ingest sets from one export file, logging any other keys it drops."""
    with open(filename) as infile:
        data = json.load(infile)
    member_releases[filename] = data.get('alliance_member_release_version')
    header_keys = ('linkml_version', 'alliance_member_release_version')
    for key in data:
        if key not in header_keys and key not in wanted_sets:
            log.info(f'{os.path.basename(filename)}: dropping "{key}" ({len(data[key])} entries).')
    sets = {}
    for set_name in wanted_sets:
        if not data.get(set_name):
            raise ValueError(f'The "{set_name}" is missing or empty in {filename}.')
        sets[set_name] = data[set_name]
    return sets


def check_no_construct_components(constructs):
    """Fail if any construct carries construct_component_dtos, which v2.19.0 removed."""
    with_components = [c.get('primary_external_id') for c in constructs if c.get('construct_component_dtos')]
    if with_components:
        raise ValueError(f'{len(with_components)} constructs have construct_component_dtos (e.g. '
                         f'{with_components[:5]}); rerun the construct export with ADD_CASS_TO_CONSTRUCT=YES.')


def main():
    """Combine the construct-family exports into one entity file and one association file."""
    parser = argparse.ArgumentParser(
        description='Merge construct, cassette and transgenic tool exports into two combined JSON files.',
    )
    parser.add_argument('-l', '--linkml_release', required=True,
                        help='The "agr_curation_schema" LinkML release number.')
    parser.add_argument('-i', '--input_dir', help='Directory holding the six retrieval JSON files.', required=True)
    parser.add_argument('-o', '--output_dir', help='Directory for the combined files (default: input_dir).')
    args = parser.parse_args()
    output_dir = args.output_dir or args.input_dir
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

    member_releases = {}
    entity_sets = {}
    for stem, set_name in ENTITY_SETS.items():
        entity_sets.update(load_sets(find_input_file(args.input_dir, stem), [set_name], member_releases))
    association_sets = {}
    for stem, set_names in ASSOC_SETS.items():
        association_sets.update(load_sets(find_input_file(args.input_dir, stem), set_names, member_releases))

    check_no_construct_components(entity_sets['construct_ingest_set'])
    releases = set(member_releases.values())
    if len(releases) != 1:
        raise ValueError(f'Input files come from different FlyBase releases: {member_releases}')
    header = {
        'linkml_version': args.linkml_release,
        'alliance_member_release_version': releases.pop(),
    }

    # Name the outputs after the construct file so they keep its release/db suffix.
    construct_basename = os.path.basename(find_input_file(args.input_dir, 'construct_curation_'))
    outputs = (
        (entity_sets, construct_basename.replace('construct_curation_', 'construct_family_curation_')),
        (association_sets, construct_basename.replace('construct_curation_', 'construct_family_association_curation_')),
    )
    for sets, basename in outputs:
        for set_name, entries in sets.items():
            log.info(f'{basename}: {set_name} has {len(entries)} entries.')
        generate_export_file({**header, **sets}, log, os.path.join(output_dir, basename))
        log.info(f'Wrote {os.path.join(output_dir, basename)}')


if __name__ == "__main__":
    main()
