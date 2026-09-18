#!/usr/bin/env python3
"""Report how every cross-reference prefix FlyBase emits resolves at the Alliance.

Synopsis:
    A standalone check for the resource descriptor plumbing added in FTA-263. It fetches
    the Alliance resource descriptors the same way an export run does - Cognito
    client_credentials exchange, falling back to ALLIANCETOKEN - and prints what would
    happen to every prefix in DataHandler.fb_agr_db_dict.

    Two uses:
      1. Verifying COGNITO_ADMIN_CLIENT_ID, COGNITO_ADMIN_CLIENT_SECRET and
         COGNITO_TOKEN_URL are set correctly, without running a full export. Exit code 0
         means descriptors were fetched; 1 means they were not, and the log says why.
      2. Seeing, before a submission, which prefixes are respelled, which fall back to
         "default", and which will be dropped because the Alliance does not know them.

Usage:
    python report_resource_descriptors.py [datatype]

    datatype defaults to "gene"; pass "allele", "construct", "strain" etc. to see how the
    page areas for that export would resolve.

Author(s):
    Ian Longden ianlongden@morgan.harvard.edu

"""

import logging
import sys

from handler import DataHandler
from resource_descriptors import (
    DEFAULT_PAGE_AREA, PageAreaResolver, fetch_resource_descriptors
)


def build_report(resolver: PageAreaResolver, prefixes: list, datatype: str):
    """Group prefixes by what the resolver does with them for this datatype."""
    report = {
        f'declare the "{datatype}" page, unchanged': [],
        f'fall back to "{DEFAULT_PAGE_AREA}"': [],
        'respelled to match the Alliance': [],
        'no valid page area - cross-references dropped': [],
        'unrecognized prefix - cross-references dropped': [],
    }
    for prefix in prefixes:
        canonical = resolver.canonical_prefix(prefix)
        if canonical is None:
            report['unrecognized prefix - cross-references dropped'].append(prefix)
            continue
        page_area = resolver.resolve(canonical, datatype)
        if page_area is None:
            report['no valid page area - cross-references dropped'].append(prefix)
        elif canonical != prefix:
            report['respelled to match the Alliance'].append(f'{prefix} -> {canonical} ({page_area})')
        elif page_area == datatype:
            report[f'declare the "{datatype}" page, unchanged'].append(prefix)
        else:
            report[f'fall back to "{DEFAULT_PAGE_AREA}"'].append(prefix)
    return report


def main():
    """Fetch descriptors and print the prefix report."""
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s', stream=sys.stdout)
    log = logging.getLogger('resource_descriptors')
    datatype = sys.argv[1] if len(sys.argv) > 1 else 'gene'

    descriptors = fetch_resource_descriptors(log)
    if descriptors is None:
        log.error('No resource descriptors fetched. An export run would fall back to keeping '
                  'FlyBase page areas and sending every other prefix to '
                  f'"{DEFAULT_PAGE_AREA}", without detecting unrecognized prefixes.')
        return 1

    resolver = PageAreaResolver(log, descriptors)
    prefixes = sorted(set(DataHandler.fb_agr_db_dict.values()))
    print(f'\n{len(descriptors)} Alliance resource descriptors; '
          f'{len(prefixes)} prefixes emitted by fb_agr_db_dict; datatype "{datatype}".')
    for heading, entries in build_report(resolver, prefixes, datatype).items():
        print(f'\n{heading} ({len(entries)}):')
        print('   ' + (', '.join(entries) if entries else 'none'))

    fb_pages = sorted(resolver.pages_by_prefix.get('FB', set()))
    print(f'\nFlyBase resource pages the Alliance declares ({len(fb_pages)}):')
    print('   ' + ', '.join(fb_pages))
    return 0


if __name__ == '__main__':
    sys.exit(main())
