# !/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Run the "retrieve_genotypes.py" script in docker.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    run_retrieve_genotypes.py [-h]
    [-p PUB] [-i GENOTYPE_INPUT] [-f GENOTYPES_FILE]

Example (-i and -f are mutually exclusive):
    python run_retrieve_genotypes.py
    -p FBrf0123456
    -i "&agr;Tub84D[1]/Df(2L)xd Scer\GAL4[wg-Gal4] &agr;Tub84D[UAS.HA]"
    -f ./genotypes_to_make.txt

Notes:
    Give this script a single genotype name or a file with a list of names
    to get FlyBase genotype (FBgo) IDs for them. If the genotype does not
    yet exist in chado, it will be created on the spot. If the genotype is
    not yet present at the Alliance, it will be exported there immediately.
    A genotype name should be a string of feature SGML symbols, as would be
    given in a proforma record: e.g., &agrTub67C[3]. Two alleles at the same
    locus should be separated by a "/" character; use spaces only to separate
    components at different loci. Allowed features: alleles, aberrations,
    balancers, and internal "bogus symbols" (e.g., wg[+]). Insertions (FBti)
    and constructs (FBtp) are also allowed, but discouraged.
    If any quality checks fail, processing stops and the errors are reported.

"""

import argparse
import configparser
import os
import subprocess
import sys

# Process input parameters.
parser = argparse.ArgumentParser(
    description='The "run_retrieve_genotypes.py" script gets IDs for genotypes during curation. It creates them as needed.',
    epilog="""
    Tips for specifying a genotype:
    Enclose the genotype string, a list of allele, aberration or insertion symbols, in double quotes.
    For alleles, use of "SGML" for Greek chars is preferred: e.g., "&agrTub67C[1]". However, plain text works too: e.g., "alphaTub67C[1]".
    Also allowed are "bogus symbol" features that represent an internal wild-type allele: e.g., "wg[+]".
    Alleles sharing a locus should be separated by a "/" character: e.g., "wg[1]/wg[l-17]".
    An allele and deficiency can share a locus: e.g., "wg[1]/Df(2L)xyz".
    Use spaces to separate components of different loci: e.g., "wg[1]/wg[+] sona[18]/sona[18]".
    By convention, transgenic alleles for a gene are listed in their own loci, separate from classical allele loci: e.g., "wg[1]/wg[1] UAS-wg[HA]".
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('-p', '--pub', help='The FBrf ID for the publication.', required=True)
run_mode = parser.add_mutually_exclusive_group(required=True)
run_mode.add_argument('-i', '--genotype_input', help='The genotype name to get or create.', required=False)
run_mode.add_argument('-f', '--genotypes_file', help='A file of genotype names to get or create.', required=False)

try:
    args = parser.parse_args()
    genotype_input = args.genotype_input
    genotype_input_file = args.genotypes_file
    fbrf_pub_id = args.pub
except SystemExit as e:
    print('ERROR: Must supply two arguments: -p/--pub (FBrf ID), and one of -i/--genotype_input or -f/--genotypes_file.')
    sys.exit(e.code)

# Open config for chado communication.
config = configparser.ConfigParser()
config.read('/data/credentials/production/config.cfg')
server = config['chiacur']['Server']
database = config['chiacur']['Database']
user = config['chiacur']['User']
pg_pwd = config['chiacur']['PGPassword']
agr_token = config['chiacur']['AllianceCurationAPIToken']

# Construct command for running script in docker.
command = 'rm -f ./genotypes_retrieved*.report && '
command += 'rm -f ./genotypes_retrieved*.log && '
command += 'docker run --user $(id -u):$(id -g) '
command += '-v ./:/src/input/ '
command += '-v ./:/src/output/ '
command += '-v ./:/src/logs/ '
command += f'--network=\"host\" -e SERVER={server} '
command += f'-e DATABASE={database} '
command += f'-e USER={user} '
command += f'-e PGPASSWORD={pg_pwd} '
command += '-e RELEASE=production '
command += '--entrypoint /usr/bin/python3 test_export_to_linkml '    # BOB: change bulk_update docker image for prod
command += f'/src/retrieve_genotypes.py -p {fbrf_pub_id} '
if genotype_input_file:
    command += f'-f /src/input/{genotype_input_file} '
else:
    command += f'-i \"{genotype_input}\" '

# For debugging.
print(command)
subprocess.run(["bash", "-c", command])

# Using devnull to suppress confusing matplotlib and bioservices warnings that I can't resolve.
# with open(os.devnull, 'w') as devnull:
#     subprocess.run(["bash", "-c", command], stdout=devnull, stderr=devnull)

report = open('./genotypes_retrieved.report', 'r')
for line in report:
    print(line.rstrip())
