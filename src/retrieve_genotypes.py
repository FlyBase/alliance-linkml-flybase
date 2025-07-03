# !/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Retrieve (get or create) a genotype for curation.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    retrieve_genotypes.py [-h] [-v VERBOSE] [-t TESTING] [-c CONFIG]
    [-p PUB] [-i GENOTYPE_INPUT] [-f GENOTYPES_FILE]

Example:
    python retrieve_genotypes.py -v -t
    -c /foo/bar/config.cfg
    -p FBrf0123456
    -i "&agr;Tub84D[1]/Df(2L)x Scer\GAL4[wg-Gal4] &agr;Tub84D[UAS.HA]"
    -f ./genotypes_to_make.txt

Notes:
    Given a genotype name (a string of feature SGML symbols) in the command
    line (enclosed in double quotes), or, a file of genotype names (one per
    line), this script will get the existing genotype from chado, or, create a
    new genotype in chado. Allowed features: alleles, aberrations, balancers,
    and internal "bogus symbols" (e.g., wg[+]). Insertions (FBti) and
    constructs (FBtp) are allowed but discouraged. No genotype is created if
    any quality checks fail.
    The script returns genotype IDs and descriptions for all input genotypes.
    Writing to chado entails writing to genotype and feature_genotype tables,
    as well as assigning an FBgo ID (genotype_dbxref), and a current synonym
    (genotype_synonym). For each new genotype, a unique string for the
    genotype is built from the IDs of its component features. This is saved in
    the genotype.description. Note that feature_genotype cgroup and rank
    values reflect the alphabetical sorting of the genotype component SGML
    symbols.
    This script is designed to be used by curators during curation of genotype
    associated data (starting with disease ontology annotations).
    FUTURE: when a new genotype is created in chado, it is then immediately
    created at the Alliance as well (in the persistent store), to support
    curation using their interface.

"""

import argparse
# import datetime
import json
import os
import requests
import sys
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound
from harvdev_utils.production import (
    Feature, FeaturePub, Pub
)
from harvdev_utils.genotype_utilities import GenotypeAnnotation
from harvdev_utils.char_conversions import (
    sgml_to_plain_text, sub_sup_sgml_to_plain_text
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Important label for output files.
report_label = 'genotypes_retrieved'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
output_dir = set_up_dict['output_dir']
genotype_report_filename = '/src/logs/genotypes_retrieved.report'
log = set_up_dict['log']
TESTING = set_up_dict['testing']
AGR_TOKEN = os.environ['ALLIANCETOKEN']
if not AGR_TOKEN:
    raise ValueError("ALLIANCETOKEN environment variable is required")
AGR_BASE_URL = os.environ.get('AGR_BASE_URL', 'https://curation.alliancegenome.org')

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
ENGINE = create_engine(engine_var_rep)
insp = inspect(ENGINE)

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
run_mode = parser.add_mutually_exclusive_group(required=True)
run_mode.add_argument('-i', '--genotype_input', help='The genotype name to get or create.', required=False)
run_mode.add_argument('-f', '--genotypes_file', help='A file of genotype names to get or create.', required=False)
parser.add_argument('-p', '--pub_id', help='The FBrf ID for the publication.', required=True)
parser.add_argument('--relax', action='store_true', help='Relax stringency to allow for processing of genotype input with warnings.', required=False)

# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
try:
    args, extra_args = parser.parse_known_args()
    log.info(f'Parsing args specific to this script: {args}.')
    GENOTYPE_INPUT = args.genotype_input
    GENOTYPE_FILE = args.genotypes_file
    FBRF_PUB_ID = args.pub_id
    RELAX = args.relax
except SystemExit as e:
    log.error('ERROR: Must supply two arguments: -p/--pub (FBrf ID), and one of -i/--genotype_input or -f/--genotypes_file.')
    sys.exit(e.code)


# The main process.
def main():
    """Run the steps for checking and updating genotypes."""
    global GENOTYPE_INPUT
    global GENOTYPE_FILE
    global FBRF_PUB_ID
    global AGR_TOKEN
    log.info('Running script "{}"'.format(__file__))
    log.info(f'Print genotype report to this file location: {genotype_report_filename}')
    log.info('STARTED MAIN FUNCTION.\n')
    if GENOTYPE_INPUT:
        genotype_input_list = [GENOTYPE_INPUT]
    else:
        try:
            genotype_file_contents = open(GENOTYPE_FILE, 'r')
            genotype_input_list = [i.strip() for i in genotype_file_contents if i.strip() != '']
        except FileNotFoundError:
            log.error(f'Cannot open "{GENOTYPE_FILE}". Make sure the file is in directory mounted to docker /src/input/')
            raise FileNotFoundError
    genotype_handler_instance = GenotypeHandler(genotype_input_list, FBRF_PUB_ID, AGR_TOKEN, RELAX)
    db_transaction(genotype_handler_instance)
    log.info('ENDED MAIN FUNCTION.\n')


class GenotypeHandler(object):
    """This object processes genotype name inputs and gets or creates chado genotypes."""
    def __init__(self, genotype_input_list, fbrf_pub_id, agr_token, relaxed_stringency):
        """Create the GenotypeHandler object.

        Args:
            genotype_input_list (list): A list of genotype names (allele SGML symbols).
            fbrf_pub_id (str): The FBrf ID for the genotype.
            agr_token (str): The Alliance API token required for interacting with the persistent store.
            relaxed_stringency (bool): If True, GenotypeAnnotations with warnings (but not errors) will be processed. If False, not processed.

        Returns:
            A GenotypeHandler object.

        """
        self.genotype_input_list = genotype_input_list
        self.fbrf_pub_id = fbrf_pub_id
        self.pub_id = None
        self.pub_associated_feature_ids = []    # List of current and public features having FB IDs associated with this pub.
        self.genotype_annotations = []          # List of GenotypeAnnotation objects generated from the input list.
        self.uname_genotype_annotations = {}    # uniquename-keyed GenotypeAnnotations (for grouping redundant entries).
        self.agr_token = agr_token
        self.relaxed_stringency = relaxed_stringency

    def get_pub_id(self, session):
        """Get pub.pub_id for given FBrf ID."""
        log.info(f'Look up chado pub for "{self.fbrf_pub_id}".')
        filters = (
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(r'^FBrf[0-9]{7}$'),
            Pub.uniquename == self.fbrf_pub_id,
        )
        try:
            result = session.query(Pub).filter(*filters).one()
            self.pub_id = result.pub_id
            log.info(f'Found ONE current pub in chado for "{self.fbrf_pub_id}".')
        except NoResultFound as ne:
            log.error(f'Found ZERO current pubs in chado for "{self.fbrf_pub_id}": {ne}.')
            raise NoResultFound
        except MultipleResultsFound as me:
            log.error(f'Found MANY current pubs in chado for "{self.fbrf_pub_id}": {me}.')
            raise MultipleResultsFound
        return

    def get_pub_associated_features(self, session):
        """Get features associated with the input pub."""
        log.info('Get features associated with the input pub.')
        filters = (
            FeaturePub.pub_id == self.pub_id,
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
        )
        results = session.query(Feature).\
            select_from(Feature).\
            join(FeaturePub, (FeaturePub.feature_id == Feature.feature_id)).\
            filter(*filters).\
            distinct()
        for result in results:
            self.pub_associated_feature_ids.append(result.feature_id)
            # log.debug(f'Found this feature_id={result.feature_id}, {result.name} ({result.uniquename})')
        log.info(f'Found {len(self.pub_associated_feature_ids)} current public features associated with "{self.fbrf_pub_id}".')
        return

    def parse_genotype_input_list(self, session):
        """Parse input genotype names into GenotypeAnnotation objects."""
        log.info('Parse input genotype names into GenotypeAnnotation objects.\n\n\n')
        log.info(f'Have {len(self.genotype_input_list)} genotypes to process.')
        for name in self.genotype_input_list:
            geno_anno = GenotypeAnnotation(name, session, log, self.pub_id)
            self.genotype_annotations.append(geno_anno)
        return

    def check_pub_associations(self):
        """Check that genotype features are already associated with the input pub FBrf ID."""
        log.info('Check that genotype features are already associated with the input pub FBrf ID.')
        no_counter = 0
        for geno_anno in self.genotype_annotations:
            geno_specific_features = geno_anno.features.keys() - set(self.pub_associated_feature_ids)
            for feature_id in geno_specific_features:
                input_symbol = geno_anno.features[feature_id]['input_symbol']
                # Relaxed feature-pub constraint.
                geno_anno.warnings.append(f'"{input_symbol}" is not associated with {self.fbrf_pub_id}')
                log.warning(f'"{input_symbol}" is not associated with {self.fbrf_pub_id}')
                # Stringent feature-pub constraint.
                # geno_anno.errors.append(f'"{input_symbol}" is not associated with {self.fbrf_pub_id}')
                # log.error(f'"{input_symbol}" is not associated with {self.fbrf_pub_id}')
                no_counter += 1
        if no_counter > 0:
            log.error(f'Found {no_counter} listed features NOT associated with {self.fbrf_pub_id}')
        return

    def find_redundant_genotype_entries(self):
        """Find redundant genotype name inputs."""
        log.info('Find redundant genotype name inputs.')
        for geno_anno in self.genotype_annotations:
            if geno_anno.errors:
                continue
            if geno_anno.uniquename in self.uname_genotype_annotations.keys():
                modified_geno_uname = geno_anno.uniquename.replace('<up>', '[').replace('</up>', ']')
                geno_anno.errors.append(f'{geno_anno} is redundant with another input genotype, sharing this uniquename: {modified_geno_uname}')
                log.error(f'{geno_anno} is redundant with another input genotype, sharing this uniquename: {geno_anno.uniquename}')
            else:
                self.uname_genotype_annotations[geno_anno.uniquename] = geno_anno
        return

    def report_errors(self):
        """Log genotype errors."""
        log.info('Reporting all input genotypes having errors.')
        for geno_anno in self.genotype_annotations:
            if geno_anno.errors:
                log.error(f'STOP processing "{geno_anno.input_genotype_name}" due to these errors: {";".join(geno_anno.errors)}.')
            elif geno_anno.warnings and self.relaxed_stringency is False:
                log.error(f'STOP processing "{geno_anno.input_genotype_name}" due to these warnings: {";".join(geno_anno.errors)}.')
        return

    def get_or_create_genotypes(self, session):
        """Find genotypes that already exist in chado, or make them if new."""
        log.info('Find genotypes that already exist in chado, or make them if new.')
        known_counter = 0
        new_counter = 0
        unclassified_counter = 0
        newly_created_counter = 0
        for geno_anno in self.uname_genotype_annotations.values():
            if geno_anno.errors:
                continue
            elif geno_anno.warnings and self.relaxed_stringency is False:
                continue
            geno_anno.get_known_or_create_new_genotype(session)
            if geno_anno.is_new is True:
                new_counter += 1
                if geno_anno.curie is not None:
                    newly_created_counter += 1
            elif geno_anno.is_new is False:
                known_counter += 1
            else:
                unclassified_counter += 1
        input_count = len(self.genotype_annotations)
        error_count = input_count - len(self.uname_genotype_annotations.keys())
        log.info(f'Had {input_count} input genotypes.')
        log.info(f'Found {known_counter} known genotypes.')
        log.info(f'Found {new_counter} new genotypes: created {newly_created_counter} new genotypes in chado.')
        log.info(f'Ignored {error_count} genotypes having errors detected.')
        if unclassified_counter > 0:
            log.error(f'Found {unclassified_counter} UNCLASSIFIED genotypes.')
        return

    def sync_with_alliance(self):
        """Synchronize genotypes at the Alliance."""
        log.info('Synchronize genotypes at the Alliance.')
        for geno_anno in self.uname_genotype_annotations.values():
            if geno_anno.errors:
                continue
            agr_curie = f'FB:{geno_anno.curie}'
            log.debug(f'Check Alliance for {agr_curie}: {geno_anno}')
            genotype_at_alliance = False
            get_url = f'{AGR_BASE_URL}/api/agm/{agr_curie}'
            headers = {
                'accept': 'application/json',
                'Authorization': f'Bearer {self.agr_token}',
            }
            get_response = requests.get(get_url, headers=headers)
            log.debug(f'Got this raw response looking for {agr_curie} at the Alliance:\n{get_response.text}')
            if get_response.status_code == 200:
                try:
                    data = get_response.json()
                    if 'primaryExternalId' in data['entity']:
                        mod_entity_id = get_response.json()['entity']['primaryExternalId']
                        log.debug(f'SUCCESS: Found {mod_entity_id} at the Alliance.')
                        genotype_at_alliance = True
                    elif 'modEntityId' in data['entity']:
                        mod_entity_id = get_response.json()['entity']['modEntityId']
                        log.debug(f'SUCCESS: Found {mod_entity_id} at the Alliance.')
                        genotype_at_alliance = True
                    else:
                        log.error('FAILURE: Got a response but could not find ID attribute.')
                        raise
                except KeyError:
                    log.debug(f'FAILURE: Could not find {agr_curie} at the Alliance.')
            else:
                log.error(f'FAILURE: Lookup of {agr_curie} did not return any response from the Alliance API.')
                raise
            if genotype_at_alliance is False:
                genotype_display_name = sub_sup_sgml_to_plain_text(geno_anno.uniquename)
                genotype_display_name = sgml_to_plain_text(genotype_display_name)
                log.debug(f'Load {geno_anno} into the Alliance using this name: {genotype_display_name}')
                linkml_genotype = {
                    'type': 'AffectedGenomicModel',
                    'subtype': {
                        'obsolete': False,
                        'internal': False,
                        'name': 'genotype'
                    },
                    'taxon': {
                        'obsolete': False,
                        'internal': False,
                        'curie': 'NCBITaxon:7227'
                    },
                    'primaryExternalId': agr_curie,
                    'name': genotype_display_name,
                    'createdBy': {
                        'obsolete': False,
                        'internal': False,
                        'uniqueId': 'FB:FB_curator'
                    },
                    'obsolete': False,
                    'internal': True,
                    'dataProvider': {
                        'obsolete': False,
                        'internal': False,
                        'abbreviation': 'FB'
                    }
                }
                json_data = json.dumps(linkml_genotype)
                log.debug(f'Have this LinkML AGM genotype JSON:\n{json_data}')
                post_url = f'{AGR_BASE_URL}/api/agm/'
                post_headers = {
                    'Content-Type': 'application/json',
                    'accept': 'application/json',
                    'Authorization': f'Bearer {self.agr_token}',
                }
                post_response = requests.post(post_url, headers=post_headers, data=json_data)
                log.debug(f'Got this raw response posting {agr_curie} at the Alliance:\n{post_response.text}')
                if post_response.status_code == 200:
                    log.debug('SUCCESS IN POSTING AGM.')
                else:
                    log.debug(f'Status code = {post_response.status_code}')
                    log.error('FAILURE TO POST AGM.')
        return

    def print_curator_genotype_report(self):
        """Print out genotype report."""
        log.info('Print out genotype report.')
        report = open(genotype_report_filename, 'w')
        lines_to_write = []
        for geno_anno in self.genotype_annotations:
            lines_to_write.append(f'\nINPUT GENOTYPE NAME: {geno_anno.input_genotype_name}')
            if geno_anno.curie:
                lines_to_write.append(f'\tGENOTYPE_ID: {geno_anno.genotype_id}')
                lines_to_write.append(f'\tCURIE: FB:{geno_anno.curie}')
            else:
                lines_to_write.append('\tCURIE:')
            if geno_anno.errors:
                status = 'ERRORS FOUND'
            elif geno_anno.is_new is True:
                status = 'NEW GENOTYPE CREATED'
            elif geno_anno.is_new is False:
                status = 'KNOWN CHADO GENOTYPE'
            lines_to_write.append(f'\tSTATUS: {status}')
            if geno_anno.uniquename:
                uniquename_output = geno_anno.uniquename.replace('<up>', '[').replace('</up>', ']')
            else:
                uniquename_output = ''
            lines_to_write.append(f'\tUNIQUENAME: {uniquename_output}')
            lines_to_write.append(f'\tDESCRIPTION: {geno_anno.description}')
            for note in geno_anno.notes:
                lines_to_write.append(f'\tNOTE: {note}')
            for warning in geno_anno.warnings:
                lines_to_write.append(f'\tWARNING: {warning}')
            for error in geno_anno.errors:
                lines_to_write.append(f'\tERROR: {error}')
            lines_to_write.append('')
        for line in lines_to_write:
            report.write(f'{line}\n')
        return

    def run(self, session):
        """Run all GenotypeHandler methods in the correct order."""
        self.get_pub_id(session)
        self.get_pub_associated_features(session)
        self.parse_genotype_input_list(session)
        self.check_pub_associations()
        self.find_redundant_genotype_entries()
        self.report_errors()
        self.get_or_create_genotypes(session)
        self.sync_with_alliance()
        self.print_curator_genotype_report()
        return


def db_transaction(object_to_execute):
    """Interact with the chado database given an object that has a "run()" method.

    Args:
        object_to_execute: Some object that has an SQL ORM "run()" method.

    Returns:
        An sqlalchemy session that can be used later on for updating the queried objects.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    global ENGINE
    global TESTING
    Session = sessionmaker(bind=ENGINE)
    session = Session()
    try:
        object_to_execute.run(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    if TESTING is True:
        log.info('Since "testing" is True, rolling back all transactions.')
        session.rollback()
    else:
        log.info('Since "testing" is False, committing transactions.')
        session.commit()
    return


if __name__ == "__main__":
    main()
