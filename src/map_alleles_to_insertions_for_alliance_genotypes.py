# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Map alleles to more representative insertions for Alliance reporting in genotypes.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    map_alleles_to_insertions_for_alliance_genotypes.py [-h] [--commit]
    [-v VERBOSE] [-c CONFIG] [-t TESTING]

Example:
    python map_alleles_to_insertions_for_alliance_genotypes.py --commit
    -v -t -c /path/to/config.cfg

Notes:
    For each FBal allele, this script determines if it is better represented by
    a single at-locus FBti insertion, or, an FBab aberration. If so, the script
    makes a feature_relationship (FBal-is_represented_at_alliance_as-FBti).
    There is no mapping if there are many insertions, or if the allele has other
    lesions.
    Notes on usage:
    This script should be run after each epicycle proforma load.
    This script flushes/replaces of represented_at_alliance_as relationships.
    One must specify "--commit" in the command line to commit the changes!
"""

# Notes:
# feature_lookup dicts have these keys: feature_id, uniquename, curie, is_obsolete, type, organism_id, name, symbol, exported.

import argparse
import re
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker, aliased
from harvdev_utils.psycopg_functions import set_up_db_reading
# from harvdev_utils.char_conversions import sgml_to_plain_text
from harvdev_utils.chado_functions import get_or_create
from harvdev_utils.production import (
    Cvterm, Feature, FeatureRelationship, Featureprop
)
from allele_handlers import AlleleHandler
import fb_datatypes

# Data types handled by this script.
REPORT_LABEL = 'map_alleles_to_insertions_for_alliance_genotypes'

# Now proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
log = set_up_dict['log']
TESTING = set_up_dict['testing']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('--commit', action='store_true', help="Commit writing of new feature_relationships")
args, extra_args = parser.parse_known_args()
COMMIT = args.commit
log.info(f'These args are handled by this specific script: {args}')
log.info(f'These args are handled by modules: {extra_args}')

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
ENGINE = create_engine(engine_var_rep)
inspect(ENGINE)


# The main process.
def main():
    """Run the main function."""
    log.info(f'Running script "{__file__}"')
    log.info('Started main function.')
    allele_mapper = AlleleMapper(log, TESTING)
    run_mapper(allele_mapper)
    log.info('Ended main function.\n')


class AlleleMapper(AlleleHandler):
    """An object that maps alleles to insertions for Alliance genotype reporting."""

    def __init__(self, log, testing):
        """Create the AlleleMapper object."""
        super().__init__(log, testing)
        self.allele_name_lookup = {}    # Allele-name-keyed dict of feature dicts.

    # Method to initially flush all "is_represented_at_alliance_as" feature_relationships.
    def initial_flush(self, session):
        """Flush existing "is_represented_at_alliance_as" feature_relationships to create a blank slate."""
        self.log.info('Flush existing "is_represented_at_alliance_as" feature_relationships to create a blank slate.')
        filters = (Cvterm.name == 'is_represented_at_alliance_as', )
        fr_type_cvterm_id = session.query(Cvterm).filter(*filters).one().cvterm_id
        results = session.query(FeatureRelationship).filter(FeatureRelationship.type_id == fr_type_cvterm_id).distinct()
        counter = 0
        for result in results:
            counter += 1
        self.log.info(f'There are currently {counter} "is_represented_at_alliance_as" feature_relationships (before updating).')
        session.query(FeatureRelationship).filter(FeatureRelationship.type_id == fr_type_cvterm_id).delete()
        self.log.info('Flushed all "is_represented_at_alliance_as" feature_relationships before updating.')
        return

    # Add methods to be run by get_general_data() below.
    def build_allele_name_lookup(self):
        """Build name-keyed dict of alleles."""
        self.log.info('Build name-keyed dict of alleles.')
        for feature in self.feature_lookup.values():
            if feature['type'] == 'allele' and feature['is_obsolete'] is False:
                self.allele_name_lookup[feature['name']] = feature
        self.log.info(f'Have {len(self.allele_name_lookup)} current alleles in the allele-by-name lookup.')
        return

    # Define get_general_data() for the AlleleMapper.
    def get_general_data(self, session):
        """Extend the method for the AGMDiseaseHandler."""
        self.log.info('GET GENERAL FLYBASE DATA FROM CHADO.')
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.get_key_cvterm_sets(session)
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'gene', 'insertion', 'construct', 'variation'])
        self.build_allele_name_lookup()
        self.build_allele_gene_lookup(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_indirect_progenitor_insertion_relationships(self, session):
        """Find FBti insertions associated_with the progenitor allele of an allele."""
        self.log.info('Find FBti insertions associated_with the progenitor allele of an allele.')
        allele = aliased(Feature, name='allele')
        progenitor_allele = aliased(Feature, name='progenitor_allele')
        insertion = aliased(Feature, name='insertion')
        allele_rel_type = aliased(Cvterm, name='allele_rel_type')
        ins_rel_type = aliased(Cvterm, name='ins_rel_type')
        al_al = aliased(FeatureRelationship, name='al_al')
        al_ti = aliased(FeatureRelationship, name='al_ti')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            progenitor_allele.is_obsolete.is_(False),
            progenitor_allele.uniquename.op('~')(self.regex['allele']),
            insertion.is_obsolete.is_(False),
            insertion.uniquename.op('~')(self.regex['insertion']),
            allele_rel_type.name == 'progenitor',
            ins_rel_type.name == 'associated_with',
        )
        results = session.query(allele, progenitor_allele, al_ti, insertion).\
            select_from(allele).\
            join(al_al, (al_al.subject_id == allele.feature_id)).\
            join(progenitor_allele, (progenitor_allele.feature_id == al_al.object_id)).\
            join(allele_rel_type, (allele_rel_type.cvterm_id == al_al.type_id)).\
            join(al_ti, (al_ti.subject_id == progenitor_allele.feature_id)).\
            join(insertion, (insertion.feature_id == al_ti.object_id)).\
            join(ins_rel_type, (ins_rel_type.cvterm_id == al_ti.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            result_str = f'Allele {result.allele.name} ({result.allele.uniquename}) has progenitor '
            result_str += f'allele {result.progenitor_allele.name} ({result.progenitor_allele.uniquename}) '
            result_str += f'which is associated_with {result.insertion.name} ({result.insertion.uniquename})'
            # self.log.debug(f'{result_str}')
            rel = fb_datatypes.FBRelationship(result.al_ti, 'feature_relationship')
            rel_type = 'indirect_progenitor_insertion_rels'
            rel_id = result.al_ti.feature_relationship_id
            # Index the feature_relationship by id.
            self.fb_data_entities[result.allele.feature_id].rels_by_id[rel_id] = rel
            # Index the feature_relationship_id by type.
            try:
                self.fb_data_entities[result.allele.feature_id].sbj_rel_ids_by_type[rel_type].append(rel_id)
            except KeyError:
                self.fb_data_entities[result.allele.feature_id].sbj_rel_ids_by_type[rel_type] = [rel_id]
            counter += 1
        self.log.info(f'Found {counter} indirect FBal-FBti associations via progenitor alleles.')
        return

    def flag_alleles_for_which_name_check_should_be_ignored(self, session):
        """Flag alleles for which name check should be ignored."""
        self.log.info('Flag alleles for which name check should be ignored.')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['allele']),
            Cvterm.name == 'internal_notes',
            Featureprop.value == 'FTA: genotype conversion - ignore unconventional name warning'
        )
        results = session.query(Feature).\
            select_from(Feature).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for allele in results:
            self.fb_data_entities[allele.feature_id].ignore_atypical_name = True
            counter += 1
        self.log.info(f'Flagged {counter} alleles for which name check should be ignored.')
        return

    # Define get_datatype_data() for the AlleleMapper.
    def get_datatype_data(self, session):
        """Get FlyBase data from chado."""
        self.log.info('GET RELEVANT ALLELE DATA FROM CHADO.')
        self.get_entities(session)
        self.get_indirect_progenitor_insertion_relationships(session)
        self.get_entity_relationships(session, 'object', rel_type='partof', entity_type='variation')
        self.get_entity_relationships(session, 'object', rel_type='associated_with', entity_type='aberration', entity_regex=self.regex['aberration'])
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='construct', entity_regex=self.regex['construct'])
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='insertion', entity_regex=self.regex['insertion'])
        self.get_entity_relationships(session, 'subject', rel_type='progenitor', entity_type='insertion', entity_regex=self.regex['insertion'])
        self.flag_alleles_for_which_name_check_should_be_ignored(session)
        return

    # Additional methods for making the allele-insertion calls.
    def extract_allele_suffix_from_insertion_name(self, allele_feature_id, insertion_feature_id):
        """Check that the allele name is conventional.

        This function checks that the name of the insertion-associated allele
        is conventional. Unconventional names are likely to represent complex
        alleles that should not be reported as insertion at the Alliance.
        This function requires that the handler has already built a feature
        lookup for alleles, genes and insertions, as well as an allele-gene
        lookup.

        Args:
            allele_feature_id (int): The feature_id for the FBal allele.
            insertion_feature_id (int): The feature_id for the FBti insertion.

        Returns:
            conventional_name (bool): True if the allele name is conventional.

        """
        conventional_name = True
        # Gather allele and insertion names and name parts.
        allele = self.feature_lookup[allele_feature_id]
        allele_name = allele['name']
        insertion = self.feature_lookup[insertion_feature_id]
        insertion_name = insertion['name']
        insertion_suffix = insertion_name.split('}')[1]
        # Getting the allele superscript is tricky.
        # The allele superscript is stuff in the square brackets after the gene name.
        # So, we remove the first instance of the gene name in the string.
        # Note that the gene name can itself have square brackets, hence the approach here.
        gene = self.feature_lookup[self.allele_gene_lookup[allele_feature_id]]
        gene_name = gene['name']
        allele_parts = allele_name.split(gene_name)
        allele_superscript = gene_name.join(allele_parts[1:]).lstrip('[').rstrip(']')    # Remove flanking brackets.
        initial_msg = 'Have these parts to assess: '
        initial_msg += f'allele_name="{allele_name}", '
        initial_msg += f'allele_superscript="{allele_superscript}", '
        initial_msg += f'insertion_name="{insertion_name}", '
        initial_msg += f'insertion_suffix="{insertion_suffix}", '
        initial_msg += f'gene_name="{gene_name}".'
        self.log.debug(f'Assess allele name. {initial_msg}')
        # Record reasons for unconventional name.
        notes = []
        double_curly_rgx = r'}.*}'
        curly_rgx = r'{.*}'
        superscript_rgx = r'(^.+)(\[[^\[]+\]$)'
        # Check 1. Return None for FBti insertion names having many curly bracket sets.
        if re.search(double_curly_rgx, insertion_name):
            conventional_name = False
            notes.append('Double curly bracket in insertion name')
        # Check 2. Return None for FBti insertion names missing any curly bracket set.
        if not re.search(curly_rgx, insertion_name):
            conventional_name = False
            notes.append('No curly bracket in insertion name')
        # Check 3. Start Gillian's checks here.
        if insertion_suffix in self.allele_name_lookup.keys():
            insertion_suffix_matches_allele_symbol = True
            # Mask double square brackets to make it easier to draw out the allele superscript.
            masked_insertion_suffix = insertion_suffix.replace('[[', 'SUBSCRIPT_START').replace(']]', 'SUBSCRIPT_END')
            masked_insertion_superscript = re.search(superscript_rgx, masked_insertion_suffix).group(2).lstrip('[').rstrip(']')    # Remove flanking brackets.
            insertion_superscript = masked_insertion_superscript.replace('SUBSCRIPT_START', '[[').replace('SUBSCRIPT_END', ']]')
            # self.log.debug(f'{insertion_name} has this superscript: {insertion_superscript}')
            if allele_superscript == insertion_superscript:
                # self.log.debug(f'PASS1: insertion_superscript == allele_superscript: allele_name={allele_name}, ins_name={insertion_name}')
                pass
            elif allele_superscript.endswith(f'-{insertion_superscript}') or allele_superscript.endswith(f'{insertion_superscript}-X'):
                # self.log.debug(f'PASS2: allele_superscript endswith insertion_superscript: allele_name={allele_name}, ins_name={insertion_name}')
                pass
            else:
                conventional_name = False
                self.log.debug(f'FAIL1: allele_superscript does not contain insertion_superscript: allele_name={allele_name}, ins_name={insertion_name}')
                notes.append('insertion-allele name mismatch A')
        else:
            insertion_suffix_matches_allele_symbol = False
            if allele_superscript == insertion_suffix:
                # self.log.debug(f'PASS3: insertion_suffix == allele_superscript: allele_name={allele_name}, ins_name={insertion_name}')
                pass
            elif allele_superscript.endswith(f'-{insertion_suffix}'):
                # self.log.debug(f'PASS4: allele_superscript endswith insertion_suffix: allele_name={allele_name}, ins_name={insertion_name}')
                pass
            else:
                conventional_name = False
                self.log.debug(f'FAIL2: allele_superscript does not contain insertion_superscript: allele_name={allele_name}, ins_name={insertion_name}')
                notes.append('insertion-allele name mismatch B')
        # Once all checks are done, print out unconventional names for debug and curator review.
        if conventional_name:
            msg = f'Conventional name for {allele_name} ({allele["uniquename"]}) '
            msg += f'associated with {insertion_name} {insertion["uniquename"]}'
            msg2 = None
            # self.log.debug(msg)
        else:
            msg = f'UNCONVENTIONAL name for {allele_name} ({allele["uniquename"]}) '
            msg += f'associated with {insertion_name} {insertion["uniquename"]}. Reasons: {";".join(notes)}'
            self.log.warning(msg)
            msg2 = f'UNCONVENTIONAL NAME: {allele_name}\t{allele["uniquename"]}'
            msg2 += f'\t{insertion_name}\t{insertion["uniquename"]}\t{insertion_suffix_matches_allele_symbol}'
        return conventional_name, msg2

    def map_alleles_to_insertions(self):
        """Map alleles to insertions, if applicable."""
        self.log.info('Map alleles to insertions, if applicable.')
        input_counter = 0
        mapped_counter = 0
        for allele in self.fb_data_entities.values():
            if allele.chado_obj.is_obsolete is True:
                continue
            self.log.debug(f'Assessing "{allele.chado_obj.name}" ({allele.chado_obj.uniquename}).')
            input_counter += 1
            # Gather feature_relationship info.
            fbab_rels = allele.recall_relationships(self.log, entity_role='object', rel_types='associated_with',
                                                    rel_entity_types=self.feature_subtypes['aberration'])
            arg_rels = allele.recall_relationships(self.log, entity_role='object', rel_types='partof',
                                                   rel_entity_types=self.feature_subtypes['variation'])
            fbtp_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='associated_with',
                                                    rel_entity_types=self.feature_subtypes['construct'])
            fbti_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='associated_with',
                                                    rel_entity_types=self.feature_subtypes['insertion'])
            prog_rel_types = ['progenitor', 'indirect_progenitor_insertion_rels']
            prog_fbti_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types=prog_rel_types)
            self.log.debug(f'Found {len(arg_rels)} ARG relationships.')
            self.log.debug(f'Found {len(fbtp_rels)} FBtp relationships.')
            self.log.debug(f'Found {len(fbti_rels)} FBti relationships.')
            self.log.debug(f'Found {len(prog_fbti_rels)} progenitor FBti relationships.')
            # Start assessment
            fbti_mappable = True
            allele_name_is_ok = True
            notes = []
            if arg_rels:
                fbti_mappable = False
                notes.append('Has ARGS')
            if fbtp_rels:
                fbti_mappable = False
                notes.append('Has FBtp')
            distinct_fbti_feature_ids = list(set([i.chado_obj.object_id for i in fbti_rels]))
            self.log.debug(f'{allele} has these "associated_with" FBti feature_ids: {distinct_fbti_feature_ids}')
            distinct_fbti_progenitor_feature_ids = set([i.chado_obj.object_id for i in prog_fbti_rels])
            self.log.debug(f'{allele} has these "progenitor" FBti feature_ids: {distinct_fbti_progenitor_feature_ids}')
            if len(distinct_fbti_feature_ids) == 0:
                fbti_mappable = False
                notes.append('Has ZERO FBti')
            elif len(distinct_fbti_feature_ids) > 1:
                fbti_mappable = False
                notes.append('Has MANY FBti(s)')
            # If there is a single FBti, apply two more tests.
            else:
                # 1. Ensure FBti is not also a progenitor.
                if distinct_fbti_feature_ids[0] in distinct_fbti_progenitor_feature_ids:
                    fbti_mappable = False
                    notes.append('Associated FBti is also a progenitor FBti')
                # 2. Look at allele name to detect complex allele.
                allele_name_is_ok, name_check_msg = self.extract_allele_suffix_from_insertion_name(allele.chado_obj.feature_id, distinct_fbti_feature_ids[0])
                if allele_name_is_ok is False and allele.ignore_atypical_name is False:
                    self.log.debug(f'{name_check_msg}\t{"; ".join(notes)}')
                    notes.append('Unconventional allele name')
                elif allele_name_is_ok is False and allele.ignore_atypical_name is True:
                    notes.append('Unconventional allele name is ok')
                    allele_name_is_ok = True
            if fbti_mappable is True and allele_name_is_ok is True:
                insertion = self.feature_lookup[distinct_fbti_feature_ids[0]]
                # Check for cases where an FBab is appropriate.
                aberration = None
                if len(fbab_rels) == 1:
                    aberration = self.feature_lookup[fbab_rels[0].subject_id]
                if aberration and aberration['name'].startswith('Df(') and insertion['name'].startswith('TI{'):
                    allele.maps_to_feature_id = aberration['feature_id']
                    mapped_counter += 1
                    mapping_str = f'\t{allele.chado_obj.uniquename}\t{allele.chado_obj.name}\t{aberration["uniquename"]}\t{aberration["name"]}'
                    self.log.debug(f'MAPPING: {mapping_str})')
                else:
                    allele.maps_to_feature_id = insertion['feature_id']
                    mapped_counter += 1
                    mapping_str = f'\t{allele.chado_obj.uniquename}\t{allele.chado_obj.name}\t{insertion["uniquename"]}\t{insertion["name"]}'
                    self.log.debug(f'MAPPING: {mapping_str})')
            else:
                self.log.debug(f'NO MAPPING: {allele} could not be mapped to an associated insertion: {"; ".join(notes)}')
        self.log.info(f'Mapped {mapped_counter}/{input_counter} current alleles to a single FBti insertion unambiguously.')
        return

    def write_new_feature_relationships(self, session):
        """Create new "is_represented_at_alliance_as" feature_relationships."""
        self.log.info('Create new "is_represented_at_alliance_as" feature_relationships.')
        new_counter = 0
        filters = (Cvterm.name == 'is_represented_at_alliance_as', )
        fr_type_cvterm_id = session.query(Cvterm).filter(*filters).one().cvterm_id
        self.log.info(f'The "is_represented_at_alliance_as" CV term corresponds to cvterm.cvterm_id={fr_type_cvterm_id}')
        for allele in self.fb_data_entities.values():
            if allele.maps_to_feature_id:
                _, _ = get_or_create(session, FeatureRelationship, subject_id=allele.chado_obj.feature_id,
                                     object_id=allele.maps_to_feature_id, type_id=fr_type_cvterm_id)
                new_counter += 1
        self.log.info(f'Created {new_counter} new "is_represented_at_alliance_as" feature_relationships.')
        return

    def run(self, session):
        """Run all methods in sequence."""
        self.log.info('Run all methods in sequence.')
        self.initial_flush(session)
        self.get_general_data(session)
        self.get_datatype_data(session)
        self.map_alleles_to_insertions()
        self.write_new_feature_relationships(session)
        return


def run_mapper(object_to_execute):
    """Run the handler.

    Args:
        object_to_execute (DataHandler): An object having a query_chado_and_export() method.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    global log
    global ENGINE
    global COMMIT
    Session = sessionmaker(bind=ENGINE)
    inspect(ENGINE)
    session = Session()
    try:
        object_to_execute.run(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during main chado query; rolling back and exiting.')
        raise
    if COMMIT is True:
        log.info('Since "commit" is True, committing transactions.')
        session.commit()
    else:
        log.info('Since "commit" is False, rolling back all transactions.')
        session.rollback()
    return


if __name__ == "__main__":
    main()
