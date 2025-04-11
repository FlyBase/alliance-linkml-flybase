# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Map alleles to more representative insertions for Alliance reporting in genotypes.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    map_alleles_to_insertions_for_alliance_genotypes.py [-h] [-v VERBOSE] [-c CONFIG] [-t TESTING]

Example:
    python map_alleles_to_insertions_for_alliance_genotypes.py -v -t
    -c /path/to/config.cfg

Notes:
    This script reviews alleles and determines if each one is better represented
    in a genotype by an associated FBti insertion. If so, it makes a
    feature_relationships (FBal - TBD - FBti). This script should be run after
    each epicycle proforma load.

"""

# Notes:
# feature_lookup dicts have these keys: feature_id, uniquename, curie, is_obsolete, type, organism_id, name, symbol, exported.

import argparse
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker, aliased
from harvdev_utils.psycopg_functions import set_up_db_reading
# from harvdev_utils.char_conversions import sgml_to_plain_text
from harvdev_utils.production import (
    Cvterm, Feature, FeatureRelationship
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
args, extra_args = parser.parse_known_args()
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

    # Add methods to be run by get_general_data() below.
    # Placeholder

    # Define get_general_data() for the AlleleMapper.
    def get_general_data(self, session):
        """Extend the method for the AGMDiseaseHandler."""
        self.log.info('GET GENERAL FLYBASE DATA FROM CHADO.')
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.get_key_cvterm_sets(session)
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'gene', 'insertion', 'construct', 'variation'])
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
        results = session.query(allele, al_ti).\
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

    # Define get_datatype_data() for the AlleleMapper.
    def get_datatype_data(self, session):
        """Get FlyBase data from chado."""
        self.log.info('GET RELEVANT ALLELE DATA FROM CHADO.')
        self.get_entities(session)
        self.get_indirect_progenitor_insertion_relationships(session)
        self.get_entity_relationships(session, 'object', rel_type='partof', entity_type='variation')
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='construct', entity_regex=self.regex['construct'])
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='insertion', entity_regex=self.regex['insertion'])
        self.get_entity_relationships(session, 'subject', rel_type='progenitor', entity_type='insertion', entity_regex=self.regex['insertion'])
        return

    # Add methods to be run by synthesize_data() below.
    def map_allele_to_insertion(self):
        """Map alleles to insertions, if applicable."""
        self.log.info('Map alleles to insertions, if applicable.')
        sample_alleles = [
            'Nmnat[W129G.EGFP]',    # FBti is both associated_with and progenitor.
            # 'Mkp3[5]',              # FBti + ARG.
            # 'mei-P26[fs1]',         # FBti + ARG.
            # 'TrpA1[-ACD-G4]',       # FBti is both associated_with and progenitor (indirectly via "TrpA1[-CD-G4]").
            # 'chrb[180]',            # Many FBti.
            # 'sd[ETX81]',            # Many FBti.
            # 'neb[k06334]',          # Many FBti.
            # 'Gpdh1[AKO107]',        # Many FBti.
            # 'gt[1]',                # Many FBti.
            # 'ac[Hw-BS]',            # Many FBti.
            # 'ovo[yct]',             # Many FBti.
            # 'sn[w]',                # Many FBti.
            # 'eyg[P20MD1]',          # Many FBti.
            # 'twin[KG00877]',        # Many FBti.
            # 'sd[+58b]',             # Many FBti.
            # 'bmm[EY06577]',         # Simple at-locus FBti.
            # 'Scer_GAL4[sLNvs]',     # Simple trap FBti.
            # 'Arf6[EP2612]',         # Shared trap FBti.
            # 'CG8155[EP2612]',       # Shared trap FBti.
            # 'lbe[UAS.cJa]',         # Allele should be mapped to construct-insertion.
            # 'wg[l-12]',             # Classical mutation.
        ]
        input_counter = 0
        mapped_counter = 0
        for allele in self.fb_data_entities.values():
            # BOB - just for testing.
            if allele.chado_obj.name not in sample_alleles:
                continue
            input_counter += 1
            arg_rels = allele.recall_relationships(self.log, entity_role='object', rel_types='partof',
                                                   rel_entity_types=self.feature_subtypes['variation'])
            fbtp_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='associated_with',
                                                    rel_entity_types=self.feature_subtypes['construct'])
            fbti_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='associated_with',
                                                    rel_entity_types=self.feature_subtypes['insertion'])
            prog_rel_types = ['progenitor', 'indirect_progenitor_insertion_rels']
            prog_fbti_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types=prog_rel_types,
                                                         rel_entity_types=self.feature_subtypes['insertion'])
            single_fbti_feature_id = None
            notes = []
            if arg_rels:
                notes.append('Has ARG(s)')
            if fbtp_rels:
                notes.append('Has FBtp(s)')
            if fbti_rels:
                distinct_fbti_feature_ids = list(set([i.chado_obj.object_id for i in fbti_rels]))
                self.log.debug(f'{allele} has these "associated_with" FBti feature_ids: {distinct_fbti_feature_ids}')
                distinct_fbti_progenitor_feature_ids = set([i.chado_obj.object_id for i in prog_fbti_rels])
                self.log.debug(f'{allele} has these "progenitor" FBti feature_ids: {distinct_fbti_progenitor_feature_ids}')
                if len(distinct_fbti_feature_ids) == 0:
                    notes.append('ERROR: Has ZERO FBti(s)')
                elif len(distinct_fbti_feature_ids) > 1:
                    notes.append('Has MANY FBti(s)')
                else:
                    single_fbti_feature_id = distinct_fbti_feature_ids[0]
                    self.log.debug(f'GILLY: Have this feature_id: {single_fbti_feature_id}')
                    self.log.debug(f'GILLY: Have these prog ids: {distinct_fbti_progenitor_feature_ids}')
                    eval = single_fbti_feature_id in distinct_fbti_progenitor_feature_ids
                    self.log.debug(f'GILLY: feature_id in prog ids? {eval}')
                    if eval is True:
                        notes.append('Associated FBti is also a progenitor FBti')
                    else:
                        notes.append(f'Found one FBti: {single_fbti_feature_id}')
            else:
                notes.append('No FBti')
            self.log.debug(f'GILLY: {allele} has these notes: {notes}')
            allele.single_fbti_feature_id = single_fbti_feature_id
            if single_fbti_feature_id:
                mapped_counter += 1
                insertion = self.feature_lookup[allele.single_fbti_feature_id]
                mapping_str = f'\t{allele.chado_obj.uniquename}\t{allele.chado_obj.name}\t{insertion["uniquename"]}\t{insertion["name"]}'
            self.log.debug(f'Assessing "{allele.chado_obj.name}" ({allele.chado_obj.uniquename}).')
            self.log.debug(f'Found {len(arg_rels)} ARG relationships.')
            self.log.debug(f'Found {len(fbtp_rels)} FBtp relationships.')
            self.log.debug(f'Found {len(fbti_rels)} FBti relationships.')
            self.log.debug(f'Found {len(prog_fbti_rels)} progenitor FBti relationships.')
            if allele.single_fbti_feature_id:
                self.log.debug(f'BOB MAPPING: {mapping_str})')
            else:
                self.log.debug(f'BOB NOPE: {allele} could not be mapped to an associated insertion: {"; ".join(notes)})')
        self.log.info(f'Mapped {mapped_counter}/{input_counter} alleles to a single FBti insertion unambiguously.')
        return

    # Define synthesize_data() for the AlleleMapper.
    def synthesize_data(self):
        """Synthesize data."""
        self.log.info('SYNTHESIZE DATA.')
        self.map_allele_to_insertion()
        return

    def run(self, session):
        """Run all methods in sequence."""
        self.log.info('Run all methods in sequence.')
        self.get_general_data(session)
        self.get_datatype_data(session)
        self.synthesize_data()
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
    load_testing = True    # Hard-coding this because "testing" means something different for Alliance data handlers.
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
    if load_testing is True:
        log.info('Since "testing" is True, rolling back all transactions.')
        session.rollback()
    else:
        log.info('Since "testing" is False, committing transactions.')
        session.commit()
    return


if __name__ == "__main__":
    main()
