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
import csv
import re
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker, aliased
from harvdev_utils.psycopg_functions import set_up_db_reading
from harvdev_utils.char_conversions import sgml_to_plain_text
from harvdev_utils.production import (
    Cv, Cvterm, Cvtermsynonym, Db, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureDbxref, FeaturePub, FeatureRelationship, FeatureSynonym, Pub, Synonym
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

    def test_query(self, session):
        counter = 0
        filters = (
            Feature.uniquename.op('~')(self.regex['gene']),
            Feature.uniquename == 'FBgn0284084',
            Feature.name == 'wg',
        )
        results = session.query(Feature).\
            filter(*filters).\
            distinct()
        for result in results:
            self.log.info(f'Found this gene: name={result.name}, uniquename={result.uniquename}')
            counter += 1
        self.log.info(f'Found {counter} results.')
        return

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
        self.get_transgenic_allele_ids(session)
        self.get_in_vitro_allele_ids(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_indirect_progenitor_insertions(self, session):
        """Find FBti insertions associated_with the progenitor allele of an allele."""
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
                self.fb_data_entities[result.allele.feature_id]['sbj_rel_ids_by_type'][rel_type].append(rel_id)
            except:
                self.fb_data_entities[result.allele.feature_id]['sbj_rel_ids_by_type'][rel_type] = [rel_id]
            counter += 1
            self.log.info(f'Found {counter} indirect FBal-FBti associations via progenitor alleles.')
        return

    # Define get_datatype_data() for the AGMDiseaseHandler.
    def get_datatype_data(self, session):
        self.log.info(f'GET FLYBASE {self.datatype.upper()} DATA FROM CHADO.')
        self.get_entities(session)
        self.get_indirect_progenitor_insertions(session)
        self.get_entity_relationships(session, 'object', rel_type='partof', entity_type='variation')
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='construct', entity_regex=self.regex['construct'])
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='insertion', entity_regex=self.regex['insertion'])
        self.get_entity_relationships(session, 'subject', rel_type='progenitor', entity_type='insertion', entity_regex=self.regex['insertion'])
        return

    def run(self, session):
        self.test_query(session)
        self.get_general_data(session)
        self.get_datatype_data(session)
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
