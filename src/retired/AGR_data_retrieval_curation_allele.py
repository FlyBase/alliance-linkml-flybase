# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase alleles for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_allele.py [-h] [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_curation_allele.py -v -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script makes a JSON file conforming to LinkML specs for the curation
    (i.e., "persistent") database; distinct from ALLELE file specs for the
    original Neo4j drop-and-reload database. As the "audit_chado" table may be
    required for proper determination of "date_updated" values, use a reporting
    database with full audit_chado table present.

"""

import argparse
import datetime
import json
import strict_rfc3339
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cv, Cvterm, CvtermRelationship, Db, Dbxref, Feature, FeatureCvterm, FeatureDbxref, FeatureGenotype,
    FeaturePub, FeatureRelationship, FeatureSynonym, Featureprop, Genotype,
    Library, LibraryFeature, LibraryFeatureprop, Organism, Organismprop,
    OrganismDbxref, Phenotype, PhenotypeCvterm, Phenstatement, Pub, PubDbxref,
    Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Now proceed with generic setup.
report_label = 'allele_curation'
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
reference_assembly = set_up_dict['assembly']
input_dir = set_up_dict['input_dir']
output_filename = set_up_dict['output_filename'].replace('tsv', 'json')
log = set_up_dict['log']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-l', '--linkml_release', help='The "agr_curation_schema" LinkML release number.', required=True)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase Alleles."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')
    log.info('Output JSON file corresponds to "agr_curation_schema" release: {}'.format(linkml_release))

    # Instantiate the object, get the data, synthesize it, export it.
    allele_handler = AlleleHandler()
    db_query_transaction(allele_handler)
    allele_handler.generate_export_file()
    log.info('Ended main function.\n')


class AllianceAllele(object):
    """An allele with it's associated FlyBase data and Alliance LinkML properties."""
    def __init__(self, feature):
        """Create a base AllianceAllele object.

        Args:
            arg1 (feature): (Feature) The Feature object corresponding to the allele.

        Returns:
            An object of the AllianceAllele class.

        """
        # Attributes representing unprocessed FlyBase data.
        # Note: use attribute names that do not match an Alliance LinkML slot name.
        self.feature = feature                             # The Feature object corresponding to the FlyBase allele, insertion or aberration/balancer.
        self.organism_abbr = None                          # Will be the organism.abbreviation for the allele's species of origin.
        self.adj_organism_abbr = 'Dmel'                    # Assume allele is Dmel (classical/transgenic) unless allele is of classical type in another insect.
        self.in_vitro = False                              # Change to True if allele associated with "in vitro%" cvterm.
        self.constructs = []                               # Will be a list of FBtp IDs for this allele's constructs.
        self.dmel_insertions = []                          # Will be a list of FBti IDs for this allele's Dmel insertions.
        self.non_dmel_insertions = []                      # Will be a list of FBti IDs for this allele's non-Dmel insertions.
        self.args = []                                     # Will be a list of ARG Features.
        self.parent_gene = None                            # Will be the FBgn ID of the allele's gene.
        self.allele_of_internal_gene = False               # Will change to True if is allele of Dmel internal-type gene (e.g., origin_of_replication).
        self.curr_symbol_name = None                       # Will be the current symbol synonym.synonym_sgml, processed by sub_sup_sgml_to_html().
        self.curr_fullname = None                          # Will be the current fullname synonym.synonym_sgml, processed by sub_sup_sgml_to_html().
        self.feature_synonyms = []                         # Will be list of all FeatureSynonym objects.
        self.dbxrefs = []                                  # Will be list of dbxrefs as sql result groupings: Db, Dbxref, FeatureDbxref.
        self.alt_fb_ids = []                               # Will be list of Dbxrefs for 2o FlyBase IDs.
        self.timestamps = []                               # Add all timestamps here.
        self.fb_references = []                            # Will be list of pub_ids from feature_pub, feature_synonym.
        self.featureprops = {}                             # A CVterm-keyed dict of Featureprop lists.
        self.phenstatements = []                           # Will be a list of SQLAlchemy (Feature, Genotype, Phenotype, Cvterm, Pub) from Phenstatements.
        self.direct_libraries = []                         # Will be a list of Library objects directly related to the allele.
        self.ins_libraries = []                            # Will be a list of Library objects related to the allele via insertion (FBti).
        self.cons_libraries = []                           # Will be a list of Library objects related to the allele via construct (FBtp).
        self.sf_libraries = []                             # Will be a list of Library objects related to the allele via seq. feature (FBsf).
        # Attributes for the Alliance AuditedObjectDTO.
        self.obsolete = feature.is_obsolete                    # Will be the FlyBase value here.
        self.internal = False                                  # Change to true if not public on FlyBase.
        self.created_by_curie = 'FB:FB_curator'                # Use placeholder value since no Person object at FlyBase.
        self.updated_by_curie = 'FB:FB_curator'                # Use placeholder value since no Person object at FlyBase.
        self.date_created = None                               # Earliest timestamp.
        self.date_updated = None                               # Latest timestamp.
        # Attributes for the Alliance SubmittedObjectDTO.
        self.mod_entity_id = 'FB:{}'.format(feature.uniquename)
        self.mod_internal_id = str(self.feature.feature_id)
        self.data_provider_dto = None                          # Will be DataProviderDTO object.
        # Attributes for the Alliance BiologicalEntityDTO.
        self.taxon_curie = None                                # A string representing the NCBI taxon ID. We have no NCBI taxonID for 223 alleles.
        # Attributes for the Alliance GenomicEntityDTO.
        self.cross_reference_dtos = []                         # Report only select dbs, using AGR-accepted db_prefix.
        # Attributes for the Alliance AlleleDTO.
        self.allele_symbol_dto = None                          # Will be a single SymbolSlotAnnotationDTO.
        self.allele_full_name_dto = None                       # Will be a single FullNameSlotAnnotation.
        self.reference_curies = []                             # Will be a list of reference curies (directly or indirectly related).
        self.in_collection_name = None                         # Will be library.name.
        self.laboratory_of_origin_curie = None                 # N/A (WB).
        self.is_extinct = None                                 # Make True if extinction reported; make False is stock exists; leave as None otherwise.
        self.is_extrachromosomal = None                        # N/A (WB).
        self.is_integrated = None                              # N/A (WB).
        self.transgene_chromosome_location_curie = None        # ToDo - get chr via FBtp from FBti floc, derived_chromosome_location featureprop, or dock site.
        self.note_dtos = []                                    # ToDo - Waiting on "Allele Note Type" CV. Get from featureprop.
        self.allele_mutation_type_dtos = []                    # Will be list of slot annotations.
        self.allele_inheritance_mode_dtos = []                 # Will be list of slot annotations. TEMPORARY: Suppress phenotype_curie_term.
        self.allele_germline_transmission_status_dto = None    # N/A (MGI).
        self.allele_functional_impact_dtos = []                # ToDo - Waiting on "Functional Impact" CV. Get feature_cvterm, child of "allele class" term.
        self.allele_database_status_dto = None                 # Use term from "allele database status" CV.
        self.allele_secondary_id_dtos = []                     # Only 2o FlyBase IDs
        self.allele_nomenclature_event_dtos = []               # N/A.
        self.allele_synonym_dtos = []                          # Will be list of NameSlotAnnotationDTO objects.
        # Future ToDo:
        # Possibly relevant
        # cvterms from "FlyBase miscellaneous CV" in the "allele class" branch.
        # cvterms from "FlyBase miscellaneous CV" having "tool_uses" feature_cvtermprop type.
        # cvterms from "SO" having "transgenic_product_class" feature_cvtermprop type.
        # Notes associated with the object.
        self.for_alliance_export = True                        # Change to False if object should be excluded from export.
        self.internal_reasons = []                             # Reasons for marking an object as internal in the export file.
        self.export_warnings = []                              # Reasons for suppressing an object from the export file.

    def __str__(self):
        """Succinct text string describing the AllianceAllele object."""
        desc = '{} ({})'.format(self.feature.name, self.feature.uniquename)
        return desc


class AlleleHandler(object):
    """This object gets alleles and related info, synthesizes/filters the data, then exports it as LinkML JSON."""
    def __init__(self):
        """Create the AlleleHandler object."""
        self.allele_dict = {}         # An FBalID-keyed dict of AllianceAllele objects.
        self.all_pubs_dict = {}       # A pub_id-keyed dict of pub curies (PMID or FBrf).
        self.all_synonyms_dict = {}   # A synonym_id-keyed dict of Synonym objects.
        self.drosophilid_list = []    # A list of organism_ids for "Drosophilid" species in chado.
        self.cvterm_dict = {}         # Will be cvterm_id-keyed dicts: {'name': 'cvterm.name', 'curie': db.name:dbx.accession}
        self.total_feat_cnt = 0       # Count of all alleles found in starting query.
        self.export_feat_cnt = 0      # Count of all alleles exported to file.
        self.internal_feat_cnt = 0    # Count of all alleles marked as internal=True in export file.
    # Key chado CV term sets.
    allele_class_terms = []           # Will be cvterm_ids for child terms of "allele_class" (FBcv:0000286).
    allele_mutant_type_terms = []     # Will be cvterm_ids: child of chromosome_structure_variation or sequence_alteration.

    # Generic objects with which to build Alliance DTOs.
    generic_audited_object = {
        'internal': False,
        'obsolete': False,
        'created_by_curie': 'FB:FB_curator',
        'updated_by_curie': 'FB:FB_curator'
    }
    generic_data_provider_dto = generic_audited_object.copy()
    generic_data_provider_dto['source_organization_abbreviation'] = 'FB'
    generic_cross_reference_dto = {'prefix': 'FB', 'page_area': 'allele', 'internal': False}
    # Regexes.
    gene_regex = r'^FBgn[0-9]{7}$'
    allele_regex = r'^FBal[0-9]{7}$'
    insertion_regex = r'^FBti[0-9]{7}$'
    aberration_regex = r'^FB(ab|ba)[0-9]{7}$'
    construct_regex = r'^FBtp[0-9]{7}$'
    seqfeat_regex = r'^FBsf[0-9]{10}$'
    feature_regex = r'^FB(tp|ti)[0-9]{7}$'
    lib_regex = r'^FBlc[0-9]{7}$'
    pub_regex = r'^(FBrf[0-9]{7}|unattributed)$'
    # Sample set.
    test_alleles = []
    # Export fields.
    required_fields = [
        'allele_symbol_dto',
        'mod_entity_id',
        'data_provider_dto',
        'internal',
        'taxon_curie',
    ]
    output_fields = [
        'allele_database_status_dto',
        'allele_full_name_dto',
        # 'allele_functional_impact_dtos',
        'allele_inheritance_mode_dtos',
        'allele_mutation_type_dtos',
        # 'note_dtos',
        'allele_secondary_id_dtos',
        'allele_symbol_dto',
        'allele_synonym_dtos',
        'created_by_curie',
        'cross_reference_dtos',
        'mod_entity_id',
        'mod_internal_id',
        'data_provider_dto',
        'date_created',
        'date_updated',
        'in_collection_name',
        'internal',
        'is_extinct',
        'obsolete',
        # 'reference_curies',    # TEMPORARY until AGR load is faster
        'taxon_curie',
        # 'transgene_chromosome_location_curie',
        'updated_by_curie',
    ]
    fb_agr_db_dict = {
        'FlyBase': 'FB'
    }

    def get_cvterm_info(self, session):
        """Get all CV terms in chado and their curies."""
        results = session.query(Cvterm, Db, Dbxref).\
            select_from(Cvterm).\
            join(Dbxref, (Dbxref.dbxref_id == Cvterm.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            distinct()
        counter = 0
        for result in results:
            cvterm = {
                'name': result.Cvterm.name,
                'curie': f'{result.Db.name}:{result.Dbxref.accession}'
            }
            self.cvterm_dict[result.Cvterm.cvterm_id] = cvterm
            counter += 1
        log.info(f'Found {counter} CV terms in chado.')
        return

    def __get_child_cvterms(self, session, starting_cvterm_name, starting_cvterm_cv_name):
        """Get all cvterm_ids for some branch of an ontology defined by the parent term.

        Args:
            arg1 (self): (AlleleHandler) The object that runs the method.
            arg2 (session): (SQLAlchemy session)
            arg3 (starting_cvterm_name): (str) The name of the parent CV term.
            arg4 (starting_cvterm_cv_name): (str) The name of the parent CV term's CV.

        Returns:
            List of cvterm.cvterm_ids, including that of the starting parent CV term.

        Raises:
            Raise NoResultFound if starting CV term cannot be found in chado.

        """
        log.info(f'Get all child terms of "{starting_cvterm_name}" ("{starting_cvterm_cv_name}")')
        # 1. Get the parent CV term.
        filters = (
            Cvterm.name == starting_cvterm_name,
            Cv.name == starting_cvterm_cv_name)
        try:
            starting_cvterm = session.query(Cvterm).\
                join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
                filter(*filters).\
                one()
            log.debug(f'Found chado CV term found for "{starting_cvterm_name}" ("{starting_cvterm_cv_name}")')
        except NoResultFound:
            log.error(f'No chado CV term found for "{starting_cvterm_name}" ("{starting_cvterm_cv_name}")')
            raise NoResultFound
        # 2. Define the start of the recursive query for child terms of the starting term.
        # Recursive query built up using the sorta helpful instructions here:
        # https://docs.sqlalchemy.org/en/13/orm/query.html
        # https://sanjayasubedi.com.np/python/sqlalchemy/recursive-query-in-postgresql-with-sqlalchemy/
        cvterm_subject1, cvterm_subject2, cvterm_object = aliased(Cvterm), aliased(Cvterm), aliased(Cvterm)
        cvterm_relationship1, cvterm_relationship2 = aliased(CvtermRelationship), aliased(CvtermRelationship)
        cv1, cv2 = aliased(Cv), aliased(Cv)
        # 2a. Define the starting point of the query (output of which to be used recursively to get more results).
        # Basically, looking for all child terms (subject) of the initial CV term (object) in cvterm_relationship.
        # The "recursive_query_start" below is not actually results, but a 'sqlalchemy.sql.base.ImmutableColumnCollection' class.
        filters = (cvterm_object.cvterm_id == starting_cvterm.cvterm_id,
                   cv1.name == starting_cvterm_cv_name)
        recursive_query_start = session.query(cvterm_subject1).\
            join(cvterm_relationship1, (cvterm_relationship1.subject_id == cvterm_subject1.cvterm_id)).\
            join(cvterm_object, (cvterm_object.cvterm_id == cvterm_relationship1.object_id)).\
            join(cv1, (cv1.cv_id == cvterm_subject1.cv_id)).\
            filter(*filters).\
            cte(recursive=True)    # Important bit here.
        # 3. Define the second part of the recursive query.
        # Essentially, we want child terms of child terms in the starting query, recursively.
        # So, subject_ids in "recursive_query_start" results will be the objects in our query of cvterm_relationship.
        # Importantly, the columns of the initial "recursive_query_start" query are accessed by the ".c" attribute.
        filters2 = (cv2.name == starting_cvterm_cv_name,)
        recursive_query_repeat = session.query(cvterm_subject2).\
            join(cvterm_relationship2, (cvterm_relationship2.subject_id == cvterm_subject2.cvterm_id)).\
            join(recursive_query_start, (recursive_query_start.c.cvterm_id == cvterm_relationship2.object_id)).\
            join(cv2, (cv2.cv_id == cvterm_subject2.cv_id)).\
            filter(*filters2)
        # 4. Define a query that takes the union of the starting and recursive queries.
        recursive_query_total = recursive_query_start.union(recursive_query_repeat)
        # 5. And finally, get the result for the start and recursive query parts). Piece of cake ;)
        recursive_query_total_results = session.query(recursive_query_total)
        # Build the list from the results.
        cvterm_id_list = [i[0] for i in recursive_query_total_results]
        log.info(f'Found {len(cvterm_id_list)} terms under "{starting_cvterm_name}" ("{starting_cvterm_cv_name}")')
        cvterm_id_list.append(starting_cvterm.cvterm_id)
        return cvterm_id_list

    def get_key_cvterm_sets(self, session):
        """Get key CV term sets from chado."""
        log.info('Get key CV term sets from chado.')
        self.allele_class_terms.extend(self.__get_child_cvterms(session, 'allele class', 'FlyBase miscellaneous CV'))
        self.allele_mutant_type_terms.extend(self.__get_child_cvterms(session, 'chromosome_structure_variation', 'SO'))
        self.allele_mutant_type_terms.extend(self.__get_child_cvterms(session, 'sequence_alteration', 'SO'))
        return

    def get_all_references(self, session):
        """Get all references."""
        log.info('Get all references.')
        # First get all current pubs having an FBrf uniquename.
        filters = (
            Pub.uniquename.op('~')(self.pub_regex),
            Pub.is_obsolete.is_(False)
        )
        results = session.query(Pub).\
            filter(*filters).\
            distinct()
        pub_counter = 0
        for pub in results:
            self.all_pubs_dict[pub.pub_id] = f'FB:{pub.uniquename}'
            pub_counter += 1
        # Next find PMIDs if available and replace the curie in the all_pubs_dict.
        filters = (
            Pub.uniquename.op('~')(self.pub_regex),
            Pub.is_obsolete.is_(False),
            Db.name == 'pubmed',
            PubDbxref.is_current.is_(True)
        )
        pmid_xrefs = session.query(Pub, Dbxref).\
            join(PubDbxref, (PubDbxref.pub_id == Pub.pub_id)).\
            join(Dbxref, (Dbxref.dbxref_id == PubDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        pmid_counter = 0
        for xref in pmid_xrefs:
            self.all_pubs_dict[xref.Pub.pub_id] = f'PMID:{xref.Dbxref.accession}'
            pmid_counter += 1
        log.info(f'Found {pmid_counter} PMID IDs for {pub_counter} current FB publications.')
        return

    def get_alleles(self, session):
        """Get all alleles."""
        log.info('Querying chado for alleles.')

        # First get all allele features from chado.
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'allele'
        )
        allele_results = session.query(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        self.total_feat_cnt = 0
        for result in allele_results:
            self.total_feat_cnt += 1
            self.allele_dict[result.uniquename] = AllianceAllele(result)
            self.allele_dict[result.uniquename].organism_abbr = result.organism.abbreviation
        log.info('Found {} alleles.'.format(self.total_feat_cnt))

    def get_allele_gene(self, session):
        """For current alleles, get the FBgn ID of allele's current gene."""
        gene = aliased(Feature, name='gene')
        allele = aliased(Feature, name='allele')
        filters = (
            gene.is_obsolete.is_(False),
            gene.is_analysis.is_(False),
            gene.uniquename.op('~')(self.gene_regex),
            allele.is_obsolete.is_(False),
            allele.is_analysis.is_(False),
            allele.uniquename.op('~')(self.allele_regex),
            Cvterm.name == 'alleleof'
        )
        allele_gene_results = session.query(allele, gene).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_gene_results:
            counter += 1
            self.allele_dict[result.allele.uniquename].parent_gene = result.gene.uniquename
        log.info('Found {} allele-gene relationships.'.format(counter))
        return

    def flag_alleles_of_internal_genes(self, session):
        """Flag alleles belonging to internal-type Dmel genes like origin_of_replication."""
        log.info('Flagging alleles of internal-type Dmel genes.')
        # First build list of internal Dmel genes.
        internal_gene_types = [
            'engineered_fusion_gene',
            'engineered_region',
            'gene_group',
            'gene_with_polycistronic_transcript',
            'insulator',
            'mitochondrial_sequence',
            'origin_of_replication',
            'region',
            'regulatory_region',
            'repeat_region',
            'satellite_DNA',
            'transposable_element_gene'
        ]
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='promoted_gene_type')
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            Organism.abbreviation == 'Dmel',
            prop_type.name == 'promoted_gene_type'
        )
        gene_type_results = session.query(Featureprop).\
            join(Feature, (Feature.feature_id == Featureprop.feature_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        gene_counter = 0
        internal_dmel_genes = []
        for result in gene_type_results:
            gene_type_name = result.value[11:-1]
            if gene_type_name in internal_gene_types:
                gene_counter += 1
                internal_dmel_genes.append(result.feature.uniquename)
        log.info('Found {} internal type Dmel genes.'.format(gene_counter))
        # Now check alleles against this internal Dmel gene list.
        allele_counter = 0
        for allele in self.allele_dict.values():
            if allele.parent_gene in internal_dmel_genes:
                allele_counter += 1
                allele.allele_of_internal_gene = True
        log.info('Found {} alleles of internal type Dmel genes.'.format(allele_counter))
        return

    def flag_in_vitro_alleles(self, session):
        """Flag alleles associated with "in vitro" type CV term."""
        log.info('Flag in vitro alleles.')
        cvterm_name_regex = '^in vitro construct'
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Cvterm.name.op('~')(cvterm_name_regex)
        )
        ivt_alleles = session.query(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for allele in ivt_alleles:
            counter += 1
            self.allele_dict[allele.uniquename].in_vitro = True
        log.info(f'Flagged {counter} alleles as "in vitro"')
        return

    def get_allele_constructs(self, session):
        """Find FBtp constructs associated with the allele."""
        log.info('Find constructs related to alleles.')
        allele = aliased(Feature, name='allele')
        construct = aliased(Feature, name='construct')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.allele_regex),
            construct.is_obsolete.is_(False),
            construct.uniquename.op('~')(self.construct_regex),
            Cvterm.name == 'derived_tp_assoc_alleles'
        )
        construct_results = session.query(allele, construct).\
            select_from(construct).\
            join(FeatureRelationship, (FeatureRelationship.object_id == construct.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.subject_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in construct_results:
            counter += 1
            self.allele_dict[result.allele.uniquename].constructs.append(result.construct)
        log.info('Found {} constructs related to alleles.'.format(counter))
        return

    def get_allele_insertions(self, session):
        """Find FBti insertions associated with the allele."""
        log.info('Find insertions related to alleles.')
        allele = aliased(Feature, name='allele')
        insertion = aliased(Feature, name='insertion')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.allele_regex),
            insertion.is_obsolete.is_(False),
            insertion.is_analysis.is_(False),
            insertion.uniquename.op('~')(self.insertion_regex),
            Cvterm.name == 'associated_with'
        )
        insertion_results = session.query(Organism, allele, insertion).\
            join(FeatureRelationship, (FeatureRelationship.object_id == insertion.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.subject_id)).\
            join(Organism, (Organism.organism_id == insertion.organism_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        dmel_counter = 0
        non_dmel_counter = 0
        for result in insertion_results:
            if result.Organism.abbreviation == 'Dmel':
                self.allele_dict[result.allele.uniquename].dmel_insertions.append(result.insertion)
                dmel_counter += 1
            else:
                self.allele_dict[result.allele.uniquename].non_dmel_insertions.append(result.insertion)
                non_dmel_counter += 1
        log.info(f'Found {dmel_counter} Dmel and {non_dmel_counter} non-Dmel insertions related to alleles.')
        return

    def get_drosophilid_organisms(self, session):
        """Find organisms for Drosophilid species."""
        filters = (
            Cvterm.name == 'taxgroup',
            Organismprop.value == 'drosophilid'
        )
        drosophilid_results = session.query(Organism).\
            join(Organismprop, (Organismprop.organism_id == Organism.organism_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Organismprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for drosophilid in drosophilid_results:
            self.drosophilid_list.append(drosophilid.organism_id)
            counter += 1
        log.info(f'Found {counter} Drosophilid organisms in chado.')
        return

    def adjust_allele_org(self, session):
        """Find classical alleles in non-Dmel species and adjust organism accordingly."""
        log.info('Adjust organism for classical non-Dmel alleles.')
        counter = 0
        for allele in self.allele_dict.values():
            # Skip alleles that unambiguously occur in Dmel.
            if allele.organism_abbr == 'Dmel':
                continue
            if allele.dmel_insertions:
                continue
            # How to classify alleles like FBal0048225? They look like non-Dmel TE transgenic alleles carried in Dmel.
            # FBal0048225 has mutagen=natural population, but still seems to be transgenic.
            if allele.constructs:
                continue
            # Find clear evidence that allele is non-Dmel classical allele.
            is_non_dmel_classical = False
            if allele.non_dmel_insertions:
                is_non_dmel_classical = True
            elif allele.feature.organism_id in self.drosophilid_list and allele.in_vitro is False:
                is_non_dmel_classical = True
            if is_non_dmel_classical is True:
                allele.adj_organism_abbr = allele.organism_abbr
                log.debug(f'Non-Dmel allele: id={allele.mod_entity_id}, name={allele.feature.name}, org_abbr={allele.organism_abbr}')
                counter += 1
        log.info('Adjusted organism to be "non-Dmel" for {} alleles.'.format(counter))
        return

    def get_allele_taxons(self, session):
        """Get taxon IDs for alleles. Depends on all organisms for features having an abbreviation."""
        log.info('Getting allele taxon IDs.')
        # First make a dict of organism abbr to NCBI taxon IDs.
        filters = (
            OrganismDbxref.is_current.is_(True),
            Db.name == 'NCBITaxon'
        )
        organism_dbxref_results = session.query(Organism, Dbxref).\
            join(OrganismDbxref, (OrganismDbxref.organism_id == Organism.organism_id)).\
            join(Dbxref, (Dbxref.dbxref_id == OrganismDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        organism_taxon_dict = {}
        for result in organism_dbxref_results:
            organism_taxon_dict[result.Organism.abbreviation] = result.Dbxref.accession
        # Now fill in the info for alleles.
        for allele in self.allele_dict.values():
            try:
                allele.taxon_curie = f'NCBITaxon:{organism_taxon_dict[allele.adj_organism_abbr]}'
            except KeyError:
                log.warning(f'Use "unidentified" NCBITaxon ID for {allele}')
                allele.taxon_curie = 'NCBITaxon:32644'
        return

    def get_synonyms(self, session):
        """Get current and non-current symbols and full names for alleles."""
        log.info('Get current and non-current symbols and full names for alleles.')
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'allele'
        )
        results = session.query(Feature, FeatureSynonym, Synonym).\
            join(FeatureSynonym, (FeatureSynonym.synonym_id == Synonym.synonym_id)).\
            join(Feature, (Feature.feature_id == FeatureSynonym.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            # Skip any references to non-current pubs.
            if result.FeatureSynonym.pub_id not in self.all_pubs_dict.keys():
                continue
            # First, build the all_synonyms_dict.
            self.all_synonyms_dict[result.Synonym.synonym_id] = result.Synonym
            # Second, collect FeatureSynonyms for each allele.
            self.allele_dict[result.Feature.uniquename].feature_synonyms.append(result.FeatureSynonym)
            # Third, capture pub_ids.
            self.allele_dict[result.Feature.uniquename].fb_references.append(result.FeatureSynonym.pub_id)

            # Finally, catch current symbol and fullname strings.
            if result.FeatureSynonym.is_current is True and result.Synonym.type.name == 'symbol':
                self.allele_dict[result.Feature.uniquename].curr_symbol_name = sub_sup_sgml_to_html(result.Synonym.synonym_sgml)
            elif result.FeatureSynonym.is_current is True and result.Synonym.type.name == 'fullname':
                self.allele_dict[result.Feature.uniquename].curr_fullname = sub_sup_sgml_to_html(result.Synonym.synonym_sgml)
            counter += 1
        log.info(f'Found {counter} feature_synonyms (current pubs) for alleles.')
        return

    def get_allele_timestamps(self, session):
        """Get timestamps for alleles."""
        log.info('Getting allele timestamps.')
        # Start with timestamps in feature itself.
        for allele in self.allele_dict.values():
            allele.timestamps.append(allele.feature.timeaccessioned)
            allele.timestamps.append(allele.feature.timelastmodified)
        # We may want to get timestamps for other allele attibutes as well (for "date_updated").
        # Placeholder: current query of "feature" table in audit_chado.
        # Querying audit_chado with SQLAlchemy has not been user-friendly, so go with standard SQL.
        ########################################################################
        # audit_chado_query = """
        #     SELECT DISTINCT f.uniquename,
        #                     ac.transaction_timestamp
        #     FROM feature
        #     JOIN audit_chado ac ON (ac.record_pkey = f.feature_id AND ac.audited_table = 'feature');
        # """
        # audit_results = session.execute(audit_chado_query).fetchall()
        # log.info('Got {} audit_chado results. Will parse them out now.'.format(len(audit_results)))
        # UNIQUE_KEY = 0
        # TIMESTAMP = 1
        # for row in audit_results:
        #     try:
        #         # log.debug('For unique_key={}, have timestamp={}'.format(row[UNIQUE_KEY], row[TIMESTAMP]))
        #         self.allele_dict[row[UNIQUE_KEY]].timestamps.append(row[TIMESTAMP])
        #     except KeyError:
        #         # log.debug('Could not put this in anno dict: {}'.format(row))
        #         pass
        ########################################################################
        return

    def get_allele_dbxrefs(self, session):
        """Get all dbxrefs for alleles."""
        log.info('Getting allele dbxrefs.')
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'allele',
            Db.name.in_((self.fb_agr_db_dict.keys()))
        )
        allele_dbxref_results = session.query(Feature, FeatureDbxref, Dbxref, Db).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_dbxref_results:
            # Skip current FlyBase accessions (i.e., equiavlent to feature.uniquename)
            #     because not all features have them (e.g., FBal0137236).
            # Better to use current uniquename as source for current FlyBase accession.
            counter += 1
            if result.FeatureDbxref.is_current is True and result.Db.name == 'FlyBase':
                pass
            elif result.FeatureDbxref.is_current is False and result.Db.name == 'FlyBase':
                self.allele_dict[result.Feature.uniquename].dbxrefs.append(result)
                self.allele_dict[result.Feature.uniquename].alt_fb_ids.append(result.Dbxref)
            else:
                self.allele_dict[result.Feature.uniquename].dbxrefs.append(result)
        log.info('Found {} allele crossreferences.'.format(counter))
        return

    def get_allele_references(self, session):
        """Get references for alleles."""
        log.info('Get allele references.')
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Pub.uniquename.op('~')(self.pub_regex),
            Pub.is_obsolete.is_(False)
        )
        allele_pubs = session.query(Feature, Pub).\
            join(FeaturePub, (FeaturePub.feature_id == Feature.feature_id)).\
            join(Pub, (Pub.pub_id == FeaturePub.pub_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_pubs:
            self.allele_dict[result.Feature.uniquename].fb_references.append(result.Pub.pub_id)
            counter += 1
        log.info(f'Found {counter} allele-pub relationships.')
        return

    def get_allele_featureprops(self, session):
        """Get all allele featureprops."""
        log.info('Get allele featureprops.')
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Cvterm.is_obsolete == 0
        )
        allele_fprops = session.query(Feature, Cvterm, Featureprop).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_fprops:
            counter += 1
            try:
                self.allele_dict[result.Feature.uniquename].\
                    featureprops[result.Cvterm.name].append(result.Featureprop)
            except KeyError:
                self.allele_dict[result.Feature.uniquename].\
                    featureprops[result.Cvterm.name] = [result.Featureprop]
        log.info(f'Found {counter} allele featureprops.')
        return

    def get_args(self, session):
        """Get ARGs related to alleles."""
        log.info('Get allele ARGs.')
        arg_types = [
            'MNV',
            'complex_substitution',
            'deletion',
            'delins',
            'insertion',
            'point_mutation',
            'sequence_alteration',
            'sequence_variant',
            'rescue_region'
        ]
        allele = aliased(Feature, name='allele')
        arg = aliased(Feature, name='arg')
        argtype = aliased(Cvterm, name='argtype')
        reltype = aliased(Cvterm, name='reltype')
        filters = (
            allele.uniquename.op('~')(self.allele_regex),
            arg.is_obsolete.is_(False),
            argtype.name.in_((arg_types)),
            reltype.name == 'partof'
        )
        arg_results = session.query(arg, allele).\
            join(argtype, (argtype.cvterm_id == arg.type_id)).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == arg.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.object_id)).\
            join(reltype, (reltype.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in arg_results:
            self.allele_dict[result.allele.uniquename].args.append(result.arg)
            counter += 1
        log.info(f'Found {counter} ARG-allele relationships.')
        return

    def get_phenotypes(self, session):
        """Get phenotypes related to alleles."""
        log.info('Get phenotypes related to alleles.')
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Genotype.is_obsolete.is_(False),
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(self.pub_regex)
        )
        results = session.query(Feature, Genotype, Phenotype, Cvterm, Pub).\
            join(FeatureGenotype, (FeatureGenotype.feature_id == Feature.feature_id)).\
            join(Genotype, (Genotype.genotype_id == FeatureGenotype.genotype_id)).\
            join(Phenstatement, (Phenstatement.genotype_id == Genotype.genotype_id)).\
            join(Phenotype, (Phenotype.phenotype_id == Phenstatement.phenotype_id)).\
            join(PhenotypeCvterm, (PhenotypeCvterm.phenotype_id == Phenotype.phenotype_id)).\
            join(Cvterm, (Cvterm.cvterm_id == PhenotypeCvterm.cvterm_id)).\
            join(Pub, (Pub.pub_id == Phenstatement.pub_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.allele_dict[result.Feature.uniquename].phenstatements.append(result)
            counter += 1
        log.info(f'Found {counter} allele-phenotype associations.')
        return

    def get_direct_collections(self, session):
        """Find library collections directly related to alleles."""
        log.info('Get directly-related allele collections.')
        counter = 0
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        # First, look for direct FBal-FBlc associations.
        filters = (
            Feature.uniquename.op('~')(self.allele_regex),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection'
        )
        libraries = session.query(Feature, Library).\
            join(LibraryFeature, (LibraryFeature.feature_id == Feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in libraries:
            self.allele_dict[result.Feature.uniquename].direct_libraries.append(result.Library)
            counter += 1
        log.info(f'Found {counter} direct allele-library associations.')
        return

    def get_indirect_collections(self, session):
        """Find library collections indirectly related to alleles via insertion or construct."""
        log.info('Get indirectly-related allele collections (via insertion or construct).')
        allele = aliased(Feature, name='allele')
        feature = aliased(Feature, name='feature')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        featreltype = aliased(Cvterm, name='featreltype')
        filters = (
            allele.uniquename.op('~')(self.allele_regex),
            feature.uniquename.op('~')(self.feature_regex),
            feature.is_obsolete.is_(False),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection',
            featreltype.name == 'associated_with'
        )
        indirect_libraries = session.query(allele, feature, Library).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(featreltype, (featreltype.cvterm_id == FeatureRelationship.type_id)).\
            join(feature, (feature.feature_id == FeatureRelationship.object_id)).\
            join(LibraryFeature, (LibraryFeature.feature_id == feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        fbti_counter = 0
        fbtp_counter = 0
        for result in indirect_libraries:
            if result.feature.uniquename.startswith('FBti'):
                self.allele_dict[result.allele.uniquename].ins_libraries.append(result.Library)
                fbti_counter += 1
            elif result.feature.uniquename.startswith('FBtp'):
                self.allele_dict[result.allele.uniquename].cons_libraries.append(result.Library)
                fbtp_counter += 1
        log.info(f'Found {fbti_counter} insertion-mediated allele-library associations.')
        log.info(f'Found {fbtp_counter} construct-mediated allele-library associations.')
        return

    def get_sf_collections(self, session):
        """Find library collections indirectly related to alleles via sequence feature."""
        log.info('Get indirectly-related allele collections (via equence feature).')
        allele = aliased(Feature, name='allele')
        construct = aliased(Feature, name='construct')
        seqfeat = aliased(Feature, name='seqfeat')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        featreltype = aliased(Cvterm, name='featreltype')
        allele_construct = aliased(FeatureRelationship, name='allele_construct')
        seqfeat_construct = aliased(FeatureRelationship, name='seqfeat_construct')
        filters = (
            allele.uniquename.op('~')(self.allele_regex),
            construct.uniquename.op('~')(self.construct_regex),
            seqfeat.uniquename.op('~')(self.seqfeat_regex),
            construct.is_obsolete.is_(False),
            seqfeat.is_obsolete.is_(False),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.lib_regex),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection',
            featreltype.name == 'associated_with'
        )
        sf_libraries = session.query(allele, Library).\
            join(allele_construct, (allele_construct.subject_id == allele.feature_id)).\
            join(construct, (construct.feature_id == allele_construct.object_id)).\
            join(featreltype, (featreltype.cvterm_id == allele_construct.type_id)).\
            join(seqfeat_construct, (seqfeat_construct.object_id == construct.feature_id)).\
            join(seqfeat, (seqfeat.feature_id == seqfeat_construct.subject_id)).\
            join(LibraryFeature, (LibraryFeature.feature_id == seqfeat.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in sf_libraries:
            self.allele_dict[result.allele.uniquename].sf_libraries.append(result.Library)
            counter += 1
        log.info(f'Found {counter} sequence feature-mediated allele-library associations.')
        return

    # Synthesis of initial db info.
    def synthesize_timestamps(self, allele):
        """Process timestamps for an allele."""
        if allele.timestamps:
            allele.date_created = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(allele.timestamps)))
            allele.date_updated = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(allele.timestamps)))
        return

    def synthesize_secondary_ids(self, allele):
        """Process 2o IDs."""
        unique_fb_id_list = list(set(allele.alt_fb_ids))
        for fb_id in unique_fb_id_list:
            secondary_id_dict = self.generic_audited_object.copy()
            secondary_id_dict['secondary_id'] = f'FB:{fb_id.accession}'
            allele.allele_secondary_id_dtos.append(secondary_id_dict)
        return

    def synthesize_xrefs(self, allele):
        """Process xrefs."""
        # Start by adding allele uniquename as an xref.
        xref_dict = self.generic_audited_object.copy()
        xref_dict['referenced_curie'] = f'FB:{allele.feature.uniquename}'
        xref_dict['display_name'] = f'FB:{allele.feature.uniquename}'
        xref_dict['prefix'] = 'FB'
        xref_dict['page_area'] = 'allele'
        allele.cross_reference_dtos.append(xref_dict)
        # Add other xrefs: code below assumes xrefs are all 'FB' at the moment.
        for result in allele.dbxrefs:
            if result.Db.name not in self.fb_agr_db_dict.keys():
                continue
            xref_dict = self.generic_audited_object.copy()
            xref_dict['referenced_curie'] = f'{self.fb_agr_db_dict[result.Db.name]}:{result.Dbxref.accession}'
            xref_dict['display_name'] = f'{self.fb_agr_db_dict[result.Db.name]}:{result.Dbxref.accession}'
            xref_dict['prefix'] = self.fb_agr_db_dict[result.Db.name]
            xref_dict['page_area'] = 'allele'
            if result.FeatureDbxref.is_current is False:
                xref_dict['internal'] = True
            allele.cross_reference_dtos.append(xref_dict)
        return

    def synthesize_references(self, allele):
        """Process pubs for allele."""
        allele.fb_references = list(set(allele.fb_references))
        allele.reference_curies = [self.all_pubs_dict[i] for i in allele.fb_references if self.all_pubs_dict[i] != 'FB:unattributed']
        return

    def synthesize_insertions(self, allele):
        """Process insertions."""
        for insertion in allele.dmel_insertions:
            xref_dict = self.generic_audited_object.copy()
            xref_dict['referenced_curie'] = f'FB:{insertion.uniquename}'
            xref_dict['display_name'] = insertion.name
            xref_dict['prefix'] = 'FB'
            xref_dict['page_area'] = 'allele'
            allele.cross_reference_dtos.append(xref_dict)
        return

    def synthesize_extinction(self, allele):
        """Determine if allele definitively exists or is extinct."""
        has_stocks = False
        reported_extinct = False
        # First find evidence of extinction.
        try:
            for fprop in allele.featureprops['availability']:
                if fprop.value == 'Stated to be lost.':
                    reported_extinct = True
        except KeyError:
            pass
        # Second find evidence for existence.
        for fprop_type in allele.featureprops.keys():
            if fprop_type.startswith('derived_stock_'):
                has_stocks = True
        # Synthesize these two pieces of info.
        if has_stocks is True:
            allele.is_extinct = False
        elif reported_extinct is True:
            allele.is_extinct = True
        return

    def synthesize_inheritance_mode(self, allele):
        """Determine inheritance mode for the allele."""
        # We convert the FlyBase terms to Alliance terms as follows.
        inheritance_mode_terms = {
            'recessive': 'recessive',
            'dominant': 'dominant',
            'semidominant': 'semi-dominant',
            'codominant': 'codominant'
        }
        # Gather phenotype data relevant to inheritance mode as a dict.
        # Keys will be tuple of (inheritance_mode, phenotype_curie, phenotype_statement)
        # Value for each key will be a list of pub_ids in support of each key.
        # Map the phenotype_cvterm.cvterm_id to "inheritance_mode_name".
        # Map the phenotype.cvalue_id to phenotype_term_curie.
        # Map the phenotype.description to phenotype_statement.
        inheritance_data = {}
        for phenstmt in allele.phenstatements:
            # First, filter out phenotype entries that are not relevant to inheritance mode.
            cvterm = phenstmt.Cvterm.name
            if cvterm not in inheritance_mode_terms.keys():
                continue
            # Second, filter out complex genotypes. Start by assuming that genotype is complex.
            single_allele_genotype = False
            # Weed out multi-locus genotypes.
            if ' ' in phenstmt.Genotype.uniquename:
                single_allele_genotype = False
            elif '_' in phenstmt.Genotype.description:
                single_allele_genotype = False
            # Assess single locus having 2 features: could be homo- or heterozygous.
            elif '|' in phenstmt.Genotype.description:
                features = phenstmt.Genotype.description.split('|')
                if features[0] == features[1]:
                    single_allele_genotype = True
                else:
                    for feature in features:
                        if feature == '+':
                            single_allele_genotype = True
                        elif feature.endswith('[+]'):
                            single_allele_genotype = True
            # If no spaces, underscores or pipes in description, it's a single-feature genotype.
            else:
                single_allele_genotype = True
            if single_allele_genotype is False:
                continue
            # For relevant phenotypes of simple genotypes, we parse out the data.
            inheritance_mode_name = inheritance_mode_terms[cvterm]
            phenotype_term_curie = f'{phenstmt.Phenotype.cvalue.dbxref.db.name}:{phenstmt.Phenotype.cvalue.dbxref.accession}'
            phenotype_statement = phenstmt.Phenotype.uniquename
            pheno_key = (inheritance_mode_name, phenotype_term_curie, phenotype_statement)
            pub_id = phenstmt.Pub.pub_id
            pub_curie = self.all_pubs_dict[pub_id]
            try:
                inheritance_data[pheno_key].append(pub_curie)
            except KeyError:
                inheritance_data[pheno_key] = [pub_curie]
        # Convert data into Alliance slot annotations.
        INHERITANCE_MODE_NAME = 0
        # PHENOTYPE_CURIE_NAME = 1    # TEMPORARY: Suppress until AGR has FBcv
        PHENOTYPE_STATEMENT = 2
        for pheno_key, pub_curie_list in inheritance_data.items():
            for pub_curie in pub_curie_list:
                if pub_curie == 'FB:unattributed':
                    pub_curie_list.remove('FB:unattributed')
            allele_inheritance_mode_slot_annotation_dto = self.generic_audited_object.copy()
            allele_inheritance_mode_slot_annotation_dto['inheritance_mode_name'] = pheno_key[INHERITANCE_MODE_NAME]
            # allele_inheritance_mode_slot_annotation_dto['phenotype_term_curie'] = pheno_key[PHENOTYPE_CURIE_NAME]    # TEMPORARY: Suppress until AGR has FBcv
            allele_inheritance_mode_slot_annotation_dto['phenotype_statement'] = pheno_key[PHENOTYPE_STATEMENT]
            # allele_inheritance_mode_slot_annotation_dto['evidence_curies'] = list(set(pub_curie_list))    # TEMPORARY: Suppress until AGR loads are faster
            allele.allele_inheritance_mode_dtos.append(allele_inheritance_mode_slot_annotation_dto)
        return

    def synthesize_collections(self, allele):
        """Get names of collections to which allele belong (directly or indirectly)."""
        collection_names = None
        if allele.direct_libraries:
            collection_names = allele.direct_libraries
        elif allele.ins_libraries:
            collection_names = allele.ins_libraries
        elif allele.cons_libraries:
            collection_names = allele.cons_libraries
        elif allele.sf_libraries:
            collection_names = allele.sf_libraries
        if collection_names:
            collection_names = list(set(collection_names))
            allele.in_collection_name = collection_names[0].name
            if len(collection_names) > 1:
                log.warning(f'\tFound {len(collection_names)} collection(s) for {allele.mod_entity_id}: {allele.in_collection_name}')
        return

    def synthesize_synonyms(self, feature):
        """Generate name/synonym DTOs for a feature that has a list of FeatureSynonym objects."""
        # Dict for converting FB to AGR synonym types.
        synonym_type_conversion = {
            'symbol': 'nomenclature_symbol',
            'fullname': 'full_name',
            'nickname': 'nomenclature_symbol',
            'synonym': 'nomenclature_symbol'
        }
        default_name_dto = self.generic_audited_object.copy()
        default_name_dto['name_type_name'] = 'unspecified'
        default_name_dto['format_text'] = 'unspecified'
        default_name_dto['display_text'] = 'unspecified'
        default_name_dto['synonym_scope_name'] = 'exact'
        # Create a dict of all distinct name/synonym_sgml combinations: for each, capture synonym type(s) an pub_ids.
        # Keys are (synonym.name, synonym.synonym_sgml) tuples.
        # Values are dicts too where keys are chado synonym types and values are lists of pub_ids.
        # Value dict also has an "internal" key that stores list of FeatureSynonym.is_internal values.
        feature_synonym_dict = {}
        for f_s in feature.feature_synonyms:
            synonym = self.all_synonyms_dict[f_s.synonym_id]
            distinct_synonym_name = (synonym.name, synonym.synonym_sgml)
            if distinct_synonym_name in feature_synonym_dict.keys():
                feature_synonym_dict[distinct_synonym_name]['internal'].append(f_s.is_internal)
                if synonym.type.name in feature_synonym_dict[distinct_synonym_name].keys():
                    feature_synonym_dict[distinct_synonym_name][synonym.type.name].append(f_s.pub_id)
                else:
                    feature_synonym_dict[distinct_synonym_name][synonym.type.name] = [f_s.pub_id]
            else:
                feature_synonym_dict[distinct_synonym_name] = {synonym.type.name: [f_s.pub_id], 'internal': [f_s.is_internal]}
        # Convert to AGR name DTO objects.
        name_dto_list = []
        FORMAT_TEXT = 0
        DISPLAY_TEXT = 1
        for syno_name, syno_attributes in feature_synonym_dict.items():
            # Determine internal status. False trumps True.
            if False in set(syno_attributes['internal']):
                syno_internal = False
            else:
                syno_internal = True
            # Collect all pubs.
            pub_id_list = []
            for syno_type, syno_type_pub_list in syno_attributes.items():
                # Skip over the 'internal' attribute, which is not actually a synonym type.
                if syno_type == 'internal':
                    continue
                pub_id_list.extend(syno_type_pub_list)
            pub_id_list = list(set(pub_id_list))
            # Pick the best synonym type.
            type_tally = {}
            for syno_type, syno_type_pub_list in syno_attributes.items():
                if syno_type == 'internal':
                    continue
                type_tally[len(set(syno_type_pub_list))] = syno_type
            name_type_to_use = synonym_type_conversion[type_tally[max(type_tally.keys())]]
            output_synonym_dto = self.generic_audited_object.copy()
            output_synonym_dto['name_type_name'] = name_type_to_use
            output_synonym_dto['format_text'] = sub_sup_sgml_to_html(syno_name[FORMAT_TEXT])
            output_synonym_dto['display_text'] = sub_sup_sgml_to_html(syno_name[DISPLAY_TEXT])
            output_synonym_dto['synonym_scope_name'] = 'exact'
            # output_synonym_dto['evidence_curies'] = [self.all_pubs_dict[i] for i in pub_id_list if self.all_pubs_dict[i] != 'FB:unattributed']    # TEMP
            output_synonym_dto['internal'] = syno_internal
            name_dto_list.append(output_synonym_dto)
        # Sift through name DTOs for symbol, fullname, systematic_name, etc.
        for name_dto in name_dto_list:
            if name_dto['display_text'] == feature.curr_symbol_name:
                if name_dto['name_type_name'] != 'nomenclature_symbol':
                    log.warning(f"{feature}: Found mistyped curr symbol: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
                    name_dto['name_type_name'] = 'nomenclature_symbol'
                feature.allele_symbol_dto = name_dto
            elif name_dto['display_text'] == feature.curr_fullname:
                if name_dto['name_type_name'] != 'full_name':
                    log.warning(f"{feature}: Found mistyped curr full_name: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
                    name_dto['name_type_name'] = 'full_name'
                feature.allele_full_name_dto = name_dto
            else:
                feature.allele_synonym_dtos.append(name_dto)
        # Symbol is required. If none, fill it in.
        if feature.allele_symbol_dto is None:
            placeholder_symbol_dto = default_name_dto.copy()
            placeholder_symbol_dto['name_type_name'] = 'nomenclature_symbol'
            placeholder_symbol_dto['format_text'] = feature.feature.name
            placeholder_symbol_dto['display_text'] = feature.feature.name
            # placeholder_symbol_dto['evidence_curies'] = []
            feature.allele_symbol_dto = placeholder_symbol_dto
        return

    def synthesize_mutation_type(self, session, allele):
        """Determine mutation type."""
        # Convert term for insertions.
        insertion_conversion = {
            'transposable_element_insertion_site': 'SO:0001218',    # transgenic_insertion
            'transposable_element': 'SO:0001837',                   # mobile_element_insertion
            'insertion': 'SO:0000667'                               # insertion
        }
        relevant_features = []
        relevant_features.extend(allele.args)
        relevant_features.extend(allele.dmel_insertions)
        relevant_features.extend(allele.non_dmel_insertions)
        mutation_types = {}
        for relevant_feature in relevant_features:
            # Check insertions
            if relevant_feature.type.name in insertion_conversion.keys():
                mutation_type = insertion_conversion[relevant_feature.type.name]
            # Check ARGs.
            elif relevant_feature.type_id in self.allele_mutant_type_terms:
                mutation_type = self.cvterm_dict[relevant_feature.type_id]['curie']
            else:
                continue
            # Get pubs.
            filters = (
                Feature.feature_id == relevant_feature.feature_id,
                Pub.uniquename.op('~')(self.pub_regex),
                Pub.is_obsolete.is_(False)
            )
            pub_results = session.query(Pub).\
                select_from(Feature).\
                join(FeaturePub, (FeaturePub.feature_id == Feature.feature_id)).\
                join(Pub, (Pub.pub_id == FeaturePub.pub_id)).\
                filter(*filters).\
                distinct()
            this_pub_curie_list = [self.all_pubs_dict[i.pub_id] for i in pub_results]
            try:
                mutation_types[mutation_type].extend(this_pub_curie_list)
            except KeyError:
                mutation_types[mutation_type] = this_pub_curie_list
        for mutation_type, full_pub_curie_list in mutation_types.items():
            for pub_curie in full_pub_curie_list:
                if pub_curie == 'FB:unattributed':
                    full_pub_curie_list.remove('FB:unattributed')
            mutant_type_annotation = self.generic_audited_object.copy()
            mutant_type_annotation['mutation_type_curies'] = [mutation_type]
            # mutant_type_annotation['evidence_curies'] = list(set(full_pub_curie_list))    # TEMPORARY until AGR loads are faster
            allele.allele_mutation_type_dtos.append(mutant_type_annotation)
        return

    def add_data_provider_info(self, allele):
        """Add data_provider info."""
        allele.data_provider_dto = self.generic_data_provider_dto.copy()
        allele.data_provider_dto['cross_reference_dto'] = self.generic_cross_reference_dto.copy()
        allele.data_provider_dto['cross_reference_dto']['referenced_curie'] = f'FB:{allele.feature.uniquename}'
        allele.data_provider_dto['cross_reference_dto']['display_name'] = allele.allele_symbol_dto['display_text']
        return

    def flag_internal_alleles(self, allele):
        """Flag alleles as internal."""
        allele.allele_database_status_dto = self.generic_audited_object.copy()
        allele.allele_database_status_dto['database_status_name'] = 'approved'
        # Allele database status may change depending on checks below.
        if allele.obsolete is True:
            allele.internal = True
            allele.internal_reasons.append('Obsolete')
            allele.allele_database_status_dto['database_status_name'] = 'deleted'
        # if allele.organism_abbr != 'Dmel':
        #     allele.internal = True
        #     allele.internal_reasons.append('Non-Dmel')
        if allele.allele_of_internal_gene is True:
            allele.internal = True
            allele.internal_reasons.append('Allele of internal Dmel gene type.')
        return

    def flag_unexportable_alleles(self, allele):
        """Flag alleles missing data required for export."""
        # TEMPORARY: Suppress non-Dmel alleles from export.
        if allele.adj_organism_abbr != 'Dmel':
            allele.for_alliance_export = False
            allele.export_warnings.append(f'Suppress non-Dmel allele from export: ORG={allele.organism_abbr}')
        # Suppress objects missing required information from export.
        for attr in self.required_fields:
            if attr not in allele.__dict__.keys():
                allele.for_alliance_export = False
                allele.export_warnings.append('Missing "{}" attribute'.format(attr))
            elif getattr(allele, attr) is None:
                allele.for_alliance_export = False
                allele.export_warnings.append('Missing value for "{}" attribute'.format(attr))
        if allele.for_alliance_export is True:
            log.debug('EXPORT {}'.format(allele.mod_entity_id))
        return

    def synthesize_info(self, session):
        """Convert FlyBase allele data into an AllianceAllele representation."""
        log.info('Synthesizing allele info.')
        for allele in self.allele_dict.values():
            log.debug('Evaluating annotation: {}'.format(allele))
            self.synthesize_collections(allele)
            self.synthesize_timestamps(allele)
            self.synthesize_synonyms(allele)
            self.synthesize_secondary_ids(allele)
            self.synthesize_xrefs(allele)
            self.synthesize_references(allele)
            self.synthesize_insertions(allele)
            self.synthesize_extinction(allele)
            self.synthesize_inheritance_mode(allele)
            self.synthesize_mutation_type(session, allele)
            self.add_data_provider_info(allele)
            self.flag_internal_alleles(allele)
            self.flag_unexportable_alleles(allele)
        log.info('Done synthesizing allele info.')
        return

    def query_chado_and_export(self, session):
        """A wrapper method that runs initial db queries."""
        self.get_cvterm_info(session)
        self.get_key_cvterm_sets(session)
        self.get_all_references(session)
        self.get_alleles(session)
        self.get_direct_collections(session)
        self.get_indirect_collections(session)
        self.get_sf_collections(session)
        self.get_allele_gene(session)
        self.flag_alleles_of_internal_genes(session)
        self.flag_in_vitro_alleles(session)
        self.get_allele_constructs(session)
        self.get_allele_insertions(session)
        self.get_drosophilid_organisms(session)
        self.adjust_allele_org(session)
        self.get_allele_taxons(session)
        self.get_synonyms(session)
        self.get_allele_timestamps(session)
        self.get_allele_dbxrefs(session)
        self.get_allele_references(session)
        self.get_allele_featureprops(session)
        self.get_args(session)
        self.get_phenotypes(session)
        self.synthesize_info(session)
        return

    def generate_export_file(self):
        """Process alleles and print to a LinkML-compliant JSON file."""
        log.info('Generating output JSON file of alleles.')
        output_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
            'allele_ingest_set': []
        }
        for allele in self.allele_dict.values():
            if allele.for_alliance_export is False:
                log.debug('Suppress allele from export: {}. Reasons: {}'.format(allele, '; '.join(allele.export_warnings)))
                continue
            self.export_feat_cnt += 1
            if allele.internal is True:
                self.internal_feat_cnt += 1
                log.debug('Mark allele as internal: {}. Reasons: {}'.format(allele, '; '.join(allele.internal_reasons)))
            output_allele = {}
            for attr in self.output_fields:
                if getattr(allele, attr) is not None and getattr(allele, attr) != []:
                    output_allele[attr] = getattr(allele, attr)
            output_dict['allele_ingest_set'].append(output_allele)
        log.info('Writing data to output file.')
        with open(output_filename, 'w') as outfile:
            json.dump(output_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
            outfile.close()
        log.info('Done writing data to output file.')
        total_public_feat_cnt = self.export_feat_cnt - self.internal_feat_cnt
        log.info('Exported {} of {} alleles ({} are public).'.
                 format(self.export_feat_cnt, self.total_feat_cnt, total_public_feat_cnt))
        log.info('Suppressed {} alleles from export.'.format(self.total_feat_cnt - self.export_feat_cnt))
        return


def db_query_transaction(object_to_execute):
    """Query the chado database given an object that has a "query_chado_and_export()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 (object_to_execute): Some object that has an SQL ORM "query_chado_and_export()" method.

    Returns:
        None.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        object_to_execute.query_chado_and_export(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    return


if __name__ == "__main__":
    main()
