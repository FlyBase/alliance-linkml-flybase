"""Module:: allele_handler.

Synopsis:
    A data handler that exports FlyBase data for alleles to Alliance Gene LinkML
    objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
from agr_datatypes import AlleleDTO, AlleleMutationTypeSlotAnnotationDTO
from fb_datatypes import FBAllele
from feature_handler import FeatureHandler
from harvdev_utils.production import (
    Cvterm, Feature, FeatureCvterm, FeatureRelationship, Organism
)


class AlleleHandler(FeatureHandler):
    """This object gets, synthesizes and filters allele data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AlleleHandler object."""
        super().__init__(log, testing)
        self.datatype = 'allele'
        self.fb_export_type = FBAllele
        self.agr_export_type = AlleleDTO
        self.primary_export_set = 'allele_ingest_set'

    test_set = {
        'FBal0137236': 'gukh[142]',    # P{hsneo}Xrp1142 insertion allele.
        'FBal0018482': 'wg[1]',        # X-ray mutation.
    }

    # Additional reference info.
    allele_class_terms = []          # A list of cvterm_ids for child terms of "allele_class" (FBcv:0000286).
    allele_mutant_type_terms = []    # A list of cvterm_ids for child terms of chromosome_structure_variation or sequence_alteration.

    # Additional sub-methods for get_general_data().
    def get_key_cvterm_sets(self, session):
        """Get key CV term sets from chado."""
        self.log.info('Get key CV term sets from chado.')
        self.allele_class_terms.extend(self.get_child_cvterms(session, 'allele class', 'FlyBase miscellaneous CV'))
        self.allele_mutant_type_terms.extend(self.get_child_cvterms(session, 'chromosome_structure_variation', 'SO'))
        self.allele_mutant_type_terms.extend(self.get_child_cvterms(session, 'sequence_alteration', 'SO'))
        return

    # Elaborate on get_general_data() for the AlleleHandler.
    def get_general_data(self, session):
        """Extend the method for the AlleleHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.get_key_cvterm_sets(session)
        self.build_ncbi_taxon_lookup(session)
        self.get_drosophilid_organisms(session)
        self.build_feature_lookup(session)
        self.find_internal_genes(session)
        self.build_featureprop_evidence_lookup(session)
        self.build_feature_relationship_evidence_lookup(session)
        return

    def get_related_features(self, session):
        """Get allele-associated features."""
        self.log.info('Get allele-associated features.')
        self.get_entity_sbj_feat_rel_by_type(session, 'parent_gene_rels', rel_type='alleleof', obj_type='gene', obj_regex=self.regex['gene'])
        self.get_entity_sbj_feat_rel_by_type(session, 'constructs', rel_type='derived_tp_assoc_alleles', obj_regex=self.regex['construct'])
        arg_types = [
            'MNV',
            'complex_substitution',
            'deletion',
            'delins',
            'insertion',
            'point_mutation',
            'sequence_alteration',
            'sequence_variant',
            # 'rescue_region',    # This type of ARG is not relevant to determining allele type (better ways).
        ]
        self.get_entity_obj_feat_rel_by_type(session, 'args', rel_type='partof', sbj_type=arg_types)
        return

    def get_associated_insertions(self, session):
        """Get allele-associated insertions."""
        self.log.info('Get allele-associated insertions')
        allele = aliased(Feature, name='allele')
        insertion = aliased(Feature, name='insertion')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            insertion.is_obsolete.is_(False),
            insertion.is_analysis.is_(False),
            insertion.uniquename.op('~')(self.regex['insertion']),
            Cvterm.name == 'associated_with'
        )
        insertion_results = session.query(Organism, allele, insertion).\
            select_from(insertion).\
            join(Organism, (Organism.organism_id == insertion.organism_id)).\
            join(FeatureRelationship, (FeatureRelationship.object_id == insertion.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.subject_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        dmel_counter = 0
        non_dmel_counter = 0
        for result in insertion_results:
            if result.Organism.abbreviation == 'Dmel':
                self.fb_data_entities[result.allele.feature_id].dmel_insertions.append(result.insertion)
                dmel_counter += 1
            else:
                self.fb_data_entities[result.allele.feature_id].non_dmel_insertions.append(result.insertion)
                non_dmel_counter += 1
        self.log.info(f'Found {dmel_counter} Dmel and {non_dmel_counter} non-Dmel insertions related to alleles.')
        return

    def find_in_vitro_alleles(self, session):
        """Find in vitro alleles."""
        self.log.info('Find in vitro alleles.')
        cvterm_name_regex = '^in vitro construct'
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            Cvterm.name.op('~')(cvterm_name_regex)
        )
        ivt_alleles = session.query(Feature).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in ivt_alleles:
            self.fb_data_entities[result.feature_id].in_vitro = True
            counter += 1
        self.log.info(f'Flagged {counter} alleles as "in vitro"')
        return

    def get_datatype_data(self, session, datatype, fb_export_type, agr_export_type):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session, datatype, fb_export_type, agr_export_type)
        self.get_entities(session, self.datatype, self.fb_export_type)
        self.get_entityprops(session, self.datatype)
        self.get_entity_pubs(session, self.datatype)
        self.get_entity_synonyms(session, self.datatype)
        self.get_entity_fb_xrefs(session, self.datatype)
        self.get_entity_xrefs(session, self.datatype)
        self.get_entity_timestamps(session, self.datatype)
        self.get_related_features(session)
        self.get_associated_insertions(session)
        self.find_in_vitro_alleles(session)
        return

    # Add sub-methods to be run by synthesize_info() below.
    def synthesize_parent_genes(self):
        """Get allele parent gene IDs."""
        self.log.info('Get allele parent gene IDs.')
        allele_counter = 0
        for allele in self.fb_data_entities.values():
            parent_gene_ids = []
            for feat_rel in allele.parent_gene_rels:
                parent_gene = self.feature_lookup[feat_rel.object_id]
                if parent_gene['is_obsolete'] is False:
                    parent_gene_ids.append(parent_gene['uniquename'])
            if len(parent_gene_ids) == 1:
                allele.parent_gene_id = parent_gene_ids[0]
                allele_counter += 1
            elif len(parent_gene_ids) == 0:
                self.log.warning(f'{allele} has no parent gene!')
            else:
                self.log.warning(f'{allele} has many parent genes!')
        self.log.info(f'Found parental gene for {allele_counter} alleles.')
        return

    def flag_alleles_of_internal_genes(self):
        """Flag alleles of internal genes."""
        self.log.info('Flag alleles of internal genes.')
        for allele in self.fb_data_entities.values():
            if allele.parent_gene_id in self.internal_gene_ids:
                allele.allele_of_internal_gene = True
        return

    def adjust_allele_organism(self):
        """Adjust organism for classical non-Dmel alleles."""
        self.log.info('Adjust organism for classical non-Dmel alleles.')
        counter = 0
        for allele in self.fb_data_entities.values():
            # Skip alleles that unambiguously occur in Dmel.
            if allele.organism_abbr == 'Dmel':
                continue
            if allele.dmel_insertions:
                continue
            if allele.constructs:
                continue
            # Find clear evidence that allele is non-Dmel classical allele.
            is_non_dmel_classical = False
            if allele.non_dmel_insertions:
                is_non_dmel_classical = True
            elif allele.feature.organism_id in self.drosophilid_list and allele.in_vitro is False:
                is_non_dmel_classical = True
            # Make the adjustment.
            if is_non_dmel_classical is True:
                allele.adj_org_abbr = allele.org_abbr
                self.log.debug(f'Non-Dmel allele: id={allele.uniquename}, name={allele.name}, org_abbr={allele.org_abbr}')
                counter += 1
        self.log.info(f'Adjusted organism to be "non-Dmel" for {counter} alleles.')
        return

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_info(self, datatype, fb_export_type, agr_export_type):
        """Extend the method for the GeneHandler."""
        super().synthesize_info(datatype, fb_export_type, agr_export_type)
        self.synthesize_props(datatype)
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_parent_genes()
        self.flag_alleles_of_internal_genes()
        self.adjust_allele_organism()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_allele_basic(self, agr_export_type):
        """Map basic FlyBase allele data to the Alliance LinkML object."""
        self.log.info('Map basic allele info to Alliance object.')
        for allele in self.fb_data_entities.values():
            agr_allele = agr_export_type()
            agr_allele.obsolete = allele.chado_obj.is_obsolete
            agr_allele.mod_entity_id = f'FB:{allele.uniquename}'
            agr_allele.mod_internal_id = str(allele.chado_obj.feature_id)
            agr_allele.taxon_curie = allele.ncbi_taxon_id
            if allele.allele_of_internal_gene is True:
                agr_allele.internal = True
                allele.internal_reasons.append('Allele of internal type FB gene.')
            allele.linkmldto = agr_allele
        return

    def map_mutation_types(self):
        """Map mutation types."""
        self.log.info('Map mutation types.')
        insertion_conversion = {
            'transposable_element_insertion_site': 'SO:0001218',    # transgenic_insertion
            'transposable_element': 'SO:0001837',                   # mobile_element_insertion
            'insertion': 'SO:0000667'                               # insertion
        }
        mutation_types = {}    # Will be a dict of mutation type curies and supporting pub curies.
        for allele in self.fb_data_entities.values():
            relevant_feat_rels = []
            relevant_feat_rels.extend(allele.args)
            relevant_feat_rels.extend(allele.dmel_insertions)
            relevant_feat_rels.extend(allele.non_dmel_insertions)
            for feat_rel in relevant_feat_rels:
                type_curie = None
                fb_feat_type_id = feat_rel.object.type_id
                fb_feat_type_name = self.cvterm_lookup[fb_feat_type_id]['name']
                if fb_feat_type_name in insertion_conversion.keys():
                    type_curie = insertion_conversion[fb_feat_type_name]
                elif fb_feat_type_id in self.allele_mutant_type_terms:
                    type_curie = self.cvterm_lookup[fb_feat_type_id]['curie']
                if type_curie is None:
                    continue
                pub_ids = self.feat_rel_pub_lookup[feat_rel.feature_relationship_id]
                pub_curies = self.lookup_pub_curies(pub_ids)
                try:
                    mutation_types[type_curie].extend(pub_curies)
                except KeyError:
                    mutation_types[type_curie] = pub_curies
            for mutation_type_curie, full_pub_curie_list in mutation_types.items():
                for pub_curie in full_pub_curie_list:
                    if pub_curie == 'FB:unattributed':
                        full_pub_curie_list.remove('FB:unattributed')
                mutant_type_annotation = AlleleMutationTypeSlotAnnotationDTO(mutation_type_curie, full_pub_curie_list)
                allele.allele_mutation_type_dtos.append(mutant_type_annotation)
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_fb_data_to_alliance(self, datatype, fb_export_type, agr_export_type):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance(datatype, fb_export_type, agr_export_type)
        self.map_allele_basic(agr_export_type)
        self.map_synonyms(datatype, agr_export_type)
        self.map_data_provider_dto(datatype)
        # self.map_pubs()    # TEMPORARILY SUPPRESS UNTIL LOAD SPEED IMPROVES
        self.map_xrefs(datatype)
        self.map_timestamps()
        self.map_mutation_types()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Elaborate on query_chado_and_export() for the GeneHandler.
    def query_chado_and_export(self, session, datatype, fb_export_type, agr_export_type):
        """Elaborate on query_chado_and_export method for the GeneHandler."""
        super().query_chado_and_export(session, datatype, fb_export_type, agr_export_type)
        return
