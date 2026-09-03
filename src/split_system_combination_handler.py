"""Module:: split_system_combination_handler.

Synopsis:
    A data handler that exports FlyBase data for split system combinations
    (FBco features) to Alliance AffectedGenomicModel (AGM) LinkML objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class SplitSystemCombinationHandler(FeatureHandler):
    """This object gets, synthesizes and filters split system combination data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the SplitSystemCombinationHandler object."""
        super().__init__(log, testing)
        self.datatype = 'split system combination'
        self.fb_export_type = fb_datatypes.FBSplitSystemCombination
        self.agr_export_type = agr_datatypes.AffectedGenomicModelDTO
        self.primary_export_set = 'agm_ingest_set'

    test_set = {
        'FBco0000001': 'Scer\\GAL4[DBD.R22B12]INTERSECTIONHsap\\RELA[AD.R14E06]',    # The first FBco; has a non-current "MB012B" symbol.
        'FBco0000010': 'Scer\\GAL4[DBD.Tdc2]INTERSECTIONHsap\\RELA[AD.R24E06]',      # DBD half is a named driver, not an R-line.
        'FBco0000054': 'Scer\\GAL4[DBD.R15B01]INTERSECTIONHsap\\RELA[AD.R14C08]',    # Has many pubs.
        'FBco0000113': 'Scer\\GAL4[DBD.R19F09]INTERSECTIONHsap\\RELA[AD.R25D01]',    # Has the most pubs (9).
        'FBco0000758': 'Scer\\GAL4[DBD.R82C10]INTERSECTIONHsap\\RELA[AD.R61H01]',    # The only obsolete FBco; its AD component allele is obsolete too.
    }

    # All FBco features are of the "synthetic construct" organism in chado, but the
    # split systems themselves are only ever used in Dmel, so report them as such.
    # This also keeps them out of the "non-Dmel is internal" bucket used for other AGMs.
    SSC_TAXON_ID = 'NCBITaxon:7227'

    # The AGM subtype reported to the Alliance for these entities.
    SSC_SUBTYPE_NAME = 'split system combination'

    # Additional export sets.
    agm_component_associations = []    # A list of AgmAlleleAssociationDTOs.

    # Elaborate on get_general_data() for the SplitSystemCombinationHandler.
    def get_general_data(self, session):
        """Extend the method for the SplitSystemCombinationHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['allele'])
        return

    # Elaborate on get_datatype_data() for the SplitSystemCombinationHandler.
    def get_datatype_data(self, session):
        """Extend the method for the SplitSystemCombinationHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject', rel_type='partially_produced_by',
                                      entity_type='allele', entity_regex=self.regex['allele'])
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        return

    # Additional sub-methods for synthesize_info().
    def synthesize_ssc_components(self):
        """Determine the allele components of each split system combination."""
        self.log.info('Determine the allele components of each split system combination.')
        component_counter = 0
        obsolete_component_counter = 0
        unknown_component_counter = 0
        for ssc in self.fb_data_entities.values():
            rels = ssc.recall_relationships(self.log, entity_role='subject', rel_types='partially_produced_by',
                                            rel_entity_types='allele')
            for rel in rels:
                feature_id = rel.chado_obj.object_id
                # Obsolete alleles are not exported to the Alliance, so an association to one would dangle there.
                if feature_id not in self.feature_lookup.keys():
                    self.log.warning(f'{ssc} has a component (feature_id={feature_id}) missing from the feature_lookup; skipping it.')
                    unknown_component_counter += 1
                    continue
                if self.feature_lookup[feature_id]['is_obsolete'] is True:
                    obs_curie = self.feature_lookup[feature_id]['curie']
                    self.log.warning(f'{ssc} has an obsolete component, {obs_curie}; skipping it.')
                    obsolete_component_counter += 1
                    continue
                ssc.component_features.append(feature_id)
                component_counter += 1
        self.log.info(f'Found {component_counter} allele components for {len(self.fb_data_entities)} split system combinations.')
        self.log.info(f'Skipped {obsolete_component_counter} obsolete allele components.')
        self.log.info(f'Skipped {unknown_component_counter} allele components not found in the feature_lookup.')
        return

    # Elaborate on synthesize_info() for the SplitSystemCombinationHandler.
    def synthesize_info(self):
        """Extend the method for the SplitSystemCombinationHandler."""
        super().synthesize_info()
        self.flag_new_additions_and_obsoletes()
        self.synthesize_ssc_components()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        return

    # Additional sub-methods for map_fb_data_to_alliance().
    def map_ssc_basic(self):
        """Map basic FlyBase split system combination data to the Alliance object."""
        self.log.info('Map basic split system combination info.')
        for ssc in self.fb_data_entities.values():
            agr_ssc = self.agr_export_type()
            agr_ssc.obsolete = ssc.chado_obj.is_obsolete
            agr_ssc.primary_external_id = f'FB:{ssc.uniquename}'
            agr_ssc.taxon_curie = self.SSC_TAXON_ID
            agr_ssc.subtype_name = self.SSC_SUBTYPE_NAME
            ssc.linkmldto = agr_ssc
        return

    def map_ssc_synonyms(self):
        """Generate split system combination name/synonym DTOs for an entity."""
        self.log.info('Map split system combination synonyms.')
        for ssc in self.fb_data_entities.values():
            if ssc.linkmldto is None:
                continue
            linkml_synonym_bins = {
                'symbol_bin': None,
                'synonym_bin': []
            }
            # FBco features only ever have "symbol" synonyms in chado.
            for syno_dict in ssc.synonym_dict.values():
                name_dto = agr_datatypes.NameSlotAnnotationDTO(syno_dict['name_type_name'], syno_dict['format_text'],
                                                               syno_dict['display_text'], syno_dict['pub_curies']).dict_export()
                name_dto['internal'] = syno_dict['is_internal']
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] == 'nomenclature_symbol':
                    linkml_synonym_bins['symbol_bin'] = name_dto
                else:
                    linkml_synonym_bins['synonym_bin'].append(name_dto)
            # The AGM LinkML class has no symbol slot, so report the current FB symbol as the AGR "full_name".
            if linkml_synonym_bins['symbol_bin']:
                linkml_synonym_bins['symbol_bin']['name_type_name'] = 'full_name'
                setattr(ssc.linkmldto, 'agm_full_name_dto', linkml_synonym_bins['symbol_bin'])
            setattr(ssc.linkmldto, 'agm_synonym_dtos', linkml_synonym_bins['synonym_bin'])
        return

    def map_ssc_components(self):
        """Map split system combination components."""
        self.log.info('Map split system combination components.')
        counter = 0
        for ssc in self.fb_data_entities.values():
            if ssc.linkmldto is None:
                continue
            for feature_id in ssc.component_features:
                ssc_allele_rel = fb_datatypes.FBExportEntity()
                component_curie = self.feature_lookup[feature_id]['curie']
                ssc_allele_rel.linkmldto = agr_datatypes.AgmAlleleAssociationDTO(ssc.linkmldto.primary_external_id,
                                                                                 component_curie, 'unspecified zygosity')
                self.agm_component_associations.append(ssc_allele_rel)
                counter += 1
        self.log.info(f'Mapped {counter} split system combination component associations.')
        return

    # Elaborate on map_fb_data_to_alliance() for the SplitSystemCombinationHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the SplitSystemCombinationHandler."""
        super().map_fb_data_to_alliance()
        self.map_ssc_basic()
        self.map_ssc_synonyms()
        self.map_ssc_components()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_pubs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        self.flag_internal_fb_entities('agm_component_associations')
        return

    # Elaborate on query_chado_and_export() for the SplitSystemCombinationHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the SplitSystemCombinationHandler."""
        super().query_chado_and_export(session)
        self.flag_unexportable_entities(self.agm_component_associations, 'agm_allele_association_ingest_set')
        self.generate_export_dict(self.agm_component_associations, 'agm_allele_association_ingest_set')
        return
