"""Module:: tool_handler.

Synopsis:
    A data handler that exports FlyBase data for experimental tools to Alliance
    TransgenicTool LinkML objects.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

# import csv
# import re
from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class ExperimentalToolHandler(FeatureHandler):
    """This object gets, synthesizes and filters exp tool data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the ExperimentalToolHandler object."""
        super().__init__(log, testing)
        self.datatype = 'tool'
        self.fb_export_type = fb_datatypes.FBTool
        self.agr_export_type = agr_datatypes.TransgenicToolDTO
        self.primary_export_set = 'transgenic_tool_ingest_set'

    test_set = {
        'FBto0000001': 'C-Cerulean',  # First one
        'FBto0000027': 'EGFP',
        'FBto0000417': 'sgGFP',
        'FBto0000921': 'Sapphire',
        'FBto0000606': 'AflIII',    # Has UniProtKB:E3VX96
    }

    transgenic_tool_prop_to_note_mapping = {
        'description': ('summary', 'note_dtos'),
        'misc': ('comment', 'note_dtos'),
        # At the moment, just for code development. (line below)
        # 'internal_notes': ('internal_note', 'note_dtos'),
    }
    tool_associations = []
    tool_tool_rels = {}

    def get_general_data(self, session):
        """Extend the method for the AlleleHandler."""
        super().get_general_data(session)
        self.build_feature_lookup(session, feature_types=['tool'])
        return

    def get_datatype_data(self, session):
        """Extend the method for the ExperimentalToolHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_relationships(session, 'object', rel_type='compatible_tool',
                                      entity_type='engineered_region', entity_regex=self.regex['tool'])

        # self.build_feature_lookup(session)

    def map_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for a FlyBase entity."""
        self.log.info('Map secondary IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.linkmldto is None:
                continue
            secondary_id_dtos = []
            for secondary_id in fb_data_entity.alt_fb_ids:
                secondary_id_dtos.append(secondary_id)
            sec_id_list = getattr(fb_data_entity.linkmldto, slot_name)
            sec_id_list.extend(secondary_id_dtos)
        return

    # Elaborate on map_fb_data_to_alliance() for the ExpToolHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_tool_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_entity_props_to_notes('transgenic_tool_prop_to_note_mapping')
        self.map_secondary_ids('secondary_identifiers')
        self.map_tool_associations()

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_tool_basic(self):
        """Map basic FlyBase transgenic tool data to the Alliance LinkML object."""
        self.log.info('Map basic transgenic tool info to Alliance object.')
        for tool in self.fb_data_entities.values():
            agr_tool = self.agr_export_type()
            agr_tool.obsolete = tool.chado_obj.is_obsolete
            agr_tool.primary_external_id = f'FB:{tool.uniquename}'
            tool.linkmldto = agr_tool
        return

    def synthesize_tool_associations(self):
        """Get tool relationships."""
        self.log.info('Synthesize transgenic tool.')
        sub_tool_counter = 0
        obj_tool_counter = 0
        for tool in self.fb_data_entities.values():
            relevant_tool_rels = tool.recall_relationships(self.log, entity_role='object', rel_types='compatible_tool')
            # rel_entity_types='engineered_region')
            if relevant_tool_rels:
                sub_tool_counter += 1
            # self.log.debug(f'For {gene}, found {len(relevant_tool_rels)} tool rels to review.')
            for tool_rel in relevant_tool_rels:
                try:
                    tool_tool_key = (tool_rel.chado_obj.object_id, tool_rel.chado_obj.subject_id)
                except AttributeError:
                    self.log.error(f"problem {tool} {tool_rel}")
                    raise
                try:
                    self.tool_tool_rels[tool_tool_key].append(tool_rel)
                except KeyError:
                    self.tool_tool_rels[tool_tool_key] = [tool_rel]
                    obj_tool_counter += 1
        self.log.info(f'Found {obj_tool_counter} tools for {sub_tool_counter} tools.')
        return

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_tool_associations()

    def map_tool_associations(self):
        """Map transgenic tool associations to Alliance object."""
        self.log.info('Map tool associations to Alliance object.')
        OBJECT = 0
        SUBJECT = 1
        counter = 0
        # First, find alleles that have many associated genes (especially FBti insertions that hit many genes).
        # This will affect the allele-gene relation_type term used.
        tool_tool_counter = {}
        for tool_tool_key in self.tool_tool_rels.keys():
            self.log.debug(f'Mapping {tool_tool_key} to Alliance object. {self.tool_tool_rels[tool_tool_key]}')
            try:
                tool_tool_counter[tool_tool_key[OBJECT]] += 1
            except KeyError:
                tool_tool_counter[tool_tool_key[OBJECT]] = 1
        # Now, go through alleles and make the allele-gene associations.
        for tool_tool_key, tool_tool_rels in self.tool_tool_rels.items():
            object_feature_id = tool_tool_key[OBJECT]
            f_object = self.fb_data_entities[object_feature_id]
            object_curie = f'FB:{f_object.uniquename}'
            subject = self.feature_lookup[tool_tool_key[SUBJECT]]
            subject_curie = f'FB:{subject["uniquename"]}'
            first_feat_rel = tool_tool_rels[0]
            all_pub_ids = []
            for tool_tool_rel in tool_tool_rels:
                all_pub_ids.extend(tool_tool_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            # NOTE: pub 383755 | FlyBase Experimental Tool information Is the only one used
            # for tools. But not in lookup pub curies!
            # So pub_curies will be empty.
            pub_curies = self.lookup_pub_curies(all_pub_ids)

            # Adjust allele-gene relation_type as needed.
            rel_type_name = 'compatible_tool'
            rel_dto = agr_datatypes.TransgenicToolAssociationDTO(subject_curie, object_curie,
                 pub_curies, False, rel_type_name)
            if f_object.is_obsolete is True or subject['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            first_feat_rel.linkmldto = rel_dto
            self.tool_associations.append(first_feat_rel)
            counter += 1
        self.log.info(f'Generated {counter} tool-tool unique associations.')
        return

    # Elaborate on query_chado_and_export() for the TransgenicToolHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the TransgenicToolHandler."""
        super().query_chado_and_export(session)
        self.generate_export_dict(self.tool_associations, 'tool_association_ingest_set')

        return
