"""Module:: tool_handler.

Synopsis:
    A data handler that exports FlyBase data for experimental tools to Alliance TransgenicTool LinkML
    objects.

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
    }

    transgenic_tool_prop_to_note_mapping = {
        'misc': ('comment', 'note_dtos'),
    }

    def get_datatype_data(self, session):
        """Extend the method for the ExperimentalToolHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_xrefs(session)
        self.get_entity_fb_xrefs(session)
        # self.build_feature_lookup(session)
        return

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

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_tool_basic(self):
        """Map basic FlyBase transgenic tool data to the Alliance LinkML object."""
        self.log.info('Map basic transgenic tool info to Alliance object.')
        for tool in self.fb_data_entities.values():
            agr_tool = self.agr_export_type()
            agr_tool.obsolete = tool.chado_obj.is_obsolete
            agr_tool.primary_external_id = f'FB:{tool.uniquename}'
            agr_tool.taxon_curie = tool.ncbi_taxon_id
            tool.linkmldto = agr_tool
        return

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()

    # Elaborate on query_chado_and_export() for the TransgenicToolHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the TransgenicToolHandler."""
        super().query_chado_and_export(session)
        return
