"""Module:: tool_handler.

Synopsis:
    A data handler that exports FlyBase data for experimental tools to Alliance Gene LinkML
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
        self.datatype = 'engineered_region'
        self.fb_export_type = fb_datatypes.FBTool
        self.agr_export_type = agr_datatypes.TransgenicToolDTO
        self.primary_export_set = 'transgenic_tool_ingest_set'
        self.ncbi_taxon_id = 'NCBITaxon:7227'    # Default Dmel, adjusted later if needed.

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
        self.build_feature_lookup(session)
        return

    # Elaborate on map_fb_data_to_alliance() for the ExpToolHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_tool_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_tool_basic(self):
        """Map basic FlyBase transgenic tool data to the Alliance LinkML object."""
        self.log.info('Map basic transgenic tool info to Alliance object.')
        for tool in self.fb_data_entities.values():
            agr_tool = self.agr_export_type()
            agr_tool.obsolete = tool.chado_obj.is_obsolete
            agr_tool.primary_external_id = f'FB:{tool.uniquename}'
            # agr_gene.mod_internal_id = f'FB.feature_id={gene.db_primary_id}'
            agr_tool.taxon_curie = tool.ncbi_taxon_id
            tool.linkmldto = agr_tool
        return
