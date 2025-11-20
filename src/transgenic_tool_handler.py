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

    test_set = {
        'FBto0000001': 'C-Cerulean',  # First one
    }

    exptool_prop_to_note_mapping = {
        'misc': ('comment', 'note_dtos'),
    }

    def get_datatype_data(self, session):
        """Extend the method for the ExperimentalToolHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        return
