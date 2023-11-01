"""Module:: utils.

Synopsis:
    Data retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import logging
log = logging.getLogger(__name__)


# Classes
class DataHandler(object):
    """A generic data handler that gets some data type and exports it, and
       reports the extent of the export."""
    def __init__(self, data_type):
        """Create the generic DataHandler object."""
        self.data_type = data_type     # Label for the data type being handled: e.g., gene.
        self.total_input_count = 0     # Count of entities found in FlyBase chado database.
        self.total_export_count = 0    # Count of exported entities.
        self.internal_count = 0        # Count of exported entities marked as internal.
        return

    def report_label(self):
        log.info(f'BOB:DataHandler for this data_type: {self.data_type}')
        return
