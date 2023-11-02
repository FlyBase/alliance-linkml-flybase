"""Module:: utils.

Synopsis:
    Data retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import json
from logging import Logger
from sqlalchemy.orm import Session


# Classes
class DataHandler(object):
    """A generic data handler that gets FlyBase data and maps it to the Alliance LinkML model."""
    def __init__(self, log: Logger, fb_data_type: str, agr_data_type: str):
        """Create the generic DataHandler object.

        Args:
            log (Logger): The global Logger object in the script using the DataHandler.
            fb_data_type (str): The FlyBase data class being handled.
            agr_data_type (str): The Alliance ingest_set to which FlyBase data is being mapped: e.g., allele_ingest_set.

        """
        self.log = log
        self.fb_data_type = fb_data_type
        self.agr_data_type = agr_data_type

        # Trackers and general data collectors.
        self.total_input_count = 0            # Count of entities found in FlyBase chado database.
        self.total_export_count = 0           # Count of exported Alliance entities.
        self.internal_count = 0               # Count of exported entities marked as internal.
        self.export_data = []                 # List of data objects for export (as Alliance ingest set).
        self.all_pubs_dict = {}               # A pub_id-keyed dict of pub curies (PMID or FBrf).

        # Generic information.
        self.pub_regex = r'^(FBrf[0-9]{7}|unattributed)$'
        self.generic_audited_object = {
            'internal': False,
            'obsolete': False,
            'created_by_curie': 'FB:FB_curator',
            'updated_by_curie': 'FB:FB_curator'
        }
        self.generic_data_provider_dto = self.generic_audited_object.copy()
        self.generic_data_provider_dto['source_organization_abbreviation'] = 'FB'
        self.generic_cross_reference_dto = {'prefix': 'FB', 'internal': False}

    # Methods
    def query_chado(self, session):
        """Test."""
        self.log.info(f'This DataHandler is mapping FlyBase "{self.fb_data_type}" to Alliance "{self.agr_data_type}".')
        return


# Functions
def db_query_transaction(session: Session, log: Logger, object_to_execute: DataHandler):
    """Query the chado database given an object that has a "query_chado()" method.

    Args:
        session (Session): SQLAlchemy session for db queries.
        log (Logger): The global Logger object in the script using the DataHandler.
        object_to_execute (DataHandler): An object having a query_chado() method.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    try:
        object_to_execute.query_chado(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    return


def generate_export_file(export_dict: dict, log: Logger, output_filename: str):
    """Print Alliance LinkML data to JSON file.

    Args:
        export_dict (dict): A LinkML dict including some "ingest" list of data elements.
        log (Logger): The global Logger object in the script calling this function.
        output_filename (str): The global output_filename in the script calling this function.

    """
    log.info('Writing output Alliance LinkML data dict to JSON file.')
    with open(output_filename, 'w') as outfile:
        json.dump(export_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
        outfile.close()
    log.info('Done writing data to output JSON file.')
    return
