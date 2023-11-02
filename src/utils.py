"""Module:: utils.

Synopsis:
    Data retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

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
            agr_data_type (str): The Alliance data class to which FlyBase data is being mapped.

        """
        self.log = log
        self.fb_data_type = fb_data_type
        self.agr_data_type = agr_data_type
        self.total_input_count = 0            # Count of entities found in FlyBase chado database.
        self.total_export_count = 0           # Count of exported Alliance entities.
        self.internal_count = 0               # Count of exported entities marked as internal.
        return

    def query_chado(self, session):
        """Test."""
        self.log.info(f'This DataHandler is mapping FlyBase {self.fb_data_type} to {self.agr_data_type}.')
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
