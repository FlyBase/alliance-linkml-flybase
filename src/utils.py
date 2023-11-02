"""Module:: utils.

Synopsis:
    Data retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

# Classes
class DataHandler(object):
    """A generic data handler that gets some data type and exports it, and
       reports the extent of the export."""
    def __init__(self, data_type):
        """Create the generic DataHandler object.
        
        Args:
            data_type (arg1): (str) The label for the data type to be handled: e.g., gene.
            log (arg2): (Logger) A global Logger object from the script using the DataHandler.
        
        """
        self.data_type = data_type     # Label for the data type being handled: e.g., gene.
        self.total_input_count = 0     # Count of entities found in FlyBase chado database.
        self.total_export_count = 0    # Count of exported entities.
        self.internal_count = 0        # Count of exported entities marked as internal.
        return

    def report_label(self, log):
        log.info(f'BOB:DataHandler for this data_type: {self.data_type}')
        return


# Functions
def db_query_transaction(session, log, object_to_execute):
    """Query the chado database given an object that has a "query_chado()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 ()
        arg2 (log): (Logger) A global Logger object from the script using the DataHandler.
        arg3 (object_to_execute): Some object that has an SQL ORM "query_chado()" method.

    Returns:
        None.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        object_to_execute.query_chado(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    return