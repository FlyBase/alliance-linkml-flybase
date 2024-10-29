"""Module:: utils.

Synopsis:
    Functions that select and run data handlers, and export FlyBase chado data
    in Alliance LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import json
from logging import Logger
from sqlalchemy.orm import Session
from handler import DataHandler


def db_query_transaction(session: Session, log: Logger, object_to_execute: DataHandler):
    """Query the chado database given an object that has a "query_chado_and_export()" method.

    Args:
        session (Session): SQLAlchemy session for db queries.
        log (Logger): The global Logger object in the script using the DataHandler.
        object_to_execute (DataHandler): An object having a query_chado_and_export() method.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    try:
        object_to_execute.query_chado_and_export(session, object_to_execute.datatype, object_to_execute.fb_export_type, object_to_execute.agr_export_type)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    return


def generate_export_file(export_dict: dict, log: Logger, output_filename: str):
    """Print Alliance LinkML data to JSON file.

    Args:
        export_dict (dict): A dict of LinkML dicts for some "agr_ingest" set.
        log (Logger): The global Logger object in the script calling this function.
        output_filename (str): The global output_filename in the script calling this function.

    """
    log.info('Writing output Alliance LinkML data dict to JSON file.'.upper())
    with open(output_filename, 'w') as outfile:
        json.dump(export_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
        outfile.close()
    return
