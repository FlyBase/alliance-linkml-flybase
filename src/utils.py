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


def export_chado_data(session: Session, log: Logger, object_to_execute: DataHandler, **kwargs):
    """Query the chado database given an object that has a "query_chado_and_export()" method.

    Args:
        session (Session): SQLAlchemy session for the database from which to query and export.
        log (Logger): The global Logger object in the script using the DataHandler.
        object_to_execute (DataHandler): An object having a query_chado_and_export() method.

    Kwargs:
        reference_session (Session|None): SQLAlchemy session for an earlier reference database (for incremental export).

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    global TESTING
    if 'reference_session' in kwargs.keys():
        try:
            object_to_execute.get_entities(kwargs['reference_session'], reference=True)
        except RuntimeError:
            kwargs['reference_session'].rollback()
            log.critical('Critical transaction error occurred during reference db chado query; rolling back and exiting.')
            raise
    try:
        object_to_execute.query_chado_and_export(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during main chado query; rolling back and exiting.')
        raise
    if TESTING is True:
        log.info('Since "testing" is True, rolling back all transactions.')
        session.rollback()
    else:
        log.info('Since "testing" is False, committing transactions.')
        session.commit()
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
