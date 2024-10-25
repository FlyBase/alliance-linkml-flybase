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
from agm_handlers import StrainHandler
from allele_handler import AlleleHandler
from construct_handler import ConstructHandler
from disease_handlers import AlleleDiseaseHandler
from gene_handler import GeneHandler


def get_handler(log: Logger, fb_data_type: str, testing: bool):
    """Return the appropriate type of data handler.

    Args:
        log (Logger): The global Logger object in the script using the DataHandler.
        fb_data_type (str): The FB data type to export: e.g., strain, genotype.
        testing (bool): Whether the handler is being run in test mode or not.

    Returns:
        A data handler of the appropriate type for the FB data type.

    Raises:
        Raises a KeyError if the FB data type is not recognized.

    """
    log.info(f'Get handler for {fb_data_type}.')
    handler_dict = {
        'gene': GeneHandler,
        'allele': AlleleHandler,
        'construct': ConstructHandler,
        # 'variation': VariationHandler,
        'strain': StrainHandler,
        # 'genotype': GenotypeHandler,
        'disease': AlleleDiseaseHandler,
    }
    try:
        data_handler = handler_dict[fb_data_type](log, fb_data_type, testing)
        log.info(f'Returning: {data_handler}')
    except KeyError:
        log.error(f'Unrecognized FB data type and/or Alliance ingest set: {fb_data_type}.')
        raise
    return data_handler


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
        object_to_execute.query_chado_and_export(session)
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
