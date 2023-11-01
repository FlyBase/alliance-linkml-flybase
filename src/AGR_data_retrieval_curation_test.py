# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TEST: Data retrieval of FlyBase data for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_test.py [-h] [-r FLYBASE_RELEASE][-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_test.py -v -r FB2023_06_EP7 -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script makes a JSON file conforming to LinkML specs for the curation
    (i.e., "persistent") database.

"""

import argparse
import datetime
import json
import strict_rfc3339
from tqdm import tqdm
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, OrganismDbxref, Pub, PubDbxref, Strain, StrainDbxref,
    StrainPub, StrainSynonym, Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading
from utils import DataHandler

# Now proceed with generic setup.
report_label = 'test'
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
reference_assembly = set_up_dict['assembly']
input_dir = set_up_dict['input_dir']
output_filename = set_up_dict['output_filename'].replace('tsv', 'json')
log = set_up_dict['log']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-l', '--linkml_release', help='The "agr_curation_schema" LinkML release number.', required=True)
parser.add_argument('-r', '--fb_release', help='The FlyBase data release from which data was obtained.', required=True)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release
fb_release = args.fb_release

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase AGMs."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')
    log.info('Output JSON file corresponds to "agr_curation_schema" release: {}'.format(linkml_release))

    # Instantiate the object, get the data, synthesize it, export it.
    agm_handler = DataHandler('billy')
    agm_handler.report_label()
    log.info('Ended main function.\n')


def db_query_transaction(object_to_execute):
    """Query the chado database given an object that has a "query_chado()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 (object_to_execute): Some object that has an SQL ORM "query_chado()" method.

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


if __name__ == "__main__":
    main()
