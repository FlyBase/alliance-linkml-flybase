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
from harvdev_utils.production import (
    Cv, Cvterm, Db, Dbxref, Organism, OrganismDbxref, Pub, PubDbxref
)


# Classes
class DataHandler(object):
    """A generic data handler that gets FlyBase data and maps it to the Alliance LinkML model."""
    def __init__(self, log: Logger, fb_data_type: str, agr_data_type: str):
        """Create the generic DataHandler object.

        Args:
            log (Logger): The global Logger object in the script using the DataHandler.
            fb_data_type (str): The FlyBase data class being handled.
            agr_data_type (str): The Alliance ingest_set to which FlyBase data is mapped: e.g., allele_ingest_set.

        """
        self.log = log
        self.fb_data_type = fb_data_type
        self.agr_data_type = agr_data_type

        # Trackers and general data collectors.
        self.input_count = 0            # Count of entities found in FlyBase chado database.
        self.export_count = 0           # Count of exported Alliance entities.
        self.internal_count = 0         # Count of exported entities marked as internal.
        self.export_data = []           # List of data objects for export (as Alliance ingest set).
        self.all_pubs_dict = {}         # A pub_id-keyed dict of pub curies (PMID or FBrf).

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

        # LinkML attributes to export: lists will be expanded for each data class.
        self.required_fields = []
        self.output_fields = []
        self.fb_agr_db_dict = {'FlyBase': 'FB'}

    def get_all_references(self, session):
        """Get all references."""
        self.log.info('Get all references.')
        # First get all current pubs having an FBrf uniquename.
        filters = (
            Pub.uniquename.op('~')(self.pub_regex),
            Pub.is_obsolete.is_(False)
        )
        results = session.query(Pub).\
            filter(*filters).\
            distinct()
        pub_counter = 0
        for pub in results:
            self.all_pubs_dict[pub.pub_id] = f'FB:{pub.uniquename}'
            pub_counter += 1
        # Next find PMIDs if available and replace the curie in the all_pubs_dict.
        filters = (
            Pub.uniquename.op('~')(self.pub_regex),
            Pub.is_obsolete.is_(False),
            Db.name == 'pubmed',
            PubDbxref.is_current.is_(True)
        )
        pmid_xrefs = session.query(Pub, Dbxref).\
            join(PubDbxref, (PubDbxref.pub_id == Pub.pub_id)).\
            join(Dbxref, (Dbxref.dbxref_id == PubDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        pmid_counter = 0
        for xref in pmid_xrefs:
            self.all_pubs_dict[xref.Pub.pub_id] = f'PMID:{xref.Dbxref.accession}'
            pmid_counter += 1
        self.log.info(f'Found {pmid_counter} PMID IDs for {pub_counter} current FB publications.')
        return

    # Methods
    def query_chado(self, session):
        """A wrapper method that runs db queries."""
        self.log.info(f'This DataHandler will map FlyBase "{self.fb_data_type}" to Alliance "{self.agr_data_type}".')
        self.get_all_references(session)
        return

    def generate_export_dict(self, input_data: list):
        """Generate LinkML export dict from list of data objects.

        Args:
            input_data (list): A list of data objects representing FlyBase data, each object having an agr_mapping
                               attribute representing a corresponding LinkML dict of transformed data. The LinkML
                               dict will be filtered to remove un-exportable objects, and to suppress empty fields.
                               The filtered dicts will be added to the handler's self.export_data list.

        """
        for i in input_data:
            self.total_input_count += 1
            if i.for_alliance_export is False:
                self.log.debug(f'Suppress {i} from export: {i.export_warnings}')
                continue
            self.total_export_count += 1
            if i.agr_mapping.internal is True:
                self.internal_count += 1
                self.log.debug(f'Export {i} but keep internal at the Alliance: {i.internal_reasons}')
            export_agr_dict = {}
            for attr in self.output_fields:
                if getattr(i.agr_mapping, attr) is not None and getattr(i.agr_mapping, attr) != []:
                    export_agr_dict[attr] = getattr(i.agr_mapping, attr)
            self.export_data.append(export_agr_dict)
        public_count = self.export_count - self.internal_count
        self.log.info(f'Exported {self.export_count} of {self.input_count} {self.fb_data_type} entities.')
        self.log.info(f'{public_count} of {self.export_count} exported {self.fb_data_type}s are PUBLIC.')
        self.log.info(f'{self.internal_count} of {self.export_count} exported {self.fb_data_type}s are INTERNAL.')
        return


class StrainHandler(DataHandler):
    """This object gets, synthesizes and filters strain data for export."""
    # def __init__(self, log, fb_data_type, agr_data_type):
    def __init__(self, log, fb_data_type, agr_data_type, extra):
        """Create the StrainHandler object."""
        super().__init__(log, fb_data_type, agr_data_type)
        self.extra = extra       # BOB: Any extra params you want here.
        self.strain_dict = {}    # A curie-keyed dict of AllianceStrainAGM objects.
        self.strain_regex = r'^FBsn[0-9]{7}$'

    def query_chado(self, session):
        """Extend the query_chado method of the StrainHandler object here."""
        # Run the inherited portion of the method first.
        super().query_chado(session)
        # Then run the extra bits.
        self.log.info(f'BOB: Running query_chado modified for FlyBase {self.extra} data.')
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
