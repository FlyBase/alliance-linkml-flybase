"""Module:: cassette_handler.

Synopsis:
    A data handler that exports FlyBase data for cassette to Alliance
    TransgenicTool LinkML objects.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

# import csv
# import re
from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler
from harvdev_utils.reporting import (
    Cvterm)

class CassetteHandler(FeatureHandler):
    """This object gets, synthesizes and filters cassette data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the CassetteHandler object."""
        super().__init__(log, testing)
        self.datatype = 'cassette'
        self.fb_export_type = fb_datatypes.FBCassette
        self.agr_export_type = agr_datatypes.CassetteDTO
        self.primary_export_set = 'cassette_ingest_set'

    test_set = {
        'FBal0028611': 'C-Cerulean',  # ivc only
        'FBto0000027': 'EGFP',
        'FBto0000417': 'sgGFP',
        'FBto0000921': 'Sapphire',
        'FBto0000606': 'AflIII',    # Has UniProtKB:E3VX96
    }

    def get_general_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_feature_lookup(session, feature_types=['cassette'])
        return

    def get_entities(self, session, **kwargs):
        """Get primary FlyBase data entities.

        Args:
            session (Session): SQLAlchemy session for the query.

        Kwargs:
            reference (bool): If True, retrieves only non-obsolete objects from
                              a previous reference database; for incremental
                              updates.

        """
        reference_set = False
        if 'reference' in kwargs.keys() and kwargs['reference'] is True:
            reference_set = True
            self.incremental_update = True
        if self.datatype in self.feat_type_export.keys():
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        if reference_set is True:
            self.log.info(f'Get {self.datatype} data entities from {chado_type} table (previous reference db for incremental update).')
        else:
            self.log.info(f'Get {self.datatype} data entities from {chado_type} table.')
        chado_table = self.chado_tables['main_table'][chado_type]
        filters = ()
        if self.datatype in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.datatype]}')
            filters += (chado_table.uniquename.op('~')(self.regex[self.datatype]), )
        if self.datatype in self.feature_subtypes.keys():
            sub = self.feature_subtypes[self.datatype]
            self.log.info(f'Filter main table by these feature_subtypes: {sub}')
            filters += (Cvterm.name.in_((self.feature_subtypes[self.datatype])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (chado_table.uniquename.in_((self.test_set.keys())), )
        if filters == ():
            self.log.warning('Have no filters for the main FlyBase entity driver query.')
            raise

        if self.datatype in self.feature_subtypes.keys():
            results = session.query(chado_table).\
                select_from(chado_table).\
                join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
                filter(*filters).\
                distinct()
        else:
            results = session.query(chado_table).\
                filter(*filters).\
                distinct()
        pkey_name = f'{chado_type}_id'
        self.log.info(f'Have this primary_key name: {pkey_name}')
        counter = 0
        for result in results:
            pkey_id = getattr(result, pkey_name)
            if reference_set is True:
                self.fb_reference_entity_ids.append(pkey_id)
            else:
                self.fb_data_entities[pkey_id] = self.fb_export_type(result)
            counter += 1
        if reference_set is True:
            self.log.info(f'Found {counter} FlyBase {self.datatype} entities in reference chado instance.')
        else:
            self.log.info(f'Found {counter} FlyBase {self.datatype} entities in chado.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_relationships(session, 'object', rel_type='compatible_tool',
                                      entity_type='engineered_region', entity_regex=self.regex['tool'])
