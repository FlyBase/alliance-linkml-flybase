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
from harvdev_utils.reporting import Cvterm, Feature, FeatureRelationship, FeatureCvterm
from sqlalchemy.orm import aliased


class CassetteHandler(FeatureHandler):
    """This object gets, synthesizes and filters cassette data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the CassetteHandler object."""
        super().__init__(log, testing)
        self.datatype = 'cassette'
        self.fb_export_type = fb_datatypes.FBCassette
        self.agr_export_type = agr_datatypes.CassetteDTO
        self.primary_export_set = 'cassette_ingest_set'
        self.incremental_update = False

    test_set = {
        'FBal0322755': 'Mcm3[+tBa]',
        'FBal0322754': 'flfl[DeltaRanBD.UAS.Venus]',
        'FBal0296109': 'sSemp1[R41G.UAS]',
        'FBal0193766': 'Gr63a[UAS.cJa]',
        'FBal0239883': 'sd[RNAi.N.UAS]',
        'FBal0000531': 'Amy-p[IX]',  # in vitro
        'FBal0028766': 'Adh[DeltaIVS1]',  # in vitro
        'FBal0410565': 'Cdkl[KD.UAS]',  # in vitro
    }

    def get_general_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_feature_lookup(session, feature_types=['cassette'])

    def get_entities(self, session, **kwargs):
        """Extend the method for the CassetteHandler."""
        reference_set = False
        if 'reference' in kwargs.keys() and kwargs['reference'] is True:
            reference_set = True
            self.incremental_update = True

        # Get the main set of cassettes
        self.get_main_entities(session, reference_set)
        # Get in vitro set of cassettes
        self.add_in_vitro_allele_entries(session, reference_set)

    def add_in_vitro_allele_entries(self, session, reference_set):
        """Extend list of entities."""
        self.log.info('Add entities for alleles having "in vitro construct" annotations.')

        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['allele']),
            Cvterm.name == 'in vitro construct',
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Feature). \
            select_from(Feature). \
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)). \
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)). \
            filter(*filters). \
            distinct()
        for result in results:
            pkey_id = getattr(result, 'feature_id')
            if reference_set is True:
                self.fb_reference_entity_ids.append(pkey_id)
            else:
                self.fb_data_entities[pkey_id] = self.fb_export_type(result)

    def get_main_entities(self, session, reference_set):
        """Get simple FlyBase cassette/allele data entities.

        So Subject FBal linked to
           Object FBtp  linked via
           Cvterm 'associated_with'
        Args:
            session (Session): SQLAlchemy session for the query.
            reference_set (bool): If True, retrieves only non-obsolete objects from
                              a previous reference database; for incremental
                              updates.

        """
        if self.datatype in self.feat_type_export.keys():
            chado_type = 'feature'
        else:
            chado_type = self.datatype

        if reference_set is True:
            mess = f'Get {self.datatype} data entities from {chado_type}'
            mess += ' table (previous reference db for incremental update).'
            self.log.info(mess)
        else:
            self.log.info(f'Get {self.datatype} data entities from {chado_type} table.')

        feat_object = aliased(Feature)
        rel_type = aliased(Cvterm)
        chado_table = self.chado_tables['main_table'][chado_type]

        filters = ()
        # Subject is the allele and has s specific regex.
        filters += (chado_table.uniquename.op('~')(self.regex[self.datatype]), )
        # Allele is not obsolete
        filters += (chado_table.is_obsolete.is_(False), )
        # Object should be a FBtp object (construct)
        filters += (feat_object.uniquename.op('~')(self.regex['construct']), )
        # Construct is not obsolete
        filters += (feat_object.is_obsolete.is_(False),)
        # 'associated_with' in relationship table
        filters += (rel_type.name == 'associated_with', )

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

        results = session.query(chado_table).\
            select_from(chado_table).\
            join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == chado_table.feature_id)).\
            join(feat_object, (feat_object.feature_id == FeatureRelationship.object_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
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

    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_cassette_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        # self.map_xrefs()

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_cassette_basic(self):
        """Map basic FlyBase transgenic tool data to the Alliance LinkML object."""
        self.log.info('Map basic cassette info to Alliance object.')
        for cass in self.fb_data_entities.values():
            agr_cass = self.agr_export_type()
            agr_cass.obsolete = cass.chado_obj.is_obsolete
            agr_cass.primary_external_id = f'FB:{cass.uniquename}'
            cass.linkmldto = agr_cass

    def get_datatype_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        # self.get_entity_fb_xrefs(session)
        # self.get_entity_xrefs(session)
        # self.get_entity_relationships(session, 'object', rel_type='associated_with',
        #                              entity_type='engineered_region', entity_regex=self.regex['tool'])
