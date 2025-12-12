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
        'FBal0000531': 'Amy-p[IX]',
        'FBal0028742': 'Act88F[E334K]',
        'FBal0212171': r'Avic\GFP[UAS.FRT1]',  # in vitro only
        'FBal0290956': 'Csas[21]',  # in vitro only
        'FBal0392043': r'Avic\GFP[EYFP.3xP3.cUa]',  # in vitro only
    }

    cassette_prop_to_note_mapping = {
        'description': ('summary', 'note_dtos'),
        'misc': ('comment', 'note_dtos'),
        'aminoacid_rep': ('comment', 'note_dtos'),
        'molecular_info': ('comment', 'note_dtos'),
        'nucleotide_sub': ('comment', 'note_dtos'),
        # At the moment, just for code development. (line below)
        # 'internal_notes': ('internal_note', 'note_dtos'),
    }
    cassette_associations = []
    cassette_cassette_rels = {}

    def get_general_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_feature_lookup(session, feature_types=['cassette', 'construct', 'allele', 'tool', 'gene'])

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
        if self.testing:
            self.log.debug("BOB: print list")
            for bob in self.fb_data_entities:
                self.log.debug(f"BOB: {bob}")

    def add_in_vitro_allele_entries(self, session, reference_set):
        """Extend list of entities."""
        self.log.info('Add entities for alleles having "in vitro construct" annotations.')
        counter = 0
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
                if pkey_id not in self.fb_data_entities:
                    self.fb_data_entities[pkey_id] = self.fb_export_type(result)
                    counter += 1
        if reference_set is True:
            self.log.info(f'Found {counter} FlyBase {self.datatype} in vitro entities in reference chado instance not in main set.')
        else:
            self.log.info(f'Found {counter} FlyBase {self.datatype} in vitro entities in chado not in main set.')

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
            self.log.info(f'TESTING: limit to these entities: {self.test_set.keys()}')
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
            if self.testing:
                self.log.debug(f"BOB: {pkey_id}: {result}")
            if reference_set is True:
                self.fb_reference_entity_ids.append(pkey_id)
            else:
                self.fb_data_entities[pkey_id] = self.fb_export_type(result)
            counter += 1
        if reference_set is True:
            self.log.info(f'Found {counter} FlyBase {self.datatype} main entities in reference chado instance.')
        else:
            self.log.info(f'Found {counter} FlyBase {self.datatype} main entities in chado.')

    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_cassette_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_entity_props_to_notes('cassette_prop_to_note_mapping')
        # self.map_xrefs()
        self.map_secondary_ids('secondary_identifiers')
        self.map_cassette_associations()

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_cassette_basic(self):
        """Map basic FlyBase transgenic cassette data to the Alliance LinkML object."""
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
        self.get_entity_fb_xrefs(session)
        # self.get_entity_xrefs(session)
        # Text form From FTA-134 next 3 lines:-
        # feature_relationship, type 'has_reg_region', subject = FBal, object in (FBto, FBsf, FBgn) # from GA30e field
        # feature_relationship, type 'tagged_with', subject = FBal, object in (FBto, FBsf) # from GA30a field
        # feature_relationship, type 'carries_tool' subject = FBal, object in (FBto, FBsf) # from GA30b field

        self.get_entity_relationships(session, 'subject')
        # , rel_type)='has_reg_region',
        #                             entity_type='engineered_region', entity_regex=self.regex['tool'])

    # Elaborate on query_chado_and_export() for the CassetteHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the CassetteHandler."""
        super().query_chado_and_export(session)
        self.generate_export_dict(self.cassette_associations, 'cassette_association_ingest_set')

        return

    def map_cassette_associations(self):
        """Map transgenic cassette associations to Alliance object."""
        self.log.info('Map cassette associations to Alliance object.')
        OBJECT = 1
        SUBJECT = 0
        counter = 0

        cassette_cassette_counter = {}
        for cassette_cassette_key in self.cassette_cassette_rels.keys():
            self.log.debug(f'Mapping {cassette_cassette_key} to Alliance object. {self.cassette_cassette_rels[cassette_cassette_key]}')
            try:
                cassette_cassette_counter[cassette_cassette_key[OBJECT]] += 1
            except KeyError:
                cassette_cassette_counter[cassette_cassette_key[OBJECT]] = 1

        for cassette_cassette_key, cassette_cassette_rels in self.cassette_cassette_rels.items():
            object_feature_id = cassette_cassette_key[OBJECT]
            f_object = self.fb_data_entities[object_feature_id]
            object_curie = f'FB:{f_object.uniquename}'
            subject = self.feature_lookup[cassette_cassette_key[SUBJECT]]
            subject_curie = f'FB:{subject["uniquename"]}'
            first_feat_rel = cassette_cassette_rels[0]
            all_pub_ids = []
            for cassette_cassette_rel in cassette_cassette_rels:
                all_pub_ids.extend(cassette_cassette_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            # NOTE: pub 383755 | FlyBase Experimental Tool information Is the only one used
            # for cassettes. But not in lookup pub curies!
            # So pub_curies will be empty.
            pub_curies = self.lookup_pub_curies(all_pub_ids)

            # Adjust allele-gene relation_type as needed.
            rel_type_name = 'compatible_tool'
            self.log.debug(f"BOB: cas cass rel: {cassette_cassette_rels}")
            self.log.debug(f"BOB: chado_obj: {cassette_cassette_rels[0].chado_obj}")
            self.log.debug(f"BOB: cassette_cassette_rels[0]: {cassette_cassette_rels[0]}")
            rel_dto = agr_datatypes.CassetteAssociationDTO(
                subject_curie, object_curie,
                pub_curies, False, rel_type_name)
            if f_object.is_obsolete is True or subject['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            first_feat_rel.linkmldto = rel_dto
            self.cassette_associations.append(first_feat_rel)
            counter += 1
        self.log.info(f'Generated {counter} cassette-cassette unique associations.')
        return

    def synthesize_cassette_associations(self):
        """Get cassette relationships."""
        self.log.info('Synthesize cassette.')
        sub_cassette_counter = 0
        obj_cassette_counter = 0
        for cassette in self.fb_data_entities.values():
            relevant_cassette_rels = cassette.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['has_reg_region', 'tagged_with', 'carries_tool'])
            if relevant_cassette_rels:
                sub_cassette_counter += 1
            for cassette_rel in relevant_cassette_rels:
                if self.testing:
                    self.log.debug(f"BOB: cassette_rel:{cassette_rel}")
                try:
                    cassette_cassette_key = (cassette_rel.chado_obj.object_id, cassette_rel.chado_obj.subject_id)
                except AttributeError:
                    self.log.error(f"problem {cassette} {cassette_rel}")
                    raise
                try:
                    self.cassette_cassette_rels[cassette_cassette_key].append(cassette_rel)
                except KeyError:
                    self.cassette_cassette_rels[cassette_cassette_key] = [cassette_rel]
                    obj_cassette_counter += 1
        self.log.info(f'Found {obj_cassette_counter} cassettes for {sub_cassette_counter} cassettes.')
        return

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_cassette_associations()
