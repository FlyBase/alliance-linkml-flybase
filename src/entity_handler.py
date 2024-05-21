"""Module:: entity_handler.

Synopsis:
    A generic data handler for FlyBase first class entities that have a FlyBase
    curie and, typically, a web report; this handler is not designed to handle
    complex annotations (e.g., disease, phenotype).

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import datetime
import re
import strict_rfc3339
from logging import Logger
from sqlalchemy import or_
from sqlalchemy.orm import aliased
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop, FeatureDbxref,
    Featureprop, FeaturepropPub, FeaturePub, FeatureRelationship,
    FeatureRelationshipPub, FeatureSynonym, Organism, OrganismDbxref, Pub,
    PubDbxref, Strain, StrainCvterm, StrainCvtermprop, StrainDbxref, Strainprop,
    StrainpropPub, StrainPub, StrainSynonym, Synonym
)
import agr_datatypes
from handler import DataHandler


class PrimaryEntityHandler(DataHandler):
    """A generic data handler for that gets FlyBase data for first class entities.

    This object handles only primary FlyBase entities, things that have a FB
    curie and, typically, web reports. These entities typically have props and
    CV term associations, and they are the subjects of relationships and
    annotations (e.g., feature_relationship, phenstatement).

    """
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the generic PrimaryEntityHandler object."""
        super().__init__(log, fb_data_type, testing)

    # Mappings of input fb_data_type to various data handling objects and chado tables.
    main_chado_entity_types = {
        'allele': 'feature',
        'construct': 'feature',
        'gene': 'feature',
        'strain': 'strain',
        'variation': 'feature',
    }
    # Mappings of main data types to chado tables with associated data.
    chado_tables = {
        'primary_key': {'feature': 'feature_id', 'strain': 'strain_id'},
        'main_table': {'feature': Feature, 'strain': Strain},
        'pubs': {'feature': FeaturePub, 'strain': StrainPub},
        'synonyms': {'feature': FeatureSynonym, 'strain': StrainSynonym},
        'dbxrefs': {'feature': FeatureDbxref, 'strain': StrainDbxref},
        'props': {'feature': Featureprop, 'strain': Strainprop},
        'prop_pubs': {'feature': FeaturepropPub, 'strain': StrainpropPub},
        'cvterms': {'feature': FeatureCvterm, 'strain': StrainCvterm},
        'cvtermprops': {'feature': FeatureCvtermprop, 'strain': StrainCvtermprop}
    }
    # CVterms used to define a fb_data_type within a larger chado table.
    subtypes = {
        'allele': ['allele'],
        'construct': ['engineered_transposable_element', 'engineered_region', 'transgenic_transposable_element'],
        'gene': ['gene'],
        'variation': ['MNV', 'complex_substitution', 'deletion', 'delins', 'insertion', 'point_mutation', 'sequence_alteration', 'sequence_variant']
    }

    # Add methods to be run by get_general_data() below.

    # Add methods to be run by get_datatype_data() below.
    def get_entities(self, session):
        """Get primary FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        self.log.info(f'Get {self.fb_data_type} data entities from {chado_type} table.')
        chado_table = self.chado_tables['main_table'][chado_type]
        datatype_object = self.datatype_objects[self.fb_data_type]
        filters = ()
        if self.fb_data_type in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.fb_data_type]}')
            filters += (chado_table.uniquename.op('~')(self.regex[self.fb_data_type]), )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (Cvterm.name.in_((self.subtypes[self.fb_data_type])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (chado_table.uniquename.in_((self.test_set.keys())), )
        if filters == ():
            self.log.warning('Have no filters for the main FlyBase entity driver query.')
            raise
        if self.fb_data_type in self.subtypes.keys():
            results = session.query(chado_table).\
                select_from(chado_table).\
                join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
                filter(*filters).\
                distinct()
        else:
            results = session.query(chado_table).\
                filter(*filters).\
                distinct()
        pkey_name = self.chado_tables['primary_key'][chado_type]
        self.log.info(f'Have this primary_key name: {pkey_name}')
        counter = 0
        for result in results:
            pkey_id = getattr(result, pkey_name)
            self.fb_data_entities[pkey_id] = datatype_object(result)
            counter += 1
        self.log.info(f'Found {counter} FlyBase {self.fb_data_type} entities in chado.')
        return

    def get_entity_pubs(self, session):
        """Get pubs directly associated with FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        asso_chado_table = self.chado_tables['pubs'][chado_type]
        self.log.info(f'Get pubs for {self.fb_data_type} data entities from {asso_chado_table}.')
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
        fkey_col = self.get_foreign_key_column(asso_chado_table, main_pkey_name)
        filters = (
            fkey_col.in_((self.fb_data_entities.keys())),
        )
        results = session.query(asso_chado_table).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = getattr(result, main_pkey_name)
            try:
                self.fb_data_entities[entity_pkey_id].pubs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} pubs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} pubs for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_synonyms(self, session):
        """Get synonyms for the FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        asso_chado_table = self.chado_tables['synonyms'][chado_type]
        self.log.info(f'Get synonyms for {self.fb_data_type} data entities from {asso_chado_table}.')
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
        fkey_col = self.get_foreign_key_column(asso_chado_table, main_pkey_name)
        filters = (
            fkey_col.in_((self.fb_data_entities.keys())),
        )
        results = session.query(asso_chado_table).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = getattr(result, main_pkey_name)
            try:
                self.fb_data_entities[entity_pkey_id].synonyms.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} synonyms for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} synonyms for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_fb_xrefs(self, session):
        """Get secondary FB xrefs for the FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        asso_chado_table = self.chado_tables['dbxrefs'][chado_type]
        self.log.info(f'Get non-current FlyBase xrefs for {self.fb_data_type} data entities from {asso_chado_table}.')
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
        fkey_col = self.get_foreign_key_column(asso_chado_table, main_pkey_name)
        filters = (
            fkey_col.in_((self.fb_data_entities.keys())),
            asso_chado_table.is_current.is_(False),
            Db.name == 'FlyBase',
        )
        results = session.query(asso_chado_table).\
            join(Dbxref, (Dbxref.dbxref_id == asso_chado_table.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = getattr(result, main_pkey_name)
            try:
                self.fb_data_entities[entity_pkey_id].fb_sec_dbxrefs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} 2o FB xrefs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} 2o FB xrefs for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_xrefs(self, session):
        """Get all other xrefs for the FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        asso_chado_table = self.chado_tables['dbxrefs'][chado_type]
        self.log.info(f'Get non-FlyBase xrefs for {self.fb_data_type} data entities from {asso_chado_table}.')
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
        fkey_col = self.get_foreign_key_column(asso_chado_table, main_pkey_name)
        filters = (
            fkey_col.in_((self.fb_data_entities.keys())),
            asso_chado_table.is_current.is_(True),
            Db.name.in_((self.fb_agr_db_dict.keys()))
        )
        results = session.query(asso_chado_table).\
            join(Dbxref, (Dbxref.dbxref_id == asso_chado_table.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = getattr(result, main_pkey_name)
            try:
                self.fb_data_entities[entity_pkey_id].dbxrefs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} xrefs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} xrefs for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_props(self, session):
        """Placeholder."""
        return

    def get_entity_cvterms(self, session):
        """Placeholder."""
        return

    def get_entity_associated_data(self, session):
        """Get data associated with primary FlyBase data entities."""
    #     associated_data_types = ['pubs', 'synonyms', 'dbxrefs', 'props', 'cvterms']
    #     chado_type = self.main_chado_entity_types[self.fb_data_type]
    #     self.log.info(f'Get associated data for {self.fb_data_type} data entities from {chado_type}-related chado tables.')
    #     main_pkey_name = self.chado_tables['primary_key'][chado_type]
    #     for i in associated_data_types:
    #         asso_chado_table = self.chado_tables[i][chado_type]
    #         self.log.info(f'Get {i} for {self.fb_data_type} from {asso_chado_table}')
    #         fkey_col = self.get_foreign_key_column(asso_chado_table, main_pkey_name)
    #         filters = (
    #             fkey_col.in_((self.fb_data_entities.keys())),
    #         )
    #         results = session.query(asso_chado_table).\
    #             filter(*filters).\
    #             distinct()
    #         counter = 0
    #         pass_counter = 0
    #         for result in results:
    #             entity_pkey_id = getattr(result, main_pkey_name)
    #             try:
    #                 self.fb_data_entities[entity_pkey_id].__dict__[i].append(result)
    #                 counter += 1
    #             except KeyError:
    #                 pass_counter += 1
    #         self.log.info(f'Found {counter} {i} for {self.fb_data_type} entities.')
    #         self.log.info(f'Ignored {pass_counter} {i} for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_prop_pubs(self, session):
        """Get prop pubs for primary FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        # pkey_name = self.chado_tables['primary_key'][chado_type]
        main_chado_table = aliased(self.chado_tables['main_table'][chado_type], name='main_chado_table')
        prop_chado_table = aliased(self.chado_tables['props'][chado_type], name='prop_chado_table')
        prop_pub_chado_table = aliased(self.chado_tables['prop_pubs'][chado_type], name='prop_pub_chado_table')
        self.log.info(f'Get prop pubs for {self.fb_data_type} data entities from {prop_pub_chado_table} chado table.')
        # fkey_col = self.get_foreign_key_column(prop_chado_table, pkey_name)
        main_pkey_col = self.get_primary_key_column(main_chado_table)
        filters = (
            main_pkey_col.in_((self.fb_data_entities.keys())),
        )
        # To do - need a way to add a uniquename filter for testing.
        if self.testing:
            pass
        results = session.query(prop_chado_table, prop_pub_chado_table).\
            select_from(prop_pub_chado_table).\
            join(prop_chado_table).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = getattr(result.prop_chado_table, main_pkey_col.name)
            try:
                self.fb_data_entities[entity_pkey_id].prop_pubs.append(result.prop_pub_chado_table)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} props for {self.fb_data_type}.')
        self.log.info(f'Ignored {pass_counter} props for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_cvtermprops(self, session):
        """Get CV term props for primary FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
        cvterm_table = aliased(self.chado_tables['cvterms'][chado_type], name='cvterm_table')
        cvtermprop_table = aliased(self.chado_tables['cvtermprops'][chado_type], name='cvtermprop_table')
        self.log.info(f'Get CV term props for {self.fb_data_type} data entities from {cvtermprop_table} chado table.')
        fkey_col = self.get_foreign_key_column(cvterm_table, main_pkey_name)
        filters = (
            fkey_col.in_((self.fb_data_entities.keys())),
        )
        # BOB - need a way to add a uniquename filter for testing.
        if self.testing:
            pass
        results = session.query(cvterm_table, cvtermprop_table).\
            select_from(cvterm_table).\
            join(cvtermprop_table).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = getattr(result.cvterm_table, fkey_col.name)
            try:
                self.fb_data_entities[entity_pkey_id].cvtermprops.append(result.cvtermprop_table)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} CV term props for {self.fb_data_type}.')
        self.log.info(f'Ignored {pass_counter} CV term props for irrelevant {self.fb_data_type} entities.')
        return

    def get_entity_timestamps(self, session):
        """Get timestamps for data entities."""
        self.log.info(f'Get timestamps for FlyBase {self.fb_data_type} entities.')
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        entity_table_counter = 0
        audit_chado_counter = 0
        # Get distinct timestamps for each entity (do not distinguish by action, etc).
        for i in self.fb_data_entities.values():
            try:
                i.timestamps.append(i.timeaccessioned)
                i.timestamps.append(i.timelastmodified)
                entity_table_counter += 1
            except AttributeError:
                audit_query = f"""
                SELECT DISTINCT record_pkey, transaction_timestamp
                FROM audit_chado
                WHERE audited_table = '{chado_type}'
                AND record_pkey = {i.db_primary_id};
                """
                TIMESTAMP = 1
                audit_results = session.execute(audit_query).fetchall()
                for row in audit_results:
                    i.timestamps.append(row[TIMESTAMP])
                audit_chado_counter += 1
        self.log.info(f'Obtained timestamps for {entity_table_counter} entities directly from the {chado_type} table.')
        self.log.info(f'Obtained timestamps for {audit_chado_counter} entities directly from the audit_chado table.')
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_ncbi_taxon_id(self):
        """Determine the NCBITaxon ID for FB entities."""
        self.log.info('Determine the NCBITaxon ID for FB entities.')
        for fb_data_entity in self.fb_data_entities.values():
            # Catch cases where the FB data entity has no organism_id.
            try:
                organism_id = fb_data_entity.chado_obj.organism_id
            except AttributeError:
                self.log.warning(f'No organism_id for {fb_data_entity}.')
                return
            # Catch cases where the FB data entity has no corresponding NCBITaxon ID.
            try:
                ncbi_taxon_id = self.ncbi_taxon_lookup[organism_id]
            except KeyError:
                self.log.warning(f'Use "unidentified" NCBITaxon ID for {fb_data_entity}')
                ncbi_taxon_id = 'NCBITaxon:32644'
            fb_data_entity.ncbi_taxon_id = ncbi_taxon_id
        return

    def synthesize_secondary_ids(self):
        """Process secondary IDs and into a list of old FB uniquename curies."""
        self.log.info('Process secondary IDs and return a list of old FB uniquenames.')
        for fb_data_entity in self.fb_data_entities.values():
            secondary_ids = []
            for xref in fb_data_entity.fb_sec_dbxrefs:
                secondary_ids.append(f'FB:{xref.dbxref.accession}')
            fb_data_entity.alt_fb_ids = list(set(secondary_ids))
        return

    def synthesize_synonyms(self):
        """Synthesize synonyms from Synonym association objects."""
        self.log.info('Synthesize synonyms.')
        # Dict for converting FB synonym types to AGR synonym types.
        synonym_type_conversion = {
            'symbol': 'nomenclature_symbol',
            'fullname': 'full_name',
            'nickname': 'unspecified',
            'synonym': 'unspecified',
        }
        for fb_data_entity in self.fb_data_entities.values():
            # For each entity, gather synonym_id-keyed dict of synonym info.
            for feat_syno in fb_data_entity.synonyms:
                try:
                    fb_data_entity.synonym_dict[feat_syno.synonym_id]['is_current'].append(feat_syno.is_current)
                    fb_data_entity.synonym_dict[feat_syno.synonym_id]['is_internal'].append(feat_syno.is_internal)
                    fb_data_entity.synonym_dict[feat_syno.synonym_id]['pub_ids'].append(feat_syno.pub_id)
                except KeyError:
                    syno_dict = {
                        'name_type_name': synonym_type_conversion[feat_syno.synonym.type.name],
                        'format_text': feat_syno.synonym.name,
                        'display_text': sub_sup_sgml_to_html(feat_syno.synonym.synonym_sgml),
                        'is_current': [feat_syno.is_current],
                        'is_internal': [feat_syno.is_internal],
                        'pub_ids': [feat_syno.pub_id],
                        'pub_curies': []
                    }
                    fb_data_entity.synonym_dict[feat_syno.synonym_id] = syno_dict
            # Go back over each synonym and refine each
            for syno_dict in fb_data_entity.synonym_dict.values():
                # Then modify attributes as needed.
                # Identify systematic names.
                if re.match(self.regex['systematic_name'], syno_dict['format_text']) and syno_dict['name_type_name'] == 'nomenclature_symbol':
                    syno_dict['name_type_name'] = 'systematic_name'
                # Classify is_current (convert list of booleans into a single boolean).
                if True in syno_dict['is_current']:
                    syno_dict['is_current'] = True
                else:
                    syno_dict['is_current'] = False
                # Classify is_internal (convert list of booleans into a single boolean).
                if False in syno_dict['is_internal']:
                    syno_dict['is_internal'] = False
                else:
                    syno_dict['is_internal'] = True
                # Convert pub_ids into pub_curies.
                syno_dict['pub_curies'] = self.lookup_pub_curies(syno_dict['pub_ids'])
                # Finally, pick out current symbol for the entity.
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] in ['systematic_name', 'nomenclature_symbol']:
                    fb_data_entity.curr_fb_symbol = syno_dict['display_text']
        return

    def synthesize_props(self):
        """Process props and pubs for FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        prop_chado_table = self.chado_tables['props'][chado_type]
        prop_pkey_col = self.get_primary_key_column(prop_chado_table)
        self.log.info(f'Synthesize {prop_chado_table} data.')
        prop_types = []
        # Build prop dict.
        for fb_data_entity in self.fb_data_entities.values():
            # Build prop_type-keyed lists of props.
            for prop in fb_data_entity.props:
                prop_type = prop.type.name
                prop_types.append(prop_type)
                try:
                    fb_data_entity.prop_dict[prop_type].append(prop)
                except KeyError:
                    fb_data_entity.prop_dict[prop_type] = [prop]
            # Build prop_id-keyed lists of pub_ids.
            for prop_pub in fb_data_entity.prop_pubs:
                prop_id = getattr(prop_pub, prop_pkey_col.name)
                try:
                    fb_data_entity.prop_pub_dict[prop_id].append(prop_pub.pub_id)
                except KeyError:
                    fb_data_entity.prop_pub_dict[prop_id] = [prop_pub.pub_id]
        # Review prop synthesis.
        self.log.info(f'Found these types of {chado_type}props: {set(prop_types)}')
        # Optional debug.
        # for i in self.fb_data_entities.values():
        #     for k, v in i.prop_pub_dict.items():
        #         self.log.debug(f'{chado_type}_id={i.db_primary_id}, {chado_type}prop_id={k}, pub_ids={v}')
        return

    def synthesize_pubs(self):
        """Collect pub_ids associated directly or indirectly with the entity."""
        self.log.info('Collect pub_ids associated directly or indirectly with the entity.')
        pub_sources = ['pubs', 'synonyms', 'cvterms', 'prop_pubs']
        for fb_data_entity in self.fb_data_entities.values():
            for pub_source in pub_sources:
                fb_data_entity.all_pub_ids.extend([i.pub_id for i in getattr(fb_data_entity, pub_source)])
            fb_data_entity.all_pub_ids = list(set(fb_data_entity.all_pub_ids))
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_data_provider_dto(self):
        """Return the DataProviderDTO for the FB data entity."""
        # Note - this method is depends on previous determination of fb_data_entity.curr_fb_symbol by map_synonyms(), if applicable.
        self.log.info('Map data provider to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.curr_fb_symbol:
                display_name = fb_data_entity.curr_fb_symbol
            else:
                display_name = fb_data_entity.name
            dp_xref = agr_datatypes.CrossReferenceDTO('FB', f'FB:{fb_data_entity.uniquename}', self.fb_data_type, display_name).dict_export()
            fb_data_entity.linkmldto.data_provider_dto = agr_datatypes.DataProviderDTO(dp_xref).dict_export()
        return

    def map_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for a FlyBase entity."""
        self.log.info('Map secondary IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            secondary_id_dtos = []
            for secondary_id in fb_data_entity.alt_fb_ids:
                sec_dto = agr_datatypes.SecondaryIdSlotAnnotationDTO(secondary_id, []).dict_export()
                secondary_id_dtos.append(sec_dto)
            sec_id_list = getattr(fb_data_entity.linkmldto, slot_name)
            sec_id_list.extend(secondary_id_dtos)
        return

    def map_pubs(self):
        """Add pub curies to a FlyBase entity."""
        self.log.info('Map pubs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            for pub_id in fb_data_entity.all_pub_ids:
                try:
                    fb_data_entity.linkmldto.reference_curies.append(self.bibliography[pub_id])
                except KeyError:
                    pass
            try:
                fb_data_entity.linkmldto.reference_curies.remove('FB:unattributed')
            except ValueError:
                pass
        return

    def map_xrefs(self):
        """Add a list of Alliance CrossReferenceDTO dicts to a FlyBase entity."""
        self.log.info('Map xrefs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            cross_reference_dtos = []
            for xref in fb_data_entity.dbxrefs:
                # Build Alliance xref DTO
                prefix = self.fb_agr_db_dict[xref.dbxref.db.name]
                # The page_area assignment assumes that the self.fb_data_type has a matching value in the Alliance resourceDescriptors.yaml page.
                page_area = self.fb_data_type
                # Clean up cases where the db prefix is redundantly included at the start of the dbxref.accession.
                redundant_prefix = f'{prefix}:'
                if xref.dbxref.accession.startswith(redundant_prefix):
                    cleaned_accession = xref.dbxref.accession.replace(redundant_prefix, '', 1)
                    self.log.debug(f'Removed "{redundant_prefix}" from "{xref.dbxref.accession}" to give "{cleaned_accession}"')
                else:
                    cleaned_accession = xref.dbxref.accession
                curie = f'{prefix}:{cleaned_accession}'
                display_name = curie
                xref_dto = agr_datatypes.CrossReferenceDTO(prefix, curie, page_area, display_name).dict_export()
                cross_reference_dtos.append(xref_dto)
            fb_data_entity.linkmldto.cross_reference_dtos = cross_reference_dtos
        return

    def map_synonyms(self):
        """Generate name/synonym DTOs for an entity."""
        self.log.info('Map synonyms to Alliance object.')
        # First determine synonym slots available, if any.
        linkml_synonym_slots = {
            'symbol_bin': '_symbol_dto',
            'full_name_bin': '_full_name_dto',
            'systematic_name_bin': '_systematic_name_dto',
            'synonym_bin': '_synonym_dtos'
        }
        map_synonyms_required = False
        test_linkmldto = self.agr_linkmldto_dict[self.fb_data_type]
        for dto_key in test_linkmldto.__dict__.keys():
            for bin_type, bin_suffix in linkml_synonym_slots.items():
                if dto_key.endswith(bin_suffix):
                    linkml_synonym_slots[bin_type] = dto_key
                    self.log.debug(f'Map {bin_type} to LinkML DTO slot {dto_key} because it has suffix "{bin_suffix}".')
                    map_synonyms_required = True
        if map_synonyms_required is False:
            self.log.warning(f'The map_synonyms() method has been incorrectly called for {self.fb_data_type} objects.')
            return
        else:
            self.log.info(f'Have these linkml name dto slots to fill in: {linkml_synonym_slots.values()}')
        for fb_data_entity in self.fb_data_entities.values():
            linkml_synonym_bins = {
                'symbol_bin': [],
                'full_name_bin': [],
                'systematic_name_bin': [],
                'synonym_bin': []
            }
            # Create NameSlotAnnotationDTO objects and sort them out.
            for syno_dict in fb_data_entity.synonym_dict.values():
                # Sort into current symbol, current fullname or synonym.
                name_dto = agr_datatypes.NameSlotAnnotationDTO(syno_dict['name_type_name'], syno_dict['format_text'],
                                                               syno_dict['display_text'], syno_dict['pub_curies']).dict_export()
                name_dto['internal'] = syno_dict['is_internal']
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] in ['nomenclature_symbol', 'systematic_name']:
                    linkml_synonym_bins['symbol_bin'].append(name_dto)
                elif syno_dict['is_current'] is True and syno_dict['name_type_name'] == 'fullname':
                    linkml_synonym_bins['full_name_bin'].append(name_dto)
                else:
                    linkml_synonym_bins['synonym_bin'].append(name_dto)
                # Also add to current systematic name for current Dmel genes only.
                if syno_dict['name_type_name'] == 'systematic_name' and syno_dict['display_text'] == fb_data_entity.curr_anno_id:
                    if fb_data_entity.chado_obj.is_obsolete is False and fb_data_entity.org_abbr == 'Dmel':
                        linkml_synonym_bins['systematic_name_bin'].append(name_dto)
            # Review the linkml_synonym_bins for each fb_data_entity.
            # 1. Symbol.
            if len(linkml_synonym_bins['symbol_bin']) == 0:
                self.log.warning(f'No current symbols found for {fb_data_entity}: create a generic one.')
                generic_symbol_dto = agr_datatypes.NameSlotAnnotationDTO('nomenclature_symbol', fb_data_entity.name, fb_data_entity.name, []).dict_export()
                setattr(fb_data_entity.linkmldto, linkml_synonym_slots['symbol_bin'], generic_symbol_dto)
            else:
                setattr(fb_data_entity.linkmldto, linkml_synonym_slots['symbol_bin'], linkml_synonym_bins['symbol_bin'][0])
                if len(linkml_synonym_bins['symbol_bin']) > 1:
                    multi_symbols = ', '.join([i['format_text'] for i in linkml_synonym_bins['symbol_bin']])
                    self.log.warning(f'Found many current symbols for {fb_data_entity}: {multi_symbols}')
            # 2. Fullname.
            if len(linkml_synonym_bins['full_name_bin']) == 1:
                setattr(fb_data_entity.linkmldto, linkml_synonym_slots['full_name_bin'], linkml_synonym_bins['full_name_bin'][0])
            elif len(linkml_synonym_bins['full_name_bin']) > 1:
                multi_symbols = ', '.join([i['format_text'] for i in linkml_synonym_bins['full_name_bin']])
                self.log.warning(f'Found many current full_names for {fb_data_entity}: {multi_symbols}')
            # 3. Systematic name.
            if len(linkml_synonym_bins['systematic_name_bin']) == 0 and fb_data_entity.curr_anno_id and fb_data_entity.chado_obj.is_obsolete is False:
                self.log.warning(f'No current systematic names found for current annotated {fb_data_entity}: create a generic one.')
                sys_name_dto = agr_datatypes.NameSlotAnnotationDTO('systematic_name', fb_data_entity.curr_anno_id,
                                                                   fb_data_entity.curr_anno_id, []).dict_export()
                setattr(fb_data_entity.linkmldto, linkml_synonym_slots['systematic_name_bin'], sys_name_dto)
            elif len(linkml_synonym_bins['systematic_name_bin']) == 1:
                setattr(fb_data_entity.linkmldto, linkml_synonym_slots['systematic_name_bin'], linkml_synonym_bins['systematic_name_bin'][0])
            elif len(linkml_synonym_bins['systematic_name_bin']) > 1:
                multi_symbols = ', '.join([i['format_text'] for i in linkml_synonym_bins['systematic_name_bin']])
                self.log.warning(f'Found many current systematic_names for {fb_data_entity}: {multi_symbols}')
            # 4. Synonyms.
            setattr(fb_data_entity.linkmldto, linkml_synonym_slots['synonym_bin'], linkml_synonym_bins['synonym_bin'])
        return

    def map_timestamps(self):
        """Map timestamps to Alliance object."""
        self.log.info('Map timestamps to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.timestamps:
                fb_data_entity.linkmldto.date_created = strict_rfc3339.\
                    timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(fb_data_entity.timestamps)))
                fb_data_entity.linkmldto.date_updated_curie = strict_rfc3339.\
                    timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(fb_data_entity.timestamps)))
        return

    def flag_internal_fb_entities(self, input_list_name: str):
        """Flag obsolete FB objects in some list as internal.

        Args:
            input_list_name (str): The name of a handler list/dict with objects with LinkMLDTO objects under the linkmldto attribute.

        """
        self.log.info(f'Flag obsolete FB objects in {input_list_name}.')
        input_data = getattr(self, input_list_name)
        if type(input_data) is dict:
            input_list = list(input_data.values())
        elif type(input_data) is list:
            input_list = input_data
        for fb_data_entity in input_list:
            try:
                if fb_data_entity.linkmldto.obsolete is True:
                    fb_data_entity.linkmldto.internal = True
                    fb_data_entity.internal_reasons.append('Obsolete')
            except AttributeError:
                self.log.error('LinkMLDTO entity lacks obsolete attribute.')
        return
