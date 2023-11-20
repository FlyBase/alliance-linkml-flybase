"""Module:: utils.

Synopsis:
    Data retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import json
import datetime
import strict_rfc3339
from logging import Logger
from sqlalchemy.orm import aliased, Session
# from sqlalchemy.inspection import inspect
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, OrganismDbxref, Pub, PubDbxref,
    Strain, StrainPub, StrainSynonym, StrainDbxref, Strainprop, StrainpropPub, StrainCvterm, StrainCvtermprop,
    Feature, FeaturePub, FeatureSynonym, FeatureDbxref, Featureprop, FeaturepropPub, FeatureCvterm, FeatureCvtermprop
)
import datatypes


# Classes
class DataHandler(object):
    """A generic data handler that gets FlyBase data and maps it to a single Alliance LinkML model.

    Specific classes of DataHandler will only map a given FB data type to a single Alliance ingest
    set (i.e., a set of LinkML DTO objects, in JSON format). In some cases, multiple handlers will
    map different FB data to the same ingest set: e.g., FB strains and genotypes both go into the
    "agm_ingest_set" via distinct DataHandler classes.

    """
    def __init__(self, log: Logger, fb_data_type: str):
        """Create the generic DataHandler object.

        Args:
            log (Logger): The global Logger object in the script using the DataHandler.
            fb_data_type (str): The FlyBase data class being handled.

        """
        self.log = log
        self.fb_data_type = fb_data_type
        self.agr_ingest_type = self.agr_ingest_type_dict[fb_data_type]
        # Datatype bins.
        self.fb_data_entities = {}    # db_primary_id-keyed dict of chado objects to export.
        self.export_data = []         # List of data objects for export (as Alliance ingest set).
        # General data bins.
        self.bibliography = {}        # A pub_id-keyed dict of pub curies (PMID or FBrf).
        self.cvterm_dict = {}         # A cvterm_id-keyed dict of Cvterm objects.
        self.ncbi_taxon_dict = {}     # An organism_id-keyed dict of NCBITaxon Dbxref.accession strings.
        # Trackers.
        self.input_count = 0          # Count of entities found in FlyBase chado database.
        self.export_count = 0         # Count of exported Alliance entities.
        self.internal_count = 0       # Count of exported entities marked as internal.
        self.warnings = []            # Handler issues of note.
        self.errors = []              # Handler issues that break things.

    # Correspondence of FB data type to Alliance data transfer ingest set.
    agr_ingest_type_dict = {
        'gene': 'gene_ingest_set',
        'allele': 'allele_ingest_set',
        'construct': 'construct_ingest_set',
        'variation': 'variation_ingest_set',
        'strain': 'agm_ingest_set',
        'genotype': 'agm_ingest_set',
        'disease': 'disease_allele_ingest_set'
    }
    # Mappings of fb_data_type to a datatype Class that will be used to represent each FB entity.
    datatype_objects = {
        # 'gene': datatypes.FBGene,
        # 'allele': datatypes.FBAllele,
        # 'construct': datatypes.FBConstruct,
        # 'variation': datatypes.FBVariant,
        'strain': datatypes.FBStrain,
        # 'disease': datatypes.FBDiseaseAlleleAnnotation
    }
    # Export directions (must be filled out in detail for each specific data handler).
    # The fields in the two lists below must be present in the datatype object
    #     specified in DataHandler.datatype_objects.values().
    required_fields = []
    output_fields = []
    # A filter for xrefs to export, with dict keys as FB db.names and dict values as Alliance db names.
    # Alliance db names should correspond to the contents of this file:
    # https://github.com/alliance-genome/agr_schemas/blob/master/resourceDescriptors.yaml
    fb_agr_db_dict = {'FlyBase': 'FB'}

    # Useful regexes.
    regex = {
        'pub': r'^(FBrf[0-9]{7}|unattributed)$',
        'feature': r'^FB[a-z]{2}[0-9]{7,10}$',
        'gene': r'^FBgn[0-9]{7}$',
        'allele': r'^FBal[0-9]{7}$',
        'construct': r'^FBtp[0-9]{7}$',
        'insertion': r'^FBti[0-9]{7}$',
        'consins': r'^FB(tp|ti)[0-9]{7}$',
        'aberration': r'^FBab[0-9]{7}$',
        'balancer': r'^FBba[0-9]{7}$',
        'abbal': r'^FB(ab|ba)[0-9]{7}$',
        'seqfeat': r'^FBsf[0-9]{10}$',
        'chem': r'^FBch[0-9]{7}$',
        'genotype': r'^FBgo[0-9]{7}$',
        'strain': r'^FBsn[0-9]{7}$',
        'library': r'^FBlc[0-9]{7}$',
        'cell': r'^FBcl[0-9]{7}$'
    }

    # Methods
    def __str__(self):
        """Print out data handler description."""
        handler_description = f'A data handler that exports FB {self.fb_data_type} to Alliance LinkML {self.agr_ingest_type}.'
        return handler_description

    def build_bibliography(self, session):
        """Build FlyBase bibliography."""
        self.log.info('Build FlyBase bibliography.')
        # First get all current pubs having an FBrf uniquename.
        filters = (
            Pub.uniquename.op('~')(self.regex['pub']),
            Pub.is_obsolete.is_(False)
        )
        results = session.query(Pub).\
            filter(*filters).\
            distinct()
        pub_counter = 0
        for pub in results:
            self.bibliography[pub.pub_id] = f'FB:{pub.uniquename}'
            pub_counter += 1
        # Next find PMIDs if available and replace the curie in the bibliography.
        filters = (
            Pub.uniquename.op('~')(self.regex['pub']),
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
            self.bibliography[xref.Pub.pub_id] = f'PMID:{xref.Dbxref.accession}'
            pmid_counter += 1
        self.log.info(f'Found {pmid_counter} PMID IDs for {pub_counter} current FB publications.')
        return

    def get_cvterms(self, session):
        """Create a cvterm_id-keyed dict of Cvterms."""
        self.log.info('Create a cvterm_id-keyed dict of Cvterms.')
        # First get all current pubs having an FBrf uniquename.
        filters = (
            Cvterm.is_obsolete == 0,
        )
        results = session.query(Cvterm).filter(*filters).distinct()
        cvterm_counter = 0
        for cvterm in results:
            self.cvterm_dict[cvterm.cvterm_id] = cvterm
            cvterm_counter += 1
        self.log.info(f'Found {cvterm_counter} current CV terms in chado.')
        return

    def get_ncbi_taxon_ids(self, session):
        """Get FlyBase Organism-NCBITaxon ID correspondence."""
        self.log.info('Get FlyBase Organism-NCBITaxon ID correspondence.')
        filters = (
            OrganismDbxref.is_current.is_(True),
            Db.name == 'NCBITaxon'
        )
        results = session.query(OrganismDbxref, Dbxref).\
            join(Dbxref, (Dbxref.dbxref_id == OrganismDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for xref in results:
            self.ncbi_taxon_dict[xref.OrganismDbxref.organism_id] = xref.Dbxref.accession
            counter += 1
        self.log.info(f'Found {counter} NCBITaxon IDs for FlyBase organisms.')
        return

    def get_general_data(self, session):
        """Get general FlyBase chado data."""
        self.log.info('Get general FlyBase data from chado.')
        self.build_bibliography(session)
        self.get_cvterms(session)
        self.get_ncbi_taxon_ids(session)
        return

    def get_datatype_data(self, session):
        """Get datatype-specific FlyBase data from chado."""
        self.log.info(f'Get FlyBase {self.fb_data_type} data from chado.')
        return

    def synthesize_info(self):
        """Synthesize FB info for each data object."""
        self.log.info(f'Synthesize FlyBase "{self.fb_data_type}" data.')
        return

    def map_fb_data_to_alliance(self):
        """Map FB data to the Alliance LinkML object."""
        self.log.info(f'Map FlyBase "{self.fb_data_type}" data to the Alliance LinkML object for the "{self.agr_ingest_type}".')
        return

    def flag_unexportable_entities(self):
        """Flag entities lacking information for a required field."""
        self.log.info(f'Flag FlyBase "{self.fb_data_type}" data lacking information for a required field.')
        for i in self.fb_data_entities.values():
            for attr in self.required_fields:
                if attr not in i.linkmldto.__dict__.keys():
                    i.for_export = False
                    i.export_warnings.append(f'Missing "{attr}" attribute.')
                elif getattr(i.linkmldto, attr) is None:
                    i.for_export = False
                    i.export_warnings.append(f'Missing value "{attr}" attribute.')
            if i.for_export is False:
                self.log.debug(f'DO NOT EXPORT {i}: {i.export_warnings}')
        return

    def generate_export_dict(self):
        """Generate LinkML export dict from the primary FB data entities."""
        for i in self.fb_data_entities.values():
            self.input_count += 1
            if i.for_export is False:
                self.log.debug(f'Suppress {i} from export: {i.export_warnings}')
                continue
            self.export_count += 1
            if i.linkmldto.internal is True:
                self.internal_count += 1
                self.log.debug(f'Export {i} but keep internal at the Alliance: {i.internal_reasons}')
            export_agr_dict = {}
            for attr in self.output_fields:
                if getattr(i.linkmldto, attr) is not None and getattr(i.linkmldto, attr) != []:
                    export_agr_dict[attr] = getattr(i.linkmldto, attr)
            self.export_data.append(export_agr_dict)
        public_count = self.export_count - self.internal_count
        self.log.info(f'Exported {self.export_count} of {self.input_count} {self.fb_data_type} entities.')
        self.log.info(f'{public_count} of {self.export_count} exported {self.fb_data_type}s are PUBLIC.')
        self.log.info(f'{self.internal_count} of {self.export_count} exported {self.fb_data_type}s are INTERNAL.')
        return

    def query_chado(self, session):
        """Wrapper that runs all methods within an SQLAlchemy session."""
        self.get_general_data(session)
        self.get_datatype_data(session)
        self.synthesize_info()
        self.map_fb_data_to_alliance()
        self.flag_unexportable_entities()
        self.generate_export_dict()
        return


class PrimaryEntityHandler(DataHandler):
    """A generic handler for that gets data for FlyBase entities and maps it to the Alliance LinkML model.

    This object handles only primary FlyBase entities, things like genes, strains, etc. that typically have
    public curies and web reports, attributes like props and cvterm associations, and are the subjects
    in relationships (e.g., feature_relationship) and annotations (e.g., phenstatement). Separate handler
    classes are to be used for the export of relationships and complex annotations.

    """
    def __init__(self, log: Logger, fb_data_type: str):
        """Create the generic PrimaryEntityHandler object."""
        super().__init__(log, fb_data_type)

    # Mappings of input fb_data_type to various data handling objects and chado tables.
    main_chado_entity_types = {
        'gene': 'feature',
        'allele': 'feature',
        'construct': 'feature',
        'variation': 'feature',
        'strain': 'strain'
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
        'gene': ['gene'],
        'allele': ['allele'],
        'variation': ['MNV', 'complex_substitution', 'deletion', 'delins', 'insertion', 'point_mutation', 'sequence_alteration', 'sequence_variant']
    }

    # BOB: Quick testing of general SQLAlchemy query approaches.
    def sqlalchemy_test(self, session):
        """Test SQLAlchemy behavior."""
        self.log.info('Test SQLAlchemy behavior.')
        lbe_types = ['gene']
        filters = (
            Feature.uniquename == 'FBgn0011278',
        )
        filters += (
            Cvterm.name.in_((lbe_types)),
        )
        results = session.query(Feature).\
            select_from(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for i in results:
            self.log.info(f'Found this feature: {i.name} ({i.uniquename})')
            counter += 1
        self.log.info(f'Found {counter} test results using natural join.')
        return

    # Elaborate on get_datatype_data() sub-methods for PrimaryEntityHandler.
    def get_entity_data(self, session):
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
        if len(filters) == 0:
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
        for result in results:
            pkey_id = getattr(result, pkey_name)
            self.fb_data_entities[pkey_id] = datatype_object(result)
        self.log.info(f'Found {self.input_count} FlyBase {self.fb_data_type} entities in chado.')
        return

    def get_entity_associated_data(self, session):
        """Get data associated with primary FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        self.log.info(f'Get associated data for {self.fb_data_type} data entities from {chado_type}-related chado tables.')
        main_chado_table = self.chado_tables['main_table'][chado_type]
        pkey_name = self.chado_tables['primary_key'][chado_type]
        self.log.info(f'Use this primary key name: {pkey_name}')
        # associated_data_types = ['pubs', 'synonyms', 'dbxrefs', 'props', 'prop_pubs', 'cvterms', 'cvtermprops']
        associated_data_types = ['pubs', 'synonyms', 'dbxrefs']
        # BOB - keep filter as is, or, query one-at-a-time, setting foreign_key to allele's feature_id.
        for i in associated_data_types:
            self.log.info(f'Get {i} for {self.fb_data_type}')
            asso_chado_table = self.chado_tables[i][chado_type]
            # Get the foreign key in associated table corresponding to primary data type.
            foreign_key_column = next((column for column in asso_chado_table.__table__.c if column.foreign_keys and column.name == pkey_name), None)
            if foreign_key_column is None:
                self.log.error(f'Could not get foreign key column for {pkey_name} from {chado_type} {i} tables.')
                raise
            filters = (
                foreign_key_column.in_((self.fb_data_entities.keys())),
            )
            counter = 0
            pass_counter = 0
            results = session.query(asso_chado_table).\
                select_from(asso_chado_table).\
                join(main_chado_table).\
                filter(*filters).\
                distinct()
            for result in results:
                pkey_id = getattr(result, pkey_name)
                try:
                    self.fb_data_entities[pkey_id].__dict__[i].append(result)
                    counter += 1
                except KeyError:
                    pass_counter += 1
            self.log.info(f'Found {counter} {i} for {self.fb_data_type}.')
            self.log.info(f'Ignored {pass_counter} {i} for {self.fb_data_type}.')
        return

    def get_entity_props(self, session):
        """Get props for primary FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        main_chado_table = self.chado_tables['main_table'][chado_type]
        prop_chado_table = self.chado_tables['props'][chado_type]
        self.log.info(f'Get props for {self.fb_data_type} data entities from {chado_type}prop chado table.')
        pkey_name = self.chado_tables['primary_key'][chado_type]
        self.log.info(f'Use this primary key name: {pkey_name}')
        # Get the foreign key in associated table corresponding to primary data type.
        foreign_key_column = next((column for column in prop_chado_table.__table__.c if column.foreign_keys and column.name == pkey_name), None)
        filters = (
            foreign_key_column.in_((self.fb_data_entities.keys())),
        )
        counter = 0
        pass_counter = 0
        results = session.query(prop_chado_table).\
            select_from(prop_chado_table).\
            join(main_chado_table).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        results = session.query(prop_chado_table).\
            select_from(prop_chado_table).\
            join(main_chado_table).\
            filter(*filters).\
            distinct()
        for result in results:
            pkey_id = getattr(result, pkey_name)
            if pkey_id not in self.fb_data_entities.keys():
                pass_counter += 1
                continue
            try:
                self.fb_data_entities[pkey_id].props[result.type.name].append(result)
                counter += 1
            except KeyError:
                self.fb_data_entities[pkey_id].props[result.type.name] = [result]
                counter += 1
        self.log.info(f'Found {counter} props for {self.fb_data_type}.')
        self.log.info(f'Ignored {pass_counter} props for {self.fb_data_type}.')
        return

    def get_entity_prop_pubs(self, session):
        """Get prop pubs for FlyBase data entities."""
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        prop_chado_table = aliased(self.chado_tables['props'][chado_type], name='prop_chado_table')
        prop_pub_chado_table = aliased(self.chado_tables['prop_pubs'][chado_type], name='prop_pub_chado_table')
        self.log.info(f'Get prop pubs for {self.fb_data_type} data entities from {chado_type}prop and {chado_type}prop_pub chado tables.')
        pkey_name = self.chado_tables['primary_key'][chado_type]
        self.log.info(f'Use this primary key name: {pkey_name}')
        # Get the foreign key in associated table corresponding to primary data type.
        foreign_key_column = next((column for column in prop_chado_table.__table__.c if column.foreign_keys and column.name == pkey_name), None)
        filters = (
            foreign_key_column.in_((self.fb_data_entities.keys())),
        )
        counter = 0
        pass_counter = 0
        results = session.query(prop_chado_table, prop_pub_chado_table).\
            select_from(prop_chado_table).\
            join(prop_pub_chado_table).\
            filter(*filters).\
            distinct()
        for result in results:
            pkey_id = getattr(result.prop_chado_table, pkey_name)
            if pkey_id not in self.fb_data_entities.keys():
                pass_counter += 1
                continue
            prop_column_name = f'{chado_type}prop_id'
            prop_id = getattr(result.prop_chado_table, prop_column_name)
            pub_id = getattr(result.prop_pub_chado_table, 'pub_id')
            try:
                self.fb_data_entities[pkey_id].prop_pubs[prop_id].append(pub_id)
                counter += 1
            except KeyError:
                self.fb_data_entities[pkey_id].prop_pubs[prop_id] = [pub_id]
                counter += 1
        self.log.info(f'Found {counter} props for {self.fb_data_type}.')
        self.log.info(f'Ignored {pass_counter} props for {self.fb_data_type}.')
        return

    def get_entity_timestamps(self, session):
        """Get timestamps for data entities."""
        self.log.info(f'Get timestamps for FlyBase {self.fb_data_type} entities.')
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        counter = 0
        # Get distinct timestamps for each entity (do not distinguish by action, etc).
        for i in self.fb_data_entities.values():
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
                counter += 1
        self.log.info(f'Found {counter} timestamps.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the PrimaryEntityHandler."""
        super().get_datatype_data(session)
        self.sqlalchemy_test(session)
        self.get_entity_data(session)
        self.get_entity_associated_data(session)
        self.get_entity_props(session)
        self.get_entity_prop_pubs(session)
        # self.get_cvterms(session)
        self.get_entity_timestamps(session)
        return

    # Elaborate on synthesize_info() sub-methods for PrimaryEntityHandler.
    def synthesize_ncbi_taxon_id(self, fb_data_entity):
        """Determine the NCBITaxon ID for the entity."""
        # Catch cases where the FB data entity has no organism_id.
        try:
            organism_id = fb_data_entity.chado_obj.organism_id
        except AttributeError:
            self.log.warning(f'No organism_id for {fb_data_entity}.')
            return None
        # Catch cases where the FB data entity has no corresponding NCBITaxon ID.
        try:
            ncbi_taxon_id = f'NCBITaxon:{self.ncbi_taxon_dict[organism_id]}'
        except KeyError:
            self.log.warning(f'Use "unidentified" NCBITaxon ID for {fb_data_entity}')
            ncbi_taxon_id = 'NCBITaxon:32644'
        fb_data_entity.ncbi_taxon_id = ncbi_taxon_id
        return

    def synthesize_secondary_ids(self, fb_data_entity):
        """Process secondary IDs and return a list of old FB uniquenames."""
        secondary_ids = []
        for xref in fb_data_entity.dbxrefs:
            if xref.dbxref.db.name == 'FlyBase' and xref.is_current is False:
                secondary_ids.append(f'FB:{xref.dbxref.accession}')
        fb_data_entity.alt_fb_ids = list(set(secondary_ids))
        return

    def synthesize_pubs(self, fb_data_entity):
        """Collect pub_ids associated directly or indirectly with the entity."""
        pub_sources = ['pubs', 'synonyms', 'cvterms']
        for pub_source in pub_sources:
            fb_data_entity.all_pub_ids.extend([i.pub_id for i in getattr(fb_data_entity, pub_source)])
        for prop_pub_id_list in fb_data_entity.prop_pubs.values():
            fb_data_entity.all_pub_ids.extend(prop_pub_id_list)
        fb_data_entity.all_pub_ids = list(set(fb_data_entity.all_pub_ids))
        return

    def synthesize_info(self):
        """Extend the method for the PrimaryEntityHandler."""
        super().synthesize_info()
        for fb_data_entity in self.fb_data_entities.values():
            self.synthesize_ncbi_taxon_id(fb_data_entity)
            self.synthesize_secondary_ids(fb_data_entity)
            self.synthesize_pubs(fb_data_entity)
        return

    # Elaborate on map_fb_data_to_alliance() sub-methods for PrimaryEntityHandler.
    # However, we only create the methods here, and run them in more tailored handler Classes like StrainHandler.
    def map_secondary_ids(self, fb_data_entity):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for a FlyBase entity."""
        secondary_id_dtos = []
        for secondary_id in fb_data_entity.alt_fb_ids:
            sec_dto = datatypes.SecondaryIdSlotAnnotationDTO(secondary_id).dict_export()
            secondary_id_dtos.append(sec_dto)
        return secondary_id_dtos

    def map_pubs(self, fb_data_entity):
        """Add pub curies to a FlyBase entity."""
        for pub_id in fb_data_entity.all_pub_ids:
            fb_data_entity.linkmldto.reference_curies.append(self.bibliography[pub_id])
        try:
            fb_data_entity.reference_curies.remove('FB:unattributed')
        except ValueError:
            pass
        return

    def map_xrefs(self, fb_data_entity):
        """Add a list of Alliance CrossReferenceDTO dicts to a FlyBase entity."""
        cross_reference_dtos = []
        for xref in fb_data_entity.dbxrefs:
            # Skip xrefs from irrelevant database sources.
            if xref.dbxref.db.name not in self.fb_agr_db_dict.keys():
                continue
            # Skip FlyBase xrefs: current xref will be in data_provider_dto; others as 2o IDs.
            elif xref.dbxref.db.name == 'FlyBase':
                continue
            # Build Alliance xref DTO
            prefix = self.fb_agr_db_dict[xref.dbxref.db.name]
            curie = f'{prefix}{xref.dbxref.accession}'
            page_area = self.fb_data_type    # Must ensure that fb_data_types match Alliance resourceDescriptors.yaml page.
            display_name = xref.dbxref.description
            xref_dto = datatypes.CrossReferenceDTO(prefix, curie, page_area, display_name).dict_export()
            cross_reference_dtos.append(xref_dto)
        fb_data_entity.linkmldto.cross_reference_dtos = cross_reference_dtos
        return

    def map_timestamps(self, fb_data_entity):
        """Map timestamps to Alliance object."""
        if fb_data_entity.timestamps:
            fb_data_entity.linkmldto.date_created = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(fb_data_entity.timestamps)))
            fb_data_entity.linkmldto.date_updated_curie = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(fb_data_entity.timestamps)))

    def flag_internal_fb_entities(self, fb_data_entity):
        """Flag obsolete FB objects as internal."""
        if fb_data_entity.chado_obj.is_obsolete is True:
            fb_data_entity.linkmldto.internal = True
        return


class StrainHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters strain data for export."""
    def __init__(self, log, fb_data_type):
        """Create the StrainHandler object."""
        super().__init__(log, fb_data_type)

    # Elaborate on export filters for StrainHandler.
    required_fields = [
        'curie',
        'data_provider_dto',
        'internal',
        'subtype_name',
        'taxon_curie'
    ]
    output_fields = [
        'agm_secondary_id_dtos',
        'created_by_curie',
        'cross_reference_dtos',
        'curie',
        'data_provider_dto',
        'date_created',
        'date_updated',
        'internal',
        'updated_by_curie',
        'name',
        'obsolete',
        'reference_curies',
        'subtype_name',
        'taxon_curie'
    ]

    # Elaborate on map_fb_data_to_alliance() sub-methods for StrainHandler.
    def map_strain_basic(self, strain):
        """Map basic FlyBase strain data to the Alliance LinkML object."""
        agr_strain = datatypes.AffectedGenomicModelDTO()
        agr_strain.obsolete = strain.chado_obj.is_obsolete
        agr_strain.curie = f'FB:{strain.uniquename}'
        agr_strain.taxon_curie = strain.ncbi_taxon_id
        dp_xref = datatypes.CrossReferenceDTO('FB', strain.uniquename, 'strain', strain.chado_obj.name).dict_export()
        agr_strain.data_provider_dto = datatypes.DataProviderDTO(dp_xref).dict_export()
        agr_strain.name = strain.chado_obj.name
        agr_strain.subtype_name = 'strain'
        strain.linkmldto = agr_strain
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the StrainHandler."""
        super().map_fb_data_to_alliance()
        for strain in self.fb_data_entities.values():
            self.map_strain_basic(strain)
            self.map_pubs(strain)
            self.map_xrefs(strain)
            self.map_timestamps(strain)
            strain.linkmldto.agm_secondary_id_dtos = self.map_secondary_ids(strain)
            self.flag_internal_fb_entities(strain)
        return


# Functions
def get_handler(log: Logger, fb_data_type: str):
    """Return the appropriate type of data handler.

    Args:
        log (Logger): The global Logger object in the script using the DataHandler.
        fb_data_type (str): The FB data type to export: e.g., strain, genotype.

    Returns:
        A data handler of the appropriate type for the FB data type.

    Raises:
        Raises a KeyError if the FB data type is not recognized.

    """
    log.info(f'Get handler for {fb_data_type}.')
    handler_dict = {
        # 'gene': GeneHandler,
        # 'allele': AlleleHandler,
        # 'construct': ConstructHandler,
        # 'variation': VariationHandler,
        'strain': StrainHandler,
        # 'genotype': GenotypeHandler,
        # 'disease': DiseaseHandler,
    }
    try:
        data_handler = handler_dict[fb_data_type](log, fb_data_type)
        log.info(f'Returning: {data_handler}')
    except KeyError:
        log.error.append(f'Unrecognized FB data type and/or Alliance ingest set: {fb_data_type}.')
        raise
    return data_handler


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
