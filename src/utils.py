"""Module:: utils.

Synopsis:
    Data retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import csv
import datetime
import json
import re
import strict_rfc3339
from logging import Logger
from sqlalchemy.orm import aliased, Session
# from sqlalchemy.inspection import inspect
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Organism, OrganismDbxref, Pub, PubDbxref, Featureloc,
    Strain, StrainPub, StrainSynonym, StrainDbxref, Strainprop, StrainpropPub, StrainCvterm, StrainCvtermprop,
    Feature, FeaturePub, FeatureSynonym, FeatureDbxref, Featureprop, FeaturepropPub, FeatureCvterm, FeatureCvtermprop
)
import datatypes


# Classes
class DataHandler(object):
    """A generic, abstract data handler that gets FlyBase data and maps it to a single Alliance LinkML model.

    Specific classes of DataHandler will only map a given FB data type to a single Alliance ingest
    set (i.e., a set of LinkML DTO objects, in JSON format). In some cases, multiple handlers will
    map different FB data to the same ingest set: e.g., FB strains and genotypes both go into the
    "agm_ingest_set" via distinct DataHandler classes.

    """
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the generic DataHandler object.

        Args:
            log (Logger): The global Logger object in the script using the DataHandler.
            fb_data_type (str): The FlyBase data class being handled.
            testing (bool): Whether handler is being run in testing mode or not.

        """
        self.log = log
        self.fb_data_type = fb_data_type
        self.agr_ingest_type = self.agr_ingest_type_dict[fb_data_type]
        self.testing = testing
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

    # Sample set for faster testing: use uniquename-keyed names of objects, tailored for each handler.
    test_set = {}
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
        'gene': datatypes.FBGene,
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
        'cell': r'^FBcl[0-9]{7}$',
        'pthr': r'PTHR[0-9]{5}',
        'gene_systematic_name': r'^(D[a-z]{3}\\|)(CG|CR|G[A-Z])[0-9]{4,5}$'
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
    """A generic, abstract handler for that gets data for FlyBase entities and maps it to the Alliance LinkML model.

    This object handles only primary FlyBase entities, things like genes, strains, etc. that typically have
    public curies and web reports, attributes like props and cvterm associations, and are the subjects
    in relationships (e.g., feature_relationship) and annotations (e.g., phenstatement). Separate handler
    classes are to be used for the export of relationships and complex annotations.

    """
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the generic PrimaryEntityHandler object."""
        super().__init__(log, fb_data_type, testing)

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

    def get_primary_key_column(self, chado_table):
        """Get the primary key Column object from a specified chado table Model object."""
        primary_key_column = next((column for column in chado_table.__table__.c if column.primary_key), None)
        if primary_key_column is None:
            self.log.error(f'Could not get primary_key Column from {chado_table}')
            raise ValueError
        else:
            self.log.debug(f'Found primary_key column: {primary_key_column.name}')
        return primary_key_column

    def get_foreign_key_column(self, chado_table, column_name):
        """Get a foreign key column object, by name, from a specified chado table Model object."""
        foreign_key_column = next((column for column in chado_table.__table__.c if column.foreign_keys and column.name == column_name), None)
        if foreign_key_column is None:
            self.log.error(f'Could not get foreign_key Column {column_name} from {chado_table}')
            raise ValueError
        else:
            self.log.debug(f'Found primary_key column: {foreign_key_column.name}')
        return foreign_key_column

    # BOB: Quick testing of general SQLAlchemy query approaches.
    def sqlalchemy_test(self, session):
        """Test SQLAlchemy behavior."""
        self.log.info('Test SQLAlchemy behavior.')
        lbe_types = {'gene': 'bob'}
        table_dict = {'main': Feature}
        main_table = table_dict['main']
        pkey_col = self.get_primary_key_column(main_table)
        fkey_col = self.get_foreign_key_column(main_table, 'type_id')
        feat_ids = [3167743]
        filters = (
            main_table.uniquename == 'FBgn0011278',
            pkey_col.in_((feat_ids)),
            fkey_col == 219
        )
        filters += (
            Cvterm.name.in_((lbe_types.keys())),
        )
        results = session.query(main_table).\
            select_from(main_table).\
            join(Cvterm, (Cvterm.cvterm_id == main_table.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for i in results:
            self.log.info(f'Found this feature: {i.name} ({i.uniquename})')
            self.log.info(f'Found this feature_id: {getattr(i, pkey_col.name)}')
            counter += 1
        self.log.info(f'Found {counter} test results using natural join.')
        return

    # Elaborate on get_general_data() sub-methods for PrimaryEntityHandler.
    def get_general_data(self, session):
        """Extend the method for the PrimaryEntityHandler."""
        super().get_general_data(session)
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

    def get_entity_associated_data(self, session):
        """Get data associated with primary FlyBase data entities."""
        associated_data_types = ['pubs', 'synonyms', 'dbxrefs', 'props', 'cvterms']
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        self.log.info(f'Get associated data for {self.fb_data_type} data entities from {chado_type}-related chado tables.')
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
        # BOB - is it more efficient to do "in_" query, or, query each fb entity one-by-one?
        for i in associated_data_types:
            asso_chado_table = self.chado_tables[i][chado_type]
            self.log.info(f'Get {i} for {self.fb_data_type} from {asso_chado_table}')
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
                    self.fb_data_entities[entity_pkey_id].__dict__[i].append(result)
                    counter += 1
                except KeyError:
                    pass_counter += 1
            self.log.info(f'Found {counter} {i} for {self.fb_data_type} entities.')
            self.log.info(f'Ignored {pass_counter} {i} for irrelevant {self.fb_data_type} entities.')
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
        # BOB - need a way to add a uniquename filter for testing.
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
        # self.get_entity_prop_pubs(session)
        # self.get_entity_cvtermprops(session)
        self.get_entity_timestamps(session)
        return

    # Elaborate on synthesize_info() sub-methods for PrimaryEntityHandler.
    def synthesize_ncbi_taxon_id(self):
        """Determine the NCBITaxon ID for the entity."""
        for fb_data_entity in self.fb_data_entities.values():
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

    def synthesize_secondary_ids(self):
        """Process secondary IDs and return a list of old FB uniquenames."""
        for fb_data_entity in self.fb_data_entities.values():
            secondary_ids = []
            for xref in fb_data_entity.dbxrefs:
                if xref.dbxref.db.name == 'FlyBase' and xref.is_current is False:
                    secondary_ids.append(f'FB:{xref.dbxref.accession}')
            fb_data_entity.alt_fb_ids = list(set(secondary_ids))
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
        pub_sources = ['pubs', 'synonyms', 'cvterms', 'prop_pubs']
        for fb_data_entity in self.fb_data_entities.values():
            for pub_source in pub_sources:
                fb_data_entity.all_pub_ids.extend([i.pub_id for i in getattr(fb_data_entity, pub_source)])
            fb_data_entity.all_pub_ids = list(set(fb_data_entity.all_pub_ids))
        return

    def synthesize_info(self):
        """Extend the method for the PrimaryEntityHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_props()
        self.synthesize_pubs()
        return

    # Elaborate on map_fb_data_to_alliance() sub-methods for PrimaryEntityHandler.
    # However, as they're not useful for all data types, call them in tailored handler Classes.
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
            fb_data_entity.linkmldto.reference_curies.remove('FB:unattributed')
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

    def map_synonyms(self, fb_data_entity):
        # BOB: TO DO
        # """Generate name/synonym DTOs for a feature that has a list of FeatureSynonym objects."""
        # # Dict for converting FB to AGR synonym types.
        # synonym_type_conversion = {
        #     'symbol': 'nomenclature_symbol',
        #     'fullname': 'full_name',
        #     'nickname': 'nomenclature_symbol',
        #     'synonym': 'nomenclature_symbol'
        # }
        # default_name_dto = {
        #     'name_type_name': 'unspecified',
        #     'format_text': 'unspecified',
        #     'display_text': 'unspecified',
        #     'synonym_scope_name': 'exact',
        #     'evidence_curies': [],
        #     'internal': False,
        #     'obsolete': False
        # }
        # # Create a dict of all distinct name/synonym_sgml combinations: for each, capture synonym type(s) an pub_ids.
        # # Keys are (synonym.name, synonym.synonym_sgml) tuples.
        # # Values are dicts too where keys are chado synonym types and values are lists of pub_ids.
        # # Value dict also has an "internal" key that stores list of FeatureSynonym.is_internal values.
        # feature_synonym_dict = {}
        # for f_s in feature.feature_synonyms:
        #     synonym = self.all_synonyms_dict[f_s.synonym_id]
        #     distinct_synonym_name = (synonym.name, synonym.synonym_sgml)
        #     if distinct_synonym_name in feature_synonym_dict.keys():
        #         feature_synonym_dict[distinct_synonym_name]['internal'].append(f_s.is_internal)
        #         if synonym.type.name in feature_synonym_dict[distinct_synonym_name].keys():
        #             feature_synonym_dict[distinct_synonym_name][synonym.type.name].append(f_s.pub_id)
        #         else:
        #             feature_synonym_dict[distinct_synonym_name][synonym.type.name] = [f_s.pub_id]
        #     else:
        #         feature_synonym_dict[distinct_synonym_name] = {synonym.type.name: [f_s.pub_id], 'internal': [f_s.is_internal]}
        # # Convert to AGR name DTO objects.
        # name_dto_list = []
        # FORMAT_TEXT = 0
        # DISPLAY_TEXT = 1
        # for syno_name, syno_attributes in feature_synonym_dict.items():
        #     # Determine internal status. False trumps True.
        #     if False in set(syno_attributes['internal']):
        #         syno_internal = False
        #     else:
        #         syno_internal = True
        #     # Collect all pubs.
        #     pub_id_list = []
        #     for syno_type, syno_type_pub_list in syno_attributes.items():
        #         if syno_type == 'internal':
        #             continue
        #         pub_id_list.extend(syno_type_pub_list)
        #     pub_id_list = list(set(pub_id_list))
        #     # Out of curiosity, report cases where same synonym used as both symbol and fullname.
        #     if 'symbol' in syno_attributes.keys() and 'fullname' in syno_attributes.keys():
        #         n_symb = len(syno_attributes['symbol'])
        #         n_full = len(syno_attributes['fullname'])
        #         log.warning(f"RATIO = {round(n_symb/n_full)}, SYMBOL_USAGE: n={n_symb}, FULLNAME_USAGE: n={n_full}, GENE={feature}, SYNONYM={syno_name}.")
        #     # Pick correct name type to apply.
        #     if re.match(self.systematic_name_regex, syno_name[DISPLAY_TEXT]):
        #         name_type_to_use = 'systematic_name'
        #     else:
        #         type_tally = {}
        #         for syno_type, syno_type_pub_list in syno_attributes.items():
        #             if syno_type == 'internal':
        #                 continue
        #             type_tally[len(set(syno_type_pub_list))] = syno_type
        #         name_type_to_use = synonym_type_conversion[type_tally[max(type_tally.keys())]]
        #     output_synonym_dto = {
        #         'name_type_name': name_type_to_use,
        #         'format_text': sub_sup_sgml_to_html(syno_name[FORMAT_TEXT]),
        #         'display_text': sub_sup_sgml_to_html(syno_name[DISPLAY_TEXT]),
        #         'synonym_scope_name': 'exact',
        #         'evidence_curies': [self.all_pubs_dict[i] for i in pub_id_list if self.all_pubs_dict[i] != 'FB:unattributed'],
        #         'internal': syno_internal,
        #         'obsolete': False
        #     }
        #     name_dto_list.append(output_synonym_dto)
        # # Sift through name DTOs for symbol, fullname, systematic_name, etc.
        # for name_dto in name_dto_list:
        #     if name_dto['display_text'] == feature.curr_anno_id:
        #         if name_dto['name_type_name'] != 'systematic_name':
        #             log.warning(f"{feature}: Found mistyped curr anno ID: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
        #             name_dto['name_type_name'] = 'systematic_name'
        #         feature.gene_systematic_name_dto = name_dto
        #     if name_dto['display_text'] == feature.curr_symbol_name:
        #         if name_dto['name_type_name'] not in ['systematic_name', 'nomenclature_symbol']:
        #             log.warning(f"{feature}: Found mistyped curr symbol: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
        #             name_dto['name_type_name'] = 'nomenclature_symbol'
        #         feature.gene_symbol_dto = name_dto
        #     elif name_dto['display_text'] == feature.curr_fullname:
        #         if name_dto['name_type_name'] != 'full_name':
        #             log.warning(f"{feature}: Found mistyped curr full_name: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
        #             name_dto['name_type_name'] = 'full_name'
        #         feature.gene_full_name_dto = name_dto
        #     else:
        #         feature.gene_synonym_dtos.append(name_dto)
        # # LinkML change required: make gene_full_name_dto and gene_systematic_name_dto OPTIONAL.
        # # Symbol is required. If none, fill it in.
        # if feature.gene_symbol_dto is None:
        #     placeholder_symbol_dto = default_name_dto.copy()
        #     placeholder_symbol_dto['name_type_name'] = 'nomenclature_symbol'
        #     placeholder_symbol_dto['format_text'] = feature.feature.name
        #     placeholder_symbol_dto['display_text'] = feature.feature.name
        #     feature.gene_symbol_dto = placeholder_symbol_dto
        # # In rare cases, a gene's annotation ID has never been used as a synonym. For these, fill in the annotation ID.
        # if feature.gene_systematic_name_dto is None and feature.curr_anno_id:
        #     placeholder_systematic_name_dto = default_name_dto.copy()
        #     placeholder_systematic_name_dto['name_type_name'] = 'systematic_name'
        #     placeholder_systematic_name_dto['format_text'] = feature.curr_anno_id
        #     placeholder_systematic_name_dto['display_text'] = feature.curr_anno_id
        #     log.warning(f"{feature}: Has annoID never used as a synonym: {feature.curr_anno_id}")
        #     feature.gene_systematic_name_dto = placeholder_systematic_name_dto
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

    def map_fb_data_to_alliance(self):
        """Extend the method for the PrimaryEntityHandler."""
        super().map_fb_data_to_alliance()
        return


class FeatureHandler(PrimaryEntityHandler):
    """A generic, abstract handler for that gets data for FlyBase features and maps it to the Alliance LinkML model."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the FeatureHandler object."""
        super().__init__(log, fb_data_type, testing)
        self.chr_dict = {}    # Will be a feature_id-keyed dict of chr scaffold uniquenames.

    # Elaborate on get_general_data() sub-methods for FeatureHandler.
    # Call get_chr_info() only for more specific FeatureHandler types.
    def get_chr_info(self, session):
        """Build chr dict."""
        self.log.info('Build chr dict.')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'golden_path',
            Organism.abbreviation == 'Dmel'
        )
        chr_results = session.query(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            filter(*filters).\
            distinct()
        self.chr_dict = {}
        counter = 0
        for result in chr_results:
            self.chr_dict[result.feature_id] = result.uniquename
            counter += 1
        self.log.info(f'Got basic info for {counter} current Dmel chr scaffolds.')
        return

    def get_general_data(self, session):
        """Extend the method for the FeatureHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() sub-methods for FeatureHandler.
    # Call get_annotation_ids() only for more specific FeatureHandler types.
    def get_annotation_ids(self, session):
        """Get annotation IDs (current and non-current)."""
        self.log.info('Get annotation IDs (current and non-current).')
        filters = (
            Feature.uniquename.op('~')(self.regex[self.fb_data_type]),
            Feature.is_analysis.is_(False),
            Db.name == 'FlyBase Annotation IDs'
        )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (Cvterm.name.in_((self.subtypes[self.fb_data_type])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(FeatureDbxref).\
            select_from(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = FeatureDbxref.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].fb_anno_xrefs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} annotation IDs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} annotation IDs for {self.fb_data_type} entities.')
        return

    # Call get_chr_featurelocs() only for more specific FeatureHandler types.
    def get_chr_featurelocs(self, session):
        """Get chromosomal featureloc data."""
        self.log.info('Get chromosomal featureloc data.')
        filters = (
            Feature.is_analysis.is_(False),
            Featureloc.srcfeature_id.in_((self.chr_dict.keys()))
        )
        if self.fb_data_type in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.fb_data_type]}')
            filters += (Feature.uniquename.op('~')(self.regex[self.fb_data_type]), )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (Cvterm.name.in_((self.subtypes[self.fb_data_type])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Featureloc).\
            select_from(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Featureloc, (Featureloc.feature_id == Feature.feature_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = Featureloc.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].chr_flocs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} chromosomal featurelocs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} chromosomal featurelocs for {self.fb_data_type} entities.')
        return

    # Call get_chr_featurelocs() only for more specific FeatureHandler types.
    def get_featureprops_by_type(self, session, fprop_type):
        """Return a list of featureprops of a given type."""
        self.log.info(f'Get featureprops of type {fprop_type} for {self.fb_data_type} entities.')
        feat_type = aliased(Cvterm, name='feat_type')
        prop_type = aliased(Cvterm, name='prop_type')
        filters = (
            Feature.is_analysis.is_(False),
            prop_type.name == fprop_type
        )
        if self.fb_data_type in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.fb_data_type]}')
            filters += (Feature.uniquename.op('~')(self.regex[self.fb_data_type]), )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (feat_type.name.in_((self.subtypes[self.fb_data_type])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Featureprop).\
            select_from(Feature).\
            join(feat_type, (feat_type.cvterm_id == Feature.type_id)).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        return results

    def get_datatype_data(self, session):
        """Extend the method for the FeatureHandler."""
        super().get_datatype_data(session)
        return

    # Elaborate on synthesize_info() sub-methods for FeatureHandler.
    def synthesize_info(self):
        """Extend the method for the FeatureHandler."""
        super().synthesize_info()
        return

    # Elaborate on map_fb_data_to_alliance() sub-methods for FeatureHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the FeatureHandler."""
        super().map_fb_data_to_alliance()
        return


class GeneHandler(FeatureHandler):
    """This object gets, synthesizes and filters gene data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the GeneHandler object."""
        super().__init__(log, fb_data_type, testing)
        self.pthr_dict = {}    # Will be an 1:1 FBgn_ID-PTHR xref dict.

    # Sample set for faster testing: use uniquename-keyed names of objects, tailored for each handler.
    test_set = {
        'FBgn0284084': 'wg',             # Current annotated nuclear protein_coding gene.
        'FBgn0004009': 'wg',             # Obsolete annotated nuclear protein_coding gene.
        'FBgn0013687': 'mt:ori',         # Current localized but unannotated mitochondrial gene.
        'FBgn0013678': 'mt:Cyt-b',       # Current annotated mitochondrial protein_coding gene.
        'FBgn0019661': 'lncRNA:roX1',    # Current annotated nuclear ncRNA gene.
        'FBgn0262451': 'mir-ban',        # Current annotated nuclear miRNA gene.
        'FBgn0031087': 'CG12656',        # Current withdrawn gene.
        'FBgn0000154': 'Bar',            # Current unannotated gene.
        'FBgn0001200': 'His4',           # Current unannotated gene family.
        'FBgn0087003': 'tal',            # Current unannotated oddball.
    }
    # Elaborate on export filters for StrainHandler.
    required_fields = [
        'curie',
        'data_provider_dto',
        'gene_symbol_dto',
        'internal',
        'taxon_curie',
    ]
    output_fields = [
        'created_by_curie',
        'cross_reference_dtos',
        'curie',
        'data_provider_dto',
        'date_created',
        'date_updated',
        'gene_full_name_dto',
        'gene_symbol_dto',
        'gene_synonym_dtos',
        'gene_systematic_name_dto',
        'gene_secondary_id_dtos',
        'gene_type_curie',
        'internal',
        'obsolete',
        'taxon_curie',
        'updated_by_curie',
    ]
    fb_agr_db_dict = {
        'EntrezGene': 'NCBI_Gene',
        'FlyBase': 'FB',
        'FlyBase Annotation IDs': 'FB',
        'RNAcentral': 'RNAcentral',
        # 'UniProt/GCRP': 'UniProt/GCRP',
        'UniProt/Swiss-Prot': 'UniProtKB',
        'UniProt/TrEMBL': 'UniProtKB'
    }
    # Reference dicts.
    internal_gene_types = [
        'engineered_fusion_gene',
        'engineered_region',
        'gene_group',
        'gene_with_polycistronic_transcript',
        'insulator',
        'mitochondrial_sequence',
        'origin_of_replication',
        'region',
        'regulatory_region',
        'repeat_region',
        'satellite_DNA',
        'transposable_element_gene'
    ]

    # Elaborate on get_general_data() sub-methods for GeneHandler.
    def get_general_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_general_data(session)
        self.get_chr_info(session)
        return

    # Elaborate on get_datatype_data() sub-methods for GeneHandler.
    def get_panther_info(self, input_dir):
        """Extract panther information from file."""
        self.log.info('Extract panther information from file.')
        if input_dir == '/src/input/':
            filepath = f'{input_dir}PTHR18.0_fruit_fly'
        else:
            filepath = '/data/ortholog/panther/PTHR18.0_fruit_fly'
        tsv_file = open(filepath, "r")
        tsvin = csv.reader(tsv_file, delimiter='\t')
        FB = 0
        PTHR = 3
        counter = 0
        gene_regex = r'FBgn[0-9]{7}'
        for row in tsvin:
            fields = len(row)
            if fields:  # Ignore blank lines
                if re.search(gene_regex, row[FB]) and re.search(self.pthr_regex, row[PTHR]):
                    self.pthr_dict[re.search(gene_regex, row[FB]).group(0)] = re.search(self.pthr_regex, row[PTHR]).group(0)
                    counter += 1
        self.log.info(f'Processed {counter} lines from the panther orthology file.')
        return

    def get_gene_snapshots(self, session):
        """Get human-written gene summaries."""
        self.log.info('Get human-written gene summaries.')
        results = self.get_featureprops_by_type(session, 'gene_summary_text')
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].gene_snapshots.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} gene snapshots for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} gene snapshots for {self.fb_data_type} entities.')
        return

    def get_gene_types(self, session):
        """Get promoted_gene_type for genes."""
        self.log.info('Get promoted_gene_type for genes.')
        results = self.get_featureprops_by_type(session, 'promoted_gene_type')
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].gene_type_names.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} promoted_gene_types for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} promoted_gene_types for {self.fb_data_type} entities.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session)
        self.get_panther_info(input_dir)
        self.get_annotation_ids(session)
        self.get_chr_featurelocs(session)
        self.get_gene_snapshots(session)
        self.get_gene_types(session)
        return

    # Elaborate on synthesize_info() sub-methods for GeneHandler.
    def synthesize_info(self):
        """Extend the method for the GeneHandler."""
        super().synthesize_info()
        return

    # Elaborate on map_fb_data_to_alliance() sub-methods for GeneHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        return


class StrainHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters strain data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the StrainHandler object."""
        super().__init__(log, fb_data_type, testing)

    # Sample set for faster testing: use uniquename-keyed names of objects, tailored for each handler.
    test_set = {
        'FBsn0000001': 'Oregon-R-modENCODE',
        'FBsn0000091': 'DGRP-373',
        'FBsn0000272': 'iso-1',
        'FBsn0001072': 'DSPR-B1-019',
        'FBsn0000283': 'MV2-25',
        'FBsn0000284': 'DGRP_Flyland',
    }
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

    # Elaborate on get_general_data() sub-methods for StrainHandler.
    def get_general_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() sub-methods for StrainHandler.
    def get_datatype_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_datatype_data(session)
        return

    # Elaborate on synthesize_info() sub-methods for StrainHandler.
    def synthesize_info(self):
        """Extend the method for the StrainHandler."""
        super().synthesize_info()
        return

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
        # 'allele': AlleleHandler,
        # 'construct': ConstructHandler,
        # 'variation': VariationHandler,
        'strain': StrainHandler,
        # 'genotype': GenotypeHandler,
        # 'disease': DiseaseHandler,
    }
    try:
        data_handler = handler_dict[fb_data_type](log, fb_data_type, testing)
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
