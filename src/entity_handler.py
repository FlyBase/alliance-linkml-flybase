"""Module:: entity_handler.

Synopsis:
    A generic data handler for FlyBase first class entities that have a FlyBase
    curie and, typically, a web report; this handler is not designed to handle
    complex annotations (e.g., disease, phenotype).

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import re
from logging import Logger
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, CellLine, CellLineCvterm, CellLineCvtermprop, CellLineDbxref, CellLineprop, CellLinepropPub, CellLinePub, CellLineRelationship,
    CellLineSynonym, Feature, FeatureCvterm, FeatureCvtermprop, FeatureDbxref, Featureprop, FeaturepropPub, FeaturePub, FeatureRelationship,
    FeatureRelationshipPub, FeatureSynonym, Genotype, GenotypeCvterm, GenotypeCvtermprop, GenotypeDbxref, Genotypeprop, GenotypepropPub, GenotypePub,
    GenotypeSynonym, Grp, GrpCvterm, GrpDbxref, Grpprop, GrppropPub, GrpPub, GrpRelationship, GrpRelationshipPub, GrpSynonym, Humanhealth, HumanhealthCvterm,
    HumanhealthCvtermprop, HumanhealthDbxref, Humanhealthprop, HumanhealthpropPub, HumanhealthPub, HumanhealthRelationship, HumanhealthRelationshipPub,
    HumanhealthSynonym, Library, LibraryCvterm, LibraryCvtermprop, LibraryDbxref, Libraryprop, LibrarypropPub, LibraryPub, LibraryRelationship,
    LibraryRelationshipPub, LibrarySynonym, Strain, StrainCvterm, StrainCvtermprop, StrainDbxref, Strainprop, StrainpropPub, StrainPub, StrainRelationship,
    StrainRelationshipPub, StrainSynonym
)
import agr_datatypes
import fb_datatypes
from handler import DataHandler


class PrimaryEntityHandler(DataHandler):
    """A generic data handler for that gets FlyBase data for first class entities.

    This object handles only primary FlyBase entities, things that have a FB
    curie and, typically, web reports. These entities typically have props and
    CV term associations, and they are the subjects of relationships and
    annotations (e.g., feature_relationship, phenstatement).

    """
    def __init__(self, log: Logger, testing: bool):
        """Create the generic PrimaryEntityHandler object."""
        super().__init__(log, testing)

    # Feature sub-types that are considered their own data class.
    feature_datatypes = [
        'aberration',
        'allele',
        'balancer',
        'chemical',
        'construct',
        'gene',
        'insertion',
        'variation',
    ]

    # CVterms used to define a fb_data_type within a larger chado table.
    subtypes = {
        'aberration': ['chromosome_structure_variation'],
        'allele': ['allele'],
        'balancer': ['chromosome_structure_variation'],
        'chemical': ['chemical entity'],
        'construct': ['engineered_transposable_element', 'engineered_region', 'transgenic_transposable_element'],
        'gene': ['gene'],
        'insertion': ['insertion_site', 'transposable_element', 'transposable_element_insertion_site'],
        'variation': ['MNV', 'complex_substitution', 'deletion', 'delins', 'insertion', 'point_mutation', 'sequence_alteration', 'sequence_variant']
    }

    # Mappings of main data types to chado tables with associated data like cvterms, props, synonyms, etc.
    chado_tables = {
        'main_table': {
            'cell_line': CellLine,
            'feature': Feature,
            'genotype': Genotype,
            'grp': Grp,
            'humanhealth': Humanhealth,
            'library': Library,
            'strain': Strain,
        },
        'pubs': {
            'cell_line': CellLinePub,
            'feature': FeaturePub,
            'genotype': GenotypePub,
            'grp': GrpPub,
            'humanhealth': HumanhealthPub,
            'library': LibraryPub,
            'strain': StrainPub,
        },
        'synonyms': {
            'cell_line': CellLineSynonym,
            'feature': FeatureSynonym,
            'genotype': GenotypeSynonym,
            'grp': GrpSynonym,
            'humanhealth': HumanhealthSynonym,
            'library': LibrarySynonym,
            'strain': StrainSynonym,
        },
        'dbxrefs': {
            'cell_line': CellLineDbxref,
            'feature': FeatureDbxref,
            'genotype': GenotypeDbxref,
            'grp': GrpDbxref,
            'humanhealth': HumanhealthDbxref,
            'library': LibraryDbxref,
            'strain': StrainDbxref,
        },
        'props': {
            'cell_line': CellLineprop,
            'feature': Featureprop,
            'genotype': Genotypeprop,
            'grp': Grpprop,
            'humanhealth': Humanhealthprop,
            'library': Libraryprop,
            'strain': Strainprop,
        },
        'prop_pubs': {
            'cell_line': CellLinepropPub,
            'feature': FeaturepropPub,
            'genotype': GenotypepropPub,
            'grp': GrppropPub,
            'humanhealth': HumanhealthpropPub,
            'library': LibrarypropPub,
            'strain': StrainpropPub,
        },
        'cvterms': {
            'cell_line': CellLineCvterm,
            'feature': FeatureCvterm,
            'genotype': GenotypeCvterm,
            'grp': GrpCvterm,
            'humanhealth': HumanhealthCvterm,
            'library': LibraryCvterm,
            'strain': StrainCvterm,
        },
        'cvtermprops': {
            'cell_line': CellLineCvtermprop,
            'feature': FeatureCvtermprop,
            'genotype': GenotypeCvtermprop,
            'grp': None,
            'humanhealth': HumanhealthCvtermprop,
            'library': LibraryCvtermprop,
            'strain': StrainCvtermprop,
        },
        'relationships': {
            'cell_line': CellLineRelationship,
            'feature': FeatureRelationship,
            'genotype': None,
            'grp': GrpRelationship,
            'humanhealth': HumanhealthRelationship,
            'library': LibraryRelationship,
            'strain': StrainRelationship,
        },
        'relationship_pubs': {
            'cell_line': None,
            'feature': FeatureRelationshipPub,
            'genotype': None,
            'grp': GrpRelationshipPub,
            'humanhealth': HumanhealthRelationshipPub,
            'library': LibraryRelationshipPub,
            'strain': StrainRelationshipPub,
        },
    }

    # Add methods to be run by get_general_data() below.
    # Placeholder.

    # Add methods to be run by get_datatype_data() below.
    def get_entities(self, session):
        """Get primary FlyBase data entities."""
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        self.log.info(f'Get {self.datatype} data entities from {chado_type} table.')
        chado_table = self.chado_tables['main_table'][chado_type]
        filters = ()
        if self.datatype in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.datatype]}')
            filters += (chado_table.uniquename.op('~')(self.regex[self.datatype]), )
        if self.datatype in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.datatype]}')
            filters += (Cvterm.name.in_((self.subtypes[self.datatype])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (chado_table.uniquename.in_((self.test_set.keys())), )
        if filters == ():
            self.log.warning('Have no filters for the main FlyBase entity driver query.')
            raise
        if self.datatype in self.subtypes.keys():
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
            self.fb_data_entities[pkey_id] = self.fb_export_type(result)
            counter += 1
        self.log.info(f'Found {counter} FlyBase {self.datatype} entities in chado.')
        return

    def get_entity_relationships(self, session, role):
        """Get relationships between primary FlyBase entities and other entities in the same chado table.

        Args:
            session (SQLAlchemy session): The session.
            role (str): Filter for relationships where the primary FB entities are the 'subject' or the 'object', no other values allowed.

        """
        role_parameters = {
            'subject': 'sbj_rels_by_type',
            'object': 'obj_rels_by_type',
        }
        if role not in role_parameters.keys():
            self.log.error(f'For "get_entity_relationships()", role was specified as "{role}"; only these values are allowed: "subject", "object".')
            raise
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        self.log.info(f'Get relationships from the {chado_type}_relationship table where the {self.datatype} is the {role}.')
        chado_table = self.chado_tables['main_table'][chado_type]
        entity_key_name = f'{chado_type}_id'
        chado_rel_table = self.chado_tables['relationships'][chado_type]
        if chado_rel_table is None:
            msg = f'The get_entity_relationships() method has been called unnecessarily, because for '
            msg += f'{self.datatype}s, there is no {self.datatype}_relationship table.'
            self.log.warning(msg)
            return
        chado_rel_pub_table = self.chado_tables['relationship_pubs'][chado_type]
        # Phase 1: Get all relationships.
        filters = ()
        if self.datatype in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.datatype]}')
            filters += (chado_table.uniquename.op('~')(self.regex[self.datatype]), )
        if self.datatype in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.datatype]}')
            filters += (Cvterm.name.in_((self.subtypes[self.datatype])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (chado_table.uniquename.in_((self.test_set.keys())), )
        if filters == ():
            self.log.warning('Have no filters for the main FlyBase entity driver query.')
            raise
        if self.datatype in self.subtypes.keys():
            rel_results = session.query(chado_rel_table).\
                select_from(chado_table).\
                join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
                join(chado_rel_table, (getattr(chado_rel_table, f'{role}_id') == getattr(chado_table, entity_key_name))).\
                filter(*filters).\
                distinct()
        else:
            rel_results = session.query(chado_rel_table).\
                select_from(chado_table).\
                join(chado_rel_table, (getattr(chado_rel_table, f'{role}_id') == getattr(chado_table, entity_key_name))).\
                filter(*filters).\
                distinct()
        rel_dict = {}    # A temporary feature_relationship_id-keyed dict of FBRelationship objects.
        rel_counter = 0
        for rel_result in rel_results:
            rel_id = getattr(rel_result, f'{chado_type}_relationship_id')
            rel_dict[rel_id] = fb_datatypes.FBRelationship(rel_result)
            rel_counter += 1
        self.log.info(f'Found {rel_counter} {chado_type}_relationships where the {self.datatype} is the {role}.')
        # Phase 2. Get pubs supporting relationships.
        if chado_rel_pub_table is None:
            rel_pub_results = []
        elif self.datatype in self.subtypes.keys():
            rel_pub_results = session.query(chado_rel_pub_table).\
                select_from(chado_table).\
                join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
                join(chado_rel_table, (getattr(chado_rel_table, f'{role}_id') == getattr(chado_table, entity_key_name))).\
                join(chado_rel_pub_table).\
                filter(*filters).\
                distinct()
        else:
            rel_pub_results = session.query(chado_rel_pub_table).\
                select_from(chado_table).\
                join(chado_rel_table, (getattr(chado_rel_table, f'{role}_id') == getattr(chado_table, entity_key_name))).\
                join(chado_rel_pub_table).\
                filter(*filters).\
                distinct()
        rel_pub_counter = 0
        for rel_pub_result in rel_pub_results:
            rel_id = getattr(rel_pub_result, f'{chado_type}_relationship_id')
            try:
                rel_dict[rel_id].pubs.append(rel_pub_result.pub_id)
                rel_pub_counter += 1
            except KeyError:
                pass
        self.log.info(f'Found {rel_pub_counter} {chado_type}_relationship_pubs where the {self.datatype} is the {role}.')
        # Phase 3. Add rel info to entities.
        assignment_counter = 0
        rel_type_tally = {}
        for rel in rel_dict.values():
            # Assign the prop to the appropriate entity.
            entity_id = getattr(rel.chado_obj, f'{role}_id')
            entity_rel_dict = getattr(self.fb_data_entities[entity_id], role_parameters[role])
            try:
                entity_rel_dict[rel.chado_obj.type.name].append(rel)
                assignment_counter += 1
            except KeyError:
                entity_rel_dict[rel.chado_obj.type.name] = [rel]
                assignment_counter += 1
            # Tally
            try:
                rel_type_tally[rel.chado_obj.type.name] += 1
            except KeyError:
                rel_type_tally[rel.chado_obj.type.name] = 1
        self.log.info(f'Assigned {assignment_counter} {chado_type}_relationships where the {self.datatype} is the {role}.')
        self.log.info(f'Found these types of {chado_type}_relationship types where the {self.datatype} is the {role}:')
        ordered_rel_types = sorted(list(rel_type_tally.keys()))
        for rel_type in ordered_rel_types:
            self.log.info(f'table={chado_type}, rel_type={rel_type}, count={rel_type_tally[rel_type]}.')
        return

    def get_entityprops(self, session):
        """Get primary FlyBase data entity props."""
        self.log.info('Get primary FlyBase data entity props.')
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        self.log.info(f'Get {self.datatype} props from {chado_type}prop table.')
        chado_table = self.chado_tables['main_table'][chado_type]
        subject_key_name = f'{chado_type}_id'
        chado_prop_table = self.chado_tables['props'][chado_type]
        chado_prop_pub_table = self.chado_tables['prop_pubs'][chado_type]
        # Phase 1: Get all props.
        filters = ()
        if self.datatype in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.datatype]}')
            filters += (chado_table.uniquename.op('~')(self.regex[self.datatype]), )
        if self.datatype in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.datatype]}')
            filters += (Cvterm.name.in_((self.subtypes[self.datatype])), )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (chado_table.uniquename.in_((self.test_set.keys())), )
        if filters == ():
            self.log.warning('Have no filters for the main FlyBase entity driver query.')
            raise
        if self.datatype in self.subtypes.keys():
            prop_results = session.query(chado_prop_table).\
                select_from(chado_table).\
                join(chado_prop_table).\
                join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
                filter(*filters).\
                distinct()
        else:
            prop_results = session.query(chado_prop_table).\
                select_from(chado_table).\
                join(chado_prop_table).\
                filter(*filters).\
                distinct()
        prop_dict = {}    # A temporary prop_id-keyed dict of prop objects.
        prop_counter = 0
        for prop_result in prop_results:
            prop_id = getattr(prop_result, f'{chado_type}prop_id')
            prop_dict[prop_id] = fb_datatypes.FBProp(prop_result)
            prop_counter += 1
        self.log.info(f'Found {prop_counter} {chado_type}props for {self.datatype}s.')
        # Phase 2. Get pubs supporting props.
        if self.datatype in self.subtypes.keys():
            prop_pub_results = session.query(chado_prop_pub_table).\
                select_from(chado_table).\
                join(chado_prop_table).\
                join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
                join(chado_prop_pub_table).\
                filter(*filters).\
                distinct()
        else:
            prop_pub_results = session.query(chado_prop_pub_table).\
                select_from(chado_table).\
                join(chado_prop_table).\
                join(chado_prop_pub_table).\
                filter(*filters).\
                distinct()
        prop_pub_counter = 0
        for prop_pub_result in prop_pub_results:
            prop_id = getattr(prop_result, f'{chado_type}prop_id')
            try:
                prop_dict[prop_id].pubs.append(prop_pub_result.pub_id)
                prop_pub_counter += 1
            except KeyError:
                pass
        self.log.info(f'Found {prop_pub_counter} {chado_type}prop_pubs for {self.datatype}s.')
        # Phase 3. Add prop info to entities.
        assignment_counter = 0
        prop_type_tally = {}
        for prop in prop_dict.values():
            # Assign the prop to the appropriate entity.
            subject_id = getattr(prop.chado_obj, subject_key_name)
            try:
                self.fb_data_entities[subject_id].props_by_type[prop.chado_obj.type.name].append(prop)
                assignment_counter += 1
            except KeyError:
                self.fb_data_entities[subject_id].props_by_type[prop.chado_obj.type.name] = [prop]
                assignment_counter += 1
            # Tally
            try:
                prop_type_tally[prop.chado_obj.type.name] += 1
            except KeyError:
                prop_type_tally[prop.chado_obj.type.name] = 1
        self.log.info(f'Assigned {assignment_counter} {chado_type}props to {self.datatype}s.')
        self.log.info(f'Found these types of {chado_type}props:')
        ordered_prop_types = sorted(list(prop_type_tally.keys()))
        for prop_type in ordered_prop_types:
            self.log.info(f'table={chado_type}, prop_type={prop_type}, count={prop_type_tally[prop_type]}.')
        return

    def get_entity_pubs(self, session):
        """Get pubs directly associated with FlyBase data entities."""
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        asso_chado_table = self.chado_tables['pubs'][chado_type]
        self.log.info(f'Get pubs for {self.datatype} data entities from {asso_chado_table}.')
        main_pkey_name = f'{chado_type}_id'
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
        self.log.info(f'Found {counter} pubs for {self.datatype} entities.')
        self.log.info(f'Ignored {pass_counter} pubs for irrelevant {self.datatype} entities.')
        return

    def get_entity_synonyms(self, session):
        """Get synonyms for the FlyBase data entities."""
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        asso_chado_table = self.chado_tables['synonyms'][chado_type]
        self.log.info(f'Get synonyms for {self.datatype} data entities from {asso_chado_table}.')
        main_pkey_name = f'{chado_type}_id'
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
        self.log.info(f'Found {counter} synonyms for {self.datatype} entities.')
        self.log.info(f'Ignored {pass_counter} synonyms for irrelevant {self.datatype} entities.')
        return

    def get_entity_fb_xrefs(self, session):
        """Get secondary FB xrefs for the FlyBase data entities."""
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        asso_chado_table = self.chado_tables['dbxrefs'][chado_type]
        self.log.info(f'Get non-current FlyBase xrefs for {self.datatype} data entities from {asso_chado_table}.')
        main_pkey_name = f'{chado_type}_id'
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
        self.log.info(f'Found {counter} 2o FB xrefs for {self.datatype} entities.')
        self.log.info(f'Ignored {pass_counter} 2o FB xrefs for irrelevant {self.datatype} entities.')
        return

    def get_entity_xrefs(self, session):
        """Get all other xrefs for the FlyBase data entities."""
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        asso_chado_table = self.chado_tables['dbxrefs'][chado_type]
        self.log.info(f'Get non-FlyBase xrefs for {self.datatype} data entities from {asso_chado_table}.')
        main_pkey_name = f'{chado_type}_id'
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
        self.log.info(f'Found {counter} xrefs for {self.datatype} entities.')
        self.log.info(f'Ignored {pass_counter} xrefs for irrelevant {self.datatype} entities.')
        return

    def get_entity_timestamps(self, session):
        """Get timestamps for data entities."""
        self.log.info(f'Get timestamps for FlyBase {self.datatype} entities.')
        if self.datatype in self.feature_datatypes:
            chado_type = 'feature'
        else:
            chado_type = self.datatype
        entity_table_counter = 0
        audit_chado_counter = 0
        # Get distinct timestamps for each entity (do not distinguish by action, etc).
        for i in self.fb_data_entities.values():
            if i.timeaccessioned is not None:
                i.timestamps.append(i.timeaccessioned)
            if i.timelastmodified is not None:
                i.timestamps.append(i.timelastmodified)
                entity_table_counter += 1
            if not i.timestamps:
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
                # Pick out current symbol for the entity.
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] in ['systematic_name', 'nomenclature_symbol']:
                    fb_data_entity.curr_fb_symbol = syno_dict['display_text']
                # Pick out current full name for the entity.
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] == 'full_name':
                    fb_data_entity.curr_fb_fullname = syno_dict['display_text']

        return

    def synthesize_pubs(self):
        """Collect pub_ids associated directly or indirectly with the entity."""
        self.log.info('Collect pub_ids associated directly or indirectly with the entity.')
        pub_sources = ['pubs', 'synonyms']
        for fb_data_entity in self.fb_data_entities.values():
            for pub_source in pub_sources:
                fb_data_entity.all_pubs.extend([i.pub_id for i in getattr(fb_data_entity, pub_source)])
            for prop_list in fb_data_entity.props_by_type.values():
                for prop in prop_list:
                    fb_data_entity.all_pubs.extend(prop.pubs)
            fb_data_entity.all_pubs = list(set(fb_data_entity.all_pubs))
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_data_provider_dto(self):
        """Return the DataProviderDTO for the FB data entity."""
        # Note - this method is depends on previous determination of fb_data_entity.curr_fb_symbol by synthesize_synonyms(), if applicable.
        self.log.info('Map data provider to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.curr_fb_symbol:
                display_name = fb_data_entity.curr_fb_symbol
            else:
                display_name = fb_data_entity.name
            dp_xref = agr_datatypes.CrossReferenceDTO('FB', f'FB:{fb_data_entity.uniquename}', self.datatype, display_name).dict_export()
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
            for pub_id in fb_data_entity.all_pubs:
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
            # First, add FB xref (since FB xrefs in chado are not complete, just use the uniquename).
            if fb_data_entity.uniquename:
                curie = f'FB:{fb_data_entity.uniquename}'
                display_name = curie
                page_area = self.datatype
                fb_xref_dto = agr_datatypes.CrossReferenceDTO('FB', curie, page_area, display_name).dict_export()
                cross_reference_dtos.append(fb_xref_dto)
            # Second, add external xrefs.
            for xref in fb_data_entity.dbxrefs:
                # Build Alliance xref DTO
                prefix = self.fb_agr_db_dict[xref.dbxref.db.name]
                # The page_area assignment assumes that the self.datatype has a matching value in the Alliance resourceDescriptors.yaml page.
                try:
                    page_area = self.agr_page_area_dict[prefix]
                except KeyError:
                    page_area = self.datatype
                # Clean up cases where the db prefix is redundantly included at the start of the dbxref.accession.
                redundant_prefix = f'{prefix}:'
                if xref.dbxref.accession.startswith(redundant_prefix):
                    cleaned_accession = xref.dbxref.accession.replace(redundant_prefix, '', 1)
                    # self.log.debug(f'Removed "{redundant_prefix}" from "{xref.dbxref.accession}" to give "{cleaned_accession}"')
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
        linkml_dto_attributes = self.agr_export_type().__dict__.keys()
        for dto_key in linkml_dto_attributes:
            for bin_type, bin_suffix in linkml_synonym_slots.items():
                if dto_key.endswith(bin_suffix):
                    linkml_synonym_slots[bin_type] = dto_key
                    self.log.debug(f'Map {bin_type} to LinkML DTO slot {dto_key} because it has suffix "{bin_suffix}".')
                    map_synonyms_required = True
        if map_synonyms_required is False:
            self.log.error(f'The map_synonyms() method has been incorrectly called for {self.datatype} objects.')
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
                elif syno_dict['is_current'] is True and syno_dict['name_type_name'] == 'full_name':
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
