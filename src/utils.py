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
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Organism, OrganismDbxref, Pub, PubDbxref, Featureloc,
    Strain, StrainPub, StrainSynonym, StrainDbxref, Strainprop, StrainpropPub, StrainCvterm, StrainCvtermprop,
    Feature, FeaturePub, FeatureSynonym, FeatureDbxref, Featureprop, FeaturepropPub, FeatureCvterm, FeatureCvtermprop,
    FeatureRelationship, FeatureRelationshipPub, Synonym
)
import datatypes


# Classes
class DataHandler(object):
    """A generic, abstract data handler that gets FlyBase data and maps it to LinkML model(s).

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
        self.primary_export_set = self.primary_agr_ingest_type_dict[fb_data_type]
        self.testing = testing
        # Datatype bins.
        self.fb_data_entities = {}     # db_primary_id-keyed dict of chado objects to export.
        self.export_data = {}          # agr_ingest_set_name-keyed lists of data objects for export.
        # General data bins.
        self.bibliography = {}         # A pub_id-keyed dict of pub curies (PMID or FBrf).
        self.cvterm_lookup = {}        # A cvterm_id-keyed dict of Cvterm objects.
        self.ncbi_taxon_lookup = {}    # An organism_id-keyed dict of f'NCBITaxon:{Dbxref.accession}' strings.
        self.feature_lookup = {}       # feature_id-keyed lookup of basic feature info: uniquename, is_obsolete, name, symbol, exported, taxon_id and species.
        # Trackers.
        self.input_count = 0           # Count of entities found in FlyBase chado database.
        self.export_count = 0          # Count of exported Alliance entities.
        self.internal_count = 0        # Count of exported entities marked as internal.
        self.warnings = []             # Handler issues of note.
        self.errors = []               # Handler issues that break things.

    # Sample set for faster testing: use uniquename-keyed names of objects, tailored for each handler.
    test_set = {}
    # Correspondence of FB data type to primary Alliance LinkML object.
    agr_linkmldto_dict = {
        'gene': datatypes.GeneDTO(),
        'allele': 'TBD',
        'construct': datatypes.ConstructDTO(),
        'variation': 'TBD',
        'strain': datatypes.AffectedGenomicModelDTO(),
        'genotype': datatypes.AffectedGenomicModelDTO(),
        'disease': 'TBD'
    }
    # Correspondence of FB data type to Alliance data transfer ingest set.
    primary_agr_ingest_type_dict = {
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
        'construct': datatypes.FBConstruct,
        # 'variation': datatypes.FBVariant,
        'strain': datatypes.FBStrain,
        # 'disease': datatypes.FBDiseaseAlleleAnnotation
    }
    # Export directions (must be filled out in detail for each specific data handler), keyed by export set name.
    # For a given export set, the list of field names must be present in the datatype object specified in DataHandler.datatype_objects.values().
    required_fields = {}
    output_fields = {}
    # A filter for non-FB xrefs to export, with dict keys as FB db.names and dict values as Alliance db names.
    # Alliance db names should correspond to the contents of this file:
    # https://github.com/alliance-genome/agr_schemas/blob/master/resourceDescriptors.yaml
    fb_agr_db_dict = {}

    # Useful regexes.
    regex = {
        'pub': r'^(FBrf[0-9]{7}|unattributed)$',
        'abbal': r'^FB(ab|ba)[0-9]{7}$',
        'aberration': r'^FBab[0-9]{7}$',
        'allele': r'^FBal[0-9]{7}$',
        'balancer': r'^FBba[0-9]{7}$',
        'chem': r'^FBch[0-9]{7}$',
        'consins': r'^FB(tp|ti)[0-9]{7}$',
        'construct': r'^FBtp[0-9]{7}$',
        'fb_uniquename': r'^FB[a-z]{2}[0-9]{7,10}$',
        'gene': r'^FBgn[0-9]{7}$',
        'insertion': r'^FBti[0-9]{7}$',
        'seqfeat': r'^FBsf[0-9]{10}$',
        'tool': r'^FBto[0-9]{7}$',
        'genotype': r'^FBgo[0-9]{7}$',
        'strain': r'^FBsn[0-9]{7}$',
        'library': r'^FBlc[0-9]{7}$',
        'cell': r'^FBcl[0-9]{7}$',
        'panther': r'PTHR[0-9]{5}',
        'systematic_name': r'^(D[a-z]{3}\\|)(CG|CR|G[A-Z])[0-9]{4,5}'
    }

    # Methods
    # General utilities.
    def __str__(self):
        """Print out data handler description."""
        handler_description = f'A data handler that exports FB {self.fb_data_type} to Alliance LinkML: {list(self.output_fields.keys())}.'
        return handler_description

    def get_primary_key_column(self, chado_table):
        """Get the primary key Column object from a specified chado table Model object."""
        primary_key_column = next((column for column in chado_table.__table__.c if column.primary_key), None)
        if primary_key_column is None:
            self.log.error(f'Could not get primary_key Column from {chado_table}')
            raise ValueError
        else:
            pass
            # self.log.debug(f'Found primary_key column: {primary_key_column.name}')
        return primary_key_column

    def get_foreign_key_column(self, chado_table, column_name):
        """Get a foreign key column object, by name, from a specified chado table Model object."""
        foreign_key_column = next((column for column in chado_table.__table__.c if column.foreign_keys and column.name == column_name), None)
        if foreign_key_column is None:
            self.log.error(f'Could not get foreign_key Column {column_name} from {chado_table}')
            raise ValueError
        else:
            pass
            # self.log.debug(f'Found primary_key column: {foreign_key_column.name}')
        return foreign_key_column

    # Sub-methods for the get_general_data() wrapper.
    # Some of these sub-methods may only be called in more specific handlers, as appropriate.
    def sqlalchemy_test(self, session):
        """Test SQLAlchemy behavior."""
        self.log.info('Test SQLAlchemy behavior.')
        lbe_types = {'gene': 'gene'}
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

    def lookup_pub_curies(self, pub_id_list):
        """Return a list of curies from a list of internal chado pub_ids."""
        pub_curie_list = []
        # First, try to get curies for each pub_id.
        for pub_id in set(pub_id_list):
            try:
                pub_curie_list.append(self.bibliography[pub_id])
            except KeyError:
                pass    # This almost certainly represents an internal or obsolete publication.
        # Second, remove unattributed.
        try:
            pub_curie_list.remove('FB:unattributed')
        except ValueError:
            pass
        return pub_curie_list

    def build_cvterm_lookup(self, session):
        """Create a cvterm_id-keyed lookup of Cvterm objects."""
        self.log.info('Create a cvterm_id-keyed dict of Cvterms.')
        # First get all current pubs having an FBrf uniquename.
        filters = (
            Cvterm.is_obsolete == 0,
        )
        results = session.query(Cvterm).filter(*filters).distinct()
        cvterm_counter = 0
        for result in results:
            self.cvterm_lookup[result.cvterm_id] = result
            cvterm_counter += 1
        self.log.info(f'Found {cvterm_counter} current CV terms in chado.')
        return

    def lookup_cvterm_name(self, cvterm_id):
        """Return the name for a given chado cvterm.cvterm_id."""
        try:
            name = self.cvterm_lookup[cvterm_id].name
        except KeyError:
            self.log.error(f'The cvterm_id given is not recognized: {cvterm_id}')
            name = None
        return name

    def lookup_cvterm_cv_name(self, cvterm_id):
        """Return the CV name for a given chado cvterm.cvterm_id."""
        try:
            cv_name = self.cvterm_lookup[cvterm_id].cv.name
        except KeyError:
            self.log.error(f'The cvterm_id given is not recognized: {cvterm_id}')
            cv_name = None
        return cv_name

    def lookup_cvterm_curie(self, cvterm_id):
        """Return the curie for a given chado cvterm.cvterm_id."""
        try:
            db_name = self.cvterm_lookup[cvterm_id].dbxref.db.name
            accession = self.cvterm_lookup[cvterm_id].dbxref.accession
            curie = f'{db_name}:{accession}'
        except KeyError:
            self.log.error(f'The cvterm_id given is not recognized: {cvterm_id}')
            curie = None
        return curie

    def build_ncbi_taxon_lookup(self, session):
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
        for result in results:
            self.ncbi_taxon_lookup[result.OrganismDbxref.organism_id] = f'NCBITaxon:{result.Dbxref.accession}'
            counter += 1
        self.log.info(f'Found {counter} NCBITaxon IDs for FlyBase organisms.')
        return

    def build_feature_lookup(self, session):
        """Build a simple feature lookup."""
        # Note - depends on prior construction of self.ncbi_taxon_lookup and self.cvterm_lookup.
        self.log.info('Build a simple feature lookup.')
        # Feature types to add to the feature_lookup.
        feat_type_export = {
            'allele': True,
            'construct': True,
            'gene': True,
            'aberration': False,
            'balancer': False,
            'chem': False,
            'insertion': False,
            'seqfeat': False,
            'tool': False,
        }
        for feat_type, is_exported in feat_type_export.items():
            self.log.info(f'Looking up {feat_type} features.')
            filters = (
                Feature.uniquename.op('~')(self.regex[feat_type]),
                FeatureSynonym.is_current.is_(True),
                Cvterm.name == 'symbol',

            )
            results = session.query(Feature.feature_id, Feature.uniquename, Feature.is_obsolete,
                                      Feature.type_id, Organism.organism_id, Organism.genus,
                                      Organism.species, Feature.name, Synonym.synonym_sgml).\
                select_from(Feature).\
                join(Organism, (Organism.organism_id == Feature.organism_id)).\
                join(FeatureSynonym, (FeatureSynonym.feature_id == Feature.feature_id)).\
                join(Synonym, (Synonym.synonym_id == FeatureSynonym.synonym_id)).\
                join(Cvterm, (Cvterm.cvterm_id == Synonym.type_id)).\
                filter(*filters).\
                distinct()
            FEATURE_ID = 0
            UNIQUENAME = 1
            OBSOLETE = 2
            TYPE_ID = 3
            ORG_ID = 4
            GENUS = 5
            SPECIES = 6
            NAME = 7
            SYMBOL = 8
            counter = 0
            for result in results:
                feat_dict = {
                    'uniquename': result[UNIQUENAME],
                    'is_obsolete': result[OBSOLETE],
                    'type': self.cvterm_lookup[result[TYPE_ID]],
                    'species': f'{result[GENUS]} {result[SPECIES]}',
                    'name': result[NAME],
                    'symbol': sub_sup_sgml_to_html(result[SYMBOL]),
                    'exported': is_exported,
                }
                try:
                    feat_dict['taxon_id'] = self.ncbi_taxon_lookup[result[ORG_ID]]
                except KeyError:
                    feat_dict['taxon_id'] = 'NCBITaxon:32644'    # Unspecified taxon.
                self.feature_lookup[result[FEATURE_ID]] = feat_dict
                counter += 1
            self.log.info(f'Added {counter} {feat_type} features to the feature_lookup.')
        return

    def old_build_feature_lookup(self, session):
        """Build a simple feature lookup."""
        # Note - depends on prior construction of self.ncbi_taxon_lookup and self.cvterm_lookup.
        self.log.info('Build a simple feature lookup.')
        # Feature types to add to the feature_lookup.
        feat_type_export = {
            'allele': True,
            'construct': True,
            'gene': True,
            'aberration': False,
            'balancer': False,
            'chem': False,
            'insertion': False,
            'seqfeat': False,
            'tool': False,
        }
        for feat_type, is_exported in feat_type_export.items():
            self.log.info(f'Looking up {feat_type} features.')
            # Using conventional query rather than ORM because it's 3X faster, saving 3m.
            feature_query = f"""
                SELECT DISTINCT
                    f.feature_id,
                    f.uniquename,
                    f.is_obsolete,
                    f.type_id,
                    o.organism_id,
                    o.genus||' '||o.species,
                    f.name,
                    s.synonym_sgml
                FROM feature f
                JOIN organism o ON o.organism_id = f.organism_id
                JOIN feature_synonym fs ON fs.feature_id = f.feature_id AND fs.is_current IS TRUE
                JOIN synonym s ON s.synonym_id = fs.synonym_id
                JOIN cvterm t ON t.cvterm_id = s.type_id AND t.name = 'symbol'
                WHERE f.uniquename ~ '{self.regex[feat_type]}';
            """
            results = session.execute(feature_query).fetchall()
            FEATURE_ID = 0
            UNIQUENAME = 1
            OBSOLETE = 2
            TYPE_ID = 3
            ORG_ID = 4
            SPECIES_NAME = 5
            NAME = 6
            SYMBOL = 7
            self.log.info(f'BILLY: Process {len(results)} {feat_type} features.')
            counter = 0
            for result in results:
                feat_dict = {
                    'uniquename': result[UNIQUENAME],
                    'is_obsolete': result[OBSOLETE],
                    'type': self.cvterm_lookup[result[TYPE_ID]],
                    'species': result[SPECIES_NAME],
                    'name': result[NAME],
                    'symbol': sub_sup_sgml_to_html(result[SYMBOL]),
                    'exported': is_exported,
                }
                try:
                    feat_dict['taxon_id'] = self.ncbi_taxon_lookup[result[ORG_ID]]
                except KeyError:
                    feat_dict['taxon_id'] = 'NCBITaxon:32644'    # Unspecified taxon.
                self.feature_lookup[result[FEATURE_ID]] = feat_dict
                counter += 1
            self.log.info(f'Added {counter} {feat_type} features to the feature_lookup.')
        return

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

    def build_feature_relationship_evidence_lookup(self, session):
        """Build evidence lookup for feature_relationships."""
        self.log.info('Build evidence lookup for feature_relationships.')
        results = session.query(FeatureRelationshipPub).distinct()
        counter = 0
        fr_counter = 0
        for result in results:
            try:
                self.feat_rel_pub_lookup[result.feature_relationship_id].append(result.pub_id)
                counter += 1
            except KeyError:
                self.feat_rel_pub_lookup[result.feature_relationship_id] = [result.pub_id]
                fr_counter += 1
                counter += 1
        self.log.info(f'Found {counter} pubs supporting {fr_counter} feature_relationships.')
        return

    def lookup_feat_rel_pubs_ids(self, feature_relationship_id):
        """Return a list of pub_ids supporting a given feature_relationship."""
        try:
            pub_ids = self.feat_rel_pub_lookup[feature_relationship_id]
        except KeyError:
            pub_ids = []
        return pub_ids

    def build_allele_gene_lookup(self, session):
        """Build an allele-gene lookup dict, current features only."""
        self.log.info('Build an allele-gene lookup dict, current features only.')
        allele = aliased(Feature, name='allele')
        gene = aliased(Feature, name='gene')
        rel_type = aliased(Cvterm, name='rel_type')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(self.regex['gene']),
            rel_type.name == 'alleleof'
        )
        results = session.query(allele, gene).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.allele_gene_lookup[result.allele.feature_id] = result.gene.feature_id
            counter += 1
        self.log.info(f'Added {counter} current allele-gene relationships to the allele-gene lookup.')
        return

    def build_seqfeat_gene_lookup(self, session):
        """Build an seqfeat-gene lookup dict, current features only."""
        self.log.info('Build an seqfeat-gene lookup dict, current features only.')
        seqfeat = aliased(Feature, name='seqfeat')
        gene = aliased(Feature, name='gene')
        rel_type = aliased(Cvterm, name='rel_type')
        filters = (
            seqfeat.is_obsolete.is_(False),
            seqfeat.uniquename.op('~')(self.regex['seqfeat']),
            gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(self.regex['gene']),
            rel_type.name == 'associated_with'
        )
        results = session.query(seqfeat, gene).\
            select_from(seqfeat).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == seqfeat.feature_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.seqfeat_gene_lookup[result.seqfeat.feature_id].append(result.gene.feature_id)
                counter += 1
            except KeyError:
                self.seqfeat_gene_lookup[result.seqfeat.feature_id] = [result.gene.feature_id]
                counter += 1
        self.log.info(f'Added {counter} current seqfeat-gene relationships to the seqfeat-gene lookup.')
        return

    def build_gene_tool_lookup(self, session):
        """Build a lookup of gene-to-tool relationships for filtering reporting, current only."""
        self.log.info('Build a lookup of gene-to-tool relationships for filtering reporting.')
        gene = aliased(Feature, name='gene')
        tool = aliased(Feature, name='tool')
        filters = (
            gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(self.regex['gene']),
            tool.is_obsolete.is_(False),
            tool.uniquename.op('~')(self.regex['tool']),
            Cvterm.name == 'originates_from'
        )
        results = session.query(gene, tool).\
            select_from(gene).\
            join(FeatureRelationship, (FeatureRelationship.object_id == gene.feature_id)).\
            join(tool, (tool.feature_id == FeatureRelationship.subject_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.gene_tool_lookup[result.gene.feature_id].append(result.tool.feature_id)
                counter += 1
            except KeyError:
                self.gene_tool_lookup[result.gene.feature_id] = [result.tool.feature_id]
                counter += 1
        self.log.info(f'Found {counter} gene-to-tool relationships.')
        return

    def build_allele_class_lookup(self, session):
        """Build a lookup of allele "transgenic product class" values."""
        self.log.info('Build a lookup of allele "transgenic product class" values.')
        # Main query.
        cvterm = aliased(Cvterm, name='cvterm')
        qualifier = aliased(Cvterm, name='qualifier')
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            cvterm.is_obsolete == 0,
            qualifier.name == 'transgenic_product_class',
        )
        results = session.query(Feature, cvterm).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(cvterm, (cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(qualifier, (qualifier.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.transgenic_allele_class_lookup[result.Feature.feature_id].append(result.cvterm.name)
                counter += 1
            except KeyError:
                self.transgenic_allele_class_lookup[result.Feature.feature_id] = [result.cvterm.name]
                counter += 1
        self.log.info(f'Found {counter} transgenic product class terms for alleles.')
        # Back up. Look at mutagen terms to find old/obsolete alleles missed in the transgenic product class retrofit.
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            Cvterm.name == 'in vitro construct - RNAi',
        )
        results = session.query(Feature).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            if result.feature_id not in self.transgenic_allele_class_lookup.keys():
                self.transgenic_allele_class_lookup[result.feature_id] = ['RNAi_reagent']
                counter += 1
            elif 'RNAi_reagent' not in self.transgenic_allele_class_lookup[result.feature_id]:
                self.transgenic_allele_class_lookup[result.feature_id].append('RNAi_reagent')
                counter += 1
        self.log.info(f'Found an additional {counter} RNAi reagent alleles via mutagen terms.')
        return

    def get_general_data(self, session):
        """Get general FlyBase chado data."""
        self.log.info('GET GENERAL FLYBASE DATA FROM CHADO.')
        # self.sqlalchemy_test(session)    # For quick dev/debugging only.
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        return

    # The get_datatype_data() wrapper; sub-methods are defined and called in more specific DataHandler types.
    def get_datatype_data(self, session):
        """Get datatype-specific FlyBase data from chado."""
        self.log.info(f'GET FLYBASE {self.fb_data_type.upper()} DATA FROM CHADO.')
        return

    # The synthesize_info() wrapper; sub-methods are defined and called in more specific DataHandler types.
    def synthesize_info(self):
        """Synthesize FB info for each data object."""
        self.log.info(f'SYNTHESIZE FLYBASE {self.fb_data_type.upper()} DATA FROM CHADO.')
        return

    # The map_fb_data_to_alliance() wrapper; sub-methods are defined and called in more specific DataHandler types.
    def map_fb_data_to_alliance(self):
        """Map FB data to the Alliance LinkML object."""
        self.log.info(f'Map FlyBase "{self.fb_data_type}" data to the Alliance LinkML object for the "{self.primary_export_set}".'.upper())
        return

    # The flag_unexportable_entities() general method.
    def flag_unexportable_entities(self, input_list: list, output_set_name: str):
        """Flag entities lacking information for a required field.

        Args:
            input_list (list): A list of objects having a linkmldto attribute that is to be exported as a dict.
            output_set_name (str): The export_set_label for the list of exported dicts: e.g., 'agm_ingest_set'

        """
        self.log.info(f'Flag FlyBase data lacking information for a required field in the {output_set_name}.'.upper())
        for i in input_list:
            for attr in self.required_fields[output_set_name]:
                if attr not in i.linkmldto.__dict__.keys():
                    i.for_export = False
                    i.export_warnings.append(f'Missing "{attr}" attribute.')
                elif getattr(i.linkmldto, attr) is None:
                    i.for_export = False
                    i.export_warnings.append(f'Missing value "{attr}" attribute.')
            if i.for_export is False:
                self.log.debug(f'DO NOT EXPORT {i}: {i.export_warnings}')
        return

    # The generate_export_dict() general method.
    def generate_export_dict(self, input_list: list, output_set_name: str):
        """Generate LinkML export dict from FB data.

        Args:
            input_list (list): A list of objects with LinkMLDTO objects under the linkmldto attribute.
            output_set_name (str): The "agr_ingest_set" label for the list of exported dicts.

        """
        self.log.info(f'Generate LinkML export dict from FB data for {output_set_name}.'.upper())
        # Create the export_data_list, keyed by the agr_ingest_set name.
        self.export_data[output_set_name] = []
        for i in input_list:
            self.input_count += 1
            if i.for_export is False:
                continue
            self.export_count += 1
            if i.linkmldto.internal is True:
                self.internal_count += 1
                self.log.debug(f'Export {i} but keep internal at the Alliance: {i.internal_reasons}')
            export_agr_dict = {}
            for attr in self.output_fields[output_set_name]:
                if getattr(i.linkmldto, attr) is not None and getattr(i.linkmldto, attr) != []:
                    export_agr_dict[attr] = getattr(i.linkmldto, attr)
            self.export_data[output_set_name].append(export_agr_dict)
        public_count = self.export_count - self.internal_count
        self.log.info(f'SUMMARY FOR EXPORT OF {output_set_name}'.upper())
        self.log.info(f'Exported {self.export_count} of {self.input_count} entities.')
        self.log.info(f'{public_count} of {self.export_count} exported entities are PUBLIC.')
        self.log.info(f'{self.internal_count} of {self.export_count} exported entities are INTERNAL.')
        return

    # The query_chado() wrapper that runs sub-methods - same order of steps for every DataHandler type.
    def query_chado(self, session):
        """Wrapper that runs all methods within an SQLAlchemy session."""
        self.log.info('Run main query_chado() handler method.'.upper())
        self.get_general_data(session)
        self.get_datatype_data(session)
        self.synthesize_info()
        self.map_fb_data_to_alliance()
        self.flag_unexportable_entities(self.fb_data_entities.values(), self.primary_export_set)
        self.generate_export_dict(self.fb_data_entities.values(), self.primary_export_set)
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

    # Elaborate on get_general_data() for the PrimaryEntityHandler.
    def get_general_data(self, session):
        """Extend the method for the PrimaryEntityHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() for the PrimaryEntityHandler; sub-methods might only be used in some more specific Datahandlers.
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
        associated_data_types = ['pubs', 'synonyms', 'dbxrefs', 'props', 'cvterms']
        chado_type = self.main_chado_entity_types[self.fb_data_type]
        self.log.info(f'Get associated data for {self.fb_data_type} data entities from {chado_type}-related chado tables.')
        main_pkey_name = self.chado_tables['primary_key'][chado_type]
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
        self.get_entities(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        # self.get_entity_props(session)
        # self.get_entity_cvterms(session)
        # self.get_entity_associated_data(session)
        # self.get_entity_prop_pubs(session)
        # self.get_entity_cvtermprops(session)
        self.get_entity_timestamps(session)
        return

    # Elaborate on synthesize_info() for the PrimaryEntityHandler; sub-methods might only be used in some more specific DataHandlers.
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

    def synthesize_info(self):
        """Extend the method for the PrimaryEntityHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_props()
        self.synthesize_pubs()
        return

    # Elaborate on map_fb_data_to_alliance() for the PrimaryEntityHandler; sub-methods might only be used in some more specific DataHandlers.
    def map_data_provider_dto(self):
        """Return the DataProviderDTO for the FB data entity."""
        # Note - this method is depends on previous determination of fb_data_entity.curr_fb_symbol by map_synonyms(), if applicable.
        self.log.info('Map data provider to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.curr_fb_symbol:
                display_name = fb_data_entity.curr_fb_symbol
            else:
                display_name = fb_data_entity.name
            dp_xref = datatypes.CrossReferenceDTO('FB', f'FB:{fb_data_entity.uniquename}', self.fb_data_type, display_name).dict_export()
            fb_data_entity.linkmldto.data_provider_dto = datatypes.DataProviderDTO(dp_xref).dict_export()
        return

    def map_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for a FlyBase entity."""
        self.log.info('Map secondary IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            secondary_id_dtos = []
            for secondary_id in fb_data_entity.alt_fb_ids:
                sec_dto = datatypes.SecondaryIdSlotAnnotationDTO(secondary_id, []).dict_export()
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
                xref_dto = datatypes.CrossReferenceDTO(prefix, curie, page_area, display_name).dict_export()
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
                name_dto = datatypes.NameSlotAnnotationDTO(syno_dict['name_type_name'], syno_dict['format_text'],
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
                    if fb_data_entity.chado_obj.is_obsolete is False and fb_data_entity.organism_abbr == 'Dmel':
                        linkml_synonym_bins['systematic_name_bin'].append(name_dto)
            # Review the linkml_synonym_bins for each fb_data_entity.
            # 1. Symbol.
            if len(linkml_synonym_bins['symbol_bin']) == 0:
                self.log.warning(f'No current symbols found for {fb_data_entity}: create a generic one.')
                generic_symbol_dto = datatypes.NameSlotAnnotationDTO('nomenclature_symbol', fb_data_entity.name, fb_data_entity.name, []).dict_export()
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
                sys_name_dto = datatypes.NameSlotAnnotationDTO('systematic_name', fb_data_entity.curr_anno_id, fb_data_entity.curr_anno_id, []).dict_export()
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

    def get_general_data(self, session):
        """Extend the method for the FeatureHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() for the FeatureHandler; sub-methods might only be used in some more specific Datahandlers.
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
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].fb_anno_dbxrefs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} annotation IDs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} annotation IDs for {self.fb_data_type} entities.')
        return

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
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].chr_flocs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} chromosomal featurelocs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} chromosomal featurelocs for {self.fb_data_type} entities.')
        return

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
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Featureprop).\
            select_from(Feature).\
            join(feat_type, (feat_type.cvterm_id == Feature.type_id)).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        return results

    def get_entity_sbj_feat_rel_by_type(self, session, slot_name, **kwargs):
        """Get a list of FeatureRelationship objects for primary feature entities (as subject) by type.

        Args:
            session (Session): SQLAlchemy session for db queries.
            slot_name (str): The name of the FB data entity attribute name to which results are appended.

        Keyword Args:
            rel_type (str): The CV term name for the feature_relationship of interest. If none given, any rel_type allowed.
            obj_type (list): A list of CV terms for the object feature types. If none given, any object feature type allowed.
            obj_regex (str): The regex for the object feature uniquename. If none given, any object uniquename allowed.

        """
        self.log.info(f'Add feature_relationships to "{slot_name}" with these criteria: {kwargs}')
        subject = aliased(Feature, name='subject')
        object = aliased(Feature, name='object')
        rel_type = aliased(Cvterm, name='rel_type')
        obj_type = aliased(Cvterm, name='obj_type')
        filters = (
            subject.feature_id.in_(self.fb_data_entities.keys()),
        )
        try:
            filters += (rel_type.name == kwargs['rel_type'], )
        except KeyError:
            pass
        try:
            filters += (obj_type.name == kwargs['obj_type'], )
        except KeyError:
            pass
        try:
            filters += (object.uniquename.op('~')(kwargs['obj_regex']), )
        except KeyError:
            pass
        results = session.query(FeatureRelationship).\
            select_from(subject).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == subject.feature_id)).\
            join(object, (object.feature_id == FeatureRelationship.object_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            join(obj_type, (obj_type.cvterm_id == object.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.fb_data_entities[result.subject_id].__dict__[slot_name].append(result)
            counter += 1
        self.log.info(f'Added {counter} feature_relationship results to "{slot_name}" list.')
        return

    def get_entity_obj_feat_rel_by_type(self, session, slot_name, **kwargs):
        """Return a list of FeatureRelationship objects for handler's primary feature entities (as object) by type.

        Args:
            session (Session): SQLAlchemy session for db queries.
            slot_name (str): The name of the FB data entity attribute name to which results are appended.

        Keyword Args:
            rel_type (str): The CV term name for the feature_relationship of interest. If none given, any rel_type allowed.
            sbj_type (str): The CV term for the subject feature types. If none given, any subject feature type allowed.
            sbj_regex (str): The regex for the subject feature uniquename. If none given, any subject uniquename allowed.

        """
        self.log.info(f'Add feature_relationships to "{slot_name}" with these criteria: {kwargs}')
        subject = aliased(Feature, name='subject')
        object = aliased(Feature, name='object')
        rel_type = aliased(Cvterm, name='rel_type')
        sbj_type = aliased(Cvterm, name='sbj_type')
        filters = (
            object.feature_id.in_(self.fb_data_entities.keys()),
        )
        try:
            filters += (rel_type.name == kwargs['rel_type'], )
        except KeyError:
            pass
        try:
            filters += (sbj_type.name == kwargs['sbj_type'], )
        except KeyError:
            pass
        try:
            filters += (subject.uniquename.op('~')(kwargs['sbj_regex']), )
        except KeyError:
            pass
        results = session.query(FeatureRelationship).\
            select_from(subject).\
            join(sbj_type, (sbj_type.cvterm_id == subject.type_id)).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == subject.feature_id)).\
            join(object, (object.feature_id == FeatureRelationship.object_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.fb_data_entities[result.object_id].__dict__[slot_name].append(result)
            counter += 1
        self.log.info(f'Added {counter} feature_relationship results to "{slot_name}" list.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the FeatureHandler."""
        super().get_datatype_data(session)
        return

    # Elaborate on synthesize_info() for the FeatureHandler; sub-methods might only be used in some more specific DataHandlers.
    def synthesize_anno_ids(self):
        """Synthesize annotation IDs."""
        self.log.info('Synthesize annotation IDs.')
        for fb_data_entity in self.fb_data_entities.values():
            current_anno_ids = []
            alt_anno_ids = []
            for xref in fb_data_entity.fb_anno_dbxrefs:
                if xref.is_current is True:
                    current_anno_ids.append(xref.dbxref.accession)
                else:
                    alt_anno_ids.append(xref.dbxref.accession)
            # Get the one current annotation ID.
            if len(current_anno_ids) == 1:
                fb_data_entity.curr_anno_id = current_anno_ids[0]
            elif len(current_anno_ids) > 1:
                self.log.warning(f'{fb_data_entity} has {len(current_anno_ids)} current annotations IDs.')
                fb_data_entity.alt_anno_ids.extend(current_anno_ids)
            # Record old annotation IDs.
            fb_data_entity.alt_anno_ids.extend(alt_anno_ids)
        return

    def synthesize_info(self):
        """Extend the method for the FeatureHandler."""
        super().synthesize_info()
        return

    # Elaborate on map_fb_data_to_alliance() for the FeatureHandler; sub-methods might only be used in some more specific DataHandlers.
    def map_anno_ids_to_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for annotation IDs."""
        self.log.info('Map annotation IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            anno_ids = []
            if fb_data_entity.curr_anno_id:
                anno_ids.append(fb_data_entity.curr_anno_id)
            if fb_data_entity.alt_anno_ids:
                anno_ids.extend(fb_data_entity.alt_anno_ids)
            anno_secondary_id_dtos = []
            for anno_id in anno_ids:
                sec_dto = datatypes.SecondaryIdSlotAnnotationDTO(f'FB:{anno_id}', []).dict_export()
                anno_secondary_id_dtos.append(sec_dto)
            curr_sec_id_dtos = getattr(fb_data_entity.linkmldto, slot_name)
            curr_sec_id_dtos.extend(anno_secondary_id_dtos)
        return anno_secondary_id_dtos

    def map_fb_data_to_alliance(self):
        """Extend the method for the FeatureHandler."""
        super().map_fb_data_to_alliance()
        return


class ConstructHandler(FeatureHandler):
    """This object gets, synthesizes and filters construct data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the ConstructHandler object."""
        super().__init__(log, fb_data_type, testing)
        # Additional set for export added to the handler.
        self.construct_associations = []            # Will be a list of FBEntity objects, map to ConstructGenomicEntityAssociationDTO.
        # Lookups needed.
        self.feat_rel_pub_lookup = {}               # Will be feature_relationship_id-keyed lists of supporting pub_ids.
        self.allele_gene_lookup = {}                # Will be allele feature_id-keyed of a single gene feature_id per allele.
        self.seqfeat_gene_lookup = {}               # Will be seqfeat feature_id-keyed of a lists of gene feature_ids.
        self.transgenic_allele_class_lookup = {}    # Will be an allele feature_id-keyed list of "transgenic product class" CV terms.
        self.gene_tool_lookup = {}                  # Will be gene feature_id-keyed lists of related FBto tools.

    test_set = {
        'FBtp0008631': 'P{UAS-wg.H.T:HA1}',                       # Expresses FBgn wg, regulated by FBto UASt.
        'FBtp0010648': 'P{wg.FRT.B}',                             # Expresses FBgn wg, regulated by FBgn sev, has FBto FRT.
        'FBtp0145675': 'PBac{UAS-hHTT.ex1.Q97.S13D.mCherry}',     # Expresses FBgn Hsap\HTT, regulated by FBto UAS, tagged with FBto mCherry.
        'FBtp0000074': 'P{ftzG}',                                 # Expresses FBgn ftz, regulated by FBgn ftz.
        'FBtp0000326': 'P{SEV5}',                                 # Expresses FBgn sev, tagged with FBto MYC.
        'FBtp0161516': 'P{lush-GAL4.3}',                          # Expresses FBto GAL4, regulated by FBgn lush.
        'FBtp0057873': 'P{GMR16C10-GAL4}',                        # Expresses FBto GAL4, regulated by GMR16C10 (related to two genes, Brf and lute).
        'FBtp0032215': 'P{GD5007}',                               # Targets FBgn wg, regulated by FBto UASt.
        'FBtp0031452': 'P{GD4157}',                               # Targets FBgn lbe, regulated by FBto UASt.
        'FBtp0145396': 'P{TOE.GS00055}',                          # Targets FBgn wg, regulated by FBto UASt.
        'FBtp0145394': 'P{TKO.GS00469}',                          # Targets FBgn Alp9, Alp10, regulated by FBto UASt.
        'FBtp0000352': 'P{GawB}',                                 # Expresses FBto GAL4, FBgn Scer\GAL4. Report both?
        'FBtp0161256': 'PBac{UAS-G-CEPIA1::TM-2A-TagRFP::TM}',    # 2 FBal; expresses FBto G-CEPIA1, RFP; expresses FBgn Equa\eqFP578, GFP; regulated by UAS.
        'FBtp0051705': 'M{MtnBcDNA-MtnDcDNA.EGFP}',               # has_reg_region MtnB.
        'FBtp0080088': 'P{UAS-Brainbow}',                         # Expresses EBFP2, EGFP, mKO2, has_reg_region UAS; tagged_with HA, MYC, V5; carries lox.
        'FBtp0083738': 'P{GR}',                                   # Is regulated_by FBgn Act5C.
    }
    # Elaborate on export filters for ConstructHandler.
    required_fields = {
        'construct_ingest_set': [
            'construct_symbol_dto',
            'data_provider_dto',
            'internal',
            'mod_entity_id',
        ],
        'construct_genomic_entity_association_ingest_set': [
            'construct_identifier',
            'genomic_entity_curie',
            'genomic_entity_relation_name',
            'internal',
        ],
    }
    output_fields = {
        'construct_ingest_set': [
            'construct_component_dtos',
            'construct_full_name_dto',
            'construct_symbol_dto',
            'construct_synonym_dtos',
            'created_by_curie',
            'data_provider_dto',
            'date_created',
            'date_updated',
            'internal',
            'mod_entity_id',
            'mod_internal_id',
            'obsolete',
            'reference_curies',
            'secondary_identifiers',
            'updated_by_curie',
        ],
        'construct_genomic_entity_association_ingest_set': [
            'construct_identifier',
            'created_by_curie',
            'date_created',
            'date_updated',
            'evidence_curies',
            'genomic_entity_curie',
            'genomic_entity_relation_name',
            'internal',
            'obsolete',
            'updated_by_curie',
        ],
    }
    fb_agr_db_dict = {}

    # Elaborate on get_general_data() for the ConstructHandler.
    def get_general_data(self, session):
        """Extend the method for the ConstructHandler."""
        super().get_general_data(session)
        self.build_feature_lookup(session)
        self.build_feature_relationship_evidence_lookup(session)
        self.build_allele_class_lookup(session)
        self.build_seqfeat_gene_lookup(session)
        self.build_gene_tool_lookup(session)
        return

    # Elaborate on get_datatype_data() for the ConstructHandler.
    def get_construct_alleles(self, session):
        """Get allele(s) to which constructs belong."""
        self.log.info('Get allele(s) to which constructs belong.')
        self.get_entity_obj_feat_rel_by_type(session, 'parent_allele_rels', rel_type='associated_with', sbj_type='allele', sbj_regex=self.regex['allele'])
        return

    def get_construct_encoded_tools(self, session):
        """Get directly related encoded FBto/FBsf objects for the construct."""
        self.log.info('Get directly related encoded FBto/FBsf objects for the construct.')
        self.get_entity_sbj_feat_rel_by_type(session, 'encodes_tool_rels', rel_type='encodes_tool', obj_regex=self.regex['fb_uniquename'])
        return

    def get_construct_reg_regions(self, session):
        """Get directly related regulatory FBgn/FBto/FBsf objects for the construct."""
        self.log.info('Get directly related regulatory FBgn/FBto/FBsf objects for the construct.')
        self.get_entity_sbj_feat_rel_by_type(session, 'reg_region_rels', rel_type='has_reg_region', obj_regex=self.regex['fb_uniquename'])
        return

    def get_construct_reg_regions_old(self, session):
        """Get directly related regulatory_region FBsf objects, old type of association."""
        self.log.info('Get directly related regulatory_region FBsf objects, old type of association.')
        self.get_entity_obj_feat_rel_by_type(session, 'seqfeat_rels', rel_type='associated_with', sbj_type='regulatory_region', sbj_regex=self.regex['seqfeat'])
        return

    def get_allele_encoded_tools(self, session):
        """Get encoded FBto/FBsf objects for the constructs via alleles."""
        self.log.info('Get encoded FBto/FBsf objects for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        component = aliased(Feature, name='component')
        filters = (
            # allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            # component.is_obsolete.is_(False),
            component.uniquename.op('~')(self.regex['fb_uniquename']),
            Cvterm.name == 'encodes_tool',
        )
        results = session.query(FeatureRelationship).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(component, (component.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        # Create allele feature_id-keyed lists of allele-component "encodes_tool" FeatureRelationship objects.
        al_comp_dict = {}
        counter = 0
        for result in results:
            try:
                al_comp_dict[result.subject_id].append(result)
                counter += 1
            except KeyError:
                al_comp_dict[result.subject_id] = [result]
                counter += 1
        self.log.info(f'Found {counter} allele-to-component "encodes_tool" relationships.')
        self.log.info('Now propagate these "encodes_tool" relationships to constructs.')
        counter = 0
        for construct in self.fb_data_entities.values():
            for allele_rel in construct.parent_allele_rels:
                allele_id = allele_rel.subject_id
                try:
                    construct.al_encodes_tool_rels.extend(al_comp_dict[allele_id])
                    counter += len(construct.al_encodes_tool_rels)
                except KeyError:
                    pass
        self.log.info(f'Propagated {counter} allele-to-component "encodes_tool" relationships to related constructs.')
        return

    def get_allele_reg_regions(self, session):
        """Get regulatory FBto/FBsf/FBgn objects for the constructs via alleles."""
        self.log.info('Get regulatory FBto/FBsf/FBgn objects for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        component = aliased(Feature, name='component')
        filters = (
            # allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            # component.is_obsolete.is_(False),
            component.uniquename.op('~')(self.regex['fb_uniquename']),
            Cvterm.name == 'has_reg_region',
        )
        results = session.query(FeatureRelationship).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(component, (component.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        # Create allele feature_id-keyed lists of allele-component "has_reg_region" FeatureRelationship objects.
        al_comp_dict = {}
        counter = 0
        for result in results:
            try:
                al_comp_dict[result.subject_id].append(result)
                counter += 1
            except KeyError:
                al_comp_dict[result.subject_id] = [result]
                counter += 1
        self.log.info(f'Found {counter} allele-to-component "has_reg_region" relationships.')
        self.log.info('Now propagate "has_reg_region" relationships to constructs.')
        counter = 0
        for construct in self.fb_data_entities.values():
            for allele_rel in construct.parent_allele_rels:
                allele_id = allele_rel.subject_id
                try:
                    construct.al_reg_region_rels.extend(al_comp_dict[allele_id])
                    counter += len(construct.al_reg_region_rels)
                except KeyError:
                    pass
        self.log.info(f'Propagated {counter} allele-to-component "has_reg_region" relationships to related constructs.')
        return

    def get_allele_genes(self, session):
        """Get genes for the constructs via alleles."""
        self.log.info('Get genes for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        gene = aliased(Feature, name='gene')
        filters = (
            # allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            # gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(self.regex['gene']),
            Cvterm.name == 'alleleof',
        )
        results = session.query(FeatureRelationship).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        # Create allele feature_id-keyed lists of allele-gene "alleleof" FeatureRelationship objects.
        al_comp_dict = {}
        counter = 0
        for result in results:
            try:
                al_comp_dict[result.subject_id].append(result)
                counter += 1
            except KeyError:
                al_comp_dict[result.subject_id] = [result]
                counter += 1
        self.log.info(f'Found {counter} allele-to-gene "alleleof" relationships.')
        self.log.info('Now propagate "alleleof" relationships to constructs.')
        counter = 0
        for construct in self.fb_data_entities.values():
            for allele_rel in construct.parent_allele_rels:
                allele_id = allele_rel.subject_id
                try:
                    construct.al_genes.extend(al_comp_dict[allele_id])
                    counter += len(construct.al_genes)
                except KeyError:
                    pass
        self.log.info(f'Propagated {counter} allele-to-gene "alleleof" relationships to related constructs.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the ConstructHandler."""
        super().get_datatype_data(session)
        # self.get_construct_alleles(session)    # BOB: SUPER SLOW STEP - FIX IT
        self.get_construct_encoded_tools(session)
        self.get_construct_reg_regions(session)
        self.get_allele_encoded_tools(session)
        self.get_allele_reg_regions(session)
        self.get_construct_reg_regions_old(session)
        self.get_allele_genes(session)
        return

    # Elaborate on synthesize_info() for the ConstructHandler.
    def synthesize_encoded_tools(self):
        """Synthesize encoded components."""
        self.log.info('Synthesize encoded components.')
        counter = 0
        for construct in self.fb_data_entities.values():
            # self.log.debug(f'Assess encoded tools for {construct}.')
            # Direct encodes_tool relationships.
            for rel in construct.encodes_tool_rels:
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pubs_ids(rel.feature_relationship_id)
                try:
                    construct.expressed_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.expressed_features[component_id] = pub_ids
            # direct_count = len(construct.expressed_features.keys())
            # self.log.debug(f'For {construct}, found {direct_count} encoded tools via direct relationships.')
            # Indirect encodes_tool relationships.
            for rel in construct.al_encodes_tool_rels:
                allele_id = rel.subject_id
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pubs_ids(rel.feature_relationship_id)
                try:
                    construct.expressed_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.expressed_features[component_id] = pub_ids
                # Fold in pubs supporting the construct-allele relationship.
                for al_con_rel in construct.parent_allele_rels:
                    if al_con_rel.subject_id == allele_id:
                        al_con_pub_ids = self.lookup_feat_rel_pubs_ids(al_con_rel.feature_relationship_id)
                        construct.expressed_features[component_id].extend(al_con_pub_ids)
            # indirect_count = len(construct.expressed_features.keys()) - direct_count
            # self.log.debug(f'For {construct}, found {indirect_count} encoded tools via indirect allele relationships.')
            counter += len(construct.expressed_features.keys())
        self.log.info(f'Found {counter} encoded tools for constructs via direct and indirect allele relationships.')
        return

    def synthesize_component_genes(self):
        """Synthesize component genes."""
        self.log.info('Synthesize component genes.')
        all_expressed_gene_counter = 0
        all_targeted_gene_counter = 0
        for construct in self.fb_data_entities.values():
            this_expressed_gene_counter = 0
            this_targeted_gene_counter = 0
            # self.log.debug(f'Assess component genes for {construct}.')
            for rel in construct.al_genes:
                allele_id = rel.subject_id
                gene_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pubs_ids(rel.feature_relationship_id)
                # Slot for gene_id depends on the allele class.
                try:
                    if set(self.transgenic_allele_class_lookup[allele_id]).intersection({'RNAi_reagent', 'sgRNA', 'antisense'}):
                        gene_slot = getattr(construct, 'targeted_features')
                        this_targeted_gene_counter += 1
                    else:
                        gene_slot = getattr(construct, 'expressed_features')
                        this_expressed_gene_counter += 1
                except KeyError:
                    gene_slot = getattr(construct, 'expressed_features')
                    this_expressed_gene_counter += 1
                try:
                    gene_slot[gene_id].extend(pub_ids)
                except KeyError:
                    gene_slot[gene_id] = pub_ids
                # Fold in pubs supporting the construct-allele relationship.
                for al_con_rel in construct.parent_allele_rels:
                    if al_con_rel.subject_id == allele_id:
                        al_con_pub_ids = self.lookup_feat_rel_pubs_ids(al_con_rel.feature_relationship_id)
                        gene_slot[gene_id].extend(al_con_pub_ids)
            # self.log.debug(f'For {construct}, found {this_expressed_gene_counter} expressed genes and {this_targeted_gene_counter} targeted genes.')
            all_expressed_gene_counter += this_expressed_gene_counter
            all_targeted_gene_counter += this_targeted_gene_counter
        self.log.info(f'Found {all_expressed_gene_counter} expressed genes and {all_targeted_gene_counter} targeted genes for constructs.')
        return

    def synthesize_reg_regions(self):
        """Synthesize construct reg_region components."""
        self.log.info('Synthesize construct reg_region components.')
        counter = 0
        for construct in self.fb_data_entities.values():
            # self.log.debug(f'Assess reg_regions for {construct}.')
            # Direct has_reg_region relationships.
            for rel in construct.reg_region_rels:
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pubs_ids(rel.feature_relationship_id)
                try:
                    construct.regulating_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.regulating_features[component_id] = pub_ids
            # direct_count = len(construct.regulating_features.keys())
            # self.log.debug(f'For {construct}, found {direct_count} reg_regions via direct relationships.')
            # Direct seqfeat relationships, old_style.
            for rel in construct.seqfeat_rels:
                component_id = rel.subject_id
                pub_ids = self.lookup_feat_rel_pubs_ids(rel.feature_relationship_id)
                try:
                    construct.regulating_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.regulating_features[component_id] = pub_ids
            # direct_count_old = len(construct.regulating_features.keys()) - direct_count
            # self.log.debug(f'For {construct}, found {direct_count_old} reg_regions via direct relationships, old style.')
            # Indirect has_reg_region relationships.
            for rel in construct.al_reg_region_rels:
                allele_id = rel.subject_id
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pubs_ids(rel.feature_relationship_id)
                try:
                    construct.regulating_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.regulating_features[component_id] = pub_ids
                # Fold in pubs supporting the construct-allele relationship.
                for al_con_rel in construct.parent_allele_rels:
                    if al_con_rel.subject_id == allele_id:
                        al_con_pub_ids = self.lookup_feat_rel_pubs_ids(al_con_rel.feature_relationship_id)
                        construct.regulating_features[component_id].extend(al_con_pub_ids)
            # indirect_count = len(construct.regulating_features.keys()) - direct_count - direct_count_old
            # self.log.debug(f'For {construct}, found {indirect_count} reg_regions tools via indirect allele relationships.')
            # Indirect relationships to genes via seqfeats.
            for component_id in list(construct.regulating_features.keys()):
                pubs_ids = construct.regulating_features[component_id]
                uniquename = self.feature_lookup[component_id]['uniquename']
                feat_type = self.feature_lookup[component_id]['type']
                if uniquename.startswith('FBsf') and feat_type in ['region', 'regulatory_region']:
                    try:
                        gene_ids = self.seqfeat_gene_lookup[component_id]
                    except KeyError:
                        continue
                    for gene_id in gene_ids:
                        try:
                            construct.regulating_features[gene_id].extend(pub_ids)
                        except KeyError:
                            construct.regulating_features[gene_id] = pubs_ids
            # genes_via_seqfeat_count = len(construct.regulating_features.keys()) - direct_count - direct_count_old - indirect_count
            # self.log.debug(f'For {construct}, found an additional {genes_via_seqfeat_count} genes related to seqfeat reg_regions.')
            counter += len(construct.regulating_features.keys())
        self.log.info(f'Found {counter} reg_regions for constructs via direct and indirect allele relationships.')
        return

    def synthesize_redundant_tool_genes(self):
        """For constructs in which a gene and related tool are related, sort out redundant genes."""
        self.log.info('For constructs in which a gene and related tool are related, sort out redundant genes.')
        slot_names = {
            'expressed_features': 'expressed_tool_genes',
            'regulating_features': 'regulating_tool_genes',
        }
        for slot_name, tool_gene_slot_name in slot_names.items():
            self.log.info(f'Prune {slot_name} for constructs.')
            counter = 0
            for construct in self.fb_data_entities.values():
                pruning_list = []
                slot_bin = getattr(construct, slot_name)
                for feature_id in slot_bin.keys():
                    try:
                        all_related_tool_ids = set(self.gene_tool_lookup[feature_id])
                        tool_overlap = all_related_tool_ids.intersection(set(slot_bin.keys()))
                        # self.log.debug(f'For {construct}, {self.feature_lookup[feature_id]["name"]} has {len(tool_overlap)} redundantly associated tools.')
                        if tool_overlap:
                            pruning_list.append(feature_id)
                            pruned_gene = f'{self.feature_lookup[feature_id]["name"]} ({self.feature_lookup[feature_id]["uniquename"]})'
                            tool_overlap_str = '|'.join([f'{self.feature_lookup[i]["name"]} ({self.feature_lookup[i]["uniquename"]})' for i in tool_overlap])
                            self.log.debug(f'For {construct}, sort out {pruned_gene} since related tools are more informative: {tool_overlap_str}')
                    except KeyError:
                        pass
                for gene_id in pruning_list:
                    tool_gene_bin = getattr(construct, tool_gene_slot_name)
                    tool_gene_bin.append(gene_id)
                    counter += 1
            self.log.info(f'Pruned {counter} genes from construct {slot_name} that are better represented as tools.')
        return

    def synthesize_construct_genomic_entity_associations(self):
        """Synthesize construct-genomic entity associations."""
        self.log.info('Synthesize construct-genomic entity associations.')
        slot_rel_types = {
            'expressed_features': 'expresses',
            'targeted_features': 'targets',
            'regulating_features': 'is_regulated_by',
        }
        for feature_slot_name, rel_type in slot_rel_types.items():
            self.log.info(f'Sort out Alliance genomic entities from "{feature_slot_name}" to "{rel_type}" associations.')
            rel_type = slot_rel_types[feature_slot_name]
            counter = 0
            for construct in self.fb_data_entities.values():
                component_slot = getattr(construct, feature_slot_name)
                for feature_id, pub_ids in component_slot.items():
                    if self.feature_lookup[feature_id]['type'] != 'gene' or not self.feature_lookup[feature_id]['uniquename'].startswith('FBgn'):
                        continue
                    rel_dict = {
                        'construct_curie': f'FB:{construct.uniquename}',
                        'rel_type': rel_type,
                        'genomic_entity_curie': f'FB:{self.feature_lookup[feature_id]["uniquename"]}',
                        'pub_curies': self.lookup_pub_curies(pub_ids),
                        'obsolete': False,
                        'internal': False,
                    }
                    # If either component in the relationship is obsolete, set the relationship to obsolete.
                    if construct.chado_obj.is_obsolete is True:
                        rel_dict['obsolete'] = True
                        rel_dict['internal'] = True
                    if self.feature_lookup[feature_id]['is_obsolete'] is True:
                        rel_dict['obsolete'] = True
                        rel_dict['internal'] = True
                    feat_rel = datatypes.FBEntity()
                    feat_rel.rel_dict = rel_dict
                    feat_rel.entity_desc = f'{rel_dict["construct_curie"]}_{rel_dict["rel_type"]}_{rel_dict["genomic_entity_curie"]}'
                    self.construct_associations.append(feat_rel)
                    counter += 1
            self.log.info(f'Synthesized {counter} construct-gene associations.')
        return

    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_encoded_tools()
        self.synthesize_component_genes()
        self.synthesize_reg_regions()
        self.synthesize_redundant_tool_genes()
        self.synthesize_construct_genomic_entity_associations()
        return

    # Elaborate on map_fb_data_to_alliance() for the ConstructHandler.
    def map_construct_basic(self):
        """Map basic FlyBase construct data to the Alliance LinkML object."""
        self.log.info('Map basic construct info to Alliance object.')
        for construct in self.fb_data_entities.values():
            agr_construct = datatypes.ConstructDTO()
            agr_construct.obsolete = construct.chado_obj.is_obsolete
            agr_construct.mod_entity_id = f'FB:{construct.uniquename}'
            construct.linkmldto = agr_construct
        return

    def map_construct_components(self):
        """Map current construct components to the Alliance LinkML object."""
        self.log.info('Map current construct components to the Alliance LinkML object.')
        component_slots = {
            'expressed_features': 'expresses',
            'targeted_features': 'targets',
            'regulating_features': 'is_regulated_by',
        }
        counter = 0
        for construct in self.fb_data_entities.values():
            for slot_name, rel_type in component_slots.items():
                slot_bin = getattr(construct, slot_name)
                for feature_id, pub_ids in slot_bin.items():
                    # Do not report obsolete components.
                    if self.feature_lookup[feature_id]['is_obsolete'] is True:
                        continue
                    # Do not report genes that are better reported as tools.
                    elif slot_name == 'expressed_features' and feature_id in construct.expressed_tool_genes:
                        continue
                    elif slot_name == 'regulating_features' and feature_id in construct.regulating_tool_genes:
                        continue
                    symbol = self.feature_lookup[feature_id]['symbol']
                    pubs = self.lookup_pub_curies(pub_ids)
                    taxon_text = self.feature_lookup[feature_id]['species']
                    taxon_curie = self.feature_lookup[feature_id]['taxon_id']
                    component_dto = datatypes.ConstructComponentSlotAnnotationDTO(rel_type, symbol, taxon_curie, taxon_text, pubs).dict_export()
                    construct.linkmldto.construct_component_dtos.append(component_dto)
                    counter += 1
        self.log.info(f'Mapped construct components to {counter} ConstructcomponentDTOs.')
        return

    def map_construct_genomic_associations(self):
        """Map current construct relations to the Alliance LinkML object."""
        self.log.info('Map current construct relations to the Alliance LinkML object.')
        counter = 0
        for cons_asso in self.construct_associations:
            rel_dto = datatypes.ConstructGenomicEntityAssociationDTO(cons_asso.rel_dict['construct_curie'], cons_asso.rel_dict['rel_type'],
                                                                     cons_asso.rel_dict['genomic_entity_curie'], cons_asso.rel_dict['pub_curies'])
            rel_dto.obsolete = cons_asso.rel_dict['obsolete']
            rel_dto.internal = cons_asso.rel_dict['internal']
            cons_asso.linkmldto = rel_dto
            counter += 1
        self.log.info(f'Mapped construct_relationships to {counter} ConstructGenomicEntityAssociationDTOs.')
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the ConstructHandler."""
        super().map_fb_data_to_alliance()
        self.map_construct_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_pubs()
        self.map_timestamps()
        # Not use self.map_secondary_ids() here because for reagents, we report only strings, not SecondaryIdSlotAnnotationDTOs.
        for construct in self.fb_data_entities.values():
            construct.linkmldto.secondary_identifiers = construct.alt_fb_ids
        self.map_construct_components()
        self.flag_internal_fb_entities('fb_data_entities')
        self.map_construct_genomic_associations()
        self.flag_internal_fb_entities('construct_associations')
        return

    # Elaborate on query_chado().
    def query_chado(self, session):
        """Elaborate on query_chado method for the ConstructHandler."""
        super().query_chado(session)
        self.flag_unexportable_entities(self.construct_associations, 'construct_genomic_entity_association_ingest_set')
        self.generate_export_dict(self.construct_associations, 'construct_genomic_entity_association_ingest_set')
        return


class GeneHandler(FeatureHandler):
    """This object gets, synthesizes and filters gene data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the GeneHandler object."""
        super().__init__(log, fb_data_type, testing)
        self.pthr_dict = {}    # Will be an 1:1 FBgn_ID-PTHR xref dict.

    test_set = {
        'FBgn0284084': 'wg',                # Current annotated nuclear protein_coding gene.
        'FBgn0004009': 'wg',                # Obsolete annotated nuclear protein_coding gene.
        'FBgn0013687': 'mt:ori',            # Current localized but unannotated mitochondrial gene.
        'FBgn0013678': 'mt:Cyt-b',          # Current annotated mitochondrial protein_coding gene.
        'FBgn0019661': 'lncRNA:roX1',       # Current annotated nuclear ncRNA gene.
        'FBgn0262451': 'mir-ban',           # Current annotated nuclear miRNA gene.
        'FBgn0034365': 'CG5335',            # Current annotated gene with CG symbol.
        'FBgn0003884': 'alphaTub84B',       # Current annotated gene with non-ASCII char in symbol.
        'FBgn0263477': 'scaRNA:PsiU1-6',    # Current annotated gene needs systematic synonym dto.
        'FBgn0030179': 'CG12094',           # Obsolete unannotated gene, should not get systematic name but needs symbol.
        'FBgn0108495': 'Dere\\GG16260',     # Current unannotated non-Dmel with systematic name.
        'FBgn0031087': 'CG12656',           # Current withdrawn gene.
        'FBgn0000154': 'Bar',               # Current unannotated gene.
        'FBgn0001200': 'His4',              # Current unannotated gene family.
        'FBgn0087003': 'tal',               # Current unannotated oddball.
        'FBgn0015267': 'Mmus\\Abl1',        # Current mouse gene with MGI xref.
    }
    # Elaborate on export filters for GeneHandler.
    required_fields = {
        'gene_ingest_set': [
            'curie',
            'data_provider_dto',
            'gene_symbol_dto',
            'internal',
            'taxon_curie',
        ],
    }
    output_fields = {
        'gene_ingest_set': [
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
            # 'related_notes',    # Not present in GeneDTO.
            'taxon_curie',
            'updated_by_curie',
        ],
    }
    fb_agr_db_dict = {
        'EntrezGene': 'NCBI_Gene',
        'RNAcentral': 'RNAcentral',
        # 'UniProt/GCRP': 'UniProt/GCRP',
        'UniProt/Swiss-Prot': 'UniProtKB',
        'UniProt/TrEMBL': 'UniProtKB',
        'SGD': 'SGD',
        'WormBase': 'WB',
        'ZFIN': 'ZFIN',
        'Xenbase': 'Xenbase',
        'RGD': 'RGD',
        'MGD': 'MGI',
        'MGI': 'MGI'
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

    # Elaborate on get_general_data() for the GeneHandler.
    def get_general_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_general_data(session)
        self.get_chr_info(session)
        return

    # Elaborate on get_datatype_data() for the GeneHandler.
    def get_panther_info(self):
        """Extract panther information from file."""
        self.log.info('Extract panther information from file.')
        filepath = '/src/input/PTHR18.0_fruit_fly'
        tsv_file = open(filepath, "r")
        tsvin = csv.reader(tsv_file, delimiter='\t')
        FB = 0
        PTHR = 3
        counter = 0
        gene_regex = r'FBgn[0-9]{7}'
        for row in tsvin:
            fields = len(row)
            if fields:  # Ignore blank lines
                if re.search(gene_regex, row[FB]) and re.search(self.regex['panther'], row[PTHR]):
                    self.pthr_dict[re.search(gene_regex, row[FB]).group(0)] = re.search(self.regex['panther'], row[PTHR]).group(0)
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
        self.get_panther_info()
        self.get_annotation_ids(session)
        self.get_chr_featurelocs(session)
        self.get_gene_snapshots(session)
        self.get_gene_types(session)
        return

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_gene_type(self):
        """Synthesize gene type."""
        self.log.info('Synthesize gene type.')
        for gene in self.fb_data_entities.values():
            if len(gene.gene_type_names) == 1:
                prop_value = gene.gene_type_names[0].value
                gene.gene_type_name = prop_value[11:-1]
                gene.gene_type_id = prop_value[1:10].replace('SO', 'SO:')
            elif len(gene.gene_type_names) > 1:
                self.log.warning(f'{gene} has many promoted gene types.')
        return

    def synthesize_info(self):
        """Extend the method for the GeneHandler."""
        super().synthesize_info()
        self.synthesize_gene_type()
        self.synthesize_anno_ids()
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_gene_basic(self):
        """Map basic FlyBase gene data to the Alliance LinkML object."""
        self.log.info('Map basic gene info to Alliance object.')
        for gene in self.fb_data_entities.values():
            agr_gene = datatypes.GeneDTO()
            agr_gene.obsolete = gene.chado_obj.is_obsolete
            agr_gene.curie = f'FB:{gene.uniquename}'
            agr_gene.taxon_curie = gene.ncbi_taxon_id
            gene.linkmldto = agr_gene
        return

    def map_gene_type(self):
        """Map gene type."""
        self.log.info('Map gene type to Alliance object.')
        for gene in self.fb_data_entities.values():
            gene.linkmldto.gene_type_curie = gene.gene_type_id
        return

    def map_gene_snapshot(self):
        """Map gene snapshot."""
        self.log.info('Map gene snapshot to Alliance object.')
        for gene in self.fb_data_entities.values():
            if len(gene.gene_snapshots) == 1:
                note_type_name = 'MOD_provided_gene_description'
                free_text = gene.gene_snapshots[0].value.replace('@', '')
                pub_curies = ['FB:FBrf0232436']
                snapshot_note_dto = datatypes.NoteDTO(note_type_name, free_text, pub_curies).dict_export()
                gene.linkmldto.related_notes.append(snapshot_note_dto)
            elif len(gene.gene_snapshots) > 1:
                self.log.warning(f'{gene} has many gene snapshots.')
        return

    def map_gene_panther_xrefs(self):
        """Add panther xrefs."""
        self.log.info('Map panther xrefs to Alliance object.')
        for gene in self.fb_data_entities.values():
            if gene.uniquename not in self.pthr_dict.keys():
                return
            # Build Alliance xref DTO
            prefix = 'PANTHER'
            page_area = 'FB'
            curie = f'{prefix}:{self.pthr_dict[gene.uniquename]}'
            display_name = curie
            xref_dto = datatypes.CrossReferenceDTO(prefix, curie, page_area, display_name).dict_export()
            gene.linkmldto.cross_reference_dtos.append(xref_dto)
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_gene_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_pubs()
        self.map_xrefs()
        self.map_timestamps()
        self.map_secondary_ids('gene_secondary_id_dtos')
        self.map_gene_snapshot()
        self.map_gene_type()
        self.map_gene_panther_xrefs()
        self.map_anno_ids_to_secondary_ids('gene_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return


class StrainHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters strain data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the StrainHandler object."""
        super().__init__(log, fb_data_type, testing)

    test_set = {
        'FBsn0000001': 'Oregon-R-modENCODE',
        'FBsn0000091': 'DGRP-373',
        'FBsn0000272': 'iso-1',
        'FBsn0001072': 'DSPR-B1-019',
        'FBsn0000283': 'MV2-25',
        'FBsn0000284': 'DGRP_Flyland',
    }
    # Elaborate on export filters for StrainHandler.
    required_fields = {
        'agm_ingest_set': [
            'curie',
            'data_provider_dto',
            'internal',
            'subtype_name',
            'taxon_curie',
        ],
    }
    output_fields = {
        'agm_ingest_set': [
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
        ],
    }

    # Elaborate on get_general_data() for the StrainHandler.
    def get_general_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() for the StrainHandler.
    def get_datatype_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_datatype_data(session)
        return

    # Elaborate on synthesize_info() for the StrainHandler.
    def synthesize_info(self):
        """Extend the method for the StrainHandler."""
        super().synthesize_info()
        return

    # Elaborate on map_fb_data_to_alliance() for the StrainHandler.
    def map_strain_basic(self):
        """Map basic FlyBase strain data to the Alliance object."""
        self.log.info('Map basic strain info.')
        for strain in self.fb_data_entities.values():
            agr_strain = datatypes.AffectedGenomicModelDTO()
            agr_strain.obsolete = strain.chado_obj.is_obsolete
            agr_strain.curie = f'FB:{strain.uniquename}'
            agr_strain.taxon_curie = strain.ncbi_taxon_id
            agr_strain.name = strain.name
            agr_strain.subtype_name = 'strain'
            strain.linkmldto = agr_strain
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the StrainHandler."""
        super().map_fb_data_to_alliance()
        self.map_strain_basic()
        self.map_data_provider_dto()
        self.map_pubs()
        self.map_xrefs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
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
        'construct': ConstructHandler,
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
        export_dict (dict): A dict of LinkML dicts for some "agr_ingest" set.
        log (Logger): The global Logger object in the script calling this function.
        output_filename (str): The global output_filename in the script calling this function.

    """
    log.info('Writing output Alliance LinkML data dict to JSON file..upper()')
    with open(output_filename, 'w') as outfile:
        json.dump(export_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
        outfile.close()
    return
