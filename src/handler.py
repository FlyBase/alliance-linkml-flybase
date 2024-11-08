"""Module:: handler.

Synopsis:
    The core data handler from which datatype-specific handlers are built. These
    handlers run basic processes common to the export of all FlyBase data to the
    Alliance via LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import datetime
import strict_rfc3339
from logging import Logger
from sqlalchemy.orm import aliased
from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cv, Cvterm, CvtermRelationship, Db, Dbxref, Feature, FeatureCvterm,
    FeatureCvtermprop, FeatureRelationship,
    FeatureSynonym, Featureprop, Organism, OrganismDbxref,
    Organismprop, Pub, PubDbxref, Synonym
)


class DataHandler(object):
    """A generic data handler that gets FlyBase data and maps it to LinkML.

    This handler serves as the base for datatype-specific data handlers. In
    general, a given handler will process one datatype, or, two related
    datatypes: e.g., ConstructHandler exports constructs and construct
    associations. Scripts call on these handlers to get and transform FB data.
    In some cases, a script will combine the output of distinct data handlers
    to generate a single integrated export file:
    e.g., AGR_data_retrieval_curation_agm.py uses StrainHandler and
    GenotypeHandler to generate a single "agm_ingest_set".

    """
    def __init__(self, log: Logger, testing: bool):
        """Create the generic DataHandler object.

        Args:
            log (Logger): The global Logger object in the script using the DataHandler.
            testing (bool): Whether handler is being run in testing mode or not.

        """
        self.log = log
        self.datatype = None                   # A single word describing the datatype: e.g., 'gene'. Define for more specific handlers.
        self.fb_export_type = None             # Will be the relevant FBExportEntity object: e.g., FBGene. Define for more specific handlers.
        self.agr_export_type = None            # Will be the LinkML object to export to: e.g., GeneDTO. Define for more specific handlers.
        self.primary_agr_ingest_type = None    # Will be name of the Alliance ingest set: e.g., 'gene_ingest_set'.
        self.testing = testing
        # Datatype bins.
        self.fb_data_entities = {}          # db_primary_id-keyed dict of chado objects to export.
        self.export_data = {}               # agr_ingest_set_name-keyed lists of data objects for export.
        # General data bins.
        self.bibliography = {}              # A pub_id-keyed dict of pub curies (PMID, or, FBrf if no PMID).
        self.cvterm_lookup = {}             # A cvterm_id-keyed dict of dicts with these keys: 'name', 'cv_name', 'db_name', 'curie'.
        self.ncbi_taxon_lookup = {}         # An organism_id-keyed dict of f'NCBITaxon:{Dbxref.accession}' strings.
        self.drosophilid_list = []          # A list of organism_ids for Drosophilid species.
        self.chr_dict = {}                  # Will be a feature_id-keyed dict of chr scaffold uniquenames.
        self.feature_lookup = {}            # feature_id-keyed dicts with these keys: uniquename, is_obsolete, name, symbol, exported, taxon_id, species.
        self.internal_gene_ids = []         # A list of feature_ids for FlyBase genes that should be internal at the Alliance: e.g., 'engineered_fusion_gene'.
        self.featureprop_pub_lookup = {}    # Will be featureprop_id-keyed lists of supporting pub_ids.
        self.feat_rel_pub_lookup = {}       # Will be feature_relationship_id-keyed lists of supporting pub_ids.
        # Trackers.
        self.input_count = 0                # Count of entities found in FlyBase chado database.
        self.export_count = 0               # Count of exported Alliance entities.
        self.internal_count = 0             # Count of exported entities marked as internal.
        self.warnings = []                  # Handler issues of note.
        self.errors = []                    # Handler issues that break things.

    # Sample set for faster testing: use uniquename-keyed names of objects, tailored for each handler.
    test_set = {}

    # Alliance db names should correspond to the contents of this file:
    # https://github.com/alliance-genome/agr_schemas/blob/master/resourceDescriptors.yaml
    fb_agr_db_dict = {
        'EntrezGene': 'NCBI_Gene',
        'RNAcentral': 'RNAcentral',
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
    # Specify page_area for cross-references for specific external databases.
    # For Alliance MODs, the page_area will be the data type: e.g., gene.s
    agr_page_area_dict = {
        'NCBI_Gene': 'default',
        'RNAcentral': 'default',
        'UniProtKB': 'default',
    }
    # Useful regexes.
    regex = {
        'pub': r'^(FBrf[0-9]{7}|unattributed)$',
        'abbal': r'^FB(ab|ba)[0-9]{7}$',
        'aberration': r'^FBab[0-9]{7}$',
        'allele': r'^FBal[0-9]{7}$',
        'balancer': r'^FBba[0-9]{7}$',
        'chemical': r'^FBch[0-9]{7}$',
        'consins': r'^FB(tp|ti)[0-9]{7}$',
        'construct': r'^FBtp[0-9]{7}$',
        'fb_curie': r'^FB:FB[a-z]{2}[0-9]{7,10}$',
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
        handler_description = f'A data handler that exports FB {self.datatype} to Alliance LinkML: {type(self.agr_export_type)}.'
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
        if type(pub_id_list) is not list:
            pub_id_list = [pub_id_list]
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
            cvterm_dict = {
                'name': result.name,
                'cv_name': result.cv.name,
                'db_name': result.dbxref.db.name,
                'curie': f'{result.dbxref.db.name}:{result.dbxref.accession}'
            }
            self.cvterm_lookup[result.cvterm_id] = cvterm_dict
            cvterm_counter += 1
        self.log.info(f'Found {cvterm_counter} current CV terms in chado.')
        return

    def get_child_cvterms(self, session, starting_cvterm_name, starting_cvterm_cv_name):
        """Get all cvterm_ids for some branch of an ontology defined by the parent term.

        Args:
            session (SQLAlchemy session): The object that queries the database.
            starting_cvterm_name (str): The name of the parent CV term.
            starting_cvterm_cv_name (str): The name of the parent CV term's CV.

        Returns:
            List of cvterm.cvterm_ids, including that of the starting parent CV term.

        Raises:
            Raise NoResultFound if starting CV term cannot be found in chado.

        """
        self.log.info(f'Get all child terms of "{starting_cvterm_name}" from the "{starting_cvterm_cv_name}" CV.')
        # 1. Get the parent CV term.
        filters = (
            Cvterm.name == starting_cvterm_name,
            Cv.name == starting_cvterm_cv_name
        )
        try:
            starting_cvterm = session.query(Cvterm).\
                join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
                filter(*filters).\
                one()
            self.log.debug(f'Found chado CV term for "{starting_cvterm_name}" from the "{starting_cvterm_cv_name}" CV.')
        except NoResultFound:
            self.log.error(f'No chado CV term found for "{starting_cvterm_name}" ("{starting_cvterm_cv_name}")')
            raise NoResultFound
        # 2. Define the start of the recursive query for child terms of the starting term.
        # Recursive query built up using the sorta helpful instructions here:
        # https://docs.sqlalchemy.org/en/13/orm/query.html
        # https://sanjayasubedi.com.np/python/sqlalchemy/recursive-query-in-postgresql-with-sqlalchemy/
        cvterm_subject1, cvterm_subject2, cvterm_object = aliased(Cvterm), aliased(Cvterm), aliased(Cvterm)
        cvterm_relationship1, cvterm_relationship2 = aliased(CvtermRelationship), aliased(CvtermRelationship)
        cv1, cv2 = aliased(Cv), aliased(Cv)
        # 2a. Define the starting point of the query (output of which to be used recursively to get more results).
        # Basically, looking for all child terms (subject) of the initial CV term (object) in cvterm_relationship.
        # The "recursive_query_start" below is not actually results, but a 'sqlalchemy.sql.base.ImmutableColumnCollection' class.
        filters = (cvterm_object.cvterm_id == starting_cvterm.cvterm_id,
                   cv1.name == starting_cvterm_cv_name)
        recursive_query_start = session.query(cvterm_subject1).\
            join(cvterm_relationship1, (cvterm_relationship1.subject_id == cvterm_subject1.cvterm_id)).\
            join(cvterm_object, (cvterm_object.cvterm_id == cvterm_relationship1.object_id)).\
            join(cv1, (cv1.cv_id == cvterm_subject1.cv_id)).\
            filter(*filters).\
            cte(recursive=True)    # Important bit here.
        # 3. Define the second part of the recursive query.
        # Essentially, we want child terms of child terms in the starting query, recursively.
        # So, subject_ids in "recursive_query_start" results will be the objects in our query of cvterm_relationship.
        # Importantly, the columns of the initial "recursive_query_start" query are accessed by the ".c" attribute.
        filters2 = (cv2.name == starting_cvterm_cv_name,)
        recursive_query_repeat = session.query(cvterm_subject2).\
            join(cvterm_relationship2, (cvterm_relationship2.subject_id == cvterm_subject2.cvterm_id)).\
            join(recursive_query_start, (recursive_query_start.c.cvterm_id == cvterm_relationship2.object_id)).\
            join(cv2, (cv2.cv_id == cvterm_subject2.cv_id)).\
            filter(*filters2)
        # 4. Define a query that takes the union of the starting and recursive queries.
        recursive_query_total = recursive_query_start.union(recursive_query_repeat)
        # 5. And finally, get the result for the start and recursive query parts). Piece of cake ;)
        recursive_query_total_results = session.query(recursive_query_total)
        # Build the list from the results.
        cvterm_id_list = [i[0] for i in recursive_query_total_results]
        self.log.info(f'Found {len(cvterm_id_list)} terms under "{starting_cvterm_name}" from the "{starting_cvterm_cv_name}" CV.')
        cvterm_id_list.append(starting_cvterm.cvterm_id)
        return cvterm_id_list

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

    def get_drosophilid_organisms(self, session):
        """Find Drosophilid species."""
        self.log.info('Find Drosophilid species.')
        filters = (
            Cvterm.name == 'taxgroup',
            Organismprop.value == 'drosophilid'
        )
        drosophilid_results = session.query(Organism).\
            join(Organismprop, (Organismprop.organism_id == Organism.organism_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Organismprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for drosophilid in drosophilid_results:
            self.drosophilid_list.append(drosophilid.organism_id)
            counter += 1
        self.log.info(f'Found {counter} Drosophilid organisms in chado.')
        return

    def build_public_feature_lookup(self, session, **kwargs):
        """Build a simple feature lookup for FlyBase features having an FB ID.

        Args:
            session (SQLAlchemy session): The object that queries the database.

        Keyword Args:
            feature_types (str|list): A list of public feature types to add to the lookup.

        """
        # Note - depends on prior construction of self.ncbi_taxon_lookup and self.cvterm_lookup.
        feat_type_export = {
            'aberration': True,
            'allele': True,
            'balancer': True,
            'chemical': False,
            'construct': True,
            'gene': True,
            'insertion': True,
            'seqfeat': False,
            'tool': False,
        }
        if 'feature_types' in kwargs.keys():
            if type(kwargs['feature_types']) is str:
                kwargs['feature_types'] = [kwargs['feature_types']]
            elif type(kwargs['feature_types']) is not list or type(kwargs['feature_types'][0]) is not str:
                error_msg = 'The build_public_feature_lookup() method accepts only a string or list of strings,'
                error_msg += f'but the "feature_types" argument given was this: {kwargs['feature_types']}.'
                self.log.error(error_msg)
                raise
        else:
            kwargs['feature_types'] = list(feat_type_export.keys())
        self.log.info(f'Build a simple feature lookup for these public feature types: {kwargs["feature_types"]}.')

        # First get features.
        for feat_type in kwargs['feature_types']:
            if feat_type not in feat_type_export.keys():
                self.log.error(f'The feature type given, "{feat_type}" is not in the acceptable list: {feat_type_export.keys()}')
                raise
            self.log.info(f'Looking up {feat_type} features.')
            feat_filters = (
                Feature.uniquename.op('~')(self.regex[feat_type]),
            )
            feat_results = session.query(Feature.feature_id, Feature.uniquename, Feature.is_obsolete,
                                         Feature.type_id, Organism.organism_id, Organism.genus,
                                         Organism.species, Feature.name).\
                select_from(Feature).\
                join(Organism, (Organism.organism_id == Feature.organism_id)).\
                filter(*feat_filters).\
                distinct()
            FEATURE_ID = 0
            UNIQUENAME = 1
            OBSOLETE = 2
            TYPE_ID = 3
            ORG_ID = 4
            GENUS = 5
            SPECIES = 6
            NAME = 7
            feat_counter = 0
            for result in feat_results:
                feat_dict = {
                    'uniquename': result[UNIQUENAME],
                    'is_obsolete': result[OBSOLETE],
                    'type': self.cvterm_lookup[result[TYPE_ID]]['name'],
                    'species': f'{result[GENUS]} {result[SPECIES]}',
                    'name': result[NAME],
                    'symbol': result[NAME],
                    'exported': feat_type_export[feat_type],
                }
                try:
                    feat_dict['taxon_id'] = self.ncbi_taxon_lookup[result[ORG_ID]]
                except KeyError:
                    feat_dict['taxon_id'] = 'NCBITaxon:32644'    # Unspecified taxon.
                self.feature_lookup[result[FEATURE_ID]] = feat_dict
                feat_counter += 1
            self.log.info(f'Added {feat_counter} {feat_type} features to the feature_lookup.')
            # Second, get current symbol synonym for each feature, if it exists.
            syno_filters = (
                Feature.uniquename.op('~')(self.regex[feat_type]),
                FeatureSynonym.is_current.is_(True),
                Cvterm.name == 'symbol',
            )
            syno_results = session.query(Feature.feature_id, Synonym.name).\
                select_from(Feature).\
                join(FeatureSynonym, (FeatureSynonym.feature_id == Feature.feature_id)).\
                join(Synonym, (Synonym.synonym_id == FeatureSynonym.synonym_id)).\
                join(Cvterm, (Cvterm.cvterm_id == Synonym.type_id)).\
                filter(*syno_filters).\
                distinct()
            FEATURE_ID = 0
            SYMBOL = 1
            syno_counter = 0
            for result in syno_results:
                try:
                    self.feature_lookup[result[FEATURE_ID]]['symbol'] = sub_sup_sgml_to_html(result[SYMBOL])
                    syno_counter += 1
                except KeyError:
                    pass
            self.log.info(f'Found {syno_counter} symbol synonyms for {feat_counter} {feat_type} features.')
        return

    def get_internal_genes(self, session):
        """Find FlyBase genes that should be internal at the Alliance due to their type."""
        self.log.info('Find FlyBase genes that should be internal at the Alliance due to their type.')
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
        filters = (
            Feature.uniquename.op('~')(self.regex['gene']),
            Cvterm.name == 'promoted_gene_type'
        )
        gene_type_results = session.query(Featureprop).\
            select_from(Feature).\
            join(Featureprop, (Feature.feature_id == Featureprop.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in gene_type_results:
            gene_type_name = result.value[11:-1]
            if gene_type_name in internal_gene_types:
                counter += 1
                self.internal_gene_ids.append(result.feature_id)
        self.log.info(f'Found {counter} internal type genes.')
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

    # DETRITUS
    # def build_feature_relationship_evidence_lookup(self, session):
    #     """Build evidence lookup for feature_relationships."""
    #     self.log.info('Build evidence lookup for feature_relationships.')
    #     results = session.query(FeatureRelationshipPub.feature_relationship_id,
    #                             FeatureRelationshipPub.pub_id).distinct()
    #     FEAT_REL_ID = 0
    #     PUB_ID = 1
    #     counter = 0
    #     fr_counter = 0
    #     for result in results:
    #         try:
    #             self.feat_rel_pub_lookup[result[FEAT_REL_ID]].append(result[PUB_ID])
    #             counter += 1
    #         except KeyError:
    #             self.feat_rel_pub_lookup[result[FEAT_REL_ID]] = [result[PUB_ID]]
    #             fr_counter += 1
    #             counter += 1
    #     self.log.info(f'Found {counter} pubs supporting {fr_counter} feature_relationships.')
    #     return

    # def lookup_feat_rel_pub_ids(self, feature_relationship_id):
    #     """Return a list of pub_ids supporting a given feature_relationship."""
    #     try:
    #         pub_ids = self.feat_rel_pub_lookup[feature_relationship_id]
    #     except KeyError:
    #         pub_ids = []
    #     return pub_ids

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

    # The get_general_data() wrapper; sub-methods are defined and called in more specific DataHandler types.
    def get_general_data(self, session):
        """Get general FlyBase chado data."""
        self.log.info('GET GENERAL FLYBASE DATA FROM CHADO.')
        return

    # The get_datatype_data() wrapper; sub-methods are defined and called in more specific DataHandler types.
    def get_datatype_data(self, session):
        """Get datatype-specific FlyBase data from chado."""
        self.log.info(f'GET FLYBASE {self.datatype.upper()} DATA FROM CHADO.')
        return

    # The synthesize_info() wrapper; sub-methods are defined and called in more specific DataHandler types.
    def synthesize_info(self):
        """Synthesize FB info for each data object."""
        self.log.info(f'SYNTHESIZE FLYBASE {self.datatype.upper()} DATA FROM CHADO.')
        return

    # Sub-methods for the map_fb_data_to_alliance() wrapper.
    def map_timestamps(self):
        """Map timestamps to Alliance object."""
        self.log.info('Map timestamps to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.timestamps and fb_data_entity.linkmldto is not None:
                fb_data_entity.linkmldto.date_created = strict_rfc3339.\
                    timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(fb_data_entity.timestamps)))
                fb_data_entity.linkmldto.date_updated = strict_rfc3339.\
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

    # The map_fb_data_to_alliance() wrapper; sub-methods are called (and usually defined) in more specific DataHandler types.
    def map_fb_data_to_alliance(self):
        """Map FB data to the Alliance LinkML object."""
        self.log.info(f'MAP FLYBASE "{self.datatype}" DATA FROM {self.fb_export_type} TO {self.agr_export_type}.')
        return

    def flag_unexportable_entities(self, input_list: list, output_set_name: str):
        """Flag entities lacking information for a required field.

        Args:
            input_list (list): A list of objects having a linkmldto attribute that is to be exported as a dict.
            output_set_name (str): The export_set_label for the list of exported dicts: e.g., 'agm_ingest_set'

        """
        self.log.info(f'Flag FlyBase data lacking information for a required field in the {output_set_name}.'.upper())
        input_counter = 0
        no_linkml_mapping_counter = 0
        missing_required_info_counter = 0
        for i in input_list:
            input_counter += 1
            # for attr in self.required_fields[output_set_name]:
            if i.linkmldto is None:
                i.for_export = False
                i.export_warnings.append('Not mappable to LinkML at all.')
                no_linkml_mapping_counter += 1
                continue
            for attr in i.linkmldto.required_fields:
                if attr not in i.linkmldto.__dict__.keys():
                    i.for_export = False
                    i.export_warnings.append(f'Missing "{attr}" attribute.')
                elif getattr(i.linkmldto, attr) is None:
                    i.for_export = False
                    i.export_warnings.append(f'Missing value "{attr}" attribute.')
            if i.for_export is False:
                self.log.debug(f'DO NOT EXPORT {i}: {i.export_warnings}')
                missing_required_info_counter += 1
        self.log.info(f'Assessed {input_counter} objects for exportability.')
        self.log.info(f'Found {no_linkml_mapping_counter} objects with no LinkML mapping at all.')
        self.log.info(f'Found {missing_required_info_counter} objects missing required LinkML info.')
        return

    def generate_export_dict(self, input_list: list, output_set_name: str):
        """Generate LinkML export dict from FB data.

        Args:
            input_list (list): A list of objects with LinkMLDTO objects under the linkmldto attribute.
            output_set_name (str): The "agr_ingest_set" label for the list of exported dicts.

        """
        self.log.info(f'Generate LinkML export dict from FB data for {output_set_name}.'.upper())
        # Create the export_data_list, keyed by the agr_ingest_set name.
        self.export_data[output_set_name] = []
        self.input_count = 0
        self.export_count = 0
        self.internal_count = 0
        for i in input_list:
            self.input_count += 1
            if i.for_export is False:
                continue
            self.export_count += 1
            if i.linkmldto.internal is True:
                self.internal_count += 1
                self.log.debug(f'Export {i} but keep INTERNAL at the Alliance: {i.internal_reasons}')
            export_agr_dict = {}
            for attr in i.linkmldto.__dict__.keys():
                # self.log.debug(f'Assess this attr: {attr}')
                if attr in i.linkmldto.internal_fields:
                    # self.log.debug(f'Skip this field: {attr}')
                    continue
                elif attr in i.linkmldto.required_fields:
                    # self.log.debug(f'Export required field: {attr}')
                    export_agr_dict[attr] = getattr(i.linkmldto, attr)
                elif getattr(i.linkmldto, attr) is not None and getattr(i.linkmldto, attr) != []:
                    # self.log.debug(f'Export optional non-empty field: {attr}')
                    export_agr_dict[attr] = getattr(i.linkmldto, attr)
                else:
                    # self.log.debug(f'Empty value for this attr: {attr}')
                    pass
            self.export_data[output_set_name].append(export_agr_dict)
        public_count = self.export_count - self.internal_count
        self.log.info(f'SUMMARY FOR EXPORT OF {output_set_name}'.upper())
        self.log.info(f'Exported {self.export_count} of {self.input_count} entities.')
        self.log.info(f'{public_count} of {self.export_count} exported entities are PUBLIC.')
        self.log.info(f'{self.internal_count} of {self.export_count} exported entities are INTERNAL.')
        return

    # The query_chado_and_export() wrapper that runs sub-methods - same order of steps for every DataHandler type.
    def query_chado_and_export(self, session):
        """Wrapper that runs all methods within an SQLAlchemy session."""
        self.log.info('Run main query_chado_and_export() handler method.'.upper())
        self.get_general_data(session)
        self.get_datatype_data(session)
        self.synthesize_info()
        self.map_fb_data_to_alliance()
        self.flag_unexportable_entities(self.fb_data_entities.values(), self.primary_export_set)
        self.generate_export_dict(self.fb_data_entities.values(), self.primary_export_set)
        return
