# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase genes for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_gene_curation.py [-h] [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_gene_curation.py -v -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script makes a JSON file conforming to LinkML specs for the curation
    (i.e., "persistent") database; distinct from BGI gene file specs for the
    original Neo4j drop-and-reload database. As the "audit_chado" table may be
    required for proper determination of "date_updated" values, use a reporting
    database with full audit_chado table present.

"""

import argparse
import csv
import datetime
import json
import re
import strict_rfc3339
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Feature, FeatureDbxref, FeatureSynonym,
    Featureloc, Featureprop, OrganismDbxref, Pub, PubDbxref, Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Generic setup.
report_label = 'gene_curation'
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
reference_assembly = set_up_dict['assembly']
input_dir = set_up_dict['input_dir']
output_filename = set_up_dict['output_filename'].replace('tsv', 'json')
log = set_up_dict['log']

# Process additional input parameters not handled by the set_up_db_reading() function above.
parser = argparse.ArgumentParser(description='inputs')
parser.add_argument('-l', '--linkml_release', help='The "agr_curation_schema" LinkML release number.', required=True)
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))
linkml_release = args.linkml_release

# Create SQL Alchemy engines from environmental variables.
engine_var_rep = 'postgresql://' + username + ":" + password + '@' + server + '/' + database
engine = create_engine(engine_var_rep)
insp = inspect(engine)


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase Genes."""
    log.info('Running main() for script "{}"'.format(__file__))
    log.info('Output corresponds to "agr_curation_schema" release: {}'.format(linkml_release))
    # Instantiate the object, get the data, synthesize it, export it.
    gene_handler = GeneHandler()
    db_query_transaction(gene_handler)
    gene_handler.synthesize_info()
    gene_handler.generate_export_file()
    log.info('Ended main function.\n')


class AllianceGene(object):
    """A gene with it's associated FlyBase data and Alliance LinkML properties."""
    def __init__(self, feature):
        """Create a base AllianceGene object.

        Args:
            arg1 (feature): (Feature) The Feature object corresponding to the gene.

        Returns:
            An object of the AllianceGene class.

        """
        # Attributes representing unprocessed FlyBase data.
        # Note: use attribute names that do not match an Alliance LinkML slot name.
        # For initial load, the Alliance A-Team just needs minimum info:
        # GENE: curie, taxon, symbol, name.
        # ALLELE: curie, taxon, symbol, description.
        # Problems with Gene LinkML:
        # 1. Gene.name is requested (not required), but not all genes have a fullname.
        # 2. Gene.taxon_curie is required, but even after updating NCBITaxon info at FlyBase, not all genes will have NCBI taxon ID.
        # 3. GenomicLocation lacks strand info.
        self.feature = feature                                # The Feature object corresponding to the FlyBase gene.
        self.organism_abbr = None                             # Will be the organism.abbreviation for the gene's species of origin.
        self.featureloc = None                                # Will be Featureloc object for the gene.
        self.gene_type_name = None                            # Will be the cvterm.name for "promoted_gene_type" featureprop.
        self.gene_snapshot = None                             # Will be the "gene_summary_text" Featureprop object.
        self.curr_symbol_name = None                          # Will be the current symbol synonym.synonym_sgml, processed by sub_sup_sgml_to_html().
        self.curr_fullname = None                             # Will be the current fullname synonym.synonym_sgml, processed by sub_sup_sgml_to_html().
        self.curr_anno_id = None                              # Will be current annotation ID for the gene (str).
        self.feature_synonyms = []                            # Will be list of all FeatureSynonym objects.
        self.dbxrefs = []                                     # Will be list of dbxrefs as sql result groupings: Db, Dbxref, FeatureDbxref.
        self.alt_fb_ids = []                                  # Will be list of Dbxrefs for 2o FlyBase IDs.
        self.annotation_ids = []                              # Will be list of Dbxrefs for annotation IDs.
        self.timestamps = []                                  # Add all timestamps here.
        # Attributes for the Alliance AuditedObjectDTO.
        self.obsolete = feature.is_obsolete                   # Will be the FlyBase value here.
        self.internal = False                                 # Change to true if not public on FlyBase.
        self.created_by_curie = 'FB:FB_curator'               # Use placeholder value since no Person object at FlyBase.
        self.updated_by_curie = 'FB:FB_curator'               # Use placeholder value since no Person object at FlyBase.
        self.date_created = None                              # Not straightforward as half of relevant annotations are derived in the reporting build.
        self.date_updated = None                              # Not straightforward as half of relevant annotations are derived in the reporting build.
        # Attributes for the Alliance BiologicalEntityDTO. BiologicalEntityDTO is_a AuditedObjectDTO.
        self.curie = 'FB:{}'.format(feature.uniquename)
        self.taxon_curie = None                               # A string representing the NCBI taxon ID. We have no NCBI taxonID for 561 genes (72 species).
        # Attributes for the Alliance GenomicEntityDTO. GenomicEntityDTO is_a BiologicalEntityDTO.
        self.cross_reference_dtos = []                        # Report only select dbs, using AGR-accepted db_prefix.
        self.secondary_identifiers = []                       # Annotation IDs and 2o FlyBase IDs.
        self.genomic_location_dtos = []                       # Will need to be list of GenomicLocation objects.
        # Attributes for the Alliance GeneDTO. GeneDTO is_a GenomicEntityDTO.
        self.gene_symbol_dto = None                           # Will be a single SymbolSlotAnnotationDTO.
        self.gene_full_name_dto = None                        # Will be a single FullNameSlotAnnotation.
        self.gene_systematic_name_dto = None                  # Will be a single SystematicNameSlotAnnotation.
        self.gene_synonym_dtos = []                           # Will be list of NameSlotAnnotationDTO objects.
        self.gene_type_curie = None                           # Will be the SO term ID corresponding to the gene's promoted_gene_type.
        # Notes associated with the object.
        self.for_alliance_export = True                       # Change to False if object should be excluded from export.
        self.internal_reasons = []                            # Reasons for marking an object as internal in the export file.
        self.export_warnings = []                             # Reasons for suppressing an object from the export file.

    def __str__(self):
        """Succinct text string describing the AllianceGene object."""
        desc = '{} ({})'.format(self.feature.name, self.feature.uniquename)
        return desc


class GeneHandler(object):
    """This object gets genes and related info, synthesizes/filters the data, then exports it as LinkML JSON."""
    def __init__(self):
        """Create the GeneHandler object."""
        self.gene_dict = {}           # An FBgnID-keyed dict of AllianceGene objects.
        self.all_pubs_dict = {}       # A pub_id-keyed dict of pub curies (PMID or FBrf).
        self.all_synonyms_dict = {}   # A synonym_id-keyed dict of Synonym objects.
        self.pthr_dict = {}           # Will be an 1:1 FBgnID-PTHR xref dict.
        self.chr_dict = {}            # Will be a feature_id-keyed dict of chr scaffold uniquenames.
        self.total_feat_cnt = 0       # Count of all genes found in starting query.
        self.export_feat_cnt = 0      # Count of all genes exported to file.
        self.internal_feat_cnt = 0    # Count of all genes marked as internal=True in export file.

    # Regexes.
    gene_regex = r'^FBgn[0-9]{7}$'
    pthr_regex = r'PTHR[0-9]{5}'
    pub_regex = r'^(FBrf[0-9]{7}|unattributed)$'
    systematic_name_regex = r'^(D[a-z]{3}\\|)(CG|CR|G[A-Z])[0-9]{4,5}$'
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
    fb_agr_db_dict = {
        'EntrezGene': 'NCBI_Gene',
        'FlyBase': 'FB',
        'FlyBase Annotation IDs': 'FB',
        'RNAcentral': 'RNAcentral',
        # 'UniProt/GCRP': 'UniProt/GCRP',
        'UniProt/Swiss-Prot': 'UniProtKB',
        'UniProt/TrEMBL': 'UniProtKB'
    }
    # Sample set.
    test_genes = ['wg', 'mt:ori', 'lncRNA:roX1', 'CG12656']
    # Export fields.
    required_fields = [
        'curie',
        'gene_symbol_dto',
        'internal',
        'taxon_curie',
    ]
    output_fields = [
        'created_by_curie',
        'cross_reference_dtos',
        'curie',
        'date_created',
        'date_updated',
        'gene_full_name_dto',
        'gene_symbol_dto',
        'gene_synonym_dtos',
        'gene_systematic_name_dto',
        'gene_type_curie',
        'internal',
        'obsolete',
        'secondary_identifiers',
        'taxon_curie',
        'updated_by_curie',
    ]

    def open_panther_file(self):
        """Extract panther information from file."""
        log.info('Open PANTHER file.')
        if input_dir == '/src/input/':
            filepath = f'{input_dir}PTHR17.0_fruit_fly'
        else:
            filepath = '/data/ortholog/panther/PTHR17.0_fruit_fly'
        tsv_file = open(filepath, "r")
        tsvin = csv.reader(tsv_file, delimiter='\t')
        FB = 0
        PTHR = 3
        counter = 0
        gene_regex = r'FBgn[0-9]{7}'    # Since the FBgn ID does not represent the entire column entry, do not use self.gene_regex here.
        for row in tsvin:
            fields = len(row)
            if fields:  # Ignore blank lines
                if re.search(gene_regex, row[FB]) and re.search(self.pthr_regex, row[PTHR]):
                    self.pthr_dict[re.search(gene_regex, row[FB]).group(0)] = re.search(self.pthr_regex, row[PTHR]).group(0)
                    counter += 1
        log.info(f'Processed {counter} lines from the panther orthology file.')
        return

    def get_all_references(self, session):
        """Get all references."""
        log.info('Get all references.')
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
        log.info(f'Found {pmid_counter} PMID IDs for {pub_counter} current FB publications.')
        return

    def get_genes(self, session):
        """Get all genes."""
        log.info('Querying chado for genes.')
        # First get all gene features from chado.
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'gene'
        )
        gene_results = session.query(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        self.total_feat_cnt = 0
        for result in gene_results:
            self.total_feat_cnt += 1
            self.gene_dict[result.uniquename] = AllianceGene(result)
            self.gene_dict[result.uniquename].organism_abbr = result.organism.abbreviation
        log.info('Found {} genes.'.format(self.total_feat_cnt))

    def get_gene_taxons(self, session):
        """Get taxon IDs for genes."""
        log.info('Getting gene taxon IDs.')
        filters = (
            OrganismDbxref.is_current.is_(True),
            Db.name == 'NCBITaxon'
        )
        organism_dbxref_results = session.query(OrganismDbxref, Dbxref).\
            join(Dbxref, (Dbxref.dbxref_id == OrganismDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        organism_taxon_dict = {}
        org_counter = 0
        gene_counter = 0
        for result in organism_dbxref_results:
            organism_taxon_dict[result.OrganismDbxref.organism_id] = result.Dbxref.accession
            org_counter += 1
        for gene in self.gene_dict.values():
            try:
                gene.taxon_curie = 'NCBITaxon:{}'.format(organism_taxon_dict[gene.feature.organism_id])
                gene_counter += 1
            except KeyError:
                log.debug('No NCBI taxon ID available for: {}'.format(gene))
        log.info(f'Found {org_counter} distinct NCBITaxon IDs for {gene_counter} genes.')
        return

    def get_synonyms(self, session):
        """Get current and non-current symbols and full names for genes."""
        log.info('Get current and non-current symbols and full names for genes.')
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'gene'
        )
        results = session.query(Feature, FeatureSynonym, Synonym).\
            join(FeatureSynonym, (FeatureSynonym.synonym_id == Synonym.synonym_id)).\
            join(Feature, (Feature.feature_id == FeatureSynonym.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            # Skip any references to non-current pubs.
            if result.FeatureSynonym.pub_id not in self.all_pubs_dict.keys():
                continue
            # First, build the all_synonyms_dict.
            self.all_synonyms_dict[result.Synonym.synonym_id] = result.Synonym
            # Second, collect FeatureSynonyms for each gene.
            self.gene_dict[result.Feature.uniquename].feature_synonyms.append(result.FeatureSynonym)
            # Catch current symbol and fullname strings.
            if result.FeatureSynonym.is_current is True and result.Synonym.type.name == 'symbol':
                self.gene_dict[result.Feature.uniquename].curr_symbol_name = sub_sup_sgml_to_html(result.Synonym.synonym_sgml)
            elif result.FeatureSynonym.is_current is True and result.Synonym.type.name == 'fullname':
                self.gene_dict[result.Feature.uniquename].curr_fullname = sub_sup_sgml_to_html(result.Synonym.synonym_sgml)
            counter += 1
        log.info(f'Found {counter} feature_synonyms (current pubs) for genes.')
        return

    def get_annotation_ids(self, session):
        """Get current annotation IDs."""
        log.info('Get current annotation IDs.')
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            Feature.is_obsolete.is_(False),
            Cvterm.name == 'gene',
            FeatureDbxref.is_current.is_(True),
            Db.name == 'FlyBase Annotation IDs'
        )
        results = session.query(Feature, Dbxref).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.gene_dict[result.Feature.uniquename].curr_anno_id = result.Dbxref.accession
            counter += 1
        log.info(f'Found {counter} current annotation IDs for FlyBase genes.')
        return

    def get_gene_dbxrefs(self, session):
        """Get all dbxrefs for genes. This will take 10-15 minutes."""
        log.info('Getting gene dbxrefs.')
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'gene',
            Db.name.in_((self.fb_agr_db_dict.keys())),
        )
        results = session.query(Feature, FeatureDbxref, Dbxref, Db).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            if counter % 100000 == 0:
                log.debug('Processing xref #{}'.format(counter))
            # Skip current FlyBase accessions because these are not comprehensive.
            # When they exist, they're always equal to the feature.uniquename.
            # But they're not always present, so these dbxrefs can't be relied upon (e.g., FBgn0085177)
            if result.FeatureDbxref.is_current is True and result.Db.name == 'FlyBase':
                pass
            elif result.FeatureDbxref.is_current is False and result.Db.name == 'FlyBase':
                self.gene_dict[result.Feature.uniquename].alt_fb_ids.append(result.Dbxref)
            elif result.Db.name == 'FlyBase Annotation IDs':
                self.gene_dict[result.Feature.uniquename].annotation_ids.append(result.Dbxref)
            else:
                self.gene_dict[result.Feature.uniquename].dbxrefs.append(result)
        log.info(f'Found {counter} gene dbxrefs.')
        return

    def get_gene_snapshots(self, session):
        """Get human-written gene summaries."""
        log.info('Getting gene snapshots.')
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='gene_summary_text')
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            feature_type.name == 'gene',
            prop_type.name == 'gene_summary_text'
        )
        gene_snapshot_results = session.query(Featureprop).\
            join(Feature, (Feature.feature_id == Featureprop.feature_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in gene_snapshot_results:
            self.gene_dict[result.feature.uniquename].gene_snapshot = result
            counter += 1
        log.info(f'Found {counter} gene snapshots.')
        return

    def get_gene_types(self, session):
        """Get and parse "promoted_gene_type" featureprop for genes."""
        log.info('Getting gene types.')
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='promoted_gene_type')
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            prop_type.name == 'promoted_gene_type'
        )
        gene_type_results = session.query(Featureprop).\
            join(Feature, (Feature.feature_id == Featureprop.feature_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in gene_type_results:
            self.gene_dict[result.feature.uniquename].gene_type_curie = result.value[1:10].replace('SO', 'SO:')
            self.gene_dict[result.feature.uniquename].gene_type_name = result.value[11:-1]
            counter += 1
        log.info(f'Found {counter} gene types for genes.')
        return

    def get_gene_timestamps(self, session):
        """Get timestamps for genes."""
        log.info('Getting gene timestamps.')
        # Start with timestamps in feature itself.
        for gene in self.gene_dict.values():
            gene.timestamps.append(gene.feature.timeaccessioned)
            gene.timestamps.append(gene.feature.timelastmodified)
        # We may want to get timestamps for other gene attibutes as well (for "date_updated").
        # Placeholder: current query of "feature" table in audit_chado.
        # Querying audit_chado with SQLAlchemy has not been user-friendly, so go with standard SQL.
        ########################################################################
        # audit_chado_query = """
        #     SELECT DISTINCT f.uniquename,
        #                     ac.transaction_timestamp
        #     FROM feature
        #     JOIN audit_chado ac ON (ac.record_pkey = f.feature_id AND ac.audited_table = 'feature');
        # """
        # audit_results = session.execute(audit_chado_query).fetchall()
        # log.info('Got {} audit_chado results. Will parse them out now.'.format(len(audit_results)))
        # UNIQUE_KEY = 0
        # TIMESTAMP = 1
        # for row in audit_results:
        #     try:
        #         # log.debug('For unique_key={}, have timestamp={}'.format(row[UNIQUE_KEY], row[TIMESTAMP]))
        #         self.gene_dict[row[UNIQUE_KEY]].timestamps.append(row[TIMESTAMP])
        #     except KeyError:
        #         # log.debug('Could not put this in anno dict: {}'.format(row))
        #         pass
        ########################################################################
        return

    def get_gene_featureloc(self, session):
        """Getting gene featureloc."""
        log.info('Getting gene genomic locations.')
        # First get a simple chr scaffold dict. We'll need this later.
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'golden_path'
        )
        chr_results = session.query(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        self.chr_dict = {}
        chr_counter = 0
        for result in chr_results:
            if result.organism.abbreviation != 'Dmel':
                continue
            self.chr_dict[result.feature_id] = result.uniquename
            chr_counter += 1
        log.info(f'Got basic info for {chr_counter} current Dmel chr scaffolds.')
        # Now get gene featureloc.
        filters = (
            Feature.uniquename.op('~')(self.gene_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'gene'
        )
        gene_featureloc_results = session.query(Featureloc).\
            join(Feature, (Feature.feature_id == Featureloc.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        gene_counter = 0
        for result in gene_featureloc_results:
            self.gene_dict[result.feature.uniquename].featureloc = result
            gene_counter += 1
        log.info(f'Found {gene_counter} genomic locations for genes.')
        return

    def query_chado(self, session):
        """A wrapper method that runs initial db queries."""
        self.open_panther_file()
        self.get_all_references(session)
        self.get_genes(session)
        self.get_gene_taxons(session)
        # self.get_gene_dbxrefs(session)    # BOB - suppress for faster test iterations.
        self.get_synonyms(session)
        self.get_annotation_ids(session)
        self.get_gene_snapshots(session)
        self.get_gene_types(session)
        self.get_gene_timestamps(session)
        self.get_gene_featureloc(session)
        return

    def process_feature_synonyms(self, feature):
        """Generate name/synonym DTOs for a feature that has a list of FeatureSynonym objects."""
        # Dict for converting FB to AGR synonym types.
        synonym_type_conversion = {
            'symbol': 'nomenclature_symbol',
            'fullname': 'full_name',
            'nickname': 'nomenclature_symbol',
            'synonym': 'nomenclature_symbol'
        }
        default_name_dto = {
            'name_type_name': 'unspecified',
            'format_text': 'unspecified',
            'display_text': 'unspecified',
            'synonym_scope_name': 'exact',
            'evidence_curies': [],
            'internal': False,
            'obsolete': False
        }
        # Create a dict of all distinct name/synonym_sgml combinations: for each, capture synonym type(s) an pub_ids.
        # Keys are (synonym.name, synonym.synonym_sgml) tuples.
        # Values are dicts too where keys are chado synonym types and values are lists of pub_ids.
        # Value dict also has an "internal" key that stores list of FeatureSynonym.is_internal values.
        feature_synonym_dict = {}
        for f_s in feature.feature_synonyms:
            synonym = self.all_synonyms_dict[f_s.synonym_id]
            distinct_synonym_name = (synonym.name, synonym.synonym_sgml)
            if distinct_synonym_name in feature_synonym_dict.keys():
                feature_synonym_dict[distinct_synonym_name]['internal'].append(f_s.is_internal)
                if synonym.type.name in feature_synonym_dict[distinct_synonym_name].keys():
                    feature_synonym_dict[distinct_synonym_name][synonym.type.name].append(f_s.pub_id)
                else:
                    feature_synonym_dict[distinct_synonym_name][synonym.type.name] = [f_s.pub_id]
            else:
                feature_synonym_dict[distinct_synonym_name] = {synonym.type.name: [f_s.pub_id], 'internal': [f_s.is_internal]}
        # Convert to AGR name DTO objects.
        name_dto_list = []
        FORMAT_TEXT = 0
        DISPLAY_TEXT = 1
        for syno_name, syno_attributes in feature_synonym_dict.items():
            # Determine internal status. False trumps True.
            if False in set(syno_attributes['internal']):
                syno_internal = False
            else:
                syno_internal = True
            # Collect all pubs.
            pub_id_list = []
            for syno_type, syno_type_pub_list in syno_attributes.items():
                if syno_type == 'internal':
                    continue
                pub_id_list.extend(syno_type_pub_list)
            pub_id_list = list(set(pub_id_list))
            # Out of curiosity, report cases where same synonym used as both symbol and fullname.
            if 'symbol' in syno_attributes.keys() and 'fullname' in syno_attributes.keys():
                n_symb = len(syno_attributes['symbol'])
                n_full = len(syno_attributes['fullname'])
                log.warning(f"RATIO = {round(n_symb/n_full)}, SYMBOL_USAGE: n={n_symb}, FULLNAME_USAGE: n={n_full}, GENE={feature}, SYNONYM={syno_name}.")
            # Pick correct name type to apply.
            if re.match(self.systematic_name_regex, syno_name[DISPLAY_TEXT]):
                name_type_to_use = 'systematic_name'
            else:
                type_tally = {}
                for syno_type, syno_type_pub_list in syno_attributes.items():
                    if syno_type == 'internal':
                        continue
                    type_tally[len(set(syno_type_pub_list))] = syno_type
                name_type_to_use = synonym_type_conversion[type_tally[max(type_tally.keys())]]
            output_synonym_dto = {
                'name_type_name': name_type_to_use,
                'format_text': sub_sup_sgml_to_html(syno_name[FORMAT_TEXT]),
                'display_text': sub_sup_sgml_to_html(syno_name[DISPLAY_TEXT]),
                'synonym_scope_name': 'exact',
                'evidence_curies': [self.all_pubs_dict[i] for i in pub_id_list if self.all_pubs_dict[i] != 'FB:unattributed'],
                'internal': syno_internal,
                'obsolete': False
            }
            name_dto_list.append(output_synonym_dto)
        # Sift through name DTOs for symbol, fullname, systematic_name, etc.
        for name_dto in name_dto_list:
            if name_dto['display_text'] == feature.curr_anno_id:
                if name_dto['name_type_name'] != 'systematic_name':
                    log.warning(f"{feature}: Found mistyped curr anno ID: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
                    name_dto['name_type_name'] = 'systematic_name'
                feature.gene_systematic_name_dto = name_dto
            if name_dto['display_text'] == feature.curr_symbol_name:
                if name_dto['name_type_name'] not in ['systematic_name', 'nomenclature_symbol']:
                    log.warning(f"{feature}: Found mistyped curr symbol: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
                    name_dto['name_type_name'] = 'nomenclature_symbol'
                feature.gene_symbol_dto = name_dto
            elif name_dto['display_text'] == feature.curr_fullname:
                if name_dto['name_type_name'] != 'full_name':
                    log.warning(f"{feature}: Found mistyped curr full_name: type={name_dto['name_type_name']}, name={name_dto['display_text']}")
                    name_dto['name_type_name'] = 'full_name'
                feature.gene_full_name_dto = name_dto
            else:
                feature.gene_synonym_dtos.append(name_dto)
        # LinkML change required: make gene_full_name_dto and gene_systematic_name_dto OPTIONAL.
        # Symbol is required. If none, fill it in.
        if feature.gene_symbol_dto is None:
            placeholder_symbol_dto = default_name_dto.copy()
            placeholder_symbol_dto['name_type_name'] = 'nomenclature_symbol'
            placeholder_symbol_dto['format_text'] = feature.feature.name
            placeholder_symbol_dto['display_text'] = feature.feature.name
            feature.gene_symbol_dto = placeholder_symbol_dto
        # In rare cases, a gene's annotation ID has never been used as a synonym. For these, fill in the annotation ID.
        if feature.gene_systematic_name_dto is None and feature.curr_anno_id:
            placeholder_systematic_name_dto = default_name_dto.copy()
            placeholder_systematic_name_dto['name_type_name'] = 'systematic_name'
            placeholder_systematic_name_dto['format_text'] = feature.curr_anno_id
            placeholder_systematic_name_dto['display_text'] = feature.curr_anno_id
            log.warning(f"{feature}: Has annoID never used as a synonym: {feature.curr_anno_id}")
            feature.gene_systematic_name_dto = placeholder_systematic_name_dto
        return

    # Synthesis of initial db info.
    def synthesize_info(self):
        """Convert FlyBase gene data into an AllianceGene representation."""
        log.info('Synthesizing gene info.')
        for gene in self.gene_dict.values():
            log.debug(f'Evaluating annotation: {gene}')
            self.process_feature_synonyms(gene)
            # Get timestamps.
            if gene.timestamps:
                gene.date_created = strict_rfc3339.\
                    timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(gene.timestamps)))
                gene.date_updated = strict_rfc3339.\
                    timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(gene.timestamps)))
            # Get genomic_locations.
            if gene.featureloc:
                genomic_location_dict = {
                    'internal': False,
                    'obsolete': False,
                    'created_by_curie': 'FB:FB_curator',
                    'updated_by_curie': 'FB:FB_curator',
                    'genomic_entity_curie': gene.curie,
                    'predicate': 'localizes_to',
                    'chromosome_curie': 'FB:{}'.format(self.chr_dict[gene.featureloc.srcfeature_id]),
                    'assembly_curie': reference_assembly
                }
                if gene.featureloc.strand == -1:
                    genomic_location_dict['start'] = str(gene.featureloc.fmax)
                    genomic_location_dict['end'] = str(gene.featureloc.fmin + 1)
                else:
                    genomic_location_dict['start'] = str(gene.featureloc.fmin + 1)
                    genomic_location_dict['end'] = str(gene.featureloc.fmax)
                gene.genomic_location_dtos.append(genomic_location_dict)
            # Add gene synopsis.
            if gene.gene_snapshot:
                gene.gene_synopsis = gene.gene_snapshot.value
            # Get secondary IDs (FBgn and annotation IDs).
            for fb_id in gene.alt_fb_ids:
                gene.secondary_identifiers.append('FB:{}'.format(fb_id.accession))
            for anno_id in gene.annotation_ids:
                gene.secondary_identifiers.append('FB:{}'.format(anno_id.accession))
            # Get crossreferences.
            # Start by adding gene uniquename as an xref.
            xref_dict = {
                'referenced_curie': 'FB:{}'.format(gene.feature.uniquename),
                'display_name': 'FB:{}'.format(gene.feature.uniquename),
                'prefix': 'FB',
                'page_area': 'gene',
                'created_by_curie': 'FB:FB_curator',
                'obsolete': False,
                'internal': False
            }
            gene.cross_reference_dtos.append(xref_dict)
            # Then add PANTHER xref (from file).
            if gene.feature.uniquename in self.pthr_dict.keys():
                pthr_xref_dict = {
                    'referenced_curie': 'PANTHER:{}'.format(self.pthr_dict[gene.feature.uniquename]),
                    'display_name': 'PANTHER:{}'.format(self.pthr_dict[gene.feature.uniquename]),
                    'prefix': 'PANTHER',
                    'page_area': 'gene',
                    'obsolete': False,
                    'internal': False
                }
                gene.cross_reference_dtos.append(pthr_xref_dict)
            # Get other xrefs.
            for result in gene.dbxrefs:
                if result.Db.name in self.fb_agr_db_dict.keys():
                    xref_dict = {
                        'referenced_curie': '{}:{}'.format(self.fb_agr_db_dict[result.Db.name], result.Dbxref.accession),
                        'display_name': '{}:{}'.format(self.fb_agr_db_dict[result.Db.name], result.Dbxref.accession),
                        'prefix': self.fb_agr_db_dict[result.Db.name],
                        'page_area': 'gene',
                        'created_by_curie': 'FB:FB_curator',
                        'obsolete': False,
                        'internal': False
                    }
                    if result.FeatureDbxref.is_current is False:
                        xref_dict['internal'] = True
                    gene.cross_reference_dtos.append(xref_dict)
            # Flag internal features.
            if gene.obsolete is True:
                gene.internal = True
                gene.internal_reasons.append('Obsolete')
            # TEMPORARY: Suppress non-Dmel genes from export.
            if gene.organism_abbr != 'Dmel':
                gene.for_alliance_export = False
                gene.export_warnings.append(f'Suppress non-Dmel genes from export: ORG={gene.organism_abbr}')
            # Suppress objects missing required information from export.
            for attr in self.required_fields:
                if attr not in gene.__dict__.keys():
                    gene.for_alliance_export = False
                    gene.export_warnings.append('Missing "{}" attribute'.format(attr))
                elif getattr(gene, attr) is None:
                    gene.for_alliance_export = False
                    gene.export_warnings.append('Missing value for "{}" attribute'.format(attr))
            if gene.internal is False and gene.for_alliance_export is True:
                log.debug('EXPORT {}'.format(gene.curie))
        log.info('Done synthesizing gene info.')
        return

    def generate_export_file(self):
        """Process genes and print to a LinkML-compliant JSON file."""
        log.info('Generating output JSON file of genes.')
        output_dict = {
            'linkml_version': linkml_release,
            'gene_ingest_set': []
        }
        for gene in self.gene_dict.values():
            if gene.for_alliance_export is False:
                log.debug('Suppress gene from export: {}. Reasons: {}'.format(gene, '; '.join(gene.export_warnings)))
                continue
            self.export_feat_cnt += 1
            if gene.internal is True:
                self.internal_feat_cnt += 1
                log.debug('Mark gene as internal: {}. Reasons: {}'.format(gene, '; '.join(gene.internal_reasons)))
            output_gene = {}
            for attr in self.output_fields:
                if getattr(gene, attr) is not None and getattr(gene, attr) != []:
                    output_gene[attr] = getattr(gene, attr)
            output_dict['gene_ingest_set'].append(output_gene)
        log.info('Writing data to output file.')
        with open(output_filename, 'w') as outfile:
            json.dump(output_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
            outfile.close()
        log.info('Done writing data to output file.')
        total_public_feat_cnt = self.export_feat_cnt - self.internal_feat_cnt
        log.info('Exported {} of {} genes ({} are public).'.
                 format(self.export_feat_cnt, self.total_feat_cnt, total_public_feat_cnt))
        log.info('Suppressed {} genes from export.'.format(self.total_feat_cnt - self.export_feat_cnt))
        return


def db_query_transaction(object_to_execute):
    """Query the chado database given an object that has a "query_chado()" method.

    Function assumes a global "engine" variable for SQLAlchemy processes.

    Args:
        arg1 (object_to_execute): Some object that has an SQL ORM "query_chado()" method.

    Returns:
        None.

    Raises:
        Raises a RuntimeError if there are problems with executing the query.

    """
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        object_to_execute.query_chado(session)
        session.flush()
    except RuntimeError:
        session.rollback()
        log.critical('Critical transaction error occurred during chado query; rolling back and exiting.')
        raise
    return


if __name__ == "__main__":
    main()
