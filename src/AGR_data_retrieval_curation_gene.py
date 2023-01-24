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

# Now proceed with generic setup.
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
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')
    log.info('Output JSON file corresponds to "agr_curation_schema" release: {}'.format(linkml_release))

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
        # 2. Gene.taxon is required, but even after updating NCBITaxon info at FlyBase, not all genes will have NCBI taxon ID.
        # 3. GenomicLocation lacks strand info.
        self.feature = feature                                # The Feature object corresponding to the FlyBase gene.
        self.organism_abbr = None                             # Will be the organism.abbreviation for the gene's species of origin.
        self.taxon_dbxref = None                              # Will be the NCBITaxon (Db, Dbxref) tuple for the organism.
        self.featureloc = None                                # Will be Featureloc object for the gene.
        self.gene_type_name = None                            # Will be the cvterm.name for "promoted_gene_type" featureprop.
        self.gene_snapshot = None                             # Will be the "gene_summary_text" Featureprop object.
        self.curr_anno_id = None                              # Will be current annotation ID for the gene (str).
        self.curr_fb_symbol = []                              # Will be all FeatureSynonym objects in support of the current symbol Synonym object.
        self.curr_fb_fullname = []                            # Will be all FeatureSynonym objects in support of the current fullname Synonym object.
        self.systematic_name = []                             # Will be all FeatureSynonym objects using the systematic name of the gene.
        self.other_synonyms = []                              # Will be all FeatureSynonym objects in support of non-current synonyms.
        self.dbxrefs = []                                     # Will be list of dbxrefs as sql result groupings: Db, Dbxref, FeatureDbxref.
        self.alt_fb_ids = []                                  # Will be list of Dbxrefs for 2o FlyBase IDs.
        self.annotation_ids = []                              # Will be list of Dbxrefs for annotation IDs.
        self.timestamps = []                                  # Add all timestamps here.
        # Attributes for the Alliance AuditedObjectDTO.
        self.obsolete = False                                 # Never True. All FB annotations are deleted if no longer current.
        self.internal = False                                 # Will be internal if annotation should not be exported to Alliance for some reason.
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
        self.gene_full_name_dto = None                        # Will be a single GeneFullNameSlotAnnotation.
        self.gene_systematic_name_dto = None                  # Will be a single GeneSystematicNameSlotAnnotation.
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

    test_genes = ['wg', 'mt:ori', 'lncRNA:roX1', 'CG12656']
    required_fields = [
        'curie',
        'gene_symbol_dto',
        'gene_full_name_dto',
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
        'gene_type_curie',
        'genomic_location_dtos',
        'internal',
        'obsolete',
        'secondary_identifiers',
        'taxon_curie',
        'updated_by_curie',
    ]
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

    def open_panther_file(self):
        """Extract panther information from file."""
        log.info('Open PANTHER file.')
        if input_dir == '/src/input/':
            filepath = f'{input_dir}PTHR17.0_fruit_fly'
        else:
            filepath = '/data/ortholog/panther/PTHR17.0_fruit_fly'
        tsv_file = open(filepath, "r")
        tsvin = csv.reader(tsv_file, delimiter='\t')
        fb_regex = r'FBgn[0-9]{7}'
        pthr_regex = r'PTHR[0-9]{5}'
        FB = 0
        PTHR = 3
        for row in tsvin:
            fields = len(row)
            if fields:  # Ignore blank lines
                if re.search(fb_regex, row[FB]) and re.search(pthr_regex, row[PTHR]):
                    self.pthr_dict[re.search(fb_regex, row[FB]).group(0)] = re.search(pthr_regex, row[PTHR]).group(0)
        return

    def get_references(self, session):
        """Get all references."""
        log.info('Get all references.')
        # First get all current pubs having an FBrf uniquename.
        fbrf_regex = r'^(FBrf[0-9]{7}|unattributed)$'
        filters = (
            Pub.uniquename.op('~')(fbrf_regex),
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
            Pub.uniquename.op('~')(fbrf_regex),
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
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
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
        for result in organism_dbxref_results:
            organism_taxon_dict[result.OrganismDbxref.organism_id] = result.Dbxref.accession
        for gene in self.gene_dict.values():
            try:
                gene.taxon_curie = 'NCBITaxon:{}'.format(organism_taxon_dict[gene.feature.organism_id])
            except KeyError:
                log.debug('No NCBI taxon ID available for: {}'.format(gene))
        return

    def get_gene_dbxrefs(self, session):
        """Get all dbxrefs for genes. This will take 10-15 minutes."""
        log.info('Getting gene dbxrefs.')
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'gene',
            Db.name.in_((self.fb_agr_db_dict.keys()))
        )
        gene_dbxref_results = session.query(Feature, FeatureDbxref, Dbxref, Db).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in gene_dbxref_results:
            counter += 1
            if counter % 100000 == 0:
                log.debug('Processing xref #{}'.format(counter))
            # Skip current FlyBase accessions.
            # If present, these are same as feature.uniquename.
            # However, not present for all genes (e.g., FBgn0085177), so cannot be relied upon.
            if result.FeatureDbxref.is_current is True and result.Db.name == 'FlyBase':
                pass
            elif result.FeatureDbxref.is_current is False and result.Db.name == 'FlyBase':
                self.gene_dict[result.Feature.uniquename].alt_fb_ids.append(result.Dbxref)
            elif result.Db.name == 'FlyBase Annotation IDs':
                self.gene_dict[result.Feature.uniquename].annotation_ids.append(result.Dbxref)
                if result.FeatureDbxref.is_current is True:
                    self.gene_dict[result.Feature.uniquename].curr_anno_id = result.Dbxref.accession
            else:
                self.gene_dict[result.Feature.uniquename].dbxrefs.append(result)
        return

    def get_synonyms(self, session):
        """Get current and non-current symbols and full names for genes."""
        log.info('Getting gene synonyms.')
        feature_type = aliased(Cvterm, name='feature_type')
        synonym_type = aliased(Cvterm, name='synonym_type')
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
            Feature.is_analysis.is_(False),
            feature_type.name == 'gene'
        )
        gene_curr_symbol_results = session.query(synonym_type, Feature, FeatureSynonym, Synonym).\
            join(FeatureSynonym, (FeatureSynonym.synonym_id == Synonym.synonym_id)).\
            join(Feature, (Feature.feature_id == FeatureSynonym.feature_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(synonym_type, (synonym_type.cvterm_id == Synonym.type_id)).\
            filter(*filters).\
            distinct()
        for result in gene_curr_symbol_results:
            # First, build the all_synonyms_dict.
            self.all_synonyms_dict[result.Synonym.synonym_id] = Synonym
            # Second, collect FeatureSynonym objects by type.
            if result.FeatureSynonym.is_current is True:
                if result.synonym_type.name == 'symbol':
                    self.gene_dict[result.Feature.uniquename].curr_fb_symbol.append(result.FeatureSynonym)
                elif result.synonym_type.name == 'fullname':
                    self.gene_dict[result.Feature.uniquename].curr_fb_fullname.append(result.FeatureSynonym)
            else:
                self.gene_dict[result.Feature.uniquename].other_synonyms.append(result.FeatureSynonym)
            # Third, catch synonyms that match the annotation ID (aka, systematic_name).
            if result.Synonym.name == self.gene_dict[result.Feature.uniquename].curr_anno_id:
                self.gene_dict[result.Feature.uniquename].systematic_name.append(result.FeatureSynonym)
        return

    def get_gene_snapshots(self, session):
        """Get human-written gene summaries."""
        log.info('Getting gene snapshots.')
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='gene_summary_text')
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
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
        for result in gene_snapshot_results:
            self.gene_dict[result.feature.uniquename].gene_snapshot = result
        return

    def get_gene_types(self, session):
        """Get and parse "promoted_gene_type" featureprop for genes."""
        log.info('Getting gene types.')
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='promoted_gene_type')
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
            Feature.is_analysis.is_(False),
            prop_type.name == 'promoted_gene_type'
        )
        gene_type_results = session.query(Featureprop).\
            join(Feature, (Feature.feature_id == Featureprop.feature_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        for result in gene_type_results:
            self.gene_dict[result.feature.uniquename].gene_type_curie = result.value[1:10].replace('SO', 'SO:')
            self.gene_dict[result.feature.uniquename].gene_type_name = result.value[11:-1]
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
        for result in chr_results:
            self.chr_dict[result.feature_id] = result.uniquename
        # Now get gene featureloc.
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'gene'
        )
        gene_featureloc_results = session.query(Featureloc).\
            join(Feature, (Feature.feature_id == Featureloc.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        for result in gene_featureloc_results:
            self.gene_dict[result.feature.uniquename].featureloc = result
        return

    def query_chado(self, session):
        """A wrapper method that runs initial db queries."""
        self.open_panther_file()
        self.get_references(session)
        self.get_genes(session)
        self.get_gene_taxons(session)
        # self.get_gene_dbxrefs(session)    # BOB - suppress for faster dev.
        self.get_synonyms(session)
        self.get_gene_snapshots(session)
        self.get_gene_types(session)
        self.get_gene_timestamps(session)
        self.get_gene_featureloc(session)
        return

    # BOB: new method for synonyms.
    def process_feature_synonyms(self, input, name_type, return_single_value):
        """Convert a string or list of FeatureSynonym objects into single or many DTO objects for export.

        Args:
            arg1 (input): (str or list) A string, or, a list of FeatureSynonym objects.
            arg2 (name_type): (str) The type of name to return. If "unspecified" is given, go by Synonym type.
            arg3 (return_single_value): (bool) True if output should be a single DTO, False if a list is to be returned.

        Returns:
            A single or list of name DTO objects.

        Raises:
            Raises error if in put is not a string/list.
            Raises error if return_single_value set to True, but many synonyms found in the input list.

        """
        if type(input) is not str and type(input) is not list:
            log.error('Input must be a string or list of FeatureSynonym objects.')
            raise
        # First handle the simplest case where a string is given as the input.
        if type(input) is str:
            output_synonym_dto = {
                'name_type_name': name_type,
                'format_text': input,
                'display_text': input,
                'synonym_scope': 'exact',
                'evidence_curies': [],
                'internal': False,
                'obsolete': False
            }
            if return_single_value is False:
                output_synonym_dto = [output_synonym_dto]
            return output_synonym_dto
        # Next handle a list of FeatureSynonym objects.
        # Collect pub_ids for each synonym (keyed by synonym_id).
        feature_synonym_dict = {}
        output_synonym_dto_list = []
        for f_s in input:
            try:
                feature_synonym_dict[f_s.synonym_id].append(f_s.pub_id)
            except KeyError:
                feature_synonym_dict[f_s.synonym_id] = [f_s.pub_id]
        for synonym_id, pub_list in feature_synonym_dict.items():
            synonym = self.all_synonyms_dict[synonym_id]
            if name_type == 'unspecified':
                name_type_to_use = synonym.type.name
            else:
                name_type_to_use = name_type
            output_synonym_dto = {
                'name_type_name': name_type_to_use,
                'format_text': synonym.name,
                'display_text': sub_sup_sgml_to_html(synonym.synonym_sgml),
                'synonym_scope': 'exact',
                'evidence_curies': [f'{self.all_pubs_dict[i]}' for i in pub_list if self.all_pubs_dict[i] != 'unattributed'],
                'internal': False,
                'obsolete': False
            }
            output_synonym_dto_list.append(output_synonym_dto)
        if return_single_value is True and len(output_synonym_dto_list) != 1:
            log.error('Found many synonyms but was expecting only one.')
            raise
        elif return_single_value is True and len(output_synonym_dto_list) == 1:
            return output_synonym_dto_list[0]
        else:
            return output_synonym_dto_list

    # Synthesis of initial db info.
    def synthesize_info(self):
        """Convert FlyBase gene data into an AllianceGene representation."""
        log.info('Synthesizing gene info.')
        for gene in self.gene_dict.values():
            log.debug(f'Evaluating annotation: {gene}')
            # BOB: Handle synonyms.
            log.debug(f'Handle synonyms for {gene}')
            if gene.curr_fb_symbol:
                gene.gene_symbol_dto = self.process_feature_synonyms(gene.curr_fb_symbol, 'nomenclature_symbol', True)
            else:
                gene.gene_symbol_dto = self.process_feature_synonyms(gene.feature.name, 'nomenclature_symbol', True)
            if gene.curr_fb_fullname:
                gene.gene_full_name_dto = self.process_feature_synonyms(gene.curr_fb_fullname, 'full_name', True)
            else:
                gene.gene_full_name_dto = self.process_feature_synonyms(gene.feature.name, 'full_name', True)
            if gene.systematic_name:
                gene.gene_systematic_name_dto = self.process_feature_synonyms(gene.systematic_name, 'systematic_name', True)
            if gene.other_synonyms:
                gene.gene_synonym_dtos = self.process_feature_synonyms(gene.other_synonyms, 'unspecified', False)
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
                'curie': 'FB:{}'.format(gene.feature.uniquename),
                'display_name': 'FB:{}'.format(gene.feature.uniquename),
                'prefix': 'FB',
                'page_areas': ['gene'],
                'created_by_curie': 'FB:FB_curator',
                'obsolete': False,
                'internal': False
            }
            gene.cross_reference_dtos.append(xref_dict)
            # Then add PANTHER xref (from file).
            if gene.feature.uniquename in self.pthr_dict.keys():
                pthr_xref_dict = {
                    'curie': 'PANTHER:{}'.format(self.pthr_dict[gene.feature.uniquename]),
                    'display_name': 'PANTHER:{}'.format(self.pthr_dict[gene.feature.uniquename]),
                    'prefix': 'PANTHER',
                    'page_areas': ['gene'],
                    'obsolete': False,
                    'internal': False
                }
                gene.cross_reference_dtos.append(pthr_xref_dict)
            # Get other xrefs.
            for result in gene.dbxrefs:
                if result.Db.name in self.fb_agr_db_dict.keys():
                    xref_dict = {
                        'curie': '{}:{}'.format(self.fb_agr_db_dict[result.Db.name], result.Dbxref.accession),
                        'display_name': '{}:{}'.format(self.fb_agr_db_dict[result.Db.name], result.Dbxref.accession),
                        'prefix': self.fb_agr_db_dict[result.Db.name],
                        'page_areas': ['gene'],
                        'created_by_curie': 'FB:FB_curator',
                        'obsolete': False,
                        'internal': False
                    }
                    if result.FeatureDbxref.is_current is False:
                        xref_dict['internal'] = True
                    gene.cross_reference_dtos.append(xref_dict)
            # Flag internal features.
            if gene.organism_abbr != 'Dmel':
                gene.internal = True
                gene.internal_reasons.append('Non-Dmel')
            if gene.obsolete is True:
                gene.internal = True
                gene.internal_reasons.append('Obsolete')
            if gene.gene_type_curie is None:
                gene.internal = True
                gene.internal_reasons.append('Lacks gene type')
            if gene.gene_type_name in self.internal_gene_types:
                gene.internal = True
                gene.internal_reasons.append('Internal gene type {} ({})'.format(gene.gene_type_name, gene.gene_type_curie))
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
