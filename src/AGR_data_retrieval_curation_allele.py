# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase alleles for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_curation_allele.py [-h] [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_curation_allele.py -v -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script makes a JSON file conforming to LinkML specs for the curation
    (i.e., "persistent") database; distinct from ALLELE file specs for the
    original Neo4j drop-and-reload database. As the "audit_chado" table may be
    required for proper determination of "date_updated" values, use a reporting
    database with full audit_chado table present.

"""

import argparse
import datetime
import json
import strict_rfc3339
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Feature, FeatureCvterm, FeatureDbxref, FeatureGenotype,
    FeaturePub, FeatureRelationship, FeatureSynonym, Featureprop, Genotype,
    Library, LibraryFeature, LibraryFeatureprop, Organism, Organismprop,
    OrganismDbxref, Phenotype, PhenotypeCvterm, Phenstatement, Pub, PubDbxref,
    Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Now proceed with generic setup.
report_label = 'allele_curation'
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
    """Run the steps for exporting LinkML-compliant FlyBase Alleles."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')
    log.info('Output JSON file corresponds to "agr_curation_schema" release: {}'.format(linkml_release))

    # Instantiate the object, get the data, synthesize it, export it.
    allele_handler = AlleleHandler()
    db_query_transaction(allele_handler)
    allele_handler.synthesize_info()
    allele_handler.generate_export_file()
    log.info('Ended main function.\n')


class AllianceAllele(object):
    """An allele with it's associated FlyBase data and Alliance LinkML properties."""
    def __init__(self, feature):
        """Create a base AllianceAllele object.

        Args:
            arg1 (feature): (Feature) The Feature object corresponding to the allele.

        Returns:
            An object of the AllianceAllele class.

        """
        # Attributes representing unprocessed FlyBase data.
        # Note: use attribute names that do not match an Alliance LinkML slot name.
        # For initial load, the Alliance A-Team just needs minimum info.
        # ALLELE: curie, taxon, symbol, description, internal, obsolete.
        # Problems with Allele LinkML:
        # 1. Allele.taxon is required, but even after updating NCBITaxon info at FlyBase, not all alleles will have NCBI taxon ID.
        self.feature = feature                    # The Feature object corresponding to the FlyBase allele.
        self.organism_abbr = None                 # Will be the organism.abbreviation for the allele's species of origin.
        self.adj_organism_abbr = 'Dmel'           # Assume allele is Dmel (classical or transgenic) unless allele is of classical type in another insect.
        self.in_vitro = False                     # Change to True if allele associated with "in vitro%" cvterm.
        self.constructs = []                      # Will be a list of FBtp IDs for this allele's constructs.
        self.dmel_insertions = []                 # Will be a list of FBti IDs for this allele's Dmel insertions.
        self.non_dmel_insertions = []             # Will be a list of FBti IDs for this allele's non-Dmel insertions.
        self.args = []                            # Will be a list of ARGs Features (variants).
        self.parent_gene = None                   # Will be the FBgn ID of the allele's gene.
        self.allele_of_internal_gene = False      # Will change to True if is allele of Dmel internal-type gene (e.g., origin_of_replication).
        self.taxon_dbxref = None                  # Will be the NCBITaxon (Db, Dbxref) tuple for the organism.
        self.curr_fb_symbol = None                # Will be the current symbol Synonym object.
        self.curr_fb_fullname = None              # Will be the current fullname Synonym object.
        self.internal_synonyms = []               # Will be list of internal synonym names (and synonym_sgml if different).
        self.public_synonyms = []                 # Will be list of public synonym names (and synonym_sgml if different).
        self.dbxrefs = []                         # Will be list of dbxrefs as sql result groupings: Db, Dbxref, FeatureDbxref.
        self.alt_fb_ids = []                      # Will be list of Dbxrefs for 2o FlyBase IDs.
        self.timestamps = []                      # Add all timestamps here.
        self.fb_references = []                   # Will be list of FBrf IDs related to an allele: directly and indirectly.
        self.featureprops = {}                    # A CVterm-keyed dict of Featureprop lists.
        self.phenotypes = []                      # Will be a list of SQLAlchemy (Feature, Genotype, Phenotype, Cvterm) results.
        self.direct_libraries = []                # Will be a list of Library objects directly related to the allele.
        self.ins_libraries = []                   # Will be a list of Library objects related to the allele via insertion (FBti).
        self.cons_libraries = []                  # Will be a list of Library objects related to the allele via construct (FBtp).
        self.sf_libraries = []                    # Will be a list of Library objects related to the allele via seq. feature (FBsf).
        # Attributes for the Alliance AuditedObject.
        self.obsolete = feature.is_obsolete       # Will be the FlyBase value here.
        self.internal = False                     # Change to true if allele not intended for display at Alliance website.
        self.created_by = 'FB:FB_curator'         # Use placeholder value since no Person object at FlyBase.
        self.updated_by = 'FB:FB_curator'         # Use placeholder value since no Person object at FlyBase.
        self.date_created = None                  # Earliest timestamp.
        self.date_updated = None                  # Latest timestamp.
        # self.data_provider = 'FB'                 # The MOD abbreviation.
        # Attributes for the Alliance BiologicalEntity. BiologicalEntity is_a AuditedObject.
        self.curie = 'FB:{}'.format(feature.uniquename)
        self.taxon = None                         # A string representing the NCBI taxon ID. We have no NCBI taxonID for 223 alleles.
        # Attributes for the Alliance GenomicEntity. GenomicEntity is_a BiologicalEntity.
        self.name = None                          # Will be current fullname synonym - report ascii or utf8 (sgml) version?
        self.synonyms = []                        # All current and non-current ASCII and SGML synonyms.
        self.cross_references = []                # Report only select dbs, using AGR-accepted db_prefix.
        self.secondary_identifiers = []           # Annotation IDs and 2o FlyBase IDs.
        # Attributes for the Alliance Allele. Allele is_a GenomicEntity.
        self.symbol = None                        # Will be a string (ascii or utf8)?
        self.references = []                      # KANBAN-237: READY: Will be a list of pubs (PMID or FB:FBrf IDs) for the allele.
        self.is_extinct = None                    # KANBAN-237: READY: Change to true if extinction has been reported. Otherwise, leave blank.
        self.inheritence_mode = []                # KANBAN-237: READY: Will be a list of CV terms.
        self.in_collection = []                   # KANBAN-237: TO DO: Will be a library names.
        self.sequencing_status = None             # KANBAN-237: TO DO: Will be a CV term? TBD. Might be dropped.
        # Notes associated with the object.
        self.for_alliance_export = True           # Change to False if object should be excluded from export.
        self.internal_reasons = []                # Reasons for marking an object as internal in the export file.
        self.export_warnings = []                 # Reasons for suppressing an object from the export file.

    def __str__(self):
        """Succinct text string describing the AllianceAllele object."""
        desc = '{} ({})'.format(self.feature.name, self.feature.uniquename)
        return desc


class AlleleHandler(object):
    """This object gets alleles and related info, synthesizes/filters the data, then exports it as LinkML JSON."""
    def __init__(self):
        """Create the AlleleHandler object."""
        self.allele_dict = {}         # An FBalID-keyed dict of AllianceAllele objects.
        self.drosophilid_list = []    # A list of organism_ids for "Drosophilid" species in chado.
        self.total_feat_cnt = 0       # Count of all alleles found in starting query.
        self.export_feat_cnt = 0      # Count of all alleles exported to file.
        self.internal_feat_cnt = 0    # Count of all alleles marked as internal=True in export file.
        self.fbrf_pmid_dict = {}      # Will be a dict of FBrf-to-PMID xrefs.

    test_alleles = []
    required_fields = [
        'curie',
        'taxon',
        'symbol',
        'internal'
    ]
    output_fields = [
        'created_by',
        'cross_references',
        'curie',
        'date_created',
        'date_updated',
        'in_collection',
        'inheritence_mode',
        'internal',
        'is_extinct',
        'updated_by',
        'name',
        'obsolete',
        'references',
        'secondary_identifiers',
        # 'sequencing_status',    # KANBAN-237: Not implemented yet in LinkML v1.2.4
        'symbol',
        'synonyms',
        'taxon'
    ]
    fb_agr_db_dict = {
        'FlyBase': 'FB'
    }

    def get_alleles(self, session):
        """Get all alleles."""
        log.info('Querying chado for alleles.')

        # First get all allele features from chado.
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'allele'
        )
        allele_results = session.query(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            filter(*filters).\
            distinct()
        self.total_feat_cnt = 0
        for result in allele_results:
            self.total_feat_cnt += 1
            self.allele_dict[result.uniquename] = AllianceAllele(result)
            self.allele_dict[result.uniquename].organism_abbr = result.organism.abbreviation
        log.info('Found {} alleles.'.format(self.total_feat_cnt))

    def get_allele_gene(self, session):
        """For current alleles, get the FBgn ID of allele's current gene."""
        gene = aliased(Feature, name='gene')
        allele = aliased(Feature, name='allele')
        gene_regex = r'^FBgn[0-9]{7}$'
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            gene.is_obsolete.is_(False),
            gene.is_analysis.is_(False),
            gene.uniquename.op('~')(gene_regex),
            allele.is_obsolete.is_(False),
            allele.is_analysis.is_(False),
            allele.uniquename.op('~')(allele_regex),
            Cvterm.name == 'alleleof'
        )
        allele_gene_results = session.query(allele, gene).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_gene_results:
            counter += 1
            self.allele_dict[result.allele.uniquename].parent_gene = result.gene.uniquename
        log.info('Found {} allele-gene relationships.'.format(counter))
        return

    def flag_alleles_of_internal_genes(self, session):
        """Flag alleles belonging to internal-type Dmel genes like origin_of_replication."""
        log.info('Flagging alleles of internal-type Dmel genes.')
        # First build list of internal Dmel genes.
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
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='promoted_gene_type')
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(gene_regex),
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            Organism.abbreviation == 'Dmel',
            prop_type.name == 'promoted_gene_type'
        )
        gene_type_results = session.query(Featureprop).\
            join(Feature, (Feature.feature_id == Featureprop.feature_id)).\
            join(Organism, (Organism.organism_id == Feature.organism_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        gene_counter = 0
        internal_dmel_genes = []
        for result in gene_type_results:
            gene_type_name = result.value[11:-1]
            if gene_type_name in internal_gene_types:
                gene_counter += 1
                internal_dmel_genes.append(result.feature.uniquename)
        log.info('Found {} internal type Dmel genes.'.format(gene_counter))
        # Now check alleles against this internal Dmel gene list.
        allele_counter = 0
        for allele in self.allele_dict.values():
            if allele.parent_gene in internal_dmel_genes:
                allele_counter += 1
                allele.allele_of_internal_gene = True
        log.info('Found {} alleles of internal type Dmel genes.'.format(allele_counter))
        return

    def flag_in_vitro_alleles(self, session):
        """Flag alleles associated with "in vitro" type CV term."""
        log.info('Flag in vitro alleles.')
        allele_regex = r'^FBal[0-9]{7}$'
        cvterm_name_regex = '^in vitro construct'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Cvterm.name.op('~')(cvterm_name_regex)
        )
        ivt_alleles = session.query(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for allele in ivt_alleles:
            counter += 1
            self.allele_dict[allele.uniquename].in_vitro = True
        log.info(f'Flagged {counter} alleles as "in vitro"')
        return

    def get_allele_constructs(self, session):
        """Find FBtp constructs associated with the allele."""
        log.info('Find constructs related to alleles.')
        allele = aliased(Feature, name='allele')
        construct = aliased(Feature, name='construct')
        allele_regex = r'^FBal[0-9]{7}$'
        construct_regex = r'^FBtp[0-9]{7}$'
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(allele_regex),
            construct.is_obsolete.is_(False),
            construct.uniquename.op('~')(construct_regex)
        )
        construct_results = session.query(allele, construct).\
            join(FeatureRelationship, (FeatureRelationship.object_id == construct.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.subject_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in construct_results:
            counter += 1
            self.allele_dict[result.allele.uniquename].constructs.append(result.construct)
        log.info('Found {} constructs related to alleles.'.format(counter))
        return

    def get_allele_insertions(self, session):
        """Find FBti insertions associated with the allele."""
        log.info('Find insertions related to alleles.')
        allele = aliased(Feature, name='allele')
        insertion = aliased(Feature, name='insertion')
        allele_regex = r'^FBal[0-9]{7}$'
        insertion_regex = r'^FBti[0-9]{7}$'
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(allele_regex),
            insertion.is_obsolete.is_(False),
            insertion.uniquename.op('~')(insertion_regex)
        )
        insertion_results = session.query(Organism, allele, insertion).\
            join(FeatureRelationship, (FeatureRelationship.object_id == insertion.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.subject_id)).\
            join(Organism, (Organism.organism_id == insertion.organism_id)).\
            filter(*filters).\
            distinct()
        dmel_counter = 0
        non_dmel_counter = 0
        for result in insertion_results:
            if result.Organism.abbreviation == 'Dmel':
                self.allele_dict[result.allele.uniquename].dmel_insertions.append(result.insertion)
                dmel_counter += 1
            else:
                self.allele_dict[result.allele.uniquename].non_dmel_insertions.append(result.insertion)
                non_dmel_counter += 1
        log.info(f'Found {dmel_counter} Dmel and {non_dmel_counter} non-Dmel insertions related to alleles.')
        return

    def get_drosophilid_organisms(self, session):
        """Find organisms for Drosophilid species."""
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
        log.info(f'Found {counter} Drosophilid organisms in chado.')
        return

    def adjust_allele_org(self, session):
        """Find classical alleles in non-Dmel species and adjust organism accordingly."""
        log.info('Adjust organism for classical non-Dmel alleles.')
        counter = 0
        for allele in self.allele_dict.values():
            # Skip alleles that unambiguously occur in Dmel.
            if allele.organism_abbr == 'Dmel':
                continue
            if allele.dmel_insertions:
                continue
            # How to classify alleles like FBal0048225? They look like non-Dmel TE transgenic alleles carried in Dmel.
            # FBal0048225 has mutagen=natural population, but still seems to be transgenic.
            if allele.constructs:
                continue
            # Find clear evidence that allele is non-Dmel classical allele.
            is_non_dmel_classical = False
            if allele.non_dmel_insertions:
                is_non_dmel_classical = True
            elif allele.feature.organism_id in self.drosophilid_list and allele.in_vitro is False:
                is_non_dmel_classical = True
            if is_non_dmel_classical is True:
                allele.adj_organism_abbr = allele.organism_abbr
                log.debug(f'Non-Dmel allele: id={allele.curie}, name={allele.feature.name}, org_abbr={allele.organism_abbr}')
                counter += 1
        log.info('Adjusted organism to be "non-Dmel" for {} alleles.'.format(counter))
        return

    def get_allele_taxons(self, session):
        """Get taxon IDs for alleles. Depends on all organisms for features having an abbreviation."""
        log.info('Getting allele taxon IDs.')
        filters = (
            OrganismDbxref.is_current.is_(True),
            Db.name == 'NCBITaxon'
        )
        organism_dbxref_results = session.query(Organism, Dbxref).\
            join(OrganismDbxref, (OrganismDbxref.organism_id == Organism.organism_id)).\
            join(Dbxref, (Dbxref.dbxref_id == OrganismDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        organism_taxon_dict = {}
        for result in organism_dbxref_results:
            organism_taxon_dict[result.Organism.abbreviation] = result.Dbxref.accession
        for allele in self.allele_dict.values():
            try:
                allele.taxon = 'NCBITaxon:{}'.format(organism_taxon_dict[allele.adj_organism_abbr])
            except KeyError:
                log.debug('No NCBI taxon ID available for: {}'.format(allele))
        return

    def get_synonyms(self, session):
        """Get current and non-current symbols and full names for alleles."""
        log.info('Getting allele synonyms.')
        feature_type = aliased(Cvterm, name='feature_type')
        synonym_type = aliased(Cvterm, name='synonym_type')
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Feature.is_analysis.is_(False),
            feature_type.name == 'allele'
        )
        allele_curr_symbol_results = session.query(synonym_type, Feature, FeatureSynonym, Synonym).\
            join(FeatureSynonym, (FeatureSynonym.synonym_id == Synonym.synonym_id)).\
            join(Feature, (Feature.feature_id == FeatureSynonym.feature_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(synonym_type, (synonym_type.cvterm_id == Synonym.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_curr_symbol_results:
            counter += 1
            if result.FeatureSynonym.is_current is True:
                if result.synonym_type.name == 'symbol':
                    self.allele_dict[result.Feature.uniquename].curr_fb_symbol = result.Synonym
                elif result.synonym_type.name == 'fullname':
                    self.allele_dict[result.Feature.uniquename].curr_fb_fullname = result.Synonym
            elif result.FeatureSynonym.is_internal is True:
                self.allele_dict[result.Feature.uniquename].internal_synonyms.append(result.Synonym.name)
                self.allele_dict[result.Feature.uniquename].internal_synonyms.append(sub_sup_sgml_to_html(result.Synonym.synonym_sgml))
            else:
                self.allele_dict[result.Feature.uniquename].public_synonyms.append(result.Synonym.name)
                self.allele_dict[result.Feature.uniquename].public_synonyms.append(sub_sup_sgml_to_html(result.Synonym.synonym_sgml))
        log.info('Found {} allele synonyms.'.format(counter))
        return

    def get_allele_timestamps(self, session):
        """Get timestamps for alleles."""
        log.info('Getting allele timestamps.')
        # Start with timestamps in feature itself.
        for allele in self.allele_dict.values():
            allele.timestamps.append(allele.feature.timeaccessioned)
            allele.timestamps.append(allele.feature.timelastmodified)
        # We may want to get timestamps for other allele attibutes as well (for "date_updated").
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
        #         self.allele_dict[row[UNIQUE_KEY]].timestamps.append(row[TIMESTAMP])
        #     except KeyError:
        #         # log.debug('Could not put this in anno dict: {}'.format(row))
        #         pass
        ########################################################################
        return

    def get_allele_dbxrefs(self, session):
        """Get all dbxrefs for alleles."""
        log.info('Getting allele dbxrefs.')
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Feature.is_analysis.is_(False),
            Cvterm.name == 'allele',
            Db.name.in_((self.fb_agr_db_dict.keys()))
        )
        allele_dbxref_results = session.query(Feature, FeatureDbxref, Dbxref, Db).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_dbxref_results:
            # Skip current FlyBase accessions (i.e., equiavlent to feature.uniquename)
            #     because not all features have them (e.g., FBal0137236).
            # Better to use current uniquename as source for current FlyBase accession.
            counter += 1
            if result.FeatureDbxref.is_current is True and result.Db.name == 'FlyBase':
                pass
            elif result.FeatureDbxref.is_current is False and result.Db.name == 'FlyBase':
                self.allele_dict[result.Feature.uniquename].alt_fb_ids.append(result.Dbxref)
            else:
                self.allele_dict[result.Feature.uniquename].dbxrefs.append(result)
        log.info('Found {} allele crossreferences.'.format(counter))
        return

    def get_references(self, session):
        """Get references for alleles."""
        log.info('Get allele references.')
        allele_regex = r'^FBal[0-9]{7}$'
        fbrf_regex = r'^FBrf[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Pub.uniquename.op('~')(fbrf_regex),
            Pub.is_obsolete.is_(False)
        )
        allele_pubs = session.query(Feature, Pub).\
            join(FeaturePub, (FeaturePub.feature_id == Feature.feature_id)).\
            join(Pub, (Pub.pub_id == FeaturePub.pub_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_pubs:
            self.allele_dict[result.Feature.uniquename].fb_references.append(result.Pub.uniquename)
            counter += 1
        log.info(f'Found {counter} allele-pub relationships.')
        return

    def get_pmid_xrefs(self, session):
        """Create a dict of FBrf to PMID for publications."""
        log.info('Getting PMID IDs for FB publications.')
        filters = (
            Db.name == 'pubmed',
            Pub.is_obsolete.is_(False),
            PubDbxref.is_current.is_(True)
        )
        pmid_xrefs = session.query(Pub, Dbxref).\
            join(PubDbxref, (PubDbxref.pub_id == Pub.pub_id)).\
            join(Dbxref, (Dbxref.dbxref_id == PubDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        for xref in pmid_xrefs:
            self.fbrf_pmid_dict[xref.Pub.uniquename] = xref.Dbxref.accession
        log.info(f'Found {len(self.fbrf_pmid_dict.keys())} PMID IDs for FB publications.')
        return

    def get_allele_featureprops(self, session):
        """Get all allele featureprops."""
        log.info('Get allele featureprops.')
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Cvterm.is_obsolete == 0
        )
        allele_fprops = session.query(Feature, Cvterm, Featureprop).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in allele_fprops:
            counter += 1
            try:
                self.allele_dict[result.Feature.uniquename].\
                    featureprops[result.Cvterm.name].append(result.Featureprop)
            except KeyError:
                self.allele_dict[result.Feature.uniquename].\
                    featureprops[result.Cvterm.name] = [result.Featureprop]
        log.info(f'Found {counter} allele featureprops.')
        return

    def get_args(self, session):
        """Get ARGs related to alleles."""
        log.info('Get allele ARGs.')
        allele_regex = r'^FBal[0-9]{7}$'
        arg_types = [
            'MNV',
            'complex_substitution',
            'deletion',
            'delins',
            'insertion',
            'point_mutation',
            'sequence_alteration',
            'sequence_variant',
            'rescue_region'
        ]
        allele = aliased(Feature, name='allele')
        arg = aliased(Feature, name='arg')
        argtype = aliased(Cvterm, name='argtype')
        reltype = aliased(Cvterm, name='reltype')
        filters = (
            allele.uniquename.op('~')(allele_regex),
            arg.is_obsolete.is_(False),
            argtype.name.in_((arg_types)),
            reltype.name == 'partof'
        )
        arg_results = session.query(arg, allele).\
            join(argtype, (argtype.cvterm_id == arg.type_id)).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == arg.feature_id)).\
            join(allele, (allele.feature_id == FeatureRelationship.object_id)).\
            join(reltype, (reltype.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in arg_results:
            self.allele_dict[result.allele.uniquename].args.append(result.arg)
            counter += 1
        log.info(f'Found {counter} ARG-allele relationships.')
        return

    def get_phenotypes(self, session):
        """Get phenotypes related to alleles."""
        log.info('Get phenotypes related to alleles.')
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Genotype.is_obsolete.is_(False)
        )
        results = session.query(Feature, Genotype, Phenotype, Cvterm).\
            join(FeatureGenotype, (FeatureGenotype.feature_id == Feature.feature_id)).\
            join(Genotype, (Genotype.genotype_id == FeatureGenotype.genotype_id)).\
            join(Phenstatement, (Phenstatement.genotype_id == Genotype.genotype_id)).\
            join(Phenotype, (Phenotype.phenotype_id == Phenstatement.phenotype_id)).\
            join(PhenotypeCvterm, (PhenotypeCvterm.phenotype_id == Phenotype.phenotype_id)).\
            join(Cvterm, (Cvterm.cvterm_id == PhenotypeCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for rst in results:
            self.allele_dict[rst.Feature.uniquename].phenotypes.append(rst)
            counter += 1
        log.info(f'Found {counter} allele-phenotype associations.')
        return

    def get_direct_collections(self, session):
        """Find library collections directly related to alleles."""
        log.info('Get directly-related allele collections.')
        counter = 0
        allele_regex = r'^FBal[0-9]{7}$'
        lib_regex = r'^FBlc[0-9]{7}$'
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        # First, look for direct FBal-FBlc associations.
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(lib_regex),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection'
        )
        libraries = session.query(Feature, Library).\
            join(LibraryFeature, (LibraryFeature.feature_id == Feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in libraries:
            self.allele_dict[result.Feature.uniquename].direct_libraries.append(result.Library)
            counter += 1
        log.info(f'Found {counter} direct allele-library associations.')
        return

    def get_indirect_collections(self, session):
        """Find library collections indirectly related to alleles via insertion or construct."""
        log.info('Get indirectly-related allele collections (via insertion or construct).')
        allele_regex = r'^FBal[0-9]{7}$'
        feature_regex = r'^FB(tp|ti)[0-9]{7}$'
        lib_regex = r'^FBlc[0-9]{7}$'
        allele = aliased(Feature, name='allele')
        feature = aliased(Feature, name='feature')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        featreltype = aliased(Cvterm, name='featreltype')
        filters = (
            allele.uniquename.op('~')(allele_regex),
            feature.uniquename.op('~')(feature_regex),
            feature.is_obsolete.is_(False),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(lib_regex),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection',
            featreltype.name == 'associated_with'
        )
        indirect_libraries = session.query(allele, feature, Library).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(featreltype, (featreltype.cvterm_id == FeatureRelationship.type_id)).\
            join(feature, (feature.feature_id == FeatureRelationship.object_id)).\
            join(LibraryFeature, (LibraryFeature.feature_id == feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        fbti_counter = 0
        fbtp_counter = 0
        for result in indirect_libraries:
            if result.feature.uniquename.startswith('FBti'):
                self.allele_dict[result.allele.uniquename].ins_libraries.append(result.Library)
                fbti_counter += 1
            elif result.feature.uniquename.startswith('FBtp'):
                self.allele_dict[result.allele.uniquename].cons_libraries.append(result.Library)
                fbtp_counter += 1
        log.info(f'Found {fbti_counter} insertion-mediated allele-library associations.')
        log.info(f'Found {fbtp_counter} construct-mediated allele-library associations.')
        return

    def get_sf_collections(self, session):
        """Find library collections indirectly related to alleles via sequence feature."""
        log.info('Get indirectly-related allele collections (via equence feature).')
        allele_regex = r'^FBal[0-9]{7}$'
        cons_regex = r'^FBtp[0-9]{7}$'
        sf_regex = r'^FBsf[0-9]{10}$'
        lib_regex = r'^FBlc[0-9]{7}$'
        allele = aliased(Feature, name='allele')
        construct = aliased(Feature, name='construct')
        seqfeat = aliased(Feature, name='seqfeat')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        featreltype = aliased(Cvterm, name='featreltype')
        allele_construct = aliased(FeatureRelationship, name='allele_construct')
        seqfeat_construct = aliased(FeatureRelationship, name='seqfeat_construct')
        filters = (
            allele.uniquename.op('~')(allele_regex),
            construct.uniquename.op('~')(cons_regex),
            seqfeat.uniquename.op('~')(sf_regex),
            construct.is_obsolete.is_(False),
            seqfeat.is_obsolete.is_(False),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(lib_regex),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection',
            featreltype.name == 'associated_with'
        )
        sf_libraries = session.query(allele, Library).\
            join(allele_construct, (allele_construct.subject_id == allele.feature_id)).\
            join(construct, (construct.feature_id == allele_construct.object_id)).\
            join(featreltype, (featreltype.cvterm_id == allele_construct.type_id)).\
            join(seqfeat_construct, (seqfeat_construct.object_id == construct.feature_id)).\
            join(seqfeat, (seqfeat.feature_id == seqfeat_construct.subject_id)).\
            join(LibraryFeature, (LibraryFeature.feature_id == seqfeat.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in sf_libraries:
            self.allele_dict[result.allele.uniquename].sf_libraries.append(result.Library)
            counter += 1
        log.info(f'Found {counter} sequence feature-mediated allele-library associations.')
        return

    def query_chado(self, session):
        """A wrapper method that runs initial db queries."""
        self.get_alleles(session)
        self.get_direct_collections(session)
        self.get_indirect_collections(session)
        self.get_sf_collections(session)
        self.get_allele_gene(session)
        self.flag_alleles_of_internal_genes(session)
        self.flag_in_vitro_alleles(session)
        self.get_allele_constructs(session)
        self.get_allele_insertions(session)
        self.get_drosophilid_organisms(session)
        self.adjust_allele_org(session)
        self.get_allele_taxons(session)
        self.get_synonyms(session)
        self.get_allele_timestamps(session)
        self.get_allele_dbxrefs(session)
        self.get_references(session)
        self.get_pmid_xrefs(session)
        self.get_allele_featureprops(session)
        self.get_args(session)
        self.get_phenotypes(session)
        return

    # Synthesis of initial db info.
    def synthesize_timestamps(self, allele):
        """Process timestamps for an allele."""
        if allele.timestamps:
            allele.date_created = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(allele.timestamps)))
            allele.date_updated = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(allele.timestamps)))
        return

    def synthesize_symbol(self, allele):
        """Process symbol for an allele."""
        if allele.curr_fb_symbol:
            allele.symbol = sub_sup_sgml_to_html(allele.curr_fb_symbol.synonym_sgml)
        else:
            allele.symbol = allele.feature.name
        return

    def synthesize_fullname(self, allele):
        """Process allele fullname."""
        if allele.curr_fb_fullname:
            allele.name = sub_sup_sgml_to_html(allele.curr_fb_fullname.synonym_sgml)
        elif allele.curr_fb_symbol:
            allele.name = sub_sup_sgml_to_html(allele.curr_fb_symbol.synonym_sgml)
        else:
            allele.name = allele.feature.name
        return

    def synthesize_synonyms(self, allele):
        """Process allele synonyms."""
        internal_synonym_set = set(allele.internal_synonyms)
        for internal_synonym in internal_synonym_set:
            internal_synonym_dict = {
                'name': internal_synonym,
                'created_by': 'FB:FB_curator',
                'obsolete': False,
                'internal': True
            }
            allele.synonyms.append(internal_synonym_dict)
        public_synonym_set = set(allele.public_synonyms)
        for public_synonym in public_synonym_set:
            public_synonym_dict = {
                'name': public_synonym,
                'created_by': 'FB:FB_curator',
                'obsolete': False,
                'internal': False
            }
            allele.synonyms.append(public_synonym_dict)
        return

    def synthesize_secondary_ids(self, allele):
        """Process 2o IDs."""
        for fb_id in allele.alt_fb_ids:
            allele.secondary_identifiers.append('FB:{}'.format(fb_id.accession))
        return

    def synthesize_xrefs(self, allele):
        """Process xrefs."""
        # Start by adding allele uniquename as an xref.
        xref_dict = {
            'curie': 'FB:{}'.format(allele.feature.uniquename),
            'display_name': 'FB:{}'.format(allele.feature.uniquename),
            'prefix': 'FB',
            'page_areas': ['allele'],
            'created_by': 'FB:FB_curator',
            'obsolete': False,
            'internal': False
        }
        allele.cross_references.append(xref_dict)
        # Add other xrefs.
        for result in allele.dbxrefs:
            if result.Db.name in self.fb_agr_db_dict.keys():
                xref_dict = {
                    'curie': '{}:{}'.format(self.fb_agr_db_dict[result.Db.name], result.Dbxref.accession),
                    'display_name': '{}:{}'.format(self.fb_agr_db_dict[result.Db.name], result.Dbxref.accession),
                    'prefix': self.fb_agr_db_dict[result.Db.name],
                    'page_areas': ['allele'],
                    'created_by': 'FB:FB_curator',
                    'obsolete': False,
                    'internal': False
                }
                if result.FeatureDbxref.is_current is False:
                    xref_dict['internal'] = True
                allele.cross_references.append(xref_dict)
        return

    def synthesize_references(self, allele):
        """Process pubs for allele."""
        for fbrf_id in allele.fb_references:
            try:
                allele.references.append(f'PMID:{self.fbrf_pmid_dict[fbrf_id]}')
            except KeyError:
                allele.references.append(f'FB:{fbrf_id}')
        return

    def synthesize_insertions(self, allele):
        """Process insertions."""
        for insertion in allele.dmel_insertions:
            xref_dict = {
                'curie': '{}:{}'.format('FB', insertion.uniquename),
                'display_name': '{}:{}'.format('FB', insertion.uniquename),
                'prefix': 'FB',
                'page_areas': ['allele'],
                'created_by': 'FB:FB_curator',
                'obsolete': False,
                'internal': False
            }
            allele.cross_references.append(xref_dict)
        return

    def flag_internal_alleles(self, allele):
        """Flag alleles as internal and/or obsolete, or not."""
        if allele.organism_abbr != 'Dmel':
            allele.internal = True
            allele.internal_reasons.append('Non-Dmel')
        if allele.obsolete is True:
            allele.internal = True
            allele.internal_reasons.append('Obsolete')
        if allele.allele_of_internal_gene is True:
            allele.internal = True
            allele.internal_reasons.append('Allele of internal Dmel gene type.')
        return

    def flag_unexportable_alleles(self, allele):
        """Flag alleles missing data required for export."""
        for attr in self.required_fields:
            if attr not in allele.__dict__.keys():
                allele.for_alliance_export = False
                allele.export_warnings.append('Missing "{}" attribute'.format(attr))
            elif getattr(allele, attr) is None:
                allele.for_alliance_export = False
                allele.export_warnings.append('Missing value for "{}" attribute'.format(attr))
        if allele.internal is False and allele.for_alliance_export is True:
            log.debug('EXPORT {}'.format(allele.curie))
        return

    def synthesize_extinction(self, allele):
        """Determine if allele is definitively extinct."""
        try:
            for fprop in allele.featureprops['availability']:
                if fprop.value == 'Stated to be lost.':
                    allele.is_extinct = True
        except KeyError:
            pass
        return

    def synthesize_inheritance_mode(self, allele):
        """Determine inheritance mode for the allele."""
        inheritance_modes = {
            'recessive': 'recessive',
            'dominant': 'dominant',
            'semidominant': 'semi-dominant',
            'codominant': 'codominant'
        }
        reported_modes = []
        mode_context_list = []
        for phenotype in allele.phenotypes:
            cvterm = phenotype.Cvterm.name
            if cvterm in inheritance_modes.keys():
                # Start with assumption that it is not a single allele genotype
                single_allele_genotype = False
                # First weed out multi locus genotypes.
                if ' ' in phenotype.Genotype.uniquename:
                    single_allele_genotype = False
                elif '_' in phenotype.Genotype.description:
                    single_allele_genotype = False
                # For single locus genotype, check for hemi- or homo- state.
                else:
                    features = phenotype.Genotype.description.split('|')
                    if features[0] == features[1]:
                        single_allele_genotype = True
                    else:
                        for feature in features:
                            if feature == '+':
                                single_allele_genotype = True
                            elif feature.endswith('[+]'):
                                single_allele_genotype = True
                if single_allele_genotype is True:
                    reported_modes.append(inheritance_modes[cvterm])
                    geno = phenotype.Genotype.uniquename
                    pheno = phenotype.Phenotype.uniquename
                    mode_context = f'{allele.curie}\t{cvterm}\t{geno}\t{pheno}'
                    mode_context_list.append(mode_context)
        if reported_modes:
            reported_modes = list(set(reported_modes))
            allele.inheritence_mode = '|'.join(reported_modes)
            log.debug(f'\tFound {len(reported_modes)} inheritance mode(s): {allele.curie}: {allele.inheritence_mode}')
            # Log cases of multiple inheritance modes for curator review.
            if len(reported_modes) > 1:
                mode_context_list = list(set(mode_context_list))
                for i in mode_context_list:
                    log.warning(f'MULTIPLE_INHERITANCE_MODES:\t{i}')
        else:
            allele.inheritence_mode = 'unknown'
        return

    def synthesize_collections(self, allele):
        """Get names of collections to which allele belong (directly or indirectly)."""
        collection_names = None
        if allele.direct_libraries:
            collection_names = allele.direct_libraries
        elif allele.ins_libraries:
            collection_names = allele.ins_libraries
        elif allele.cons_libraries:
            collection_names = allele.cons_libraries
        elif allele.sf_libraries:
            collection_names = allele.sf_libraries
        if collection_names:
            collection_names = list(set(collection_names))
            # allele.in_collection = '|'.join(lib.name for lib in collection_names)
            allele.in_collection = collection_names[0].name
            log.debug(f'\tFound {len(collection_names)} collection(s): {allele.curie}: {allele.in_collection}')
        return

    def synthesize_info(self):
        """Convert FlyBase allele data into an AllianceAllele representation."""
        log.info('Synthesizing allele info.')
        for allele in self.allele_dict.values():
            log.debug('Evaluating annotation: {}'.format(allele))
            self.synthesize_collections(allele)
            self.synthesize_timestamps(allele)
            self.synthesize_symbol(allele)
            self.synthesize_fullname(allele)
            self.synthesize_synonyms(allele)
            self.synthesize_secondary_ids(allele)
            self.synthesize_xrefs(allele)
            self.synthesize_references(allele)
            self.synthesize_insertions(allele)
            self.flag_internal_alleles(allele)
            self.flag_unexportable_alleles(allele)
            self.synthesize_extinction(allele)
            self.synthesize_inheritance_mode(allele)
        log.info('Done synthesizing allele info.')
        return

    def generate_export_file(self):
        """Process alleles and print to a LinkML-compliant JSON file."""
        log.info('Generating output JSON file of alleles.')
        output_dict = {
            'linkml_version': linkml_release,
            'allele_ingest_set': []
        }
        for allele in self.allele_dict.values():
            if allele.for_alliance_export is False:
                log.debug('Suppress allele from export: {}. Reasons: {}'.format(allele, '; '.join(allele.export_warnings)))
                continue
            self.export_feat_cnt += 1
            if allele.internal is True:
                self.internal_feat_cnt += 1
                log.debug('Mark allele as internal: {}. Reasons: {}'.format(allele, '; '.join(allele.internal_reasons)))
            output_allele = {}
            for attr in self.output_fields:
                if getattr(allele, attr) is not None and getattr(allele, attr) != []:
                    output_allele[attr] = getattr(allele, attr)
            output_dict['allele_ingest_set'].append(output_allele)
        log.info('Writing data to output file.')
        with open(output_filename, 'w') as outfile:
            json.dump(output_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
            outfile.close()
        log.info('Done writing data to output file.')
        total_public_feat_cnt = self.export_feat_cnt - self.internal_feat_cnt
        log.info('Exported {} of {} alleles ({} are public).'.
                 format(self.export_feat_cnt, self.total_feat_cnt, total_public_feat_cnt))
        log.info('Suppressed {} alleles from export.'.format(self.total_feat_cnt - self.export_feat_cnt))
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
