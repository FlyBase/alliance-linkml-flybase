# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase AGMs for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_agm_curation.py [-h] [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_agm_curation.py -v -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script makes a JSON file conforming to LinkML specs for the curation
    (i.e., "persistent") database; distinct from AGM file specs for the
    original Neo4j drop-and-reload database. As the "audit_chado" table may be
    required for proper determination of "date_updated" values, use a reporting
    database with full audit_chado table present.

"""

import argparse
import datetime
import json
import strict_rfc3339
from tqdm import tqdm
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, OrganismDbxref, Pub, PubDbxref, Strain, StrainDbxref,
    StrainPub, StrainSynonym, Synonym
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Now proceed with generic setup.
report_label = 'agm_curation'
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
    """Run the steps for exporting LinkML-compliant FlyBase AGMs."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')
    log.info('Output JSON file corresponds to "agr_curation_schema" release: {}'.format(linkml_release))

    # Instantiate the object, get the data, synthesize it, export it.
    agm_handler = AGMHandler()
    db_query_transaction(agm_handler)
    agm_handler.synthesize_info()
    agm_handler.generate_export_file()
    log.info('Ended main function.\n')


class AllianceStrainAGM(object):
    """A strain with it's associated FlyBase data and Alliance LinkML properties."""
    def __init__(self, strain):
        """Create a base AllianceStrainAGM object.

        Args:
            arg1 (strain): (Strain) The chado Strain object corresponding to the strain.

        Returns:
            An object of the AllianceStrainAGM class.

        """
        # Attributes representing unprocessed FlyBase data.
        # Note: use attribute names that do not match an Alliance LinkML slot name.
        # For initial load, the Alliance A-Team just needs minimum info:
        # AGM: curie, taxon, name, subtype.
        self.strain = strain                       # The Strain object corresponding to the FlyBase strain.
        self.organism_abbr = None                  # Will be the organism.abbreviation for the strain's species of origin.
        self.curr_fb_symbol = None                 # Will be the current symbol Synonym object.
        self.curr_fb_fullname = None               # Will be the current fullname Synonym object.
        self.internal_synonyms = []                # Will be list of internal synonym names (and synonym_sgml if different).
        self.public_synonyms = []                  # Will be list of public synonym names (and synonym_sgml if different).
        self.dbxrefs = []                          # Will be list of dbxrefs as sql result groupings: Db, Dbxref, StrainDbxref.
        self.alt_fb_ids = []                       # Will be list of Dbxrefs for 2o FlyBase IDs.
        self.timestamps = []                       # Add all timestamps here.
        self.fb_references = []                    # Will be list of pub_ids from strain_pub.
        # Attributes for the Alliance AuditedObject.
        self.obsolete = strain.is_obsolete         # Will be the FlyBase value here.
        self.internal = False                      # Change to true if strain not intended for display at Alliance website.
        self.created_by_curie = 'FB:FB_curator'    # Use placeholder value since no Person object at FlyBase.
        self.updated_by_curie = 'FB:FB_curator'    # Use placeholder value since no Person object at FlyBase.
        self.date_created = None                   # Earliest timestamp.
        self.date_updated = None                   # Latest timestamp.
        # Attributes for the Alliance BiologicalEntityDTO. BiologicalEntityDTO is_a AuditedObjectDTO.
        self.curie = 'FB:{}'.format(strain.uniquename)
        self.taxon_curie = None                    # A string representing the NCBI taxon ID.
        self.data_provider_dto = None              # Will be DataProviderDTO object.
        # Attributes for the Alliance GenomicEntityDTO. GenomicEntity is_a BiologicalEntityDTO.
        self.cross_reference_dtos = []             # Report only select dbs, using AGR-accepted db_prefix.
        self.genomic_location_dtos = []            # Not applicable to strains or genotypes.
        # Attributes for the Alliance AffectedGenomicModelDTO. AffectedGenomicModelDTO is_a GenomicEntityDTO.
        self.name = None                           # Will be current fullname synonym - report ascii or utf8 (sgml) version?
        self.subtype_name = 'strain'               # Here we only get strains, for now.
        self.agm_secondary_id_dtos = []            # Annotation IDs and 2o FlyBase IDs.
        self.reference_curies = []
        # Notes associated with the object.
        self.for_alliance_export = True         # Change to False if object should be excluded from export.
        self.internal_reasons = []              # Reasons for marking an object as internal in the export file.
        self.export_warnings = []               # Reasons for suppressing an object from the export file.

    def __str__(self):
        """Succinct text string describing the AllianceStrainAGM object."""
        desc = '{} ({})'.format(self.strain.name, self.strain.uniquename)
        return desc


class AGMHandler(object):
    """This object gets strains, synthesizes/filters the data, then exports it as LinkML JSON."""
    def __init__(self):
        """Create the AGMHandler object."""
        self.strain_dict = {}        # An FBsnID-keyed dict of AllianceStrainAGM objects.
        self.all_pubs_dict = {}      # A pub_id-keyed dict of pub curies (PMID or FBrf).
        self.total_agm_cnt = 0       # Count of all strains found in starting query.
        self.export_agm_cnt = 0      # Count of all strains exported to file.
        self.internal_agm_cnt = 0    # Count of all strains marked as internal=True in export file.

    # Generic objects with which to build Alliance DTOs.
    generic_audited_object = {
        'internal': False,
        'obsolete': False,
        'created_by_curie': 'FB:FB_curator',
        'updated_by_curie': 'FB:FB_curator'
    }
    generic_data_provider_dto = generic_audited_object.copy()
    generic_data_provider_dto['source_organization_abbreviation'] = 'FB'
    generic_cross_reference_dto = {'prefix': 'FB', 'internal': False}
    # Regexes.
    strain_regex = r'^FBsn[0-9]{7}'
    pub_regex = r'^(FBrf[0-9]{7}|unattributed)$'
    # Export fields.
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
    fb_agr_db_dict = {
        'FlyBase': 'FB'
    }

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

    def get_strains(self, session):
        """Get all strains."""
        log.info('Querying chado for strains.')
        # First get all strains
        filters = (
            Strain.uniquename.op('~')(self.strain_regex),
        )
        strain_results = session.query(Strain).\
            filter(*filters).\
            distinct()
        self.total_agm_cnt = 0
        for result in strain_results:
            self.total_agm_cnt += 1
            self.strain_dict[result.uniquename] = AllianceStrainAGM(result)
            self.strain_dict[result.uniquename].organism_abbr = result.organism.abbreviation
        log.info('Found {} strains.'.format(self.total_agm_cnt))

    def get_strain_taxons(self, session):
        """Get taxon IDs for strains."""
        log.info('Getting strain taxon IDs.')
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
        for strain in self.strain_dict.values():
            try:
                strain.taxon_curie = f'NCBITaxon:{organism_taxon_dict[strain.strain.organism_id]}'
            except KeyError:
                log.debug('No NCBI taxon ID available for: {}'.format(strain))
        return

    def get_synonyms(self, session):
        """Get current and non-current symbols and full names for strains."""
        log.info('Getting strain synonyms.')
        strain_regex = r'^FBgn[0-9]{7}$'
        filters = (
            Strain.uniquename.op('~')(strain_regex),
        )
        strain_curr_symbol_results = session.query(Cvterm, Strain, StrainSynonym, Synonym).\
            join(StrainSynonym, (StrainSynonym.synonym_id == Synonym.synonym_id)).\
            join(Strain, (Strain.strain_id == StrainSynonym.strain_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Synonym.type_id)).\
            filter(*filters).\
            distinct()
        for result in strain_curr_symbol_results:
            if result.StrainSynonym.is_current is True:
                if result.Cvterm.name == 'symbol':
                    self.strain_dict[result.Strain.uniquename].curr_fb_symbol = result.Synonym
                elif result.Cvterm.name == 'fullname':
                    self.strain_dict[result.Strain.uniquename].curr_fb_fullname = result.Synonym
            elif result.StrainSynonym.is_internal is True:
                self.strain_dict[result.Strain.uniquename].internal_synonyms.append(result.Synonym.name)
                self.strain_dict[result.Strain.uniquename].internal_synonyms.append(result.Synonym.synonym_sgml)
            else:
                self.strain_dict[result.Strain.uniquename].public_synonyms.append(result.Synonym.name)
                self.strain_dict[result.Strain.uniquename].public_synonyms.append(result.Synonym.synonym_sgml)
        return

    def get_strain_timestamps(self, session):
        """Get timestamps for strains."""
        log.info('Getting strain timestamps.')
        # Querying audit_chado with SQLAlchemy has not been user-friendly, so go with standard SQL.
        strain_related_tables = [
            'cell_line_strain',
            'library_strain',
            'strain_cvterm',
            'strain_dbxref',
            'strain_feature',
            'strain_phenotype',
            'strain_pub',
            'strain_synonym',
            'strainprop'
        ]
        primary_audit_chado_query = """
            SELECT DISTINCT s.uniquename,
                            ac.transaction_timestamp
            FROM strain s
            JOIN audit_chado ac ON (ac.record_pkey = s.strain_id AND ac.audited_table = 'strain')
            ;
        """
        related_audit_chado_query_template = """
            SELECT DISTINCT s.uniquename,
                            ac.transaction_timestamp
            FROM strain s
            JOIN {related_table} ON {related_table}.strain_id = s.strain_id
            JOIN audit_chado ac ON (ac.record_pkey = {related_table}.{related_table}_id AND ac.audited_table = '{related_table}')
        ;
        """
        UNIQUE_KEY = 0
        TIMESTAMP = 1
        primary_audit_results = session.execute(primary_audit_chado_query).fetchall()
        for row in primary_audit_results:
            try:
                self.strain_dict[row[UNIQUE_KEY]].timestamps.append(row[TIMESTAMP])
            except KeyError:
                pass
        for table in strain_related_tables:
            related_audit_chado_query = related_audit_chado_query_template.format(related_table=table)
            related_audit_results = session.execute(related_audit_chado_query).fetchall()
            for row in related_audit_results:
                try:
                    self.strain_dict[row[UNIQUE_KEY]].timestamps.append(row[TIMESTAMP])
                except KeyError:
                    pass
        return

    def get_strain_dbxrefs(self, session):
        """Get all dbxrefs for strains."""
        log.info('Getting strain dbxrefs.')
        filters = (
            Db.name.in_((self.fb_agr_db_dict.keys())),
        )
        strain_dbxref_results = session.query(Strain, StrainDbxref, Dbxref, Db).\
            join(StrainDbxref, (StrainDbxref.strain_id == Strain.strain_id)).\
            join(Dbxref, (Dbxref.dbxref_id == StrainDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        tqdm_desc = 'Processing strain dbxrefs'
        for result in tqdm(strain_dbxref_results, desc=tqdm_desc):
            if result.StrainDbxref.is_current is False and result.Db.name == 'FlyBase':
                self.strain_dict[result.Strain.uniquename].alt_fb_ids.append(result.Dbxref)
            else:
                self.strain_dict[result.Strain.uniquename].dbxrefs.append(result)
        return

    def get_strain_references(self, session):
        """Get references for strains."""
        log.info('Get strain references.')
        filters = (
            Strain.uniquename.op('~')(self.strain_regex),
            Pub.uniquename.op('~')(self.pub_regex),
            Pub.is_obsolete.is_(False)
        )
        strain_pubs = session.query(Strain, Pub).\
            select_from(Strain).\
            join(StrainPub, (StrainPub.strain_id == Strain.strain_id)).\
            join(Pub, (Pub.pub_id == StrainPub.pub_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in strain_pubs:
            self.strain_dict[result.Strain.uniquename].fb_references.append(result.Pub.pub_id)
            counter += 1
        log.info(f'Found {counter} strain-pub relationships.')
        return

    # Synthesis of initial db info.
    def synthesize_timestamps(self, agm):
        """Process timestamps for an AGM."""
        if agm.timestamps:
            agm.date_created_curie = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(agm.timestamps)))
            agm.date_updated_curie = strict_rfc3339.\
                timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(agm.timestamps)))
        return

    def synthesize_synonyms(self, agm):
        """Determine name for AGM."""
        if agm.curr_fb_symbol:
            agm.name = agm.curr_fb_symbol.synonym_sgml
        else:
            agm.name = agm.strain.name
        return

    def synthesize_secondary_ids(self, agm):
        """Process 2o IDs."""
        unique_fb_id_list = list(set(agm.alt_fb_ids))
        for fb_id in unique_fb_id_list:
            secondary_id_dict = self.generic_audited_object.copy()
            secondary_id_dict['secondary_id'] = f'FB:{fb_id.accession}'
            agm.agm_secondary_id_dtos.append(secondary_id_dict)
        return

    def synthesize_xrefs(self, agm):
        """Process xrefs."""
        # Start by adding AGM uniquename as an xref.
        xref_dict = self.generic_audited_object.copy()
        xref_dict['referenced_curie'] = agm.curie
        xref_dict['display_name'] = agm.name
        xref_dict['prefix'] = 'FB'
        xref_dict['page_area'] = agm.subtype_name
        agm.cross_reference_dtos.append(xref_dict)
        # Add other xrefs: code below assumes xrefs are all 'FB' at the moment.
        for result in agm.dbxrefs:
            # Skip irrelevant xrefs.
            if result.Db.name not in self.fb_agr_db_dict.keys():
                continue
            # Skip 1o FB xref (already handled above).
            if f'FB:{result.Dbxref.accession}' == agm.curie:
                continue
            xref_dict = self.generic_audited_object.copy()
            xref_dict['referenced_curie'] = f'{self.fb_agr_db_dict[result.Db.name]}:{result.Dbxref.accession}'
            xref_dict['display_name'] = f'{self.fb_agr_db_dict[result.Db.name]}:{result.Dbxref.accession}'
            xref_dict['prefix'] = self.fb_agr_db_dict[result.Db.name]
            xref_dict['page_area'] = agm.subtype_name
            if result.StrainDbxref.is_current is False:
                xref_dict['internal'] = True
            agm.cross_reference_dtos.append(xref_dict)
        return

    def synthesize_references(self, agm):
        """Process pubs for agm."""
        agm.fb_references = list(set(agm.fb_references))
        agm.reference_curies = [self.all_pubs_dict[i] for i in agm.fb_references if self.all_pubs_dict[i] != 'FB:unattributed']
        return

    def add_data_provider_info(self, agm):
        """Add data_provider info."""
        agm.data_provider_dto = self.generic_data_provider_dto.copy()
        agm.data_provider_dto['cross_reference_dto'] = self.generic_cross_reference_dto.copy()
        agm.data_provider_dto['cross_reference_dto']['page_area'] = agm.subtype_name
        agm.data_provider_dto['cross_reference_dto']['referenced_curie'] = f'FB:{agm.curie}'
        agm.data_provider_dto['cross_reference_dto']['display_name'] = agm.name
        return

    def flag_internal_agms(self, agm):
        """Flag agms as internal."""
        if agm.obsolete is True:
            agm.internal = True
            agm.internal_reasons.append('Obsolete')
        if agm.organism_abbr != 'Dmel':
            agm.internal = True
            agm.internal_reasons.append('Non-Dmel')
        return

    def flag_unexportable_agms(self, agm):
        """Flag agms missing data required for export."""
        # TEMPORARY: Suppress non-Dmel AGMs from export.
        if agm.organism_abbr != 'Dmel':
            agm.for_alliance_export = False
            agm.export_warnings.append(f'Suppress non-Dmel AGM from export: ORG={agm.organism_abbr}')
        # Suppress objects missing required information from export.
        for attr in self.required_fields:
            if attr not in agm.__dict__.keys():
                agm.for_alliance_export = False
                agm.export_warnings.append('Missing "{}" attribute'.format(attr))
            elif getattr(agm, attr) is None:
                agm.for_alliance_export = False
                agm.export_warnings.append('Missing value for "{}" attribute'.format(attr))
        if agm.for_alliance_export is True:
            log.debug('EXPORT {}'.format(agm.curie))
        return

    def synthesize_info(self):
        """Convert FlyBase strain data into an AllianceAGM representation."""
        log.info('Synthesizing AGM info.')
        for strain in self.strain_dict.values():
            log.debug('Evaluating annotation: {}'.format(strain))
            self.synthesize_timestamps(strain)
            self.synthesize_synonyms(strain)
            self.synthesize_secondary_ids(strain)
            self.synthesize_xrefs(strain)
            self.synthesize_references(strain)
            self.add_data_provider_info(strain)
            self.flag_internal_agms(strain)
            self.flag_unexportable_agms(strain)
        log.info('Done synthesizing strain info.')
        return

    def query_chado(self, session):
        """A wrapper method that runs initial db queries."""
        self.get_all_references(session)
        self.get_strains(session)
        self.get_strain_taxons(session)
        self.get_synonyms(session)
        self.get_strain_timestamps(session)
        self.get_strain_dbxrefs(session)
        self.get_strain_references(session)
        return

    def generate_export_file(self):
        """Process AGMs and print to a LinkML-compliant JSON file."""
        log.info('Generating output JSON file of AGMs.')
        output_dict = {
            'linkml_version': linkml_release,
            'agm_ingest_set': []
        }
        for strain in self.strain_dict.values():
            if strain.for_alliance_export is False:
                log.debug('Suppress strain from export: {}. Reasons: {}'.format(strain, '; '.join(strain.export_warnings)))
                continue
            self.export_agm_cnt += 1
            if strain.internal is True:
                self.internal_agm_cnt += 1
                log.debug('Mark strain as internal: {}. Reasons: {}'.format(strain, '; '.join(strain.internal_reasons)))
            output_strain = {}
            for attr in self.output_fields:
                if getattr(strain, attr) is not None and getattr(strain, attr) != []:
                    output_strain[attr] = getattr(strain, attr)
            output_dict['agm_ingest_set'].append(output_strain)
        log.info('Writing data to output file.')
        with open(output_filename, 'w') as outfile:
            json.dump(output_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
            outfile.close()
        log.info('Done writing data to output file.')
        total_public_agm_cnt = self.export_agm_cnt - self.internal_agm_cnt
        log.info('Exported {} of {} strains ({} are public).'.
                 format(self.export_agm_cnt, self.total_agm_cnt, total_public_agm_cnt))
        log.info('Suppressed {} strains from export.'.format(self.total_agm_cnt - self.export_agm_cnt))
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
