# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Data retrieval of FlyBase disease annotations for Alliance curation database.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    AGR_data_retrieval_doid_curation.py [-h] [-l LINKML_RELEASE] [-v VERBOSE] [-c CONFIG]

Example:
    python AGR_data_retrieval_doid_curation.py -v -l v1.1.2 -c /path/to/config.cfg

Notes:
    This script exports disease annotations to a JSON file conforming to LinkML
    specs for the curation (i.e., "persistent") database; distinct from DAF
    file specs for the original Neo4j drop-and-reload database.

"""

import argparse
# import datetime
import json
import re
# import strict_rfc3339
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import aliased, sessionmaker
# from sqlalchemy.orm.exc import NoResultFound
from harvdev_utils.production import (
    Cv, Cvterm, Db, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureDbxref, FeatureRelationship, Pub, PubDbxref
)
from harvdev_utils.psycopg_functions import set_up_db_reading

# Now proceed with generic setup.
report_label = 'disease_curation'
set_up_dict = set_up_db_reading(report_label)
server = set_up_dict['server']
database = set_up_dict['database']
username = set_up_dict['username']
password = set_up_dict['password']
database_release = set_up_dict['database_release']
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
Session = sessionmaker(bind=engine)
session = Session()


# The main process.
def main():
    """Run the steps for exporting LinkML-compliant FlyBase disease annotations."""
    log.info('Running script "{}"'.format(__file__))
    log.info('Started main function.')
    log.info('Output JSON file corresponds to "agr_curation_schema" release: {}'.format(linkml_release))

    # Instantiate the object, get the data, synthesize it, export it.
    daf_maker = DAFMaker()
    daf_maker.query_chado(session)
    daf_maker.synthesize_info(session)
    daf_maker.generate_export_file()
    log.info('Ended main function.\n')


class DiseaseAnnotation(object):
    """A disease annotation."""
    def __init__(self, feature_cvterm, provenance_prop):
        """Create a disease annotation object from each distinct FeatureCvterm/FeatureCvtermprop combo.

        Args:
            arg1 (feature_cvterm): (FeatureCvterm) The FeatureCvterm object.
            arg2 (provenance_prop): (FeatureCvtermprop) The FeatureCvtermprop object of type "provenance".

        Returns:
            An object of the DiseaseAnnotation class.

        """
        # FlyBase data
        self.unique_key = '{}_{}'.format(feature_cvterm.feature_cvterm_id, provenance_prop.rank)
        self.feature_cvterm = feature_cvterm             # The FeatureCvterm object.
        self.provenance = provenance_prop                # The "provenance" FeatureCvtermprop.
        self.evidence_code = None                        # Will be the "evidence_code" FeatureCvtermprop.
        self.qualifier = None                            # Will be the "qualifier" FeatureCvtermprop.
        self.timestamps = []                             # Will be a list of audit_chado timestamp lists.
        # Derived attribures.
        self.modifier_problem = False                    # Change to true if there's a problem finding the modifier allele.
        # Attributes for the Alliance AuditedObject.
        self.obsolete = False                            # Never True. All FB annotations are deleted if no longer current.
        self.internal = False                            # Will be internal if annotation should not be exported to Alliance for some reason.
        self.created_by = 'FB:FB_curator'                # Use placeholder value since no Person object at FlyBase.
        self.updated_by = 'FB:FB_curator'               # Use placeholder value since no Person object at FlyBase.
        self.date_created = None                         # Not straightforward as half of relevant annotations are derived in the reporting build.
        self.date_updated = None                         # Not straightforward as half of relevant annotations are derived in the reporting build.
        # Attributes for the Alliance Association
        self.subject = None                              # Provide allele curie (slot usage from AlleleDiseaseAnnotation)
        self.predicate = 'is_implicated_in'              # "Allele disease relations" CV (slot usage from AlleleDiseaseAnnotation)
        self.object = None                               # Provide DOID (slot usage from DiseaseAnnotation).
        # Attributes for the Alliance DiseaseAnnotation
        self.data_provider = 'FB'
        self.negated = False                             # Change to True for "NOT" annotations.
        self.evidence_codes = []                         # Set as appropriate.
        self.single_reference = None                     # Provide FBrf ID.
        self.annotation_type = 'manually_curated'        # "Annotation types" CV.
        self.disease_genetic_modifier = None             # Gene, Allele or AGM curie.
        self.disease_genetic_modifier_relation = None    # "Disease genetic modifier relations" CV.
        self.unique_id = self.unique_key                 # Use the unique_key (internal ID is ok).
        self.mod_entity_id = None                        # N/A to FlyBase data.
        self.inferred_gene = None                        # Gene asserted by curator to be associated with the disease annotation.
        # self.with = None                               # N/A to FlyBase data.
        self.disease_qualifiers = []                     # N/A to FlyBase data. "Disease Qualifiers" CV.
        self.condition_relations = []                    # N/A to FlyBase data.
        self.genetic_sex = None                          # N/A to FlyBase data. "Genetic sexes" CV.
        self.related_notes = []                          # N/A to FlyBase data.
        self.secondary_data_provider = None              # N/A to FlyBase data.
        # Notes associated with the object.
        self.for_alliance_export = True         # Change to False if object should be excluded from export.
        self.internal_reasons = []              # Reasons for marking an object as internal. Will be exported but not necessarily displayed at Alliance.
        self.export_warnings = []               # Reasons for suppressing an object from the export file.

    def __str__(self):
        """Succinct text string describing the disease annotation."""
        desc = 'feature_cvterm_id={}\trank={}\tfeature={} ({})\tcvterm={} (DOID:{})\tpub={}\tqual={}\tevidence={}'.\
            format(self.feature_cvterm.feature_cvterm_id,
                   self.provenance.rank,
                   self.feature_cvterm.feature.name,
                   self.feature_cvterm.feature.uniquename,
                   self.feature_cvterm.cvterm.name,
                   self.feature_cvterm.cvterm.dbxref.accession,
                   self.feature_cvterm.pub.uniquename,
                   self.qualifier.value,
                   self.evidence_code.value)
        return desc


class DAFMaker(object):
    """This object gets disease annotations, synthesizes/filters the data, then exports it as LinkML JSON."""
    def __init__(self):
        """Create the DAFMaker object."""
        self.dis_anno_dict = {}       # A dict of DiseaseAnnotations keyed by feature_cvterm_id plus rank (e.g., 1234567_0).
        self.total_anno_cnt = 0       # Count of all disease annotations found in starting query.
        self.export_anno_cnt = 0      # Count of all disease annotations exported to file.
        self.internal_anno_cnt = 0    # Count of all disease annotations marked as internal=True in export file.

    relevant_qualifiers = [
        'model of',
        'DOES NOT model'
    ]

    disease_genetic_modifier_terms = {
        'is ameliorated by': 'ameliorated_by',
        'is NOT ameliorated by': 'not_ameliorated_by',
        'is exacerbated by': 'exacerbated_by',
        'is NOT exacerbated by': 'not_exacerbated_by'
    }

    evidence_code_xrefs = {
        'CEA': 'ECO:0007013',
        'CEC': 'ECO:0007014'
    }

    required_fields = [
        'data_provider',
        'evidence_codes',
        'internal',
        'object',
        'predicate',
        'single_reference'
        'subject'
    ]

    output_fields = [
        'annotation_type',
        'created_by',
        'data_provider',
        'date_created',
        'date_updated',
        'disease_genetic_modifier',
        'disease_genetic_modifier_relation',
        'evidence_codes',
        'inferred_gene',
        'internal',
        'updated_by',
        'negated',
        'object',
        'obsolete',
        'predicate',
        'single_reference',
        'subject',
        'unique_id'    # For derived annotations, feature_cvterm_id+rank changes each release. So, suppress.
    ]

    def get_disease_annotations(self, session):
        """Get all allele-level disease annotations."""
        log.info('Querying chado for feature-level disease annotations.')

        # First get feature_cvterm and related provenance props (each combo is a distinct annotation).
        # This will ignore orthology-computed disease annotations, which are at the gene level.
        allele_regex = r'^FBal[0-9]{7}$'
        pub_regex = r'^FBrf[0-9]{7}$'
        disease_term = aliased(Cvterm, name='disease_term')
        feature_type = aliased(Cvterm, name='feature_type')
        prop_type = aliased(Cvterm, name='prop_type')
        filters = (
            Feature.uniquename.op('~')(allele_regex),
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            feature_type.name == 'allele',
            disease_term.is_obsolete == 0,
            Cv.name == 'disease_ontology',
            Db.name == 'DOID',
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(pub_regex),
            prop_type.name == 'provenance'
        )
        dis_annos = session.query(FeatureCvterm, FeatureCvtermprop).\
            join(Feature, (Feature.feature_id == FeatureCvterm.feature_id)).\
            join(feature_type, (feature_type.cvterm_id == Feature.type_id)).\
            join(disease_term, (disease_term.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == disease_term.cv_id)).\
            join(Dbxref, (Dbxref.dbxref_id == disease_term.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Pub, (Pub.pub_id == FeatureCvterm.pub_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(prop_type, (prop_type.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        for result in dis_annos:
            self.total_anno_cnt += 1
            dis_anno = DiseaseAnnotation(result.FeatureCvterm, result.FeatureCvtermprop)
            self.dis_anno_dict[dis_anno.unique_key] = dis_anno
        log.info('Found {} disease annotations.'.format(self.total_anno_cnt))

        # Get qualifiers for each disease annotation.
        log.info('Getting disease annotation qualifiers.')
        filters = (
            Cvterm.name == 'qualifier',
        )
        fcvt_qualifiers = session.query(FeatureCvtermprop).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        qualifier_count = 0
        for qualifier in fcvt_qualifiers:
            unique_key = '{}_{}'.format(qualifier.feature_cvterm_id, qualifier.rank)
            try:
                self.dis_anno_dict[unique_key].qualifier = qualifier
                qualifier_count += 1
            except KeyError:
                pass
        log.info('Found {} disease annotation qualifiers.'.format(qualifier_count))

        # Get evidence_code for each disease annotation.
        log.info('Getting disease annotation evidence codes.')
        filters = (
            Cvterm.name == 'evidence_code',
        )
        fcvt_evidence_codes = session.query(FeatureCvtermprop).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        evidence_code_count = 0
        for evidence_code in fcvt_evidence_codes:
            unique_key = '{}_{}'.format(evidence_code.feature_cvterm_id, evidence_code.rank)
            try:
                self.dis_anno_dict[unique_key].evidence_code = evidence_code
                evidence_code_count += 1
            except KeyError:
                pass
        log.info('Found {} disease annotation evidence codes.'.format(evidence_code_count))

        # Print out annotations for review and development.
        for dis_anno in self.dis_anno_dict.values():
            log.debug('\nANNOTATION: {}'.format(dis_anno))
        return

    def get_dis_anno_timestamps(self, session):
        """Get timestamps for disease annotations."""
        log.info('Getting disease annotation timestamps.')
        # Note - I'm using standard SQL queries because querying t_audit_chado table with SQLAlchemy ORM is not user-friendly.
        audit_chado_query = """
            SELECT DISTINCT fcvt.feature_cvterm_id||'_'||fcvtp.rank,
                            ac.transaction_timestamp
            FROM feature_cvterm fcvt
            JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
            JOIN cv ON (cv.cv_id = cvt.cv_id AND cv.name = 'disease_ontology')
            JOIN feature_cvtermprop fcvtp ON fcvtp.feature_cvterm_id = fcvt.feature_cvterm_id
            JOIN audit_chado ac ON (ac.record_pkey = fcvtp.feature_cvtermprop_id AND ac.audited_table = 'feature_cvtermprop');
        """
        audit_results = session.execute(audit_chado_query).fetchall()
        log.info('Got {} audit_chado results. Will parse them out now.'.format(len(audit_results)))
        UNIQUE_KEY = 0
        TIMESTAMP = 1
        for row in audit_results:
            try:
                log.debug('For unique_key={}, have timestamp={}'.format(row[UNIQUE_KEY], row[TIMESTAMP]))
                self.dis_anno_dict[row[UNIQUE_KEY]].timestamps.append(row[TIMESTAMP])
            except KeyError:
                log.debug('Could not put this in anno dict: {}'.format(row))
        log.info('Timestamps have been retrieved.')
        return

    def query_chado(self, session):
        """A wrapper method that runs initial db queries."""
        self.get_disease_annotations(session)
        # self.get_dis_anno_timestamps(session)   # Suppress since half of annotations are derived during the release build.
        return

    def get_pub_xref(self, session, uniquename):
        """Get the preferred xref for a publication: PubMed or FlyBase.

        Args:
            arg1 (session): (Session) The session for the query.
            arg2 (uniquename): (str) The pub uniquename that needs an xref.

        Returns:
            A string that is either the PubMed ID (e.g., PMID:12345) or the FlyBase ID (e.g., FB:FBrf1234567).

        """
        filters = (
            Pub.uniquename == uniquename,
            PubDbxref.is_current.is_(True),
            Db.name == 'pubmed'
        )
        pubmed_xref = session.query(Dbxref).\
            join(PubDbxref, (PubDbxref.dbxref_id == Dbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Pub, (Pub.pub_id == PubDbxref.pub_id)).\
            filter(*filters).\
            one_or_none()
        if pubmed_xref:
            pub_xref = 'PMID:{}'.format(pubmed_xref.accession)
        else:
            pub_xref = 'FB:{}'.format(uniquename)
        return pub_xref

    def get_inferred_gene(self, session, allele_feature_id):
        """Get gene for the annotation's subject allele."""
        gene_regex = r'^FBgn[0-9]{7}$'
        filters = (
            FeatureRelationship.subject_id == allele_feature_id,
            Cvterm.name == 'alleleof',
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(gene_regex)
        )
        parent_gene = session.query(Feature).\
            join(FeatureRelationship, (FeatureRelationship.object_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            one()
        gene_curie = None
        if parent_gene:
            gene_curie = f'FB:{parent_gene.uniquename}'
        return gene_curie

    def confirm_current_allele_by_uniquename(self, session, uniquename):
        """Confirm that a given uniquename corresponds to a current allele.

        Args:
            arg1 (session): (Session) The session for the query.
            arg2 (uniquename): (str) The allele uniquename to be checked.

        Returns:
            A boolean value: True if uniquename corresponds to a current allele; False otherwise.

        """
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.uniquename == uniquename,
            Feature.uniquename.op('~')(allele_regex),
            Feature.is_obsolete.is_(False)
        )
        results = session.query(Feature).filter(*filters).one_or_none()
        if results:
            return True
        else:
            return False

    def get_current_id_for_allele(self, session, old_uniquename):
        """Get the ID for the current feature corresponding to some obsolete feature.

        Args:
            arg1 (session): (Session) The session for the query.
            arg2 (old_uniquename): (str) The obsolete allele uniquename to be checked.

        Returns:
            None, or the uniquename for the current feature that corresponds to the obsolete feature.

        """
        allele_regex = r'^FBal[0-9]{7}$'
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(allele_regex),
            Feature.uniquename != old_uniquename,
            FeatureDbxref.is_current.is_(False),
            Dbxref.accession == old_uniquename,
            Db.name == 'FlyBase'
        )
        curr_alleles = session.query(Feature).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        curr_uniquenames = [i.uniquename for i in curr_alleles]
        if len(curr_uniquenames) == 1:
            log.warning('For obsolete {}, found one current allele: {}'.format(old_uniquename, curr_uniquenames[0]))
            curr_allele_id = curr_uniquenames[0]
        elif len(curr_uniquenames) > 1:
            log.warning('For obsolete {}, found many current alleles: {}'.format(old_uniquename, curr_uniquenames))
            curr_allele_id = None
        else:
            log.warning('For obsolete {}, found no current alleles.'.format(old_uniquename))
            curr_allele_id = None
        return curr_allele_id

    def synthesize_info(self, session):
        """Synthesize disease annotation info and determine exportability."""
        log.info('Synthesizing disease annotation info.')
        for dis_anno in self.dis_anno_dict.values():
            log.debug('Evaluating annotation: {}'.format(dis_anno))
            # Get subject, object and pub.
            dis_anno.subject = 'FB:{}'.format(dis_anno.feature_cvterm.feature.uniquename)
            dis_anno.object = 'DOID:{}'.format(dis_anno.feature_cvterm.cvterm.dbxref.accession)
            dis_anno.single_reference = self.get_pub_xref(session, dis_anno.feature_cvterm.pub.uniquename)
            dis_anno.inferred_gene = self.get_inferred_gene(session, dis_anno.feature_cvterm.feature.feature_id)
            # Mark negative annotations.
            if dis_anno.qualifier.value == 'DOES NOT model':
                dis_anno.negated = True
            # Get timestamps. Suppressed for now.
            # if dis_anno.timestamps:
            #     dis_anno.creation_date = strict_rfc3339.\
            #         timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(min(dis_anno.timestamps)))
            #     dis_anno.date_last_modified = strict_rfc3339.\
            #         timestamp_to_rfc3339_localoffset(datetime.datetime.timestamp(max(dis_anno.timestamps)))
            # Determine evidence_code
            if dis_anno.evidence_code.value.startswith('CEC'):
                dis_anno.evidence_codes.append(self.evidence_code_xrefs['CEC'])
            else:
                dis_anno.evidence_codes.append(self.evidence_code_xrefs['CEA'])
            # Find modifiers and their relations.
            allele_regex = r'FBal[0-9]{7}'
            for fb_term in self.disease_genetic_modifier_terms.keys():
                if fb_term in dis_anno.evidence_code.value:
                    dis_anno.disease_genetic_modifier_relation = self.disease_genetic_modifier_terms[fb_term]
                    if re.search(allele_regex, dis_anno.evidence_code.value):
                        allele_id = re.search(allele_regex, dis_anno.evidence_code.value).group(0)
                        if self.confirm_current_allele_by_uniquename(session, allele_id):
                            dis_anno.disease_genetic_modifier = 'FB:{}'.format(allele_id)
                        else:
                            # Look up current allele by 2o ID. Use that.
                            curr_allele_id = self.get_current_id_for_allele(session, allele_id)
                            if curr_allele_id:
                                dis_anno.disease_genetic_modifier = 'FB:{}'.format(curr_allele_id)
                            else:
                                dis_anno.modifier_problem = True
            # Now check for conditions that prevent export.
            self.evaluate_annot(dis_anno)
        log.info('Done synthesizing disease annotation info.')
        return

    def evaluate_annot(self, dis_anno):
        """Identify conditions that prevent export of disease annotation."""
        # To Do: Need a better way of identifying internal FlyBase alleles (in chado itself, a featureprop).
        # For now, all disease annotations marked as "internal = True".
        export_checks = {
            dis_anno.qualifier.value not in self.relevant_qualifiers:
                'Only "model of|DOES NOT model" is exportable',
            ' with FLYBASE' in dis_anno.evidence_code.value:
                'Only disease annotations modeled by a single allele are exportable',
            dis_anno.modifier_problem is True: 'Cannot find current feature for disease modifier.'
        }
        for check, msg in export_checks.items():
            if check:
                dis_anno.for_alliance_export = False
                dis_anno.export_warnings.append(msg)
                log.debug(msg)
        return

    def generate_export_file(self):
        """Process disease annotations and print to a LinkML-compliant JSON file."""
        log.info('Generating output JSON file of disease annotations.')
        output_dict = {
            'linkml_version': linkml_release,
            'disease_allele_ingest_set': []
        }
        for dis_anno in self.dis_anno_dict.values():
            if dis_anno.for_alliance_export is False:
                log.debug('Suppress disease annotation from export: {}. Reasons: {}'.format(dis_anno, '; '.join(dis_anno.export_warnings)))
                continue
            self.export_anno_cnt += 1
            if dis_anno.internal is True:
                self.internal_anno_cnt += 1
                log.debug('Mark disease annotation as internal: {}. Reasons: {}'.format(dis_anno, '; '.join(dis_anno.internal_reasons)))
            output_dis_anno = {}
            for attr in self.output_fields:
                if getattr(dis_anno, attr) is not None and getattr(dis_anno, attr) != []:
                    output_dis_anno[attr] = getattr(dis_anno, attr)
            output_dict['disease_allele_ingest_set'].append(output_dis_anno)
        log.info('Writing data to output file.')
        with open(output_filename, 'w') as outfile:
            json.dump(output_dict, outfile, sort_keys=True, indent=2, separators=(',', ': '))
            outfile.close()
        log.info('Done writing data to output file.')
        total_public_anno_cnt = self.export_anno_cnt - self.internal_anno_cnt
        log.info('Exported {} of {} disease annotations ({} are public).'.
                 format(self.export_anno_cnt, self.total_anno_cnt, total_public_anno_cnt))
        log.info('Suppressed {} disease annotations from export.'.format(self.total_anno_cnt - self.export_anno_cnt))
        return


if __name__ == "__main__":
    main()
