"""Module:: disease_handlers.

Synopsis:
    Data handlers for FlyBase disease annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

# NOTES:
# 1. There are hdm_comment featureprops for alleles that might be relevant.
#    a. No explicit connection to annotations, or even to specific disease terms.
#    b. But it might be possible to find relevant annotations by matching feature and pub.
#    c. Consider adding hdm_comment as notes to disease annotations if 1:1 match.
#    d. Consider giving curators a report of non-transferred comments for manual addition.

import re
from collections import defaultdict
from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.char_conversions import sgml_to_plain_text
from harvdev_utils.reporting import (
    Cv, Cvterm, Db, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureDbxref, FeatureRelationship, Pub
)
import agr_datatypes
import fb_datatypes
from handler import DataHandler


class AGMDiseaseHandler(DataHandler):
    """A data handler that converts allele-based disease annotations into genotype (AGM) annotations."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AGMDiseaseHandler object."""
        super().__init__(log, testing)
        self.datatype = 'agm_disease'
        self.fb_export_type = fb_datatypes.FBAlleleDiseaseAnnotation
        self.agr_export_type = agr_datatypes.AGMDiseaseAnnotationDTO
        self.primary_export_set = 'disease_agm_ingest_set'
        self.allele_name_lookup = {}                 # feature.name-keyed feature dicts, current only.
        self.doid_term_lookup = {}                   # cvterm.name-keyed cvterm dicts.
        self.model_eco_lookup = defaultdict(list)    # Evidence abbreviation lookup for "model_of" annotations.
        self.uniq_dis_dict = defaultdict(list)       # Disease annotations keyed by a unique attribute concatenation.
        self.gal4_dict = defaultdict(list)           # Gal4 info to integrate, keyed by unique disease descriptor.

    relevant_fcvtp_types = [
        'evidence_code',
        'qualifier'
    ]

    permitted_qualifier_evidence_combos = {
        'ameliorates': r'^modeled by ',
        'DOES NOT ameliorate': r'^modeled by ',
        'exacerbates': r'^modeled by ',
        'DOES NOT exacerbate': r'^modeled by ',
        'model of': r'^CE(A|C)($| with )',
        'DOES NOT model': r'^CE(A|C)($| with )',
    }

    evidence_code_xrefs = {
        'CEA': 'ECO:0007013',
        'CEC': 'ECO:0007014'
    }

    disease_genetic_modifier_terms = {
        'ameliorates': 'ameliorated_by',
        'DOES NOT ameliorate': 'not_ameliorated_by',
        'exacerbates': 'exacerbated_by',
        'DOES NOT exacerbate': 'not_exacerbated_by'
    }

    # Add methods to be run by get_general_data() below.
    def build_allele_name_lookup(self):
        """Build name-keyed dict of alleles."""
        self.log.info('Build name-keyed dict of alleles.')
        for feature in self.feature_lookup.values():
            if feature['type'] == 'allele' and feature['is_obsolete'] is False:
                self.allele_name_lookup[feature['name']] = feature
        self.log.info(f'Have {len(self.allele_name_lookup)} current alleles in the allele-by-name lookup.')
        return

    def build_doid_term_lookup(self):
        """Build name-keyed dict of DOID terms."""
        self.log.info('Build name-keyed dict of DOID terms.')
        for term in self.cvterm_lookup.values():
            if term['cv_name'] == 'disease_ontology' and term['db_name'] == 'DOID':
                self.doid_term_lookup[term['name']] = term
        self.log.info(f'Have {len(self.cvterm_lookup)} DOID terms in the CV term-by-name lookup.')
        return

    # Elaborate on get_general_data() for the AGMDiseaseHandler.
    def get_general_data(self, session):
        """Extend the method for the AGMDiseaseHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'gene', 'insertion', 'construct'])
        self.build_uname_feature_lookup()
        self.get_transgenic_allele_ids(session)
        self.get_in_vitro_allele_ids(session)
        self.build_allele_gene_lookup(session)
        self.build_allele_name_lookup()
        self.build_doid_term_lookup()
        return

    # Add methods to be run by get_datatype_data() below.
    def get_allele_disease_annotations(self, session):
        """Get allele-based disease annotations from chado."""
        self.log.info('Get allele-based disease annotations from chado.')
        disease_term = aliased(Cvterm, name='disease_term')
        prop_type = aliased(Cvterm, name='prop_type')
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            Cv.name == 'disease_ontology',
            disease_term.is_obsolete == 0,
            Db.name == 'DOID',
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(self.regex['pub']),
            prop_type.name == 'provenance',
        )
        results = session.query(FeatureCvterm, FeatureCvtermprop).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(disease_term, (disease_term.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == disease_term.cv_id)).\
            join(Dbxref, (Dbxref.dbxref_id == disease_term.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Pub, (Pub.pub_id == FeatureCvterm.pub_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(prop_type, (prop_type.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            dis_anno = self.fb_export_type(result.FeatureCvterm, result.FeatureCvtermprop)
            self.fb_data_entities[dis_anno.db_primary_id] = dis_anno
            counter += 1
        self.log.info(f'Found {counter} allele-based disease annotations from chado (excludes annotations to obsolete alleles).')
        return

    def get_disease_qualifiers(self, session):
        """Get disease annotation qualifiers."""
        self.log.info('Get disease annotation qualifiers.')
        filters = (
            Cvterm.name == 'qualifier',
        )
        results = session.query(FeatureCvtermprop).\
            select_from(FeatureCvtermprop).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        qualifier_count = 0
        for qualifier in results:
            db_primary_id = '{}_{}'.format(qualifier.feature_cvterm_id, qualifier.rank)
            try:
                self.fb_data_entities[db_primary_id].qualifier = qualifier
                qualifier_count += 1
            except KeyError:
                pass
        self.log.info('Found {} disease annotation qualifiers.'.format(qualifier_count))
        return

    def get_disease_evidence_codes(self, session):
        """Get disease annotation evidence codes."""
        self.log.info('Get disease annotation evidence codes.')
        filters = (
            Cvterm.name == 'evidence_code',
        )
        results = session.query(FeatureCvtermprop).\
            select_from(FeatureCvtermprop).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        evidence_code_count = 0
        for evidence_code in results:
            db_primary_id = '{}_{}'.format(evidence_code.feature_cvterm_id, evidence_code.rank)
            try:
                self.fb_data_entities[db_primary_id].evidence_code = evidence_code
                evidence_code_count += 1
            except KeyError:
                pass
        self.log.info('Found {} disease annotation evidence codes.'.format(evidence_code_count))
        # Update annotation descriptions now that all info is in.
        for dis_anno in self.fb_data_entities.values():
            dis_anno.set_entity_desc()
        return

    def get_disease_timestamps(self, session):
        """Get timestamps for disease annotations."""
        self.log.info('Get timestamps for disease annotations.')
        # Note - using standard SQL query because SQLAlchemy ORM for audit_chado is not user-friendly.
        audit_chado_query = """
            SELECT DISTINCT fcvt.feature_cvterm_id||'_'||fcvtp.rank, ac.transaction_timestamp
            FROM feature_cvterm fcvt
            JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
            JOIN cv ON (cv.cv_id = cvt.cv_id AND cv.name = 'disease_ontology')
            JOIN feature_cvtermprop fcvtp ON fcvtp.feature_cvterm_id = fcvt.feature_cvterm_id
            JOIN audit_chado ac ON (ac.record_pkey = fcvtp.feature_cvtermprop_id AND ac.audited_table = 'feature_cvtermprop');
        """
        audit_results = session.execute(audit_chado_query).fetchall()
        self.log.info(f'Got {len(audit_results)} audit_chado results to parse.')
        DB_PRIMARY_ID = 0
        TIMESTAMP = 1
        counter = 0
        for row in audit_results:
            try:
                self.fb_data_entities[row[DB_PRIMARY_ID]].timestamps.append(row[TIMESTAMP])
                counter += 1
            except KeyError:
                # self.log.debug(f'Could not put this in anno dict: {row}')
                pass
        self.log.info(f'Obtained {counter} timestamps for annotations directly from the audit_chado table.')
        return

    def confirm_current_allele_by_uniquename(self, session, uniquename):
        """Confirm that a given uniquename corresponds to a current allele.

        Args:
            session (Session): The session for the query.
            uniquename (str): The allele uniquename to be checked.

        Returns:
            A boolean value: True if uniquename corresponds to a current allele; False otherwise.

        """
        filters = (
            Feature.uniquename == uniquename,
            Feature.uniquename.op('~')(self.regex['allele']),
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
            session (Session): The session for the query.
            old_uniquename (str): The obsolete allele uniquename to be checked.

        Returns:
            None, or the uniquename for the current feature that corresponds to the obsolete feature.

        """
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['allele']),
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
            self.log.debug('For obsolete {}, found one current allele: {}'.format(old_uniquename, curr_uniquenames[0]))
            curr_allele_id = curr_uniquenames[0]
        elif len(curr_uniquenames) > 1:
            self.log.debug('For obsolete {}, found many current alleles: {}'.format(old_uniquename, curr_uniquenames))
            curr_allele_id = None
        else:
            self.log.debug('For obsolete {}, found no current alleles.'.format(old_uniquename))
            curr_allele_id = None
        return curr_allele_id

    def extract_text_embedded_alleles(self, session):
        """Extract alleles from annotation text."""
        self.log.info('Extract alleles from annotation text.')
        embedded_allele_regex = r'FBal[0-9]{7}'
        current_allele_id_counter = 0
        updated_allele_id_counter = 0
        cannot_update_allele_id_counter = 0
        allele_id_prob_counter = 0
        for dis_anno in self.fb_data_entities.values():
            if not re.search(embedded_allele_regex, dis_anno.evidence_code.value):
                continue
            embedded_allele_ids = re.findall(embedded_allele_regex, dis_anno.evidence_code.value)
            for allele_id in embedded_allele_ids:
                if self.confirm_current_allele_by_uniquename(session, allele_id):
                    dis_anno.text_embedded_allele_ids.append(allele_id)
                    current_allele_id_counter += 1
                else:
                    curr_allele_id = self.get_current_id_for_allele(session, allele_id)
                    if curr_allele_id:
                        dis_anno.text_embedded_allele_ids.append(allele_id)
                        dis_anno.allele_id_was_updated = True
                        updated_allele_id_counter += 1
                    else:
                        dis_anno.allele_id_problem = True
                        cannot_update_allele_id_counter += 1
            if dis_anno.allele_id_problem is True:
                allele_id_prob_counter += 1
        self.log.info(f'{current_allele_id_counter} text-embedded allele IDs are current.')
        self.log.info(f'{updated_allele_id_counter} text-embedded allele IDs were not current and mapped unambiguously to a current allele.')
        self.log.info(f'{cannot_update_allele_id_counter} text-embedded allele IDs were not current and could NOT be mapped to a current allele.')
        self.log.info(f'{allele_id_prob_counter} disease annotations have one or more problematic text-embedded allele IDs.')
        return

    def extract_model_and_modifiers(self):
        """Extract model components from annotation subject, qualifier and evidence_code."""
        self.log.info('Extract model components from annotation subject, qualifier and evidence_code.')
        for dis_anno in self.fb_data_entities.values():
            dis_anno.modeled_by.extend(dis_anno.text_embedded_allele_ids)
            allele_subject_id = self.feature_lookup[dis_anno.feature_cvterm.feature_id]['uniquename']
            if dis_anno.qualifier.value in ('model of', 'DOES NOT model'):
                dis_anno.modeled_by.append(allele_subject_id)
                if dis_anno.qualifier.value == 'DOES NOT model':
                    dis_anno.is_not = True
            else:
                dis_anno.modifier_id = allele_subject_id
                dis_anno.modifier_role = self.disease_genetic_modifier_terms[dis_anno.qualifier.value]
        return

    def build_model_eco_lookup(self):
        """Build ECO lookup for model-type annotations."""
        self.log.info('Build ECO lookup for model-type annotations.')
        counter = 0
        for dis_anno in self.fb_data_entities.values():
            # Determine unique key for the model in this annotation.
            dis_anno.model_unique_key = f'{dis_anno.feature_cvterm.pub.uniquename}_'
            if dis_anno.is_not:
                dis_anno.unique_key += 'NOT_'
            dis_anno.model_unique_key = f'model={"|".join(sorted(dis_anno.modeled_by))}_'
            dis_anno.model_unique_key += f'disease_term=DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}'
            if re.match(r'^CE(A|C)', dis_anno.evidence_code.value):
                dis_anno.eco_abbr = dis_anno.evidence_code.value[0:3]
                if dis_anno.is_not is False:
                    self.model_eco_lookup[dis_anno.model_unique_key].append(dis_anno.eco_abbr)
                    counter += 1
        self.log.info(f'Have ECO abbreviations for {len(self.model_eco_lookup)} distinct disease models from {counter} annotations.')
        # Make sure the lookup is ok.
        zero_counter = 0
        one_counter = 0
        many_counter = 0
        distinct_ecos = []
        for ukey, eco_list in self.model_eco_lookup.items():
            uniqued_list = list(set(eco_list))
            self.model_eco_lookup[ukey] = uniqued_list
            distinct_ecos.extend(uniqued_list)
            if len(uniqued_list) == 0:
                self.log.debug(f'Zero ECOs for model_key={ukey}')
                zero_counter += 1
            elif len(uniqued_list) == 1:
                one_counter += 1
            else:
                self.log.debug(f'Many ECOs for model_key={ukey}: {uniqued_list}')
                many_counter += 1
        self.log.info(f'Have {zero_counter} keys in lookup with NO ECO; {one_counter} keys with ONE ECO; {many_counter} keys with MANY ECOS.')
        self.log.info(f'Have these distinct ECO codes: {set(distinct_ecos)}')
        return

    def lookup_eco_codes_for_modifier_annotations(self):
        """Lookup ECO code for modifier-type annotations."""
        self.log.info('Lookup ECO code for modifier-type annotations.')
        input_counter = 0
        skip_counter = 0
        assess_counter = 0
        match_counter = 0
        no_match_counter = 0
        for dis_anno in self.fb_data_entities.values():
            input_counter += 1
            if dis_anno.eco_abbr:
                skip_counter += 1
                continue
            assess_counter += 1
            # Choose default CEA if no info.
            if dis_anno.model_unique_key not in self.model_eco_lookup.keys():
                dis_anno.eco_abbr = 'CEA'
                no_match_counter += 1
            elif len(self.model_eco_lookup[dis_anno.model_unique_key]) == 1:
                dis_anno.eco_abbr = self.model_eco_lookup[dis_anno.model_unique_key][0]
                match_counter += 1
            # Choose default CEA if many codes.
            elif len(self.model_eco_lookup[dis_anno.model_unique_key]) > 1:
                dis_anno.eco_abbr = 'CEA'
                match_counter += 1
        self.log.info(f'Assessed {assess_counter}/{input_counter} annotations, skipped {skip_counter}/{input_counter} annotations.')
        self.log.info(f'Found ECO for {match_counter} modifier-type annotations.')
        self.log.info(f'Assigned "CEA" for {no_match_counter} modifier-type annotations with no info.')
        return

    def calculate_annotation_unique_key(self):
        """Calculate unique descriptors for disease annotations."""
        self.log.info('Calculate unique descriptors for disease annotations.')
        for dis_anno in self.fb_data_entities.values():
            dis_anno.unique_key = f'{dis_anno.feature_cvterm.pub.uniquename}_'
            if dis_anno.is_not:
                dis_anno.unique_key += 'NOT_'
            dis_anno.unique_key += f'model={"|".join(sorted(dis_anno.modeled_by))}_'
            dis_anno.unique_key += f'disease_term=DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}_'
            dis_anno.unique_key += f'eco_code={dis_anno.eco_abbr}'
            if not dis_anno.eco_abbr:
                self.log.warning(f'No ECO abbr for {dis_anno.unique_key}')
            if dis_anno.modifier_id:
                dis_anno.unique_key += f'_{dis_anno.modifier_role}={dis_anno.modifier_id}'
            # self.log.debug(f'Annotation db_primary_id={dis_anno.db_primary_id} has this unique key: {dis_anno.unique_key}')
            self.uniq_dis_dict[dis_anno.unique_key].append(dis_anno)
        return

    def integrate_gal4_info(self):
        """Integrate Gal4 driver info into annotations."""
        self.log.info('Integrate Gal4 driver info into annotations.')
        PUB_GIVEN = 0
        ALLELE_SYMBOL = 1
        QUAL = 4
        DO_TERM = 5
        EVI_CODE = 6
        ADDITIONAL_ALLELES = 7
        GAL4_INPUT = 8
        OPERATION = 12
        gal4_input = open('/src/output/gal4_driver_info.tsv')
        pub_not_found_counter = 0
        allele_not_found_counter = 0
        additional_allele_not_found_counter = 0
        do_term_not_found_counter = 0
        gal4_not_found_counter = 0
        dis_anno_not_found = 0
        line_number = 0
        matched_dis_anno_counter = 0
        unmatched_dis_anno_counter = 0
        for i in gal4_input:
            line_number += 1
            if not i.startswith('FBrf'):
                continue
            line = i.split('\t')
            gal4_info = {
                # Attributes from input file.
                'line_number': line_number,
                'pub_given': line[PUB_GIVEN],
                'allele_symbol': sgml_to_plain_text(line[ALLELE_SYMBOL]),
                'additional_alleles': line[ADDITIONAL_ALLELES].split(','),
                'qualifier': line[QUAL],
                'evi_code': line[EVI_CODE],
                'do_term': line[DO_TERM],
                'gal4_input': line[GAL4_INPUT].split(','),
                'operation': line[OPERATION].rstrip(),
                # Attributes to be obtained from chado.
                'pub': None,
                'allele_feature_id': None,
                'additional_allele_ids': [],
                'gal4_ids': [],
                'doid_term_curie': None,
                # Attributes synthesized from the above.
                'modeled_by': [],
                'is_not': False,
                'modifier_id': None,
                'modifier_role': None,
                'eco_abbr': None,
                'unique_key': None,
                'problem': False,
            }
            # Fill in info from chado.
            try:
                gal4_info['pub'] = self.fbrf_bibliography[gal4_info['pub_given']]
            except KeyError:
                self.log.error(f'Line={line_number}: could not find pub \"{gal4_info["pub_given"]}\" in chado.')
                gal4_info['problem'] = True
                pub_not_found_counter += 1
            try:
                gal4_info['allele_feature_id'] = self.allele_name_lookup[gal4_info['allele_symbol']]['uniquename']
            except KeyError:
                self.log.error(f'Line={line_number}: could not find allele \"{gal4_info["allele_symbol"]}\" in chado.')
                gal4_info['problem'] = True
                allele_not_found_counter += 1
            for allele_symbol in gal4_info['additional_alleles']:
                if allele_symbol == '':
                    continue
                converted_allele_symbol = sgml_to_plain_text(allele_symbol).strip()
                try:
                    gal4_info['additional_allele_ids'].append(self.allele_name_lookup[converted_allele_symbol]['uniquename'])
                except KeyError:
                    self.log.error(f'Line={line_number}: could not find additional allele "{allele_symbol}" in chado.')
                    gal4_info['problem'] = True
                    additional_allele_not_found_counter += 1
            try:
                gal4_info['doid_term_curie'] = self.doid_term_lookup[gal4_info['do_term']]['curie']
            except KeyError:
                self.log.error(f'Line={line_number}: could not find DO term \"{gal4_info["do_term"]}\" in chado.')
                gal4_info['problem'] = True
                do_term_not_found_counter += 1
            for gal4_symbol in gal4_info['gal4_input']:
                converted_gal4_symbol = sgml_to_plain_text(allele_symbol).strip().replace('\\', '\\\\')
                gal4_rgx = r'(GAL4|lexA|QF)'
                if not re.search(gal4_rgx, converted_gal4_symbol):
                    self.log.error(f'Line={line_number}: symbol given does not seem to represent a driver: "{gal4_symbol}".')
                    gal4_info['problem'] = True
                    gal4_not_found_counter += 1
                    continue
                try:
                    gal4_info['gal4_ids'].append(self.allele_name_lookup[converted_gal4_symbol]['uniquename'])
                except KeyError:
                    self.log.error(f'Line={line_number}: could not find driver "{gal4_symbol}" in chado.')
                    gal4_info['problem'] = True
                    gal4_not_found_counter += 1
            # Map info to annotation.
            if gal4_info['qualifier'] in self.disease_genetic_modifier_terms.keys():
                gal4_info['modifier_id'] = gal4_info['allele_feature_id']
                gal4_info['modifier_role'] = self.disease_genetic_modifier_terms[gal4_info['qualifier']]
                gal4_info['modeled_by'].extend(gal4_info['additional_allele_ids'])
                gal4_info['eco_abbr'] = 'CEA'
            else:
                if gal4_info['qualifier'] == 'DOES NOT model':
                    gal4_info['is_not'] = True
                gal4_info['modeled_by'].append(gal4_info['allele_feature_id'])
                gal4_info['modeled_by'].extend(gal4_info['additional_allele_ids'])
                gal4_info['eco_abbr'] = gal4_info['evi_code'][0:3]
            # Build an annotation descriptor.
            if gal4_info['problem'] is True:
                continue
            gal4_info['unique_key'] = f'{gal4_info["pub_given"]}_'
            if gal4_info['is_not']:
                gal4_info['unique_key'] += 'NOT_'
            gal4_info['unique_key'] += f'model={"|".join(sorted(gal4_info["modeled_by"]))}_'
            gal4_info['unique_key'] += f'disease_term={gal4_info["doid_term_curie"]}_'
            gal4_info['unique_key'] += f'eco_code={gal4_info["eco_abbr"]}'
            if gal4_info['modifier_id']:
                gal4_info['unique_key'] += f'_{gal4_info["modifier_role"]}={gal4_info["modifier_id"]}'
            self.gal4_dict[gal4_info['unique_key']].append(gal4_info)
            # Find the matching disease annotation.
            if gal4_info['unique_key'] in self.uniq_dis_dict.keys():
                matched_dis_anno_counter += 1
            else:
                self.log.info(f'Could not find dis anno for line={line_number}; unique_key={gal4_info["unique_key"]}')
                unmatched_dis_anno_counter += 1
        self.log.info(f'Could not find {pub_not_found_counter} pubs.')
        self.log.info(f'Could not find {allele_not_found_counter} alleles.')
        self.log.info(f'Could not find {additional_allele_not_found_counter} additional alleles.')
        self.log.info(f'Could not find {do_term_not_found_counter} DO terms.')
        self.log.info(f'Could not find {gal4_not_found_counter} Gal4s.')
        self.log.info(f'Could not find {dis_anno_not_found} disease annotations.')
        self.log.info(f'Found dis anno for {matched_dis_anno_counter} Gal4 lines.')
        self.log.info(f'Found NO dis anno for {unmatched_dis_anno_counter} Gal4 lines.')
        return

    # BOB: to do.
    def integrate_aberration_info(self):
        """Integrate aberration info into annotations."""
        # self.log.info('Integrate aberration info into annotations.')
        return

    # BOB: to do.
    def get_genotype(self, session):
        """Get or create the appropriate genotype for each disease annotation."""
        # self.log.info('Get or create the appropriate genotype for each disease annotation.')
        return

    # Elaborate on get_datatype_data() for the AGMDiseaseHandler.
    def get_datatype_data(self, session):
        """Extend the method for the AGMDiseaseHandler."""
        super().get_datatype_data(session)
        self.get_allele_disease_annotations(session)
        self.get_disease_qualifiers(session)
        self.get_disease_evidence_codes(session)
        self.get_disease_timestamps(session)
        self.extract_text_embedded_alleles(session)
        self.extract_model_and_modifiers()
        self.build_model_eco_lookup()
        self.lookup_eco_codes_for_modifier_annotations()
        self.calculate_annotation_unique_key()
        self.integrate_gal4_info()
        self.integrate_aberration_info()
        self.get_genotype(session)
        return

    # Add methods to be run by synthesize_info() below.
    def flag_problematic_annotations(self):
        """Flag internal annotations."""
        self.log.info('Flag internal annotations.')
        problem_counter = {
            'Obsolete modifier ID': 0,
        }
        no_export_counter = 0
        for dis_anno in self.fb_data_entities.values():
            export_checks = {
                dis_anno.modifier_problem is True: 'Obsolete modifier ID',
            }
            for check, msg in export_checks.items():
                if check:
                    dis_anno.for_export = False
                    dis_anno.export_warnings.append(msg)
                    problem_counter[msg] += 1
            if dis_anno.for_export is False:
                no_export_counter += 1
        self.log.info(f'{no_export_counter} annotations flagged as unexportable in early checking.')
        for problem, problem_count in problem_counter.items():
            self.log.info(f'Problem: "{problem}", count={problem_count}')
        return

    # Elaborate on synthesize_info() for the AGMDiseaseHandler.
    def synthesize_info(self):
        """Extend the method for the AGMDiseaseHandler."""
        super().synthesize_info()
        self.flag_problematic_annotations()
        return

    # # Add methods to be run by map_fb_data_to_alliance() below.
    # def map_allele_disease_annotation_basic(self):
    #     """Map basic FlyBase allele disease annotation to the Alliance LinkML object."""
    #     self.log.info('Map basic FlyBase allele disease annotation to the Alliance LinkML object.')
    #     for dis_anno in self.fb_data_entities.values():
    #         if dis_anno.for_export is False:
    #             continue
    #         allele_curie = f'FB:{dis_anno.feature_cvterm.feature.uniquename}'
    #         do_curie = f'DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}'
    #         pub_curie = self.lookup_single_pub_curie(dis_anno.feature_cvterm.pub_id)
    #         agr_dis_anno = self.agr_export_type(allele_curie, do_curie, pub_curie)
    #         if dis_anno.qualifier.value == 'DOES NOT model':
    #             agr_dis_anno.negated = True
    #         if dis_anno.evidence_code.value.startswith('CEC'):
    #             agr_dis_anno.evidence_code_curies.append(self.evidence_code_xrefs['CEC'])
    #         else:
    #             agr_dis_anno.evidence_code_curies.append(self.evidence_code_xrefs['CEA'])
    #         if dis_anno.fb_modifier_type in self.disease_genetic_modifier_terms.keys():
    #             agr_dis_anno.disease_genetic_modifier_relation_name = self.disease_genetic_modifier_terms[dis_anno.fb_modifier_type]
    #         if dis_anno.fb_modifier_id:
    #             agr_dis_anno.disease_genetic_modifier_identifiers = [f'FB:{dis_anno.fb_modifier_id}']
    #         agr_dis_anno.inferred_gene_identifier = dis_anno.preferred_gene_curie
    #         dis_anno.linkmldto = agr_dis_anno
    #     return

    # def map_data_provider_dto(self):
    #     """Return the DataProviderDTO for the annotation."""
    #     self.log.info('Map data provider.')
    #     for dis_anno in self.fb_data_entities.values():
    #         if dis_anno.for_export is False:
    #             continue
    #         dp_xref = agr_datatypes.CrossReferenceDTO('DOID', dis_anno.linkmldto.do_term_curie, 'disease/fb', dis_anno.linkmldto.do_term_curie).dict_export()
    #         dis_anno.linkmldto.data_provider_dto = agr_datatypes.DataProviderDTO(dp_xref).dict_export()
    #     return

    # def derive_uniq_key(self):
    #     """Derive the unique key based on defining aspects of Alliance disease annotation."""
    #     self.log.info('Derive the unique key based on defining aspects of Alliance disease annotation.')
    #     for dis_anno in self.fb_data_entities.values():
    #         if dis_anno.for_export is False:
    #             continue
    #         dis_anno.uniq_key = f'{dis_anno.linkmldto.allele_identifier}'
    #         dis_anno.uniq_key += f'||{dis_anno.linkmldto.do_term_curie}'
    #         dis_anno.uniq_key += f'||{dis_anno.linkmldto.disease_relation_name}'
    #         dis_anno.uniq_key += f'||{dis_anno.linkmldto.negated}'
    #         dis_anno.uniq_key += f'||{dis_anno.linkmldto.reference_curie}'
    #         evi_codes = sorted(list(set(dis_anno.linkmldto.evidence_code_curies)))
    #         evi_code_str = '|'.join(evi_codes)
    #         dis_anno.uniq_key += f'||{evi_code_str}'
    #         if dis_anno.linkmldto.disease_genetic_modifier_identifiers:
    #             dis_anno.uniq_key += f'||{dis_anno.linkmldto.disease_genetic_modifier_identifiers[0]}'
    #         else:
    #             dis_anno.uniq_key += f'{None}'
    #         dis_anno.uniq_key += f'||{dis_anno.linkmldto.disease_genetic_modifier_relation_name}'
    #     return

    # def group_dis_annos(self):
    #     """Group redundant disease annotations."""
    #     self.log.info('Group redundant disease annotations.')
    #     input_counter = 0
    #     redundant_counter = 0
    #     for dis_anno in self.fb_data_entities.values():
    #         if dis_anno.for_export is False:
    #             continue
    #         input_counter += 1
    #         try:
    #             self.uniq_dis_dict[dis_anno.uniq_key].append(dis_anno.db_primary_id)
    #         except KeyError:
    #             self.uniq_dis_dict[dis_anno.uniq_key] = [dis_anno.db_primary_id]
    #     grouped_counter = len(self.uniq_dis_dict.keys())
    #     self.log.info(f'Found {grouped_counter} unique keys for {input_counter} exportable disease annotations.')
    #     # Flag redundant disease annotations.
    #     for uniq_key, db_primary_ids in self.uniq_dis_dict.items():
    #         if len(db_primary_ids) > 1:
    #             self.log.warning(f'REDUNDANT: AGR_UNIQ_KEY: {uniq_key}')
    #             first_db_id = min(db_primary_ids)
    #             for db_primary_id in db_primary_ids:
    #                 self.log.debug(f'REDUNDANT: {self.fb_data_entities[db_primary_id]}')
    #                 if db_primary_id != first_db_id:
    #                     self.fb_data_entities[db_primary_id].is_redundant = True
    #                     self.fb_data_entities[db_primary_id].for_export = False
    #                     self.fb_data_entities[db_primary_id].export_warnings.append('Annotation is redundant')
    #                     redundant_counter += 1
    #     self.log.info(f'A further {redundant_counter} redundant annotations blocked from export.')
    #     return

    # Elaborate on map_fb_data_to_alliance() for the AGMDiseaseHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AGMDiseaseHandler."""
        super().map_fb_data_to_alliance()
        # self.map_allele_disease_annotation_basic()
        # self.map_data_provider_dto()
        # self.map_timestamps()
        # self.derive_uniq_key()
        # self.group_dis_annos()
        return


class AlleleDiseaseHandler(DataHandler):
    """A data handler for allele-based disease annotations."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AlleleDiseaseHandler object."""
        super().__init__(log, testing)
        self.datatype = 'disease'
        self.fb_export_type = fb_datatypes.FBAlleleDiseaseAnnotation
        self.agr_export_type = agr_datatypes.AlleleDiseaseAnnotationDTO
        self.primary_export_set = 'disease_allele_ingest_set'

    # A dict of unique disease annotation attributes.
    # This will be used after initial data pull to filter out redundant disease annotations.
    uniq_dis_dict = {}

    # Key disease annotation term sets and look ups.
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

    # Elaborate on get_general_data() for the AlleleDiseaseHandler.
    def get_general_data(self, session):
        """Extend the method for the AlleleDiseaseHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_allele_disease_annotations(self, session):
        """Get allele-based disease annotations from chado."""
        self.log.info('Get allele-based disease annotations from chado.')
        disease_term = aliased(Cvterm, name='disease_term')
        prop_type = aliased(Cvterm, name='prop_type')
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            Feature.is_obsolete.is_(False),
            Feature.is_analysis.is_(False),
            Cv.name == 'disease_ontology',
            disease_term.is_obsolete == 0,
            Db.name == 'DOID',
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(self.regex['pub']),
            prop_type.name == 'provenance',
        )
        results = session.query(FeatureCvterm, FeatureCvtermprop).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(disease_term, (disease_term.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == disease_term.cv_id)).\
            join(Dbxref, (Dbxref.dbxref_id == disease_term.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Pub, (Pub.pub_id == FeatureCvterm.pub_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(prop_type, (prop_type.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            dis_anno = self.fb_export_type(result.FeatureCvterm, result.FeatureCvtermprop)
            self.fb_data_entities[dis_anno.db_primary_id] = dis_anno
            counter += 1
        self.log.info(f'Found {counter} allele-based disease annotations from chado (excludes annotations to obsolete alleles).')
        return

    def get_disease_qualifiers(self, session):
        """Get disease annotation qualifiers."""
        self.log.info('Get disease annotation qualifiers.')
        filters = (
            Cvterm.name == 'qualifier',
        )
        results = session.query(FeatureCvtermprop).\
            select_from(FeatureCvtermprop).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        qualifier_count = 0
        for qualifier in results:
            db_primary_id = '{}_{}'.format(qualifier.feature_cvterm_id, qualifier.rank)
            try:
                self.fb_data_entities[db_primary_id].qualifier = qualifier
                qualifier_count += 1
            except KeyError:
                pass
        self.log.info('Found {} disease annotation qualifiers.'.format(qualifier_count))
        return

    def get_disease_evidence_codes(self, session):
        """Get disease annotation evidence codes."""
        self.log.info('Get disease annotation evidence codes.')
        filters = (
            Cvterm.name == 'evidence_code',
        )
        results = session.query(FeatureCvtermprop).\
            select_from(FeatureCvtermprop).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        evidence_code_count = 0
        for evidence_code in results:
            db_primary_id = '{}_{}'.format(evidence_code.feature_cvterm_id, evidence_code.rank)
            try:
                self.fb_data_entities[db_primary_id].evidence_code = evidence_code
                evidence_code_count += 1
            except KeyError:
                pass
        self.log.info('Found {} disease annotation evidence codes.'.format(evidence_code_count))
        # Update annotation descriptions now that all info is in.
        for dis_anno in self.fb_data_entities.values():
            dis_anno.set_entity_desc()
        return

    def get_disease_timestamps(self, session):
        """Get timestamps for disease annotations."""
        self.log.info('Get timestamps for disease annotations.')
        # Note - using standard SQL query because SQLAlchemy ORM for audit_chado is not user-friendly.
        audit_chado_query = """
            SELECT DISTINCT fcvt.feature_cvterm_id||'_'||fcvtp.rank, ac.transaction_timestamp
            FROM feature_cvterm fcvt
            JOIN cvterm cvt ON cvt.cvterm_id = fcvt.cvterm_id
            JOIN cv ON (cv.cv_id = cvt.cv_id AND cv.name = 'disease_ontology')
            JOIN feature_cvtermprop fcvtp ON fcvtp.feature_cvterm_id = fcvt.feature_cvterm_id
            JOIN audit_chado ac ON (ac.record_pkey = fcvtp.feature_cvtermprop_id AND ac.audited_table = 'feature_cvtermprop');
        """
        audit_results = session.execute(audit_chado_query).fetchall()
        self.log.info(f'Got {len(audit_results)} audit_chado results to parse.')
        DB_PRIMARY_ID = 0
        TIMESTAMP = 1
        counter = 0
        for row in audit_results:
            try:
                self.fb_data_entities[row[DB_PRIMARY_ID]].timestamps.append(row[TIMESTAMP])
                counter += 1
            except KeyError:
                self.log.debug(f'Could not put this in anno dict: {row}')
        self.log.info(f'Obtained {counter} timestamps for annotations directly from the audit_chado table.')
        return

    def confirm_current_allele_by_uniquename(self, session, uniquename):
        """Confirm that a given uniquename corresponds to a current allele.

        Args:
            session (Session): The session for the query.
            uniquename (str): The allele uniquename to be checked.

        Returns:
            A boolean value: True if uniquename corresponds to a current allele; False otherwise.

        """
        filters = (
            Feature.uniquename == uniquename,
            Feature.uniquename.op('~')(self.regex['allele']),
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
            session (Session): The session for the query.
            old_uniquename (str): The obsolete allele uniquename to be checked.

        Returns:
            None, or the uniquename for the current feature that corresponds to the obsolete feature.

        """
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['allele']),
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
            self.log.debug('For obsolete {}, found one current allele: {}'.format(old_uniquename, curr_uniquenames[0]))
            curr_allele_id = curr_uniquenames[0]
        elif len(curr_uniquenames) > 1:
            self.log.debug('For obsolete {}, found many current alleles: {}'.format(old_uniquename, curr_uniquenames))
            curr_allele_id = None
        else:
            self.log.debug('For obsolete {}, found no current alleles.'.format(old_uniquename))
            curr_allele_id = None
        return curr_allele_id

    def extract_modifiers(self, session):
        """Extract modifiers from annotation text."""
        self.log.info('Extract modifiers from annotation text.')
        allele_regex = r'FBal[0-9]{7}'
        current_modifier_id_counter = 0
        updated_modifier_id_counter = 0
        cannot_update_modifier_id_counter = 0
        modifier_prob_counter = 0
        for dis_anno in self.fb_data_entities.values():
            for fb_term in self.disease_genetic_modifier_terms.keys():
                if fb_term in dis_anno.evidence_code.value:
                    dis_anno.fb_modifier_type = fb_term
                    if re.search(allele_regex, dis_anno.evidence_code.value):
                        allele_id = re.search(allele_regex, dis_anno.evidence_code.value).group(0)
                    if self.confirm_current_allele_by_uniquename(session, allele_id):
                        dis_anno.fb_modifier_id = allele_id
                        current_modifier_id_counter += 1
                    else:
                        # Look up current allele by 2o ID. Use that.
                        curr_allele_id = self.get_current_id_for_allele(session, allele_id)
                        if curr_allele_id:
                            dis_anno.fb_modifier_id = ['FB:{}'.format(curr_allele_id)]
                            dis_anno.modifier_id_was_updated = True
                            updated_modifier_id_counter += 1
                        else:
                            dis_anno.modifier_problem = True
                            cannot_update_modifier_id_counter += 1
            if dis_anno.modifier_problem is True:
                modifier_prob_counter += 1
        self.log.info(f'{current_modifier_id_counter} allele modifiers mentioned use a current allele ID.')
        self.log.info(f'{updated_modifier_id_counter} allele modifiers mentioned use a non-current allele ID mapped to a current allele.')
        self.log.info(f'{cannot_update_modifier_id_counter} allele modifiers mentioned use a non-current allele ID not mappable to a current allele.')
        self.log.info(f'{modifier_prob_counter} disease annotations have one or more problematic allele modifiers.')
        return

    def get_preferred_inferred_gene_curie(self, session):
        """Get preferred curie for the allele's parental gene for each annotation."""
        self.log.info('Get preferred curie for the allele\'s parental gene for each annotation.')
        # MOD for each organism.
        mod_organisms = {
            'Scer': 'SGD',
            'Cele': 'WormBase',
            'Drer': 'ZFIN',
            'Mmus': 'MGI',
            'Rnor': 'RGD',
            'Hsap': 'HGNC',
            'Xlae': 'Xenbase',
            'Xtro': 'Xenbase',
        }
        input_counter = 0
        gene_counter = 0
        all_mod_curie_counter = 0
        for dis_anno in self.fb_data_entities.values():
            input_counter += 1
            filters = (
                FeatureRelationship.subject_id == dis_anno.feature_cvterm.feature_id,
                Cvterm.name == 'alleleof',
                Feature.is_obsolete.is_(False),
                Feature.uniquename.op('~')(self.regex['gene']),
            )
            parent_gene = session.query(Feature).\
                select_from(Feature).\
                join(FeatureRelationship, (FeatureRelationship.object_id == Feature.feature_id)).\
                join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
                filter(*filters).\
                one()
            gene_curie = None
            # By default, use FB gene curie.
            if parent_gene:
                gene_curie = f'FB:{parent_gene.uniquename}'
                gene_counter += 1
            # Look for a MOD curie and use it instead, if available and unambiguous.
            if parent_gene and parent_gene.organism.abbreviation in mod_organisms.keys():
                filters = (
                    FeatureDbxref.feature_id == parent_gene.feature_id,
                    FeatureDbxref.is_current.is_(True),
                    Db.name == mod_organisms[parent_gene.organism.abbreviation]
                )
                mod_curies = session.query(Dbxref).\
                    select_from(FeatureDbxref).\
                    join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
                    join(Db, (Db.db_id == Dbxref.db_id)).\
                    filter(*filters).\
                    distinct()
                mod_xref_counter = 0
                for mod_curie in mod_curies:
                    mod_xref_counter += 1
                if mod_xref_counter == 1:
                    gene_curie = f'{mod_organisms[parent_gene.organism.abbreviation]}:{mod_curie.accession}'
                    gene_curie = gene_curie.replace('WormBase', 'WB')
                    gene_curie = gene_curie.replace('MGI:MGI:', 'MGI:')
                    gene_curie = gene_curie.replace('GENEPAGE', 'GENE')
                    all_mod_curie_counter += 1
            dis_anno.preferred_gene_curie = gene_curie
        self.log.info(f'Found {gene_counter} parental genes for {input_counter} disease annotations.')
        self.log.info(f'Found {all_mod_curie_counter} MOD curies preferred to FB curies.')
        return

    # Elaborate on get_datatype_data() for the AlleleDiseaseHandler.
    def get_datatype_data(self, session):
        """Extend the method for the AlleleDiseaseHandler."""
        super().get_datatype_data(session)
        self.get_allele_disease_annotations(session)
        self.get_disease_qualifiers(session)
        self.get_disease_evidence_codes(session)
        self.get_disease_timestamps(session)
        self.extract_modifiers(session)
        self.get_preferred_inferred_gene_curie(session)
        return

    # Add methods to be run by synthesize_info() below.
    def flag_unexportable_annotations(self):
        """Flag internal annotations."""
        self.log.info('Flag internal annotations.')
        problem_counter = {
            'Not model annotation': 0,
            'Multi-allele model': 0,
            'Obsolete modifier ID': 0,
        }
        no_export_counter = 0
        for dis_anno in self.fb_data_entities.values():
            export_checks = {
                dis_anno.qualifier.value not in self.relevant_qualifiers: 'Not model annotation',
                ' with FLYBASE' in dis_anno.evidence_code.value: 'Multi-allele model',
                dis_anno.modifier_problem is True: 'Obsolete modifier ID',
                dis_anno.modifier_id_was_updated is True: 'Obsolete modifier ID',
            }
            for check, msg in export_checks.items():
                if check:
                    dis_anno.for_export = False
                    dis_anno.export_warnings.append(msg)
                    problem_counter[msg] += 1
            if dis_anno.for_export is False:
                no_export_counter += 1
        self.log.info(f'{no_export_counter} annotations flagged as unexportable in early checking.')
        for problem, problem_count in problem_counter.items():
            self.log.info(f'Problem: "{problem}", count={problem_count}')
        return

    # Elaborate on synthesize_info() for the AlleleDiseaseHandler.
    def synthesize_info(self):
        """Extend the method for the AlleleDiseaseHandler."""
        super().synthesize_info()
        self.flag_unexportable_annotations()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_allele_disease_annotation_basic(self):
        """Map basic FlyBase allele disease annotation to the Alliance LinkML object."""
        self.log.info('Map basic FlyBase allele disease annotation to the Alliance LinkML object.')
        for dis_anno in self.fb_data_entities.values():
            if dis_anno.for_export is False:
                continue
            allele_curie = f'FB:{dis_anno.feature_cvterm.feature.uniquename}'
            do_curie = f'DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}'
            pub_curie = self.lookup_single_pub_curie(dis_anno.feature_cvterm.pub_id)
            agr_dis_anno = self.agr_export_type(allele_curie, do_curie, pub_curie)
            if dis_anno.qualifier.value == 'DOES NOT model':
                agr_dis_anno.negated = True
            if dis_anno.evidence_code.value.startswith('CEC'):
                agr_dis_anno.evidence_code_curies.append(self.evidence_code_xrefs['CEC'])
            else:
                agr_dis_anno.evidence_code_curies.append(self.evidence_code_xrefs['CEA'])
            if dis_anno.fb_modifier_type in self.disease_genetic_modifier_terms.keys():
                agr_dis_anno.disease_genetic_modifier_relation_name = self.disease_genetic_modifier_terms[dis_anno.fb_modifier_type]
            if dis_anno.fb_modifier_id:
                agr_dis_anno.disease_genetic_modifier_identifiers = [f'FB:{dis_anno.fb_modifier_id}']
            agr_dis_anno.inferred_gene_identifier = dis_anno.preferred_gene_curie
            dis_anno.linkmldto = agr_dis_anno
        return

    def map_data_provider_dto(self):
        """Return the DataProviderDTO for the annotation."""
        self.log.info('Map data provider.')
        for dis_anno in self.fb_data_entities.values():
            if dis_anno.for_export is False:
                continue
            dp_xref = agr_datatypes.CrossReferenceDTO('DOID', dis_anno.linkmldto.do_term_curie, 'disease/fb', dis_anno.linkmldto.do_term_curie).dict_export()
            dis_anno.linkmldto.data_provider_dto = agr_datatypes.DataProviderDTO(dp_xref).dict_export()
        return

    def derive_uniq_key(self):
        """Derive the unique key based on defining aspects of Alliance disease annotation."""
        self.log.info('Derive the unique key based on defining aspects of Alliance disease annotation.')
        for dis_anno in self.fb_data_entities.values():
            if dis_anno.for_export is False:
                continue
            dis_anno.uniq_key = f'{dis_anno.linkmldto.allele_identifier}'
            dis_anno.uniq_key += f'||{dis_anno.linkmldto.do_term_curie}'
            dis_anno.uniq_key += f'||{dis_anno.linkmldto.disease_relation_name}'
            dis_anno.uniq_key += f'||{dis_anno.linkmldto.negated}'
            dis_anno.uniq_key += f'||{dis_anno.linkmldto.reference_curie}'
            evi_codes = sorted(list(set(dis_anno.linkmldto.evidence_code_curies)))
            evi_code_str = '|'.join(evi_codes)
            dis_anno.uniq_key += f'||{evi_code_str}'
            if dis_anno.linkmldto.disease_genetic_modifier_identifiers:
                dis_anno.uniq_key += f'||{dis_anno.linkmldto.disease_genetic_modifier_identifiers[0]}'
            else:
                dis_anno.uniq_key += f'{None}'
            dis_anno.uniq_key += f'||{dis_anno.linkmldto.disease_genetic_modifier_relation_name}'
        return

    def group_dis_annos(self):
        """Group redundant disease annotations."""
        self.log.info('Group redundant disease annotations.')
        input_counter = 0
        redundant_counter = 0
        for dis_anno in self.fb_data_entities.values():
            if dis_anno.for_export is False:
                continue
            input_counter += 1
            try:
                self.uniq_dis_dict[dis_anno.uniq_key].append(dis_anno.db_primary_id)
            except KeyError:
                self.uniq_dis_dict[dis_anno.uniq_key] = [dis_anno.db_primary_id]
        grouped_counter = len(self.uniq_dis_dict.keys())
        self.log.info(f'Found {grouped_counter} unique keys for {input_counter} exportable disease annotations.')
        # Flag redundant disease annotations.
        for uniq_key, db_primary_ids in self.uniq_dis_dict.items():
            if len(db_primary_ids) > 1:
                self.log.warning(f'REDUNDANT: AGR_UNIQ_KEY: {uniq_key}')
                first_db_id = min(db_primary_ids)
                for db_primary_id in db_primary_ids:
                    self.log.debug(f'REDUNDANT: {self.fb_data_entities[db_primary_id]}')
                    if db_primary_id != first_db_id:
                        self.fb_data_entities[db_primary_id].is_redundant = True
                        self.fb_data_entities[db_primary_id].for_export = False
                        self.fb_data_entities[db_primary_id].export_warnings.append('Annotation is redundant')
                        redundant_counter += 1
        self.log.info(f'A further {redundant_counter} redundant annotations blocked from export.')
        return

    # Elaborate on map_fb_data_to_alliance() for the AlleleDiseaseHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AlleleDiseaseHandler."""
        super().map_fb_data_to_alliance()
        self.map_allele_disease_annotation_basic()
        self.map_data_provider_dto()
        self.map_timestamps()
        self.derive_uniq_key()
        self.group_dis_annos()
        return
