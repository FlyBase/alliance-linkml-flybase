"""Module:: disease_handlers.

Synopsis:
    Data handlers for FlyBase disease annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import csv
import re
from collections import defaultdict
from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.char_conversions import sgml_to_plain_text
from harvdev_utils.genotype_utilities import GenotypeAnnotation
from harvdev_utils.reporting import (
    Cv, Cvterm, Cvtermsynonym, Db, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureDbxref, FeaturePub, FeatureRelationship, FeatureSynonym, Pub, Synonym
)
import agr_datatypes
import fb_datatypes
from handler import DataHandler


class AGMDiseaseHandler(DataHandler):
    """A data handler that converts allele-based disease annotations into genotype (AGM) annotations.

    Note that the path from chado annotations to final Alliance AGM-level annotations is complicated.
    1. This handler gets allele-level disease annotations from chado (self.allele_dis_annos).
       Note that there are duplicated allele-level annotations in chado.
    2. Those allele-level annotations are converted into genotype-level annotations (self.genotype_dis_annos).
       Note that distinct allele-level annotations can represent the same genotype from different perspectives.
    3. So, genotype-level annotations are grouped to deal with duplications and redundancies.
    4. This object processes driver line info (from TSV file) and folds it into genotype-level annotations.
    5. Because there can be many driver combinations for a given genotype, there can be an expansion of annotations.
       So, new annotations (for each genotype-driver combination) are created (self.fb_data_entities).
    6. Then, a separate aberration TSV file is processed to add even more genotype-level annotations to self.fb_data_entities.

    """
    def __init__(self, log: Logger, testing: bool):
        """Create the AGMDiseaseHandler object."""
        super().__init__(log, testing)
        self.datatype = 'agm_disease'
        self.fb_export_type = fb_datatypes.FBAlleleDiseaseAnnotation
        self.agr_export_type = agr_datatypes.AGMDiseaseAnnotationDTO
        self.primary_export_set = 'disease_agm_ingest_set'
        self.allele_dis_annos = {}                   # Initial allele-level disease annotations.
        self.genotype_dis_annos = {}                 # Intermediate genotype-level disease annotations.
        self.allele_name_lookup = {}                 # feature.name-keyed feature dicts for alleles and aberrations, current only.
        self.gene_name_lookup = {}                   # feature.name-keyed feature dicts for genes, current only.
        self.doid_term_lookup = {}                   # cvterm.name-keyed cvterm dicts.
        self.model_eco_lookup = defaultdict(list)    # Evidence abbreviation lookup for "model_of" annotations.
        self.driver_dict = defaultdict(list)         # Unique disease descriptors key lists of driver info to integrate.
        self.aberr_dict = defaultdict(list)          # Unique disease descriptors key lists of aberration info to integrate.
        self.rejected_driver_info = []               # List of driver dicts missing chado dis anno.
        self.misxprn_allele_ids = []                 # Allele IDs for misexpression elements.

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
            if feature['type'] in ['allele', 'chromosome_structure_variation'] and feature['is_obsolete'] is False:
                self.allele_name_lookup[feature['name']] = feature
        self.log.info(f'Have {len(self.allele_name_lookup)} current alleles in the allele-by-name lookup.')
        return

    def build_gene_name_lookup(self):
        """Build name-keyed dict of genes."""
        self.log.info('Build name-keyed dict of genes.')
        for feature in self.feature_lookup.values():
            if feature['type'] == 'gene' and feature['is_obsolete'] is False:
                self.gene_name_lookup[feature['name']] = feature
        self.log.info(f'Have {len(self.gene_name_lookup)} current genes in the gene-by-name lookup.')
        return

    def build_doid_term_lookup(self):
        """Build name-keyed dict of DOID terms."""
        self.log.info('Build name-keyed dict of DOID terms.')
        for term in self.cvterm_lookup.values():
            if term['cv_name'] == 'disease_ontology' and term['db_name'] == 'DOID':
                self.doid_term_lookup[term['name']] = term
        self.log.info(f'Have {len(self.cvterm_lookup)} DOID terms in the CV term-by-name lookup.')
        return

    def get_aberration_pub_info(self, session):
        """Get pub info for aberrations."""
        self.log.info('Get pub info for aberrations.')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['aberration']),
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(self.regex['pub']),
            Pub.uniquename != 'unattributed',
        )
        results = session.query(FeaturePub).\
            select_from(Feature).\
            join(FeaturePub, (FeaturePub.feature_id == Feature.feature_id)).\
            join(Pub, (Pub.pub_id == FeaturePub.pub_id)).\
            filter(*filters).\
            distinct()
        ab_counter = 0
        pub_counter = 0
        for result in results:
            try:
                self.feature_lookup[result.feature_id]['pub_ids'].append(result.pub_id)
                pub_counter += 1
            except KeyError:
                self.feature_lookup[result.feature_id]['pub_ids'] = [result.pub_id]
                ab_counter += 1
                pub_counter += 1
        self.log.info(f'Found {pub_counter} pubs for {ab_counter} aberrations.')
        return

    def get_aberration_gene_info(self, session):
        """Get gene info for aberrations."""
        self.log.info('Get gene info for aberrations.')
        aberration = aliased(Feature, name='aberration')
        gene = aliased(Feature, name='gene')
        fr_types = [
            'deletes',
            'duplicates',
            'molec_deletes',
            'molec_dups',
            'molec_partdeletes',
            'molec_partdups',
            'part_deletes',
            'part_duplicates',
            'useful_Df_direct',
            'useful_Df_from_cyto',
            'useful_Dp_direct',
            'useful_Dp_from_cyto',
        ]
        filters = (
            aberration.is_obsolete.is_(False),
            aberration.uniquename.op('~')(self.regex['aberration']),
            gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(self.regex['gene']),
            Cvterm.name.in_((fr_types)),
        )
        results = session.query(aberration, Cvterm, gene).\
            select_from(aberration).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == aberration.feature_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        ab_counter = 0
        gene_counter = 0
        for result in results:
            # self.log.debug(f'{result.aberration.uniquename} {result.Cvterm.name} {result.gene.name}')
            try:
                self.feature_lookup[result.aberration.feature_id]['affected_genes'].append(result.gene.feature_id)
                gene_counter += 1
            except KeyError:
                self.feature_lookup[result.aberration.feature_id]['affected_genes'] = [result.gene.feature_id]
                ab_counter += 1
                gene_counter += 1
        self.log.info(f'Found {gene_counter} genes with useful relationships to {ab_counter} aberrations.')
        return

    def get_misexpression_elements(self, session):
        """Get misexpression alleles."""
        self.log.info('Get misexpression alleles.')
        allele_feature = aliased(Feature, name='allele_feature')
        construct_feature = aliased(Feature, name='construct_feature')
        insertion_feature = aliased(Feature, name='insertion_feature')
        allele_insertion_rel = aliased(FeatureRelationship, name='allele_insertion_rel')
        insertion_construct_rel = aliased(FeatureRelationship, name='insertion_construct_rel')
        ai_rel_type = aliased(Cvterm, name='ai_rel_type')
        ic_rel_type = aliased(Cvterm, name='ic_rel_type')
        tool_type = aliased(Cvterm, name='tool_type')
        tool_rel = aliased(Cvterm, name='tool_rel')
        filters = (
            allele_feature.is_obsolete.is_(False),
            allele_feature.uniquename.op('~')(self.regex['allele']),
            construct_feature.uniquename.op('~')(self.regex['construct']),
            construct_feature.is_obsolete.is_(False),
            insertion_feature.uniquename.op('~')(self.regex['insertion']),
            insertion_feature.is_obsolete.is_(False),
            ai_rel_type.name == 'associated_with',
            ic_rel_type.name == 'producedby',
            tool_type.name == 'misexpression element',
            tool_rel.name == 'tool_uses'
        )
        results = session.query(allele_feature).\
            select_from(allele_feature).\
            join(allele_insertion_rel, (allele_insertion_rel.subject_id == allele_feature.feature_id)).\
            join(insertion_feature, (insertion_feature.feature_id == allele_insertion_rel.object_id)).\
            join(ai_rel_type, (ai_rel_type.cvterm_id == allele_insertion_rel.type_id)).\
            join(insertion_construct_rel, (insertion_construct_rel.subject_id == insertion_feature.feature_id)).\
            join(construct_feature, (construct_feature.feature_id == insertion_construct_rel.object_id)).\
            join(ic_rel_type, (ic_rel_type.cvterm_id == insertion_construct_rel.type_id)).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == construct_feature.feature_id)).\
            join(tool_type, (tool_type.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(tool_rel, (tool_rel.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.misxprn_allele_ids.append(result.feature_id)
            counter += 1
        self.log.info(f'Found {counter} misexpression alleles in chado.')
        return

    # Elaborate on get_general_data() for the AGMDiseaseHandler.
    def get_general_data(self, session):
        """Extend the method for the AGMDiseaseHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'gene', 'insertion', 'construct'])
        self.get_aberration_pub_info(session)
        self.get_aberration_gene_info(session)
        self.build_uname_feature_lookup()
        self.get_transgenic_allele_ids(session)
        self.get_in_vitro_allele_ids(session)
        self.build_allele_gene_lookup(session)
        self.build_allele_name_lookup()
        self.build_gene_name_lookup()
        self.build_doid_term_lookup()
        self.get_misexpression_elements(session)
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
            self.allele_dis_annos[dis_anno.db_primary_id] = dis_anno
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
                self.allele_dis_annos[db_primary_id].qualifier = qualifier
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
                self.allele_dis_annos[db_primary_id].evidence_code = evidence_code
                evidence_code_count += 1
            except KeyError:
                pass
        self.log.info('Found {} disease annotation evidence codes.'.format(evidence_code_count))
        # Update annotation descriptions now that all info is in.
        for dis_anno in self.allele_dis_annos.values():
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
                self.allele_dis_annos[row[DB_PRIMARY_ID]].timestamps.append(row[TIMESTAMP])
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
        for dis_anno in self.allele_dis_annos.values():
            if not re.search(embedded_allele_regex, dis_anno.evidence_code.value):
                continue
            embedded_allele_ids = re.findall(embedded_allele_regex, dis_anno.evidence_code.value)
            for allele_id in embedded_allele_ids:
                if self.confirm_current_allele_by_uniquename(session, allele_id):
                    dis_anno.text_embedded_allele_curies.append(allele_id)
                    current_allele_id_counter += 1
                else:
                    curr_allele_id = self.get_current_id_for_allele(session, allele_id)
                    if curr_allele_id:
                        dis_anno.text_embedded_allele_curies.append(curr_allele_id)
                        dis_anno.allele_id_was_updated = True
                        updated_allele_id_counter += 1
                    else:
                        self.log.error(f'Could not update allele ID {allele_id} for {dis_anno}.')
                        dis_anno.allele_id_problem = True
                        dis_anno.export_warnings.append('Obsolete allele ID could not be updated')
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
        for dis_anno in self.allele_dis_annos.values():
            if dis_anno.export_warnings:
                continue
            dis_anno.modeled_by.extend(dis_anno.text_embedded_allele_curies)
            allele_subject_curie = self.feature_lookup[dis_anno.feature_cvterm.feature_id]['uniquename']
            if dis_anno.qualifier.value in ('model of', 'DOES NOT model'):
                dis_anno.modeled_by.append(allele_subject_curie)
                if dis_anno.qualifier.value == 'DOES NOT model':
                    dis_anno.is_not = True
            else:
                dis_anno.modifier_curie = allele_subject_curie
                dis_anno.modifier_role = self.disease_genetic_modifier_terms[dis_anno.qualifier.value]
        return

    def get_parent_genes(self):
        """Get parent genes for key alleles."""
        self.log.info('Get parent genes for key alleles.')
        counter = 0
        for dis_anno in self.allele_dis_annos.values():
            if dis_anno.export_warnings:
                continue
            key_alleles = []
            key_alleles.extend(dis_anno.modeled_by)
            if dis_anno.modifier_curie:
                key_alleles.append(dis_anno.modifier_curie)
            key_alleles = set(key_alleles)
            for fbal_id in key_alleles:
                allele_feature_id = self.uname_feature_lookup[fbal_id]['feature_id']
                dis_anno.parent_gene_ids.add(self.allele_gene_lookup[allele_feature_id])
            counter += len(dis_anno.parent_gene_ids)
        self.log.info(f'Found {counter} parent genes for key alleles of disease annotations.')
        return

    def find_relevant_aberrations(self):
        """Find aberrations that might be relevant to the disease annotations."""
        self.log.info('Find aberrations that might be relevant to the disease annotations.')
        pair_counter = 0
        match_counter = 0
        for feature in self.feature_lookup.values():
            if feature['type'] != 'chromosome_structure_variation':
                continue
            elif 'pub_ids' not in feature.keys():
                continue
            for dis_anno in self.allele_dis_annos.values():
                if dis_anno.feature_cvterm.pub_id in feature['pub_ids']:
                    if 'affected_genes' not in feature:
                        continue
                    if dis_anno.parent_gene_ids.intersection(set(feature['affected_genes'])):
                        dis_anno.possible_aberrations.add(feature['feature_id'])
                        msg = f'{dis_anno.feature_cvterm.pub.uniquename}\t'
                        msg += f'{feature["uniquename"]}\t'
                        msg += f'{feature["name"]}\t'
                        msg += f'{dis_anno.feature_cvterm.feature.uniquename}\t'
                        msg += f'{dis_anno.feature_cvterm.feature.name}\t'
                        msg += f'{dis_anno.feature_cvterm.cvterm.name}\t'
                        msg += f'DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}\t'
                        msg += f'{dis_anno.qualifier.value}\t'
                        msg += f'{dis_anno.evidence_code.value}'
                        self.log.debug(f'Found FBab-dis_anno match: {msg}')
                        match_counter += 1
                pair_counter += 1
                if pair_counter % 1000000 == 0:
                    self.log.debug(f'Have checked {pair_counter} aberration-dis_anno pairs.')
        self.log.info(f'Found {match_counter} potential aberration-dis_anno matches.')
        return

    def build_model_eco_lookup(self):
        """Build ECO lookup from model-type annotations for modifier annotations."""
        self.log.info('Build ECO lookup from model-type annotations for modifier annotations.')
        counter = 0
        for dis_anno in self.allele_dis_annos.values():
            if dis_anno.export_warnings:
                continue
            # Determine unique key for the model in this annotation.
            dis_anno.model_unique_key = f'{dis_anno.feature_cvterm.pub.uniquename}_'
            if dis_anno.is_not:
                dis_anno.unique_key += 'NOT_'
            dis_anno.model_unique_key = f'model={"|".join(sorted(dis_anno.modeled_by))}_'
            dis_anno.model_unique_key += f'disease_term=DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}'
            # self.log.debug(f'{dis_anno} has model_unique_key={dis_anno.model_unique_key}')
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
                # self.log.debug(f'Many ECOs for model_key={ukey}: {uniqued_list}')
                many_counter += 1
        self.log.info(f'Have {zero_counter} keys in lookup with NO ECO; {one_counter} keys with ONE ECO; {many_counter} keys with MANY ECOS.')
        self.log.info(f'Have these distinct ECO codes: {set(distinct_ecos)}')
        return

    def lookup_eco_codes_for_modifier_annotations(self):
        """Lookup ECO code for modifier-type annotations."""
        self.log.info('Lookup ECO code for modifier-type annotations.')
        input_counter = 0
        skip_counter = 0
        bad_counter = 0
        assess_counter = 0
        match_counter = 0
        no_match_counter = 0
        for dis_anno in self.allele_dis_annos.values():
            input_counter += 1
            if dis_anno.export_warnings:
                bad_counter += 1
                continue
            if dis_anno.eco_abbr:
                skip_counter += 1
                continue
            assess_counter += 1
            if dis_anno.model_unique_key not in self.model_eco_lookup.keys():
                dis_anno.eco_abbr = 'CEC'
                no_match_counter += 1
            elif len(self.model_eco_lookup[dis_anno.model_unique_key]) == 1:
                dis_anno.eco_abbr = self.model_eco_lookup[dis_anno.model_unique_key][0]
                match_counter += 1
            # Choose default CEA if many codes.
            elif len(self.model_eco_lookup[dis_anno.model_unique_key]) > 1:
                dis_anno.eco_abbr = 'CEC'
                match_counter += 1
        self.log.info(f'Skipped {bad_counter}/{input_counter} problematic allele-level annotations.')
        self.log.info(f'Found {skip_counter}/{input_counter} annotations with explicitly curated ECO codes.')
        self.log.info(f'Filled in missing ECO codes for {assess_counter}/{input_counter} annotations.')
        self.log.info(f'Found ECO for {match_counter} modifier-type annotations.')
        self.log.info(f'Assigned "CEC" for {no_match_counter} modifier-type annotations with no info.')
        return

    def group_redundant_annotations(self):
        """Group allele-level disease annotations into genotype-level disease annotations."""
        self.log.info('Group allele-level disease annotations into genotype-level disease annotations.')
        for dis_anno in self.allele_dis_annos.values():
            if dis_anno.export_warnings:
                continue
            dis_anno.unique_key = f'{dis_anno.feature_cvterm.pub.uniquename}_'
            if dis_anno.is_not:
                dis_anno.unique_key += 'NOT_'
            dis_anno.unique_key += f'model={"|".join(sorted(dis_anno.modeled_by))}_'
            dis_anno.unique_key += f'disease_term=DOID:{dis_anno.feature_cvterm.cvterm.dbxref.accession}_'
            dis_anno.unique_key += f'eco_code={dis_anno.eco_abbr}'
            if not dis_anno.eco_abbr:
                self.log.warning(f'No ECO abbr for {dis_anno.unique_key}')
            if dis_anno.modifier_curie:
                dis_anno.unique_key += f'_{dis_anno.modifier_role}={dis_anno.modifier_curie}'
            self.log.debug(f'Annotation db_primary_id={dis_anno.db_primary_id} has this unique key: {dis_anno.unique_key}')
            if dis_anno.unique_key in self.genotype_dis_annos.keys():
                self.genotype_dis_annos[dis_anno.unique_key].allele_annotations.append(dis_anno)
            else:
                self.genotype_dis_annos[dis_anno.unique_key] = fb_datatypes.FBGenotypeDiseaseAnnotation(dis_anno.unique_key)
                self.genotype_dis_annos[dis_anno.unique_key].allele_annotations.append(dis_anno)
        return

    def find_cvterm_curie_from_name(self, session, cvterm_name):
        """Find the CV term curie from a CV term name or alias; return None if nothing found."""
        cvterm_curie = None
        if not self.doid_term_lookup:
            e = 'Must create handler.doid_term_lookup with handler.build_doid_term_lookup() '
            e += 'before calling the handler.find_cvterm_curie_from_name() method.'
            self.log.critical(e)
            raise Exception(e)
        try:
            cvterm_curie = self.doid_term_lookup[cvterm_name]['curie']
        except KeyError:
            filters = (
                Cvterm.is_obsolete == 0,
                Cv.name == 'disease_ontology',
                Db.name == 'DOID',
                Cvtermsynonym.name == cvterm_name,
            )
            results = session.query(Dbxref).\
                select_from(Cvterm).\
                join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
                join(Dbxref, (Dbxref.dbxref_id == Cvterm.dbxref_id)).\
                join(Db, (Db.db_id == Dbxref.db_id)).\
                join(Cvtermsynonym, (Cvtermsynonym.cvterm_id == Cvterm.cvterm_id)).\
                filter(*filters).\
                distinct()
            curies = []
            for result in results:
                curies.append(f'DOID:{result.accession}')
            if len(curies) == 1:
                cvterm_curie = curies[0]
            elif len(curies) > 1:
                self.log.error(f'Found MANY possible curies for "{cvterm_name}": {curies}')
        if cvterm_curie is None:
            self.log.debug(f'No result for cvterm={cvterm_name}')
        return cvterm_curie

    def find_feature_uniquename_from_name(self, session, feature_symbol, fb_id_rgx):
        """Find the uniquename for a feature from a symbol or synonym.

        Args:
            feature_symbol (str): The feature symbol (plain text or SGML), current or an alias.
            fb_id_rgx (str): The FB ID regex: e.g., r'FBal[0-9]{7}$'.

        Returns:
            The feature uniquename, or, None if there is no answer or an ambiguous answer.

        Raises:
            Raises an exception if the handler (self) does not have a name-based feature lookup dict.

        """
        converted_feature_symbol = sgml_to_plain_text(feature_symbol).strip().strip(',').strip('.')
        uniquename = None
        if 'FBal' in fb_id_rgx or 'FBab' in fb_id_rgx:
            if not self.allele_name_lookup:
                e = 'Must create handler.allele_name_lookup with handler.build_allele_name_lookup() '
                e += 'before calling the handler.find_feature_uniquename() method.'
                self.log.critical(e)
                raise Exception(e)
            try:
                uniquename = self.allele_name_lookup[converted_feature_symbol]['uniquename']
            except KeyError:
                pass
        elif 'FBgn' in fb_id_rgx:
            if not self.gene_name_lookup:
                e = 'Must create handler.gene_name_lookup with handler.build_gene_name_lookup() '
                e += 'before calling the handler.find_feature_uniquename() method.'
                self.log.critical(e)
                raise Exception(e)
            try:
                uniquename = self.gene_name_lookup[converted_feature_symbol]['uniquename']
            except KeyError:
                pass
        else:
            e = f'find_feature_uniquename_from_name() method does not support {fb_id_rgx} yet'
            self.log.critical(e)
            raise Exception(e)
        if uniquename is None:
            filters = (
                Feature.is_obsolete.is_(False),
                Feature.uniquename.op('~')(fb_id_rgx),
                Synonym.name == converted_feature_symbol,
                FeatureSynonym.is_current.is_(False),
            )
            results = session.query(Feature).\
                select_from(Feature).\
                join(FeatureSynonym, (FeatureSynonym.feature_id == Feature.feature_id)).\
                join(Synonym, (Synonym.synonym_id == FeatureSynonym.synonym_id)).\
                filter(*filters).\
                distinct()
            fb_ids = []
            for result in results:
                fb_ids.append(result.uniquename)
            if len(fb_ids) == 1:
                uniquename = fb_ids[0]
            elif len(fb_ids) > 1:
                self.log.error(f'Found MANY possible FB IDs for "{feature_symbol}": {fb_ids}')
        if uniquename is None:
            self.log.debug(f'No result for given symbol={feature_symbol}, converted={converted_feature_symbol}')
        return uniquename

    def check_disease_annotation(self, session, feat, dis, pub, qual, evi_code):
        """Check that a given disease annotations exists, or not.

        Args:
            feat (str): The FBal ID for the feature.
            dis (str): The DO term name.
            pub (str): The pub FBrf ID.
            qual (str): The qualifier value.
            evi_code (str): The evidence_code value.

        Returns:
            exists (boolean): True if any matching annotations found.

        """
        exists = False
        dis_term = aliased(Cvterm, name='dis_term')
        qual_type = aliased(Cvterm, name='qual_type')
        evi_type = aliased(Cvterm, name='evi_type')
        qualp = aliased(FeatureCvtermprop, name='qualp')
        evip = aliased(FeatureCvtermprop, name='evip')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename == feat,
            dis_term.name == dis,
            Pub.uniquename == pub,
            qual_type == 'qualifier',
            qualp.value == qual,
            evi_type == 'evidence_code',
            evip.value == evi_code,
            evip.rank == qualp.rank,
        )
        results = session.query(FeatureCvterm, qualp).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(dis_term, (dis_term.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Pub, (Pub.pub_id == FeatureCvterm.pub_id)).\
            join(qualp, (qualp.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(qual_type, (qual_type.cvterm_id == qualp.type_id)).\
            join(evip, (evip.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(evi_type, (evi_type.cvterm_id == evip.type_id)).\
            filter(*filters).\
            distinct()
        for result in results:
            msg = f'Found {result.FeatureCvterm.feature_cvterm_id}_{result.qualp.rank}'
            self.log.warning(msg)
            exists = True
        return exists

    def parse_driver_info(self, session):
        """Parse driver info to integrate."""
        self.log.info('Parse driver info to integrate.')
        PUB_GIVEN = 0
        ALLELE_SYMBOL = 1
        QUAL = 4
        DO_TERM = 5
        EVI_CODE = 6
        ADDITIONAL_ALLELES = 7
        DRIVER_INPUT = 8
        OPERATION = 12
        file_input = open('/src/output/driver_info.tsv')
        line_number = 0
        input_counter = 0
        malformed_line_counter = 0
        matched_dis_anno_counter = 0
        close_matched_dis_anno_counter = 0
        unmatched_dis_anno_counter = 0
        prob_counter = 0
        self.rejected_driver_info = []
        for i in file_input:
            # self.log.debug(f'Process this line: {i.strip()}')
            line_number += 1
            if not i.startswith('FBrf0'):
                continue
            input_counter += 1
            line = i.split('\t')
            if len(line) != 13:
                self.log.error(f'Line={line_number} is malformed: {i}')
                malformed_line_counter += 1
                continue
            driver_info = {
                # Attributes from input file.
                'line_number': line_number,
                'pub_given': line[PUB_GIVEN].strip(),
                'allele_symbol': sgml_to_plain_text(line[ALLELE_SYMBOL]).strip(),
                'additional_alleles': line[ADDITIONAL_ALLELES].split(' '),
                'qualifier': line[QUAL].strip(),
                'evi_code': line[EVI_CODE].strip(),
                'do_term': line[DO_TERM].strip(),
                'driver_input': line[DRIVER_INPUT].split(' '),
                'operation': line[OPERATION].rstrip(),
                # Attributes to be obtained from chado.
                'pub_id': None,
                'allele_feature_curie': None,
                'additional_allele_curies': [],
                'driver_curies': [],
                'doid_term_curie': None,
                # Attributes synthesized from the above.
                'modeled_by': [],
                'is_not': False,
                'modifier_curie': None,
                'modifier_role': None,
                'eco_abbr': None,
                'unique_key': None,
                'problems': [],
            }
            try:
                driver_info['pub_id'] = self.fbrf_bibliography[driver_info['pub_given']].pub_id
            except KeyError:
                self.log.error(f'Line={line_number}: could not find pub "{driver_info["pub_given"]}" in chado.')
                prob_msg = 'bad pub id'
                driver_info['problems'].append(prob_msg)

            cvterm_curie = self.find_cvterm_curie_from_name(session, driver_info['do_term'])
            if cvterm_curie:
                driver_info['doid_term_curie'] = cvterm_curie
                # self.log.debug(f'Line={line_number}: found "{cvterm_curie}" for "{driver_info["do_term"]}" in chado.')
            else:
                self.log.error(f'Line={line_number}: could not find DO term for "{driver_info["do_term"]}" in chado.')
                prob_msg = 'bad DO term name'
                driver_info['problems'].append(prob_msg)

            allele_curie = self.find_feature_uniquename_from_name(session, driver_info['allele_symbol'], self.regex['allele'])
            if allele_curie:
                driver_info['allele_feature_curie'] = allele_curie
                # self.log.error(f'Line={line_number}: found "{allele_id}" for "{driver_info["allele_symbol"]}" in chado.')
            else:
                self.log.error(f'Line={line_number}: could not find allele "{driver_info["allele_symbol"]}" in chado.')
                prob_msg = 'bad allele symbol'
                driver_info['problems'].append(prob_msg)

            for allele_symbol in driver_info['additional_alleles']:
                if allele_symbol == '' or allele_symbol == ' ' or allele_symbol == '+':
                    continue
                allele_id = self.find_feature_uniquename_from_name(session, allele_symbol, self.regex['allele'])
                if allele_id:
                    driver_info['additional_allele_curies'].append(allele_id)
                    # self.log.error(f'Line={line_number}: found "{allele_id}" for additional allele "{allele_symbol}" in chado.')
                else:
                    self.log.error(f'Line={line_number}: could not find additional allele "{allele_symbol}" in chado.')
                    prob_msg = 'bad additional allele symbol'
                    driver_info['problems'].append(prob_msg)

            for driver_symbol in driver_info['driver_input']:
                if driver_symbol == '' or driver_symbol == ' ' or driver_symbol == '+':
                    continue
                driver_rgx = r'(GAL4|GAL80|lexA|QF|FLP1|Cas|VP16|RELA|G4DBD)'
                if not re.search(driver_rgx, driver_symbol):
                    self.log.warning(f'Line={line_number}: symbol given does not seem to represent a driver: "{driver_symbol}".')
                allele_curie = self.find_feature_uniquename_from_name(session, driver_symbol, self.regex['allele'])
                if allele_curie:
                    driver_info['driver_curies'].append(allele_curie)
                    # self.log.error(f'Line={line_number}: found "{allele_id}" for driver "{driver_symbol}" in chado.')
                else:
                    self.log.error(f'Line={line_number}: could not find driver "{driver_symbol}" in chado.')
                    prob_msg = 'bad driver symbol'
                    driver_info['problems'].append(prob_msg)

            # Map info to annotation.
            if driver_info['qualifier'] in self.disease_genetic_modifier_terms.keys():
                driver_info['modifier_curie'] = driver_info['allele_feature_curie']
                driver_info['modifier_role'] = self.disease_genetic_modifier_terms[driver_info['qualifier']]
                driver_info['modeled_by'].extend(driver_info['additional_allele_curies'])
                driver_info['eco_abbr'] = 'CEC'    # The default.
            else:
                if driver_info['qualifier'] == 'DOES NOT model':
                    driver_info['is_not'] = True
                driver_info['modeled_by'].append(driver_info['allele_feature_curie'])
                driver_info['modeled_by'].extend(driver_info['additional_allele_curies'])
                driver_info['eco_abbr'] = driver_info['evi_code'][0:3]

            # Build an annotation unique key.
            if driver_info['problems']:
                prob_counter += 1
                continue
            driver_info['unique_key'] = f'{driver_info["pub_given"]}_'
            if driver_info['is_not']:
                driver_info['unique_key'] += 'NOT_'
            driver_info['unique_key'] += f'model={"|".join(sorted(driver_info["modeled_by"]))}_'
            driver_info['unique_key'] += f'disease_term={driver_info["doid_term_curie"]}_'
            driver_info['unique_key'] += f'eco_code={driver_info["eco_abbr"]}'
            if driver_info['modifier_curie']:
                driver_info['unique_key'] += f'_{driver_info["modifier_role"]}={driver_info["modifier_curie"]}'
            self.log.debug(f'Line={line_number}; ukey={driver_info["unique_key"]}')

            # Look for matches between input file annotations and chado annotations.
            if driver_info['unique_key'] in self.genotype_dis_annos.keys():
                matched_dis_anno_counter += 1
                if driver_info['driver_curies']:
                    self.driver_dict[driver_info['unique_key']].append(driver_info)
            else:
                alt_unique_key = driver_info['unique_key'].replace('eco_code=CEC', 'eco_code=CEA')
                if alt_unique_key in self.genotype_dis_annos.keys():
                    close_matched_dis_anno_counter += 1
                    if driver_info['driver_curies']:
                        driver_info['unique_key'] = alt_unique_key
                        self.driver_dict[alt_unique_key].append(driver_info)
                else:
                    line_number = driver_info['line_number']
                    # self.log.warning(f'Could not find dis anno: line={line_number}; unique_key={driver_info["unique_key"]}; line={line}; dict={driver_info}')
                    prob_msg = 'no matching dis anno'
                    driver_info['problems'].append(prob_msg)
                    self.rejected_driver_info.append(driver_info)
                    unmatched_dis_anno_counter += 1

        # Summary
        self.log.info(f'Skipped {malformed_line_counter}/{input_counter} driver info lines due to bad column formatting.')
        self.log.info(f'Processed {input_counter} driver info lines having Gal4 info.')
        fully_processed_count = input_counter - prob_counter
        self.log.info(f'Had problems finding pub/allele/term info for {prob_counter}/{input_counter} driver info lines.')
        self.log.info(f'Fully processed {fully_processed_count}/{input_counter} driver info lines having Gal4 info without issue.')
        self.log.info(f'Found dis anno for {matched_dis_anno_counter}/{fully_processed_count} fully processed driver info lines.')
        self.log.info(f'Found close dis anno for {close_matched_dis_anno_counter}/{fully_processed_count} fully processed driver info lines (ECO adjustment).')
        self.log.info(f'Could not find dis anno for {unmatched_dis_anno_counter}/{fully_processed_count} fully processed driver info lines.')
        return

    def report_unmatched_driver_lines(self, session):
        """Print out a report for driver input missing in chado."""
        self.log.info('Print out a report for driver input missing in chado.')
        # Open up the report.
        curator_report = open('/src/output/unmatched_driver_lines.tsv', 'w')
        headers = [
            'line_number',
            'pub_given',
            'allele_symbol',
            'qualifier',
            'do_term',
            'evi_code',
            'additional_alleles',
            'driver_input',
            'operation',
            'problems',
        ]
        curator_report.write('#')
        csv_writer = csv.DictWriter(curator_report, fieldnames=headers, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        csv_writer.writeheader()
        # As info is printed, we'll also confirm, independently, that the driver file input is not in chado.
        confirmed_counter = 0
        discrepancy_counter = 0
        for i in self.rejected_driver_info:
            i['problems'] = '; '.join(i['problems'])
            exists = self.check_disease_annotation(session, i['allele_feature_curie'], i['do_term'],
                                                   i['pub_given'], i['qualifier'], i['evi_code'])
            if exists:
                desc = f'line={i["line_number"]}; unique_key={i["unique_key"]}'
                self.log.error(f'This driver file line was not matched to an annotation, but appears to exist in chado: {desc}')
                discrepancy_counter += 1
            else:
                csv_writer.writerow(i)
                confirmed_counter += 1
        self.log.info(f'Confirmed that {confirmed_counter} unmatched disease annotations are not in chado.')
        self.log.info(f'Found {discrepancy_counter} discrepancies where an unmatched disease annotation is in chado.')
        return

    def integrate_driver_info(self):
        """Integrate driver info."""
        self.log.info('Integrate driver info.')
        counter = 0
        miscounter = 0
        for uniq_key, driver_info_list in self.driver_dict.items():
            # self.log.debug(f'Process {uniq_key} driver_info: {driver_info_list}')
            driver_combos = set()
            for driver_info in driver_info_list:
                # For the "and" operation, we keep the set intact - integrate the combination.
                if driver_info['operation'] == 'and':
                    driver_combo_str = '_'.join(sorted(driver_info['driver_curies']))
                    driver_combos.add(driver_combo_str)
                # For "or" and "na" operations, we integrate each driver separately.
                else:
                    for driver_curie in driver_info['driver_curies']:
                        driver_combos.add(driver_curie)
            # self.log.debug(f'Have this final set of driver combos: {driver_combos}')
            try:
                self.genotype_dis_annos[uniq_key].driver_combos = driver_combos
                counter += 1
            except KeyError:
                miscounter += 1
        self.log.info(f'Integrated driver info into {counter} genotype-level disease annotations.')
        self.log.info(f'For {miscounter} driver info objects, could not match up driver info to a genotype-level disease annotation.')
        return

    def split_out_genotype_disease_annotations(self):
        """Split out genotype annotations for each driver combo."""
        self.log.info('Split out genotype annotations for each driver combo.')
        input_counter = 0
        output_counter = 0
        for dis_anno in self.genotype_dis_annos.values():
            input_counter += 1
            if not dis_anno.driver_combos:
                dis_anno.unique_key += '_driver_curies=None'
                self.fb_data_entities[dis_anno.unique_key] = dis_anno
                output_counter += 1
            else:
                for driver_combo in dis_anno.driver_combos:
                    new_unique_key = f'{dis_anno.unique_key}_driver_curies={driver_combo}'
                    new_dis_anno = fb_datatypes.FBGenotypeDiseaseAnnotation(new_unique_key)
                    new_dis_anno.allele_annotations = dis_anno.allele_annotations
                    new_dis_anno.driver_combos = {driver_combo}
                    self.fb_data_entities[new_unique_key] = new_dis_anno
                    output_counter += 1
        self.log.info(f'Turned {input_counter} initial genotype-level disease annotations into {output_counter} genotype disease annotations.')
        self.log.info(f'Have {len(self.fb_data_entities.keys())} split out annotations.')
        return

    def propagate_allele_to_genotype_attributes(self):
        """Propagate allele annotation to genotype annotation."""
        self.log.info('Propagate allele annotation to genotype annotation.')
        counter = 0
        for dis_anno in self.fb_data_entities.values():
            if dis_anno.for_export is False:
                continue
            al_dis_anno = dis_anno.allele_annotations[0]
            dis_anno.modeled_by = al_dis_anno.modeled_by
            dis_anno.pub_curie = self.lookup_single_pub_curie(al_dis_anno.feature_cvterm.pub_id)
            dis_anno.do_term_curie = self.cvterm_lookup[al_dis_anno.feature_cvterm.cvterm_id]['curie']
            dis_anno.is_not = al_dis_anno.is_not
            dis_anno.eco_abbr = al_dis_anno.eco_abbr
            dis_anno.modifier_curie = al_dis_anno.modifier_curie
            dis_anno.modifier_role = al_dis_anno.modifier_role
            counter += 1
        self.log.info(f'Extracted key info from {counter} genotype-level annotations in chado.')
        return

    def parse_aberration_info(self, session):
        """Parse new aberration disease annotations."""
        self.log.info('Parse new aberration disease annotations.')
        DF_ACROSS_ALLELE = 0
        BANG_C = 1
        PUB_GIVEN = 2
        SUBJECT = 3
        QUAL = 4
        DO_TERM = 5
        EVI_CODE = 6
        OBJECTS = 7
        DRIVER_INPUT = 8
        DF_SBJ_GENE = 9
        DF_OBJ_GENE = 10
        GENOTYPE = 11
        file_input = open('/src/output/aberr_info.tsv')
        line_number = 0
        input_counter = 0
        malformed_line_counter = 0
        prob_counter = 0
        self.rejected_aberr_info = []
        for i in file_input:
            line_number += 1
            if 'FBrf0' not in i:
                continue
            self.log.debug(f'Process this aberr line: {i.strip()}')
            input_counter += 1
            line = i.split('\t')
            if len(line) != 12:
                self.log.error(f'Line={line_number} is malformed: {i}')
                malformed_line_counter += 1
                continue
            aberr_info = {
                # Attributes from input file.
                'line_number': line_number,
                'df_across_allele': line[DF_ACROSS_ALLELE].strip(),
                'bang_c': line[BANG_C].strip(),
                'pub_given': line[PUB_GIVEN].strip(),
                'subject': sgml_to_plain_text(line[SUBJECT]).strip(),
                'qualifier': line[QUAL].strip(),
                'do_term': line[DO_TERM].strip(),
                'evi_code': line[EVI_CODE].strip(),
                'objects': line[OBJECTS].split(' '),
                'driver_input': line[DRIVER_INPUT].split(' '),
                'df_sbj_genes': line[DF_SBJ_GENE].strip().split(' '),
                'df_obj_genes': line[DF_OBJ_GENE].strip().split(' '),
                'genotype': line[GENOTYPE].strip(),
                # Attributes to be obtained from chado.
                'pub_id': None,
                'doid_term_curie': None,
                'subject_curie': None,
                'object_curies': [],
                'driver_curies': [],
                'driver_combo_str': '',
                'asserted_gene_feature_ids': [],
                # Attributes synthesized from the above.
                'modeled_by': [],
                'is_not': False,
                'modifier_curie': None,
                'modifier_role': None,
                'eco_abbr': None,
                'unique_key': None,
                'problems': [],
            }

            if aberr_info['df_across_allele'] not in ['no', 'n/a', 'trans']:
                self.log.error(f'Line={line_number}: unexpected "DF_ACROSS_ALLELE" value: "{aberr_info["df_across_allele"]}".')
                prob_msg = 'bad df_across_allele value'
                aberr_info['problems'].append(prob_msg)

            if not aberr_info['evi_code']:
                self.log.error(f'Line={line_number}: no evidence code given.')
                prob_msg = 'no evidence code given'
                aberr_info['problems'].append(prob_msg)

            if not aberr_info['qualifier']:
                self.log.error(f'Line={line_number}: no qualifier given.')
                prob_msg = 'no qualifier given'
                aberr_info['problems'].append(prob_msg)

            try:
                aberr_info['pub_id'] = self.fbrf_bibliography[aberr_info['pub_given']].pub_id
            except KeyError:
                self.log.error(f'Line={line_number}: could not find pub "{aberr_info["pub_given"]}" in chado.')
                prob_msg = 'bad pub id'
                aberr_info['problems'].append(prob_msg)

            cvterm_curie = self.find_cvterm_curie_from_name(session, aberr_info['do_term'])
            if cvterm_curie:
                aberr_info['doid_term_curie'] = cvterm_curie
                # self.log.debug(f'Line={line_number}: found "{cvterm_curie}" for "{aberr_info["do_term"]}" in chado.')
            else:
                self.log.error(f'Line={line_number}: could not find DO term for "{aberr_info["do_term"]}" in chado.')
                prob_msg = 'bad DO term name'
                aberr_info['problems'].append(prob_msg)

            if aberr_info['qualifier'] not in ['ameliorates', 'exacerbates', 'model of', 'DOES NOT model']:
                self.log.error(f'Line={line_number}: bad qualifier "{aberr_info["qualifier"]}".')
                prob_msg = 'bad qualifier'
                aberr_info['problems'].append(prob_msg)

            subject_curie = self.find_feature_uniquename_from_name(session, aberr_info['subject'], self.regex['allele'])
            if subject_curie is None:
                subject_curie = self.find_feature_uniquename_from_name(session, aberr_info['subject'], self.regex['aberration'])
            if subject_curie:
                aberr_info['subject_curie'] = subject_curie
                # self.log.debug(f'Line={line_number}: found "{allele_id}" for "{aberr_info["allele_symbol"]}" in chado.')
            else:
                self.log.error(f'Line={line_number}: could not find subject "{aberr_info["subject"]}" in chado.')
                prob_msg = 'bad subject symbol'
                aberr_info['problems'].append(prob_msg)

            for feature_symbol in aberr_info['objects']:
                if feature_symbol == '' or feature_symbol == ' ' or feature_symbol == '+':
                    continue
                feature_id = self.find_feature_uniquename_from_name(session, feature_symbol, self.regex['allele'])
                if feature_id is None:
                    feature_id = self.find_feature_uniquename_from_name(session, feature_symbol, self.regex['aberration'])
                if feature_id:
                    aberr_info['object_curies'].append(feature_id)
                    # self.log.debug(f'Line={line_number}: found "{feature_id}" for additional allele "{feature_symbol}" in chado.')
                else:
                    self.log.error(f'Line={line_number}: could not find object feature "{feature_symbol}" in chado.')
                    prob_msg = 'bad object feature symbol'
                    aberr_info['problems'].append(prob_msg)

            for driver_symbol in aberr_info['driver_input']:
                if driver_symbol == '' or driver_symbol == ' ' or driver_symbol == '+':
                    continue
                driver_curie = self.find_feature_uniquename_from_name(session, driver_symbol, self.regex['allele'])
                if driver_curie:
                    aberr_info['driver_curies'].append(driver_curie)
                    # self.log.debug(f'Line={line_number}: found "{driver_curie}" for driver "{driver_symbol}" in chado.')
                else:
                    self.log.error(f'Line={line_number}: could not find driver "{driver_symbol}" in chado.')
                    prob_msg = 'bad driver symbol'
                    aberr_info['problems'].append(prob_msg)

            asserted_gene_symbols = []
            asserted_gene_symbols.extend(aberr_info['df_obj_genes'])
            # Only consider sbj genes if a "model of" annotation; no asserted genes for modifiers (in sbj column).
            if aberr_info['qualifier'] == 'model of':
                asserted_gene_symbols.extend(aberr_info['df_sbj_genes'])
            asserted_gene_symbols = set(asserted_gene_symbols)
            for gene_symbol in asserted_gene_symbols:
                if gene_symbol == '' or gene_symbol == '***':
                    continue
                gene_curie = self.find_feature_uniquename_from_name(session, gene_symbol, self.regex['gene'])
                if gene_curie:
                    # self.log.debug(f'Line={line_number}: found gene curie {gene_curie} for "{gene_symbol}"')
                    gene_feature_id = self.uname_feature_lookup[gene_curie]['feature_id']
                    aberr_info['asserted_gene_feature_ids'].append(gene_feature_id)
                else:
                    self.log.error(f'Line={line_number}: could not find gene "{gene_symbol}" in chado.')
                    prob_msg = 'bad gene symbol'
                    aberr_info['problems'].append(prob_msg)

            # Map info to annotation.
            if aberr_info['qualifier'] in self.disease_genetic_modifier_terms.keys():
                aberr_info['modifier_curie'] = aberr_info['subject_curie']
                aberr_info['modifier_role'] = self.disease_genetic_modifier_terms[aberr_info['qualifier']]
                aberr_info['modeled_by'].extend(aberr_info['object_curies'])
                aberr_info['eco_abbr'] = 'CEC'    # The default.
            else:
                if aberr_info['qualifier'] == 'DOES NOT model':
                    aberr_info['is_not'] = True
                aberr_info['modeled_by'].append(aberr_info['subject_curie'])
                aberr_info['modeled_by'].extend(aberr_info['object_curies'])
                aberr_info['eco_abbr'] = aberr_info['evi_code'][0:3]

            # Build an annotation unique key.
            if aberr_info['problems']:
                prob_counter += 1
                continue
            aberr_info['unique_key'] = f'{aberr_info["pub_given"]}_'
            if aberr_info['is_not']:
                aberr_info['unique_key'] += 'NOT_'
            aberr_info['unique_key'] += f'model={"|".join(sorted(aberr_info["modeled_by"]))}_'
            aberr_info['unique_key'] += f'disease_term={aberr_info["doid_term_curie"]}_'
            aberr_info['unique_key'] += f'eco_code={aberr_info["eco_abbr"]}'
            if aberr_info['modifier_curie']:
                aberr_info['unique_key'] += f'_{aberr_info["modifier_role"]}={aberr_info["modifier_curie"]}'
            if not aberr_info['driver_curies']:
                aberr_info['unique_key'] += '_driver_curies=None'
            else:
                aberr_info['driver_combo_str'] = '_'.join(sorted(aberr_info['driver_curies']))
                aberr_info['unique_key'] += f'_driver_curies={aberr_info["driver_combo_str"]}'

            # Create a new annotation.
            new_dis_anno = fb_datatypes.FBGenotypeDiseaseAnnotation(aberr_info['unique_key'])
            new_dis_anno.modeled_by = aberr_info['modeled_by']
            if aberr_info['driver_combo_str']:
                new_dis_anno.driver_combos = {aberr_info['driver_combo_str']}
            new_dis_anno.is_not = aberr_info['is_not']
            if aberr_info['df_across_allele'] == 'trans':
                new_dis_anno.aberr_trans = True
            new_dis_anno.pub_curie = self.lookup_single_pub_curie(aberr_info['pub_id'])
            new_dis_anno.do_term_curie = aberr_info['doid_term_curie']
            new_dis_anno.eco_abbr = aberr_info['eco_abbr']
            new_dis_anno.modifier_curie = aberr_info['modifier_curie']
            new_dis_anno.modifier_role = aberr_info['modifier_role']
            new_dis_anno.asserted_genes = aberr_info['asserted_gene_feature_ids']
            self.fb_data_entities[new_dis_anno.unique_key] = new_dis_anno
            self.log.debug(f'Line={line_number}; ukey={aberr_info["unique_key"]}')

        # Summary
        self.log.info(f'Skipped {malformed_line_counter}/{input_counter} aberration info lines due to bad column formatting.')
        self.log.info(f'Processed {input_counter} aberration info lines having aberration info.')
        fully_processed_count = input_counter - prob_counter
        self.log.info(f'Had problems finding pub/feature/term info for {prob_counter}/{input_counter} aberration info lines.')
        self.log.info(f'Fully processed {fully_processed_count}/{input_counter} aberration info lines having aberration info without issue.')
        return

    def derive_genotypes(self):
        """Derive genotypes for final genotype-level disease annotations."""
        self.log.info('Derive genotypes for final genotype-level disease annotations.')
        counter = 0
        for dis_anno in self.fb_data_entities.values():
            # self.log.debug(f'Derive genotype for {dis_anno}')
            # Determine if input is allele-level, or, aberration-from-spreadsheet, annotation.
            aberr_anno = True
            if dis_anno.allele_annotations:
                aberr_anno = False
            # Get all components.
            components = []
            components.extend(dis_anno.modeled_by)
            if dis_anno.driver_combos:
                driver_curies = list(dis_anno.driver_combos)[0].split('_')
                components.extend(driver_curies)
            # Flag a prespecified pair (up to one for two components that make up a model).
            prespecified_pair = []
            if aberr_anno and dis_anno.aberr_trans and len(dis_anno.modeled_by) == 2:
                prespecified_pair.extend(dis_anno.modeled_by)
                prespecified_pair.sort()
            elif aberr_anno and dis_anno.aberr_trans:
                self.log.error(f'For {dis_anno}, have "trans" specified for a model with {len(dis_anno)} components.')
            # Sort things into complementation groups.
            cgroup_dict = {}
            # Start with the prespecified pair, if present.
            if prespecified_pair:
                cgroup_dict['_'.join(prespecified_pair)] = [self.uname_feature_lookup[i] for i in prespecified_pair]
            # Add other components.
            for curie in components:
                if curie in prespecified_pair:
                    continue
                if curie == '':
                    self.log.error(f'Found empty component curie for {dis_anno}.')
                    continue
                feature = self.uname_feature_lookup[curie]
                single_cgroup = True
                # Flag transgenic alleles.
                if feature['feature_id'] in self.transgenic_allele_ids:
                    single_cgroup = False
                elif feature['feature_id'] in self.in_vitro_allele_ids:
                    single_cgroup = False
                elif feature['feature_id'] in self.misxprn_allele_ids:
                    single_cgroup = False
                # Sort transgenic alleles into their own cgroup.
                if single_cgroup is False:
                    cgroup_dict[curie] = [feature]
                # Group classical alleles by gene.
                elif curie.startswith('FBal'):
                    gene_feature_id = self.allele_gene_lookup[feature['feature_id']]
                    gene = self.feature_lookup[gene_feature_id]
                    try:
                        cgroup_dict[gene['curie']].append(feature)
                    except KeyError:
                        cgroup_dict[gene['curie']] = [feature]
                # Group aberrations by relevant gene as well.
                elif curie.startswith('FBab'):
                    cgroup_dict[curie] = [feature]
                else:
                    self.log.error(f'For {dis_anno}, could not assign {curie} to genotype cgroup.')
            cgroup_names = []
            for cgroup in cgroup_dict.values():
                cgroup_name = '/'.join([i['name'] for i in cgroup])
                cgroup_names.append(cgroup_name)
            dis_anno.genotype_name = ' '.join(cgroup_names)
            self.log.debug(f'For {dis_anno}:\n\tgenotype name: {dis_anno.genotype_name}')
            counter += 1
        self.log.info(f'Derived genotype names for {counter} disease annotations.')
        return

    def get_genotypes(self, session):
        """Get genotypes for final genotype-level disease annotations."""
        self.log.info('Get genotypes for final genotype-level disease annotations.')
        counter = 0
        prob_counter = 0
        no_counter = 0
        for dis_anno in self.fb_data_entities.values():
            if not dis_anno.genotype_name:
                self.log.error(f'No genotype_name for {dis_anno}')
                prob_counter += 1
                continue
            genotype = GenotypeAnnotation(dis_anno.genotype_name, session, self.log)
            genotype.get_known_or_create_new_genotype(session)
            self.log.debug(f'Got this curie: {genotype.curie}')
            dis_anno.genotype_curie = genotype.curie
            if genotype.curie is None:
                no_counter += 1
            else:
                counter += 1
        self.log.info(f'Skipped {prob_counter} disease annotations lacking a specified genotype name.')
        self.log.info(f'Got/created {counter} disease genotypes.')
        self.log.info(f'Could not get/create a genotype in {no_counter} cases.')
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
        self.get_parent_genes()
        # self.find_relevant_aberrations()    # Slow, run only upon request (Jira FTA-45)
        self.build_model_eco_lookup()
        self.lookup_eco_codes_for_modifier_annotations()
        self.group_redundant_annotations()
        self.parse_driver_info(session)
        self.report_unmatched_driver_lines(session)
        self.integrate_driver_info()
        self.split_out_genotype_disease_annotations()
        self.propagate_allele_to_genotype_attributes()
        self.parse_aberration_info(session)
        self.derive_genotypes()
        self.get_genotypes(session)
        return

    # Add methods to be run by synthesize_info() below.

    # Elaborate on synthesize_info() for the AGMDiseaseHandler.
    def synthesize_info(self):
        """Extend the method for the AGMDiseaseHandler."""
        super().synthesize_info()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_genotype_disease_annotation_basic(self):
        """Map basic FlyBase genotype disease annotation to the Alliance LinkML object."""
        self.log.info('Map basic FlyBase genotype disease annotation to the Alliance LinkML object.')
        for geno_dis_anno in self.fb_data_entities.values():
            if geno_dis_anno.for_export is False:
                continue
            agr_dis_anno = self.agr_export_type(f'FB:{geno_dis_anno.genotype_curie}', geno_dis_anno.do_term_curie, geno_dis_anno.pub_curie)
            if geno_dis_anno.is_not is True:
                agr_dis_anno.negated = True
            agr_dis_anno.evidence_code_curies.append(self.evidence_code_xrefs[geno_dis_anno.eco_abbr])
            if geno_dis_anno.modifier_curie:
                agr_dis_anno.disease_genetic_modifier_relation_name = geno_dis_anno.modifier_role
                agr_dis_anno.disease_genetic_modifier_identifiers = [f'FB:{geno_dis_anno.modifier_curie}']
            geno_dis_anno.linkmldto = agr_dis_anno
        return

    def add_asserted_genes_alleles(self):
        """Add asserted genes and alleles."""
        self.log.info('Add asserted genes and alleles.')
        for dis_anno in self.fb_data_entities.values():
            if dis_anno.for_export is False:
                continue
            # Determine asserted alleles.
            asserted_fbal_ids = []
            for fbal_id in dis_anno.modeled_by:
                asserted_fbal_ids.append(fbal_id)
            dis_anno.linkmldto.asserted_allele_identifier = f'FB:{asserted_fbal_ids[0]}'                # BOB - temp, want many.
            # dis_anno.linkmldto.asserted_allele_identifiers = [f'FB:{i}' for i in asserted_fbal_ids]    # BOB - what we want eventually.
            # Determine asserted genes.
            asserted_gene_feature_ids = []
            asserted_gene_feature_ids.extend(dis_anno.asserted_gene_ids)
            for fbal_id in asserted_fbal_ids:
                if not fbal_id.startswith('FBal'):
                    continue
                allele_feature_id = self.uname_feature_lookup[fbal_id]['feature_id']
                parent_gene_feature_id = self.feature_lookup[self.allele_gene_lookup[allele_feature_id]]['feature_id']
                asserted_gene_feature_ids.append(parent_gene_feature_id)
            asserted_gene_feature_ids = set(asserted_gene_feature_ids)
            for gene_feature_id in asserted_gene_feature_ids:
                gene = self.feature_lookup[gene_feature_id]
                gene_organism = self.organism_lookup[gene['organism_id']]
                if gene_organism['is_drosophilid'] is True:
                    dis_anno.linkmldto.asserted_gene_identifiers.append(gene['curie'])
                elif not gene['curie'].startswith('FB') and gene_organism['official_db']:
                    dis_anno.linkmldto.asserted_gene_identifiers.append(gene['curie'])
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

    # Elaborate on map_fb_data_to_alliance() for the AGMDiseaseHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AGMDiseaseHandler."""
        super().map_fb_data_to_alliance()
        self.map_genotype_disease_annotation_basic()
        self.add_asserted_genes_alleles()
        self.map_data_provider_dto()
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
        current_modifier_curie_counter = 0
        updated_modifier_curie_counter = 0
        cannot_update_modifier_curie_counter = 0
        modifier_prob_counter = 0
        for dis_anno in self.fb_data_entities.values():
            for fb_term in self.disease_genetic_modifier_terms.keys():
                if fb_term in dis_anno.evidence_code.value:
                    dis_anno.fb_modifier_type = fb_term
                    if re.search(allele_regex, dis_anno.evidence_code.value):
                        allele_id = re.search(allele_regex, dis_anno.evidence_code.value).group(0)
                    if self.confirm_current_allele_by_uniquename(session, allele_id):
                        dis_anno.fb_modifier_id = allele_id
                        current_modifier_curie_counter += 1
                    else:
                        # Look up current allele by 2o ID. Use that.
                        curr_allele_id = self.get_current_id_for_allele(session, allele_id)
                        if curr_allele_id:
                            dis_anno.fb_modifier_id = ['FB:{}'.format(curr_allele_id)]
                            dis_anno.modifier_curie_was_updated = True
                            updated_modifier_curie_counter += 1
                        else:
                            dis_anno.modifier_problem = True
                            cannot_update_modifier_curie_counter += 1
            if dis_anno.modifier_problem is True:
                modifier_prob_counter += 1
        self.log.info(f'{current_modifier_curie_counter} allele modifiers mentioned use a current allele ID.')
        self.log.info(f'{updated_modifier_curie_counter} allele modifiers mentioned use a non-current allele ID mapped to a current allele.')
        self.log.info(f'{cannot_update_modifier_curie_counter} allele modifiers mentioned use a non-current allele ID not mappable to a current allele.')
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
                dis_anno.modifier_curie_was_updated is True: 'Obsolete modifier ID',
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
