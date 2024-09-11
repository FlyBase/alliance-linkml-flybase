"""Module:: disease_handlers.

Synopsis:
    A data handlers for FlyBase disease annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import re
from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.production import (
    Cv, Cvterm, Db, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureDbxref, FeatureRelationship, Pub
)
import agr_datatypes
from handler import DataHandler


class AlleleDiseaseHandler(DataHandler):
    """A data handler for allele-based disease annotations."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the AlleleDiseaseHandler object."""
        super().__init__(log, fb_data_type, testing)

        self.uniq_dis_dict = {}    # A dict of unique disease annotation attributes to list of disease annotations.

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
            dis_anno = self.datatype_objects[self.fb_data_type](result.FeatureCvterm, result.FeatureCvtermprop)
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
            pub_curie = self.lookup_pub_curie(dis_anno.feature_cvterm.pub_id)
            agr_dis_anno = agr_datatypes.AlleleDiseaseAnnotationDTO(allele_curie, do_curie, pub_curie)
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
                        self.fb_data_entities[db_primary_id].export_warnings.append('Annotation is redundant.')
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
