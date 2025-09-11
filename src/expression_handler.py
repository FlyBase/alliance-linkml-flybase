"""Module:: expression_handler.

Synopsis:
    Data handler for FlyBase expression annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import re
from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.char_conversions import clean_free_text
from harvdev_utils.reporting import (
    Cv, Cvterm, Db, Dbxref, Expression, ExpressionCvterm, ExpressionCvtermprop,
    Feature, FeatureExpression, FeatureExpressionprop, FeatureRelationship
)
import fb_datatypes
from handler import DataHandler


class ExpressionHandler(DataHandler):
    """A data handler for gene expression data."""
    def __init__(self, log: Logger, testing: bool):
        """Create the ExpressionHandler object."""
        super().__init__(log, testing)
        self.datatype = 'feature_expression'
        self.agr_export_type = None
        self.slot_types = ['anatomy', 'assay', 'cellular', 'stage']
        self.placeholder = fb_datatypes.FBExpressionCvterm(None)
        self.expression_patterns = {}             # expression_id-keyed FBExpressionAnnotation objects.
        self.export_data_for_tsv = []             # List of dicts for export to TSV.
        self.isoform_gene_product_lookup = {}     # Will be feature_id-keyed feature_id that connects an isoform to an XR/XP gene product.
        self.gene_product_gene_lookup = {}        # Will be feature_id-keyed feature_id that connects an XR/XP gene product to a gene.
        self.allele_product_allele_lookup = {}    # Will be feature_id-keyed feature_id that connects an RA\PA allele product to an allele.
        self.insertion_allele_lookup = {}         # Will be FBti ID keyed lists of related allele FBal IDs (list).
        self.hemi_drivers = []                    # Will be a list of feature_ids for hemidriver alleles that have FBco parents.
        self.split_system_combos = {}             # Will be FBco ID-keyed list of FBal IDs for hemidriver components of each split system combination.
        self.split_system_combo_strs = []         # Will be a list of concatenated hemidriver FBal IDs for each split system combination feature.

    # Key info.

    # Mappings between FB stage slim terms and UBERON: embryonic (UBERON:0000068), post-embryonic, adult (UBERON:0000113).
    fb_uberon_stage_slim_map = {
        'fertilized egg stage': ['UBERON:0000068'],
        'embryonic stage': ['UBERON:0000068'],
        'larval stage': ['post embryonic, pre-adult'],
        'P-stage': ['post embryonic, pre-adult'],
        'adult stage': ['UBERON:0000113'],
        'mature adult stage': ['UBERON:0000113'],
        'oogenesis': ['UBERON:0000113'],
        'unfertilized egg stage': ['UBERON:0000113'],
        'spermatogenesis': ['UBERON:0000113'],
    }

    # Mappings between FB anatomy slim terms and UBERON.
    fb_uberon_anatomy_slim_map = {
        'adipose system': ['UBERON:0001013', 'UBERON:0002423'],
        'presumptive embryonic/larval adipose system': ['UBERON:0001013', 'UBERON:0002423'],
        'oenocyte': ['UBERON:0002423'],
        'appendage': ['UBERON:0000026'],
        'auditory system': ['UBERON:0002105'],
        'chemosensory system': ['UBERON:0005726'],
        'circulatory system': ['UBERON:0001009'],
        'presumptive embryonic/larval circulatory system': ['UBERON:0001009'],
        'digestive system': ['UBERON:0005409'],
        'presumptive embryonic/larval digestive system': ['UBERON:0005409'],
        'ectoderm': ['UBERON:0000924'],
        'endocrine system': ['UBERON:0000949'],
        'presumptive embryonic/larval endocrine system': ['UBERON:0000949'],
        'endoderm': ['UBERON:0000925'],
        'excretory system': ['UBERON:0001008'],
        'extraembryonic structure': ['UBERON:0016887'],
        'embryonic/larval lymph gland': ['UBERON:0002193', 'UBERON:0001009'],
        'hemolymph': ['UBERON:0002193', 'UBERON:0001009'],
        'imaginal tissue': ['UBERON:6005023'],
        'integumentary system': ['UBERON:0002416'],
        'presumptive embryonic/larval integumentary system': ['UBERON:0002416'],
        'mechanosensory system': ['UBERON:0007037'],
        'mesoderm': ['UBERON:0000926'],
        'muscle system': ['UBERON:0002204'],
        'presumptive embryonic/larval muscle system': ['UBERON:0002204'],
        'presumptive embryonic/larval nervous system': ['UBERON:0001016'],
        'nervous system': ['UBERON:0001016'],
        'reproductive system': ['UBERON:0000990'],
        'tracheal system': ['UBERON:0001004'],
        'presumptive embryonic/larval tracheal system': ['UBERON:0001004'],
        'sensory system': ['UBERON:0001032'],
        'embryonic sensory system': ['UBERON:0001032'],
        'visual system': ['UBERON:0002104'],
    }

    # Utility functions.
    def regex_for_anatomical_terms_in_numerical_series(self, start_term, end_term, xprn_id):
        """For two anatomical terms in a series, return a regex for all terms in that series, and the start and stop positions."""
        match = None
        start_position = None
        end_position = None
        # First get the positional part of the term.
        # Case 1. Exceptional intersegmental region terms: e.g., "embryonic/larval abdominal 1-2 intersegmental apodeme".
        intersegmental_positions = ['1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9']
        if 'intersegmental' in start_term:
            rgx = r' ([0-9]-[0-9])'
            start_match = re.search(rgx, start_term)
            end_match = re.search(rgx, end_term)
            if start_match and end_match and start_match != end_match and start_match.group(1) in intersegmental_positions:
                match = start_match
                start_position = start_match.group(1)
                end_position = end_match.group(1)
            # Filter out cases that do not meet all criteria.
            else:
                match = None
                start_position = None
                end_position = None
        # Case 2. Look for numbers at the end of the terms and confirm that they differ.
        else:
            terminal_rgx = r' ([A-Z]{0,1}[0-9]{1,2})$'
            start_terminal_match = re.search(terminal_rgx, start_term)
            end_terminal_match = re.search(terminal_rgx, end_term)
            if start_terminal_match and end_terminal_match and start_terminal_match.group(1) != end_terminal_match.group(1):
                match = start_terminal_match
                start_position = start_terminal_match.group(1)
                end_position = end_terminal_match.group(1)
            # Case 3. Look for internal numbers.
            else:
                internal_rgx = r' ([A-Z]{0,1}[0-9]{1,2})'
                start_internal_match = re.search(internal_rgx, start_term)
                end_internal_match = re.search(internal_rgx, end_term)
                if start_internal_match and end_internal_match and start_internal_match.group(1) != end_internal_match.group(1):
                    match = start_internal_match
                    start_position = start_internal_match.group(1)
                    end_position = end_internal_match.group(1)
                else:
                    match = None
                    start_position = None
                    end_position = None
        if not match:
            self.log.error(f'For xprn_id={xprn_id}, no number or letter+number pattern found for these terms: {start_term}, {end_term}')
            return None, start_position, end_position
        # Build the regex.
        before = re.escape(start_term[:match.start()])
        after = re.escape(start_term[match.end():])
        position = match.group(0)
        letter_part = re.match(r' ([A-Za-z]{0,2})', position).group(0)
        if 'intersegmental' in start_term:
            regex_number_part = r'\d-\d'
        else:
            regex_number_part = r'\d+'
        position_regex = letter_part + regex_number_part
        regex = f'^{before}({position_regex}){after}$'
        # self.log.debug(f'For {start_term} and {end_term}, found this range: {start_position}--{end_position}')
        # self.log.debug(f'For {start_term} and {end_term}, found this regex: {regex}')
        return regex, start_position, end_position

    def get_anatomical_terms_by_regex(self, session, term_regex):
        """Get anatomical terms matching a given regex."""
        # self.log.debug(f'Find terms matching this regex: {term_regex}')
        filters = (
            Cvterm.is_obsolete == 0,
            Cvterm.name.op('~')(term_regex),
            Cv.name == 'FlyBase anatomy CV',
        )
        terms = session.query(Cvterm).\
            select_from(Cvterm).\
            join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
            filter(*filters).\
            distinct()
        term_names = [i.name for i in terms]
        if term_names:
            # self.log.debug(f'For {term_regex}, found these matching FBbt terms: {term_names}')
            pass
        else:
            self.log.error(f'For {term_regex}, found NO matching FBbt terms.')
        return terms

    def select_in_range_anatomical_terms(self, terms, rgx, start, end):
        """Filter for anatomical terms within a certain range."""
        filtered_terms = []
        # For intersegmental terms, this approach just uses the first number in the \d-\d pattern for range comparison.
        num_rgx = r'[A-ZAa-z]{0,2}(\d+)'
        try:
            num_start = int(re.search(num_rgx, start).group(1))
            num_end = int(re.search(num_rgx, end).group(1))
        except TypeError:
            self.log.error('Could not find start/end numbers in input terms.')
            return filtered_terms
        # self.log.debug(f'From {start}--{end}, look for numbers between {num_start} and {num_end}.')
        for term in terms:
            # self.log.debug(f'Is "{term.name}" in range?')
            position = re.search(rgx, term.name).group(1)
            # self.log.debug(f'Found this position "{position}".')
            num_position = int(re.search(num_rgx, position).group(1))
            # self.log.debug(f'Term {term.name} is at position {position}, with number={num_position}')
            if num_position > num_start and num_position < num_end:
                filtered_terms.append(term)
        # self.log.debug(f'Found these tissue range terms: {[i.name for i in filtered_terms]}')
        return filtered_terms

    # Add methods to be run by get_general_data() below.
    def add_placeholder_cvterm(self):
        """Add an empty placeholder CV term to the CV term lookup."""
        self.log.info('Add an empty placeholder CV term to the CV term lookup.')
        # This placeholder will permit the generation of all possible term combinations when a given slot type is empty.
        placeholder_cvterm_dict = {
            'cvterm_id': '',
            'name': '',
            'cv_name': '',
            'db_name': '',
            'curie': '',
            'name_plus_curie': '',
            'slim_term_cvterm_ids': [],
        }
        self.cvterm_lookup['placeholder'] = placeholder_cvterm_dict
        return

    def get_slim_term_mappings(self, session):
        """Add anatomy and stage slim term mappings to the CV term lookup."""
        self.log.info('Add anatomy and stage slim term mappings to the CV term lookup.')
        slim_term_sets = {
            'FlyBase anatomy CV': self.fb_uberon_anatomy_slim_map.keys(),
            'FlyBase development CV': self.fb_uberon_stage_slim_map.keys(),
        }
        cv_db_correspondence = {
            'FlyBase anatomy CV': 'FBbt',
            'FlyBase development CV': 'FBdv',
        }
        for cv_name, slim_term_set in slim_term_sets.items():
            for slim_term_name in slim_term_set:
                # self.log.debug(f'Assess slim term: cv_name={cv_name}, cvterm_name={slim_term_name}')
                # First, get the slim term Cvterm from chado (need the cvterm_id).
                filters = (
                    Cvterm.is_obsolete == 0,
                    Cvterm.name == slim_term_name,
                    Cv.name == cv_name,
                    Db.name == cv_db_correspondence[cv_name],
                )
                slim_term = session.query(Cvterm).\
                    select_from(Cvterm).\
                    join(Cv, (Cv.cv_id == Cvterm.cv_id)).\
                    join(Dbxref, (Dbxref.dbxref_id == Cvterm.dbxref_id)).\
                    join(Db, (Db.db_id == Dbxref.db_id)).\
                    filter(*filters).\
                    one_or_none()
                if not slim_term:
                    self.log.error(f'Could not find FBbt term for slim term name "{slim_term_name}".')
                    continue
                else:
                    # self.log.debug(f'Slim term: cv_name={cv_name}, cvterm_name="{slim_term_name}", cvterm_id={slim_term.cvterm_id}')
                    pass
                # Second, get the child terms under each slim term.
                child_cvterm_ids = self.get_child_cvterms(session, slim_term_name, cv_name)
                for child_cvterm_id in child_cvterm_ids:
                    self.cvterm_lookup[child_cvterm_id]['slim_term_cvterm_ids'].append(slim_term.cvterm_id)
                # For evaluation of the slim-child results.
                # child_term_names = [self.cvterm_lookup[i]['name'] for i in child_cvterm_ids if i in self.cvterm_lookup.keys()]
                # child_term_names.sort()
                # child_term_name_str = '\n'.join(child_term_names)
                # self.log.debug(f'Found {len(child_term_names)} child terms for the {cv_name} slim term "{slim_term_name}":\n{child_term_name_str}')
        for cvterm in self.cvterm_lookup.values():
            if cvterm['db_name'] in ['FBbt', 'FBdv']:
                cvterm['slim_term_cvterm_ids'] = list(set(cvterm['slim_term_cvterm_ids']))
                # self.log.debug(f'For {cvterm["db_name"]} term "{cvterm["name"]}", found slim terms: {cvterm["slim_term_cvterm_ids"]}')
                # self.log.debug(f'For {cvterm["db_name"]} term "{cvterm["name"]}", found slim terms: {len(cvterm["slim_term_cvterm_ids"])}')
        return

    def get_isoform_mappings(self, session):
        """Get isoform to gene product mappings."""
        self.log.info('Get isoform to gene product mappings.')
        isoform = aliased(Feature, name='isoform')
        gene_product = aliased(Feature, name='gene_product')
        filters = (
            isoform.is_obsolete.is_(False),
            isoform.uniquename.op('~')(r'^FB(tr|pp)[0-9]{7}$'),
            gene_product.is_obsolete.is_(False),
            gene_product.uniquename.op('~')(r'^FB(tr|pp)[0-9]{7}$'),
            Cvterm.name == 'isa',
        )
        results = session.query(isoform, gene_product).\
            select_from(isoform).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == isoform.feature_id)).\
            join(gene_product, (FeatureRelationship.object_id == gene_product.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            if result.isoform.feature_id in self.isoform_gene_product_lookup.keys():
                isoform_str = f'{result.isoform.name} ({result.isoform.uniquename})'
                self.log.warning(f'Found multiple gene products for isoform {isoform_str}.')
                continue
            self.isoform_gene_product_lookup[result.isoform.feature_id] = result.gene_product.feature_id
        self.log.info(f'Found {counter} distinct isoform to gene product relationships.')
        return

    def get_gene_product_gene_mappings(self, session):
        """Get gene product to gene mappings."""
        self.log.info('Get gene product to gene mappings.')
        gene_product = aliased(Feature, name='gene_product')
        gene = aliased(Feature, name='gene')
        filters = (
            gene_product.is_obsolete.is_(False),
            gene_product.uniquename.op('~')(r'^FB(tr|pp)[0-9]{7}$'),
            gene_product.name.op('~')(r'-X(R|P)$'),
            gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(r'^FBgn[0-9]{7}$'),
            Cvterm.name == 'associated_with',
        )
        results = session.query(gene_product, gene).\
            select_from(gene_product).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == gene_product.feature_id)).\
            join(gene, (FeatureRelationship.object_id == gene.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            if result.gene_product.feature_id in self.gene_product_gene_lookup.keys():
                gene_product_str = f'{result.gene_product.name} ({result.gene_product.uniquename})'
                self.log.warning(f'Found multiple genes for gene product {gene_product_str}.')
                continue
            self.gene_product_gene_lookup[result.gene_product.feature_id] = result.gene.feature_id
        self.log.info(f'Found {counter} distinct gene_product to gene relationships.')
        return

    def get_allele_product_allele_mappings(self, session):
        """Get allele product to allele mappings."""
        self.log.info('Get allele product to allele mappings.')
        allele_product = aliased(Feature, name='allele_product')
        allele = aliased(Feature, name='allele')
        filters = (
            allele_product.is_obsolete.is_(False),
            allele_product.uniquename.op('~')(r'^FB(tr|pp)[0-9]{7}$'),
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(r'^FBal[0-9]{7}$'),
            Cvterm.name == 'associated_with',
        )
        results = session.query(allele_product, allele).\
            select_from(allele_product).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele_product.feature_id)).\
            join(allele, (FeatureRelationship.object_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            if result.allele_product.feature_id in self.allele_product_allele_lookup.keys():
                allele_product_str = f'{result.allele_product.name} ({result.allele_product.uniquename})'
                self.log.warning(f'Found multiple alleles for allele product {allele_product_str}.')
                continue
            self.allele_product_allele_lookup[result.allele_product.feature_id] = result.allele.feature_id
        self.log.info(f'Found {counter} distinct allele_product to allele relationships.')
        return

    def get_insertion_allele_mappings(self, session):
        """Get insertion to allele mappings."""
        self.log.info('Get insertion to allele mappings.')
        insertion = aliased(Feature, name='insertion')
        allele = aliased(Feature, name='allele')
        filters = (
            insertion.is_obsolete.is_(False),
            insertion.uniquename.op('~')(r'^FBti[0-9]{7}$'),
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(r'^FBal[0-9]{7}$'),
            Cvterm.name == 'associated_with',
        )
        results = session.query(insertion, allele).\
            select_from(insertion).\
            join(FeatureRelationship, (FeatureRelationship.object_id == insertion.feature_id)).\
            join(allele, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            counter += 1
            try:
                self.insertion_allele_lookup[result.insertion.uniquename].append(result.allele.uniquename)
            except KeyError:
                self.insertion_allele_lookup[result.insertion.uniquename] = [result.allele.uniquename]
        self.log.info(f'Found {counter} distinct insertion-allele relationships.')
        return

    def get_hemi_driver_info(self, session):
        """Get hemi-driver allele feature_ids."""
        self.log.info('Get hemi-driver allele feature_ids.')
        split_system = aliased(Feature, name='split_system')
        hemi_driver = aliased(Feature, name='hemi_driver')
        filters = (
            hemi_driver.is_obsolete.is_(False),
            hemi_driver.uniquename.op('~')(r'^FBal[0-9]{7}$'),
            split_system.is_obsolete.is_(False),
            split_system.uniquename.op('~')(r'^FBco[0-9]{7}$'),
            Cvterm.name == 'partially_produced_by',
        )
        results = session.query(split_system, hemi_driver).\
            select_from(hemi_driver).\
            join(FeatureRelationship, (FeatureRelationship.object_id == hemi_driver.feature_id)).\
            join(split_system, (FeatureRelationship.subject_id == split_system.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        for result in results:
            self.hemi_drivers.append(result.hemi_driver.feature_id)
            try:
                self.split_system_combos[result.split_system.uniquename].append(result.hemi_driver.uniquename)
            except KeyError:
                self.split_system_combos[result.split_system.uniquename] = [result.hemi_driver.uniquename]
        self.log.info(f'Found {len(self.hemi_drivers)} distinct hemi-driver alleles with FBco parents.')
        self.log.info(f'Found {len(self.split_system_combos.keys())} distinct split system combinations with hemi-driver components.')
        for hemidriver_list in self.split_system_combos.values():
            hemidriver_list.sort()
            self.split_system_combo_strs.append('|'.join(hemidriver_list))
        return

    # Elaborate on get_general_data() for the ExpressionHandler.
    def get_general_data(self, session):
        """Extend the method for the ExpressionHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.add_placeholder_cvterm()
        self.get_slim_term_mappings(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['transcript', 'polypeptide', 'allele', 'gene', 'insertion', 'split system combination'])
        self.get_isoform_mappings(session)
        self.get_gene_product_gene_mappings(session)
        self.get_allele_product_allele_mappings(session)
        self.get_insertion_allele_mappings(session)
        self.get_hemi_driver_info(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_expression_patterns(self, session):
        """Build a dictionary of expression patterns from the "expression" table."""
        self.log.info('Build a dictionary of expression patterns from the "expression" table.')
        # Note - get only expression patterns related to CV terms (ignore those with only a TAP statement note).
        filters = (
            Feature.is_obsolete.is_(False),
        )
        expression_patterns = session.query(Expression).\
            select_from(Expression).\
            join(FeatureExpression, (FeatureExpression.expression_id == Expression.expression_id)).\
            join(Feature, (Feature.feature_id == FeatureExpression.feature_id)).\
            join(ExpressionCvterm, (ExpressionCvterm.expression_id == FeatureExpression.expression_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for expression_pattern in expression_patterns:
            self.expression_patterns[expression_pattern.expression_id] = fb_datatypes.FBExpressionAnnotation(expression_pattern)
            counter += 1
        self.log.info(f'Found {counter} distinct expression patterns in chado.')
        return

    def get_expression_pattern_cvterms(self, session):
        """Get the cvterms for expression annotations."""
        self.log.info('Get the cvterms for expression annotations.')
        xprn_cvterm = aliased(Cvterm, name='xprn_cvterm')
        type_cvterm = aliased(Cvterm, name='type_cvterm')
        filters = (
            xprn_cvterm.is_obsolete == 0,
        )
        expression_cvterms = session.query(type_cvterm, ExpressionCvterm).\
            select_from(ExpressionCvterm).\
            join(xprn_cvterm, (xprn_cvterm.cvterm_id == ExpressionCvterm.cvterm_id)).\
            join(type_cvterm, (type_cvterm.cvterm_id == ExpressionCvterm.cvterm_type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        skip_counter = 0
        for result in expression_cvterms:
            xprn_id = result.ExpressionCvterm.expression_id
            # Skip "expression patterns" associated only with interaction or dataset tissue samples.
            if xprn_id not in self.expression_patterns.keys():
                skip_counter += 1
                continue
            # Skip reporting-build-derived interpolated stages.
            if result.ExpressionCvterm.rank == 9999:
                continue
            slot = f'{result.type_cvterm.name}_terms'
            xprn_pattern_slot = getattr(self.expression_patterns[xprn_id], slot)
            xprn_cvt_id = result.ExpressionCvterm.expression_cvterm_id
            xprn_pattern_slot[xprn_cvt_id] = fb_datatypes.FBExpressionCvterm(result.ExpressionCvterm)
            counter += 1
        self.log.info(f'Found {counter} distinct expression cvterm objects in chado.')
        self.log.info(f'Ignored {skip_counter} distinct expression cvterm objects associated only with interactions/datasets.')
        return

    def get_expression_pattern_operators(self, session):
        """Get the cvterm operators for expression annotations."""
        self.log.info('Get the cvterm operators for expression annotations.')
        xprn_cvterm = aliased(Cvterm, name='xprn_cvterm')
        type_cvterm = aliased(Cvterm, name='type_cvterm')
        operator_cvterm = aliased(Cvterm, name='operator_cvterm')
        operator_values = ['FROM', 'TO', 'OF']
        filters = (
            operator_cvterm.name == 'operator',
            Cv.name != 'FlyBase miscellaneous CV',
            xprn_cvterm.is_obsolete == 0,
            ExpressionCvtermprop.value.in_((operator_values)),
        )
        expression_cvterm_operators = session.query(type_cvterm, ExpressionCvterm, ExpressionCvtermprop).\
            select_from(ExpressionCvterm).\
            join(xprn_cvterm, (xprn_cvterm.cvterm_id == ExpressionCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == xprn_cvterm.cv_id)).\
            join(type_cvterm, (type_cvterm.cvterm_id == ExpressionCvterm.cvterm_type_id)).\
            join(ExpressionCvtermprop, (ExpressionCvtermprop.expression_cvterm_id == ExpressionCvterm.expression_cvterm_id)).\
            join(operator_cvterm, (operator_cvterm.cvterm_id == ExpressionCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in expression_cvterm_operators:
            xprn_id = result.ExpressionCvterm.expression_id
            if xprn_id not in self.expression_patterns.keys():
                continue
            # Skip reporting-build-derived interpolated stages.
            if result.ExpressionCvterm.rank == 9999:
                continue
            slot = f'{result.type_cvterm.name}_terms'
            xprn_pattern_slot = getattr(self.expression_patterns[xprn_id], slot)
            xprn_cvt_id = result.ExpressionCvterm.expression_cvterm_id
            operator_value = result.ExpressionCvtermprop.value
            xprn_pattern_slot[xprn_cvt_id].operators.append(operator_value)
            counter += 1
        self.log.info(f'Found {counter} distinct expression cvterm operators in chado.')
        return

    def get_feature_expression_data(self, session):
        """Get feature_expression data."""
        self.log.info('Get feature_expression data.')
        filters = (
            Feature.is_obsolete.is_(False),
        )
        feature_expressions = session.query(FeatureExpression).\
            select_from(Feature).\
            join(FeatureExpression, (Feature.feature_id == FeatureExpression.feature_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in feature_expressions:
            feat_xprn_id = result.feature_expression_id
            xprn_id = result.expression_id
            feat_type = result.feature.type.name
            if xprn_id not in self.expression_patterns.keys():
                continue
            feat_xprn = fb_datatypes.FBFeatureExpressionAnnotation(result)
            if 'RNA' in feat_type:
                feat_xprn.xprn_type = 'RNA'
            self.fb_data_entities[feat_xprn_id] = feat_xprn
            counter += 1
        self.log.info(f'Found {counter} distinct feature_expression annotations in chado.')
        return

    def get_xprn_notes(self, session):
        """Get feature_expression notes."""
        self.log.info('Get feature_expression notes.')
        filters = (
            Feature.is_obsolete.is_(False),
            Cvterm.name == 'comment',
        )
        feat_xprnprops = session.query(FeatureExpressionprop).\
            select_from(Feature).\
            join(FeatureExpression, (Feature.feature_id == FeatureExpression.feature_id)).\
            join(FeatureExpressionprop, (FeatureExpressionprop.feature_expression_id == FeatureExpression.feature_expression_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureExpressionprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in feat_xprnprops:
            if result.value is None:
                self.log.warning(f'feature_expression_id={result.feature_expression_id} has a NULL comment, skipping. ')
                continue
            try:
                cleaned_prop_text = clean_free_text(result.value)
                self.fb_data_entities[result.feature_expression_id].tap_stmt_notes.append(cleaned_prop_text)
                counter += 1
            except KeyError:
                continue
        self.log.info(f'Found {counter} distinct feature_expressionprop notes in chado.')
        return

    # Elaborate on get_datatype_data() for the ExpressionHandler.
    def get_datatype_data(self, session):
        """Extend the method for the ExpressionHandler."""
        super().get_datatype_data(session)
        self.get_expression_patterns(session)
        self.get_expression_pattern_cvterms(session)
        self.get_expression_pattern_operators(session)
        self.get_feature_expression_data(session)
        self.get_xprn_notes(session)
        return

    # Add methods to be run by synthesize_info() below.
    def assign_qualifiers(self):
        """Assign qualifiers to the appropriate primary CV terms."""
        self.log.info('Assign qualifiers to the appropriate primary CV terms.')
        for xprn_pattern in self.expression_patterns.values():
            for slot_type in self.slot_types:
                slot_name = f'{slot_type}_terms'
                xprn_pattern_slot = getattr(xprn_pattern, slot_name)
                if not xprn_pattern_slot:
                    continue
                # self.log.debug(f'Evaluate type={slot_type} for xprn_id={xprn_pattern.db_primary_id}')
                # First, some QC on rank values to make sure they match expectation.
                observed_rank_list = [i.chado_obj.rank for i in xprn_pattern_slot.values()]
                observed_rank_list.sort()
                expected_rank_list = list(range(0, len(observed_rank_list)))
                if observed_rank_list != expected_rank_list:
                    msg = f'Expression pattern {xprn_pattern.db_primary_id} has unexpected ranks. '
                    msg += f'For {slot_name}, has unexpected ranks: {observed_rank_list}. Expected: {expected_rank_list}.'
                    self.log.error(msg)
                    xprn_pattern.is_problematic = True
                    xprn_pattern.notes.append('Found problematic term ranks.')
                    continue
                # Sort qualifiers to their primary term, then delete them from the primary set.
                current_primary_cvt_id = None
                qualifier_xprn_cvt_ids = []
                rank_sorted_xprn_cvts = {}
                for xprn_cvt in xprn_pattern_slot.values():
                    rank_sorted_xprn_cvts[xprn_cvt.chado_obj.rank] = xprn_cvt
                for rank in expected_rank_list:
                    this_xprn_cvt = rank_sorted_xprn_cvts[rank]
                    # self.log.debug(f'Evaluate type={slot_type}, rank={this_xprn_cvt.chado_obj.rank}, term={this_xprn_cvt.cvterm_name}')
                    # Note that the check for FBcv OBO filters out the "presumptive" qualifier that has an internal ID so cannot be exported.
                    if this_xprn_cvt.cv_name == 'FlyBase miscellaneous CV':
                        qualifier_xprn_cvt_ids.append(this_xprn_cvt.db_primary_id)
                        if this_xprn_cvt.obo == 'FBcv':
                            xprn_pattern_slot[current_primary_cvt_id].qualifier_cvterm_ids.append(this_xprn_cvt.cvterm_id)
                        else:
                            self.log.debug(f'Ignoring non-FBcv qualifier: "{this_xprn_cvt.cvterm_name}".')
                    else:
                        current_primary_cvt_id = this_xprn_cvt.db_primary_id
                        # self.log.debug(f'Found primary term="{this_xprn_cvt.cvterm_name}", xprn_cvterm_id={this_xprn_cvt.db_primary_id}')
                for qualifier_xprn_cvt_id in qualifier_xprn_cvt_ids:
                    del xprn_pattern_slot[qualifier_xprn_cvt_id]
                # Final review of qualifier term assignments.
                # for xprn_cvt in xprn_pattern_slot.values():
                    # qualifiers = [self.cvterm_lookup[i]['name'] for i in xprn_cvt.qualifier_cvterm_ids]
                    # self.log.debug(f'Final term: type={slot_type}, rank={xprn_cvt.chado_obj.rank}, term={xprn_cvt.cvterm_name}, qualifiers={qualifiers}')
        return

    def identify_stage_ranges(self):
        """Identify stage ranges in expression patterns."""
        self.log.info('Identify stage ranges in expression patterns.')
        counter = 0
        prob_counter = 0
        for xprn_pattern in self.expression_patterns.values():
            start_terms = []
            end_terms = []
            for stage_term in xprn_pattern.stage_terms.values():
                if 'FROM' in stage_term.operators:
                    stage_term.is_stage_start = True
                    start_terms.append(stage_term)
                if 'TO' in stage_term.operators:
                    stage_term.is_stage_end = True
                    end_terms.append(stage_term)
            if not start_terms and not end_terms:
                continue
            elif len(start_terms) == 1 and len(end_terms) == 1:
                start_terms[0].has_stage_end = end_terms[0]
                # stage_range_string = f'{start_terms[0].cvterm_name}--{end_terms[0].cvterm_name}'
                # self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, found this stage range: {stage_range_string}')
                counter += 1
            else:
                self.log.error(f'Many/partial stage ranges found for xprn_id={xprn_pattern.db_primary_id}')
                xprn_pattern.is_problematic = True
                xprn_pattern.notes.append('Found many and/or partial stage ranges.')
                prob_counter += 1
        self.log.info(f'Found {counter} stage ranges in expression patterns.')
        self.log.info(f'Found {prob_counter} expression patterns with many/partial stage ranges that could not be processed.')
        return

    def identify_tissue_ranges(self, session):
        """Identify tissue ranges in expression patterns."""
        self.log.info('Identify tissue ranges in expression patterns.')
        counter = 0
        prob_counter = 0
        for xprn_pattern in self.expression_patterns.values():
            start_terms = []
            end_terms = []
            for anatomy_term in xprn_pattern.anatomy_terms.values():
                # Handle all anatomy terms as if they were tissue ranges.
                anatomy_term.has_anat_term_ids.append(anatomy_term.cvterm_id)
                if 'FROM' in anatomy_term.operators:
                    anatomy_term.is_anat_start = True
                    start_terms.append(anatomy_term)
                if 'TO' in anatomy_term.operators:
                    anatomy_term.is_anat_end = True
                    end_terms.append(anatomy_term)
                    if anatomy_term.qualifier_cvterm_ids:
                        # qualifiers = [self.cvterm_lookup[i]['name'] for i in anatomy_term.qualifier_cvterm_ids]
                        # self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, "{anatomy_term.cvterm_name}" has qualifiers: {qualifiers}.')
                        pass
            if not start_terms and not end_terms:
                continue
            elif len(start_terms) == 1 and len(end_terms) == 1:
                end_terms[0].has_anat_term_ids.append(start_terms[0].cvterm_id)    # Add start term to the range (under the end term object).
                end_terms[0].operators.extend(start_terms[0].operators)    # Propagate operators from start of tissue range to end.
                # Get intervening tissue range terms, then add them to the list of terms in the tissue range.
                # tissue_range_string = f'{start_terms[0].cvterm_name}--{end_terms[0].cvterm_name}'
                # self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, found this tissue range: {tissue_range_string}')
                rgx, start, end = self.regex_for_anatomical_terms_in_numerical_series(start_terms[0].cvterm_name,
                                                                                      end_terms[0].cvterm_name, xprn_pattern.db_primary_id)
                if not rgx:
                    xprn_pattern.is_problematic = True
                    xprn_pattern.notes.append('Could not process tissue range.')
                    prob_counter += 1
                    continue
                # self.log.debug(f'Look for terms between positions {start} and {end} matching this regex: {rgx}')
                anatomical_series_terms = self.get_anatomical_terms_by_regex(session, rgx)
                filtered_terms = self.select_in_range_anatomical_terms(anatomical_series_terms, rgx, start, end)
                anat_cvterm_ids = [i.cvterm_id for i in filtered_terms]
                end_terms[0].has_anat_term_ids.extend(anat_cvterm_ids)
                end_terms[0].has_anat_term_ids.sort()
                counter += 1
            else:
                self.log.error(f'Many/partial tissue ranges found for xprn_id={xprn_pattern.db_primary_id}')
                xprn_pattern.is_problematic = True
                xprn_pattern.notes.append('Found many and/or partial tissue ranges.')
                prob_counter += 1
        self.log.info(f'Found {counter} tissue ranges in expression patterns.')
        self.log.info(f'Found {prob_counter} expression patterns with many/partial tissue ranges that could not be processed.')
        return

    def identify_tissue_sub_parts(self):
        """Identify tissue sub_parts in expression patterns."""
        self.log.info('Identify tissue sub_parts in expression patterns.')
        counter = 0
        for xprn_pattern in self.expression_patterns.values():
            main_parts = []
            potential_sub_parts = []
            for anatomy_term in xprn_pattern.anatomy_terms.values():
                if 'OF' in anatomy_term.operators:
                    anatomy_term.is_main_part = True
                    main_parts.append(anatomy_term)
                else:
                    potential_sub_parts.append(anatomy_term)
            if main_parts:
                if not potential_sub_parts:
                    self.log.error(f'For xprn_id={xprn_pattern.db_primary_id}, there are "main" parts but no sub-parts.')
                else:
                    if len(potential_sub_parts) > 1:
                        self.log.warning(f'For xprn_id={xprn_pattern.db_primary_id}, there are {len(potential_sub_parts)} sub-parts.')
                    elif len(main_parts) > 1:
                        self.log.warning(f'For xprn_id={xprn_pattern.db_primary_id}, there are {len(main_parts)} main parts.')
                    # Move anatomy sub_parts into a separate dict within the expression pattern object.
                    for sub_part in potential_sub_parts:
                        xprn_pattern.sub_anatomy_terms[sub_part.db_primary_id] = sub_part
                        del xprn_pattern.anatomy_terms[sub_part.db_primary_id]
            if xprn_pattern.sub_anatomy_terms:
                counter += 1
                # n_main = len(xprn_pattern.anatomy_terms)
                # n_sub = len(xprn_pattern.sub_anatomy_terms)
                # n_combos = n_main * n_sub
                # self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, have {n_combos} main/sub_part combinations.')
            else:
                xprn_pattern.sub_anatomy_terms['placeholder'] = self.placeholder
        self.log.info(f'Found {counter} expression patterns having anatomy sub_parts.')
        return

    def generate_xprn_pattern_dict(self, xprn_id, assay_term, stage_term, anatomy_term, anatomy_sub_term, cellular_term):
        """Convert a specific combination of terms from an expression pattern into a simpler dict."""
        input_str = f'assay="{assay_term.cvterm_name}", stage="{stage_term.cvterm_name}", main_anatomy="{anatomy_term.cvterm_name}", '
        input_str += f'sub_part_anatomy="{anatomy_sub_term.cvterm_name}", cellular="{cellular_term.cvterm_name}"'
        # self.log.debug(f'Generate xprn_pattern_dict for xprn_id={xprn_id}; input: {input_str}')
        xprn_pattern_dict_list = []
        for main_part_id in anatomy_term.has_anat_term_ids:
            for sub_part_id in anatomy_sub_term.has_anat_term_ids:
                xprn_pattern_dict = {
                    'assay_cvterm_id': assay_term.cvterm_id,
                    'stage_start_cvterm_id': stage_term.cvterm_id,
                    'stage_end_cvterm_id': 'placeholder',
                    'stage_qualifier_cvterm_ids': stage_term.qualifier_cvterm_ids.copy(),
                    'stage_slim_cvterm_ids': self.cvterm_lookup[stage_term.cvterm_id]['slim_term_cvterm_ids'].copy(),
                    'anatomical_structure_cvterm_id': main_part_id,
                    'anatomical_structure_qualifier_cvterm_ids': anatomy_term.qualifier_cvterm_ids.copy(),
                    'anatomical_structure_slim_cvterm_ids': self.cvterm_lookup[main_part_id]['slim_term_cvterm_ids'].copy(),
                    'anatomical_substructure_cvterm_id': sub_part_id,
                    'anatomical_substructure_qualifier_cvterm_ids': anatomy_sub_term.qualifier_cvterm_ids.copy(),
                    'anatomical_substructure_slim_cvterm_ids': self.cvterm_lookup[sub_part_id]['slim_term_cvterm_ids'].copy(),
                    'cellular_component_cvterm_id': cellular_term.cvterm_id,
                    'cellular_component_qualifier_cvterm_ids': cellular_term.qualifier_cvterm_ids.copy(),
                }
                # Handle stage end.
                if stage_term.has_stage_end:
                    xprn_pattern_dict['stage_end_cvterm_id'] = stage_term.has_stage_end.cvterm_id
                    xprn_pattern_dict['stage_qualifier_cvterm_ids'].extend(stage_term.has_stage_end.qualifier_cvterm_ids)
                    additional_slim_term_ids = self.cvterm_lookup[stage_term.has_stage_end.cvterm_id]['slim_term_cvterm_ids'].copy()
                    xprn_pattern_dict['stage_slim_cvterm_ids'].extend(additional_slim_term_ids)
                    xprn_pattern_dict['stage_slim_cvterm_ids'] = list(set(xprn_pattern_dict['stage_slim_cvterm_ids']))
                xprn_pattern_dict_list.append(xprn_pattern_dict)
                # self.log.debug(f'For xprn_id={xprn_id}, generated xprn_pattern_dict: {xprn_pattern_dict}')
        return xprn_pattern_dict_list

    def split_out_expression_patterns(self):
        """Generate all combinations of anatomy/assay/cellular/stage terms for an expression pattern."""
        self.log.info('Generate all combinations of anatomy/assay/cellular/stage terms for an expression pattern.')
        for xprn_pattern in self.expression_patterns.values():
            xprn_id = xprn_pattern.db_primary_id
            n_combos = 0
            # Skip expression patterns already flagged as problematic.
            if xprn_pattern.is_problematic:
                continue
            # Add a placeholder to slots with zero terms.
            for slot_type in self.slot_types:
                slot_name = f'{slot_type}_terms'
                xprn_pattern_slot = getattr(xprn_pattern, slot_name)
                if not xprn_pattern_slot:
                    xprn_pattern_slot['placeholder'] = self.placeholder
            for assay_term in xprn_pattern.assay_terms.values():
                for stage_term in xprn_pattern.stage_terms.values():
                    # Skip stage range end terms, as these are folded into reporting of the stage range start term.
                    if stage_term.is_stage_end:
                        continue
                    for anatomy_term in xprn_pattern.anatomy_terms.values():
                        # Skip anatomy range start terms, as these are folded into reporting of the stage range end term.
                        if anatomy_term.is_anat_start:
                            continue
                        for anatomy_sub_term in xprn_pattern.sub_anatomy_terms.values():
                            # Skip sub_part anatomy range start terms, as these are folded into reporting of the stage range end term.
                            if anatomy_sub_term.is_anat_start:
                                continue
                            for cellular_term in xprn_pattern.cellular_terms.values():
                                xp_list = self.generate_xprn_pattern_dict(xprn_id, assay_term, stage_term, anatomy_term, anatomy_sub_term, cellular_term)
                                xprn_pattern.xprn_pattern_combos.extend(xp_list)
                                n_combos += len(xp_list)
            # self.log.debug(f'For xprn_id={xprn_id}, found {n_combos} total term combinations.')
        # Check these difficult cases in the output:
        # xprn_id=42175, <a> cell | subset &&of mesoderm | dorsal &&of parasegment 2--12
        # xprn_id=42170, <a> parasegment 3--12 &&of larval ventral nerve cord
        # xprn_id=34285, <a> parasegment 8 && parasegment 9 &&of midgut
        # xprn_id=32457, <a> neuron && glial cell &&of central nervous system
        return

    def determine_public_feature_to_report(self):
        """Determine the public feature to report for each expressed product."""
        self.log.info('Determine the public feature to report for expressed gene product.')
        split_system_counter = 0
        isoform_to_gene_counter = 0
        gene_product_to_gene_counter = 0
        allele_product_to_allele_counter = 0
        no_mapping_counter = 0
        hemidriver_counter = 0
        for feat_xprn in self.fb_data_entities.values():
            # 1. Deal with expressed split system combinations.
            if self.feature_lookup[feat_xprn.feature_id]['uniquename'].startswith('FBco'):
                feat_xprn.public_feature_id = feat_xprn.feature_id
                split_system_counter += 1
            # 2. Deal with gene product isoforms, which map to genes indirectly.
            elif feat_xprn.feature_id in self.isoform_gene_product_lookup.keys():
                gene_product_id = self.isoform_gene_product_lookup[feat_xprn.feature_id]
                try:
                    feat_xprn.public_feature_id = self.gene_product_gene_lookup[gene_product_id]
                    isoform_to_gene_counter += 1
                except KeyError:
                    feat_xprn.is_problematic = True
                    feat_xprn.notes.append('Could not find gene for gene product.')
                    self.log.error(f'Could not find gene for gene product ID {gene_product_id}.')
            # 3. Deal with gene products, which map to genes directly.
            elif feat_xprn.feature_id in self.gene_product_gene_lookup.keys():
                feat_xprn.public_feature_id = self.gene_product_gene_lookup[feat_xprn.feature_id]
                gene_product_to_gene_counter += 1
            # 4. Deal with allele products, which map to transgenic alleles.
            elif feat_xprn.feature_id in self.allele_product_allele_lookup.keys():
                allele_feature_id = self.allele_product_allele_lookup[feat_xprn.feature_id]
                insertion_rgx = r'FBti[0-9]{7}'
                partner_insertion_curies = []
                partner_allele_curies = []
                split_system_features_represented = []
                if allele_feature_id in self.hemi_drivers:
                    allele_curie = self.feature_lookup[allele_feature_id]['uniquename']
                    for note in feat_xprn.tap_stmt_notes:
                        if re.search(insertion_rgx, note):
                            partner_insertion_curies.extend(re.findall(insertion_rgx, note))
                    partner_insertion_curies = set(partner_insertion_curies)
                    for partner_insertion_curie in partner_insertion_curies:
                        partner_allele_curies.extend(self.insertion_allele_lookup[partner_insertion_curie])
                    partner_allele_curies = set(partner_allele_curies)
                    for partner_allele_curie in partner_allele_curies:
                        pair_combo = [allele_curie, partner_allele_curie]
                        pair_combo.sort()
                        pair_combo_str = '|'.join(pair_combo)
                        if pair_combo_str in self.split_system_combo_strs:
                            split_system_features_represented.append(pair_combo_str)
                if split_system_features_represented:
                    feat_xprn.is_problematic = True
                    feat_xprn.notes.append('Suppress export of hemi-driver expression having split system combination feature.')
                    hemidriver_counter += 1
                    self.log.debug(f'BOB: Suppress feat_xprn_id={feat_xprn.feature_expression_id} for split system combos {split_system_features_represented}')
                else:
                    feat_xprn.public_feature_id = allele_feature_id
                    allele_product_to_allele_counter += 1
            # 5. Deal with expressed products that cannot be mapped to a gene or allele.
            else:
                feat_xprn.is_problematic = True
                feat_xprn.notes.append('Could not map gene product to gene or allele.')
                self.log.error(f'Could not map feature ID {feat_xprn.feature_id} to a gene or allele.')
                no_mapping_counter += 1
        self.log.info(f'Found {split_system_counter} split system combination annotations.')
        self.log.info(f'Mapped isoforms indirectly to genes (via gene products) for {isoform_to_gene_counter} annotations.')
        self.log.info(f'Mapped gene products directly to genes for {gene_product_to_gene_counter} annotations.')
        self.log.info(f'Mapped allele products directly to alleles for {allele_product_to_allele_counter} annotations.')
        self.log.info(f'Could not map the expressed product to a gene or allele for {no_mapping_counter} annotations.')
        self.log.info(f'Supressing {hemidriver_counter} hemidriver annotations.')
        return

    def process_for_tsv_export(self):
        """Process expression patterns for export to TSV."""
        self.log.info('Process expression patterns for export to TSV.')
        counter = 0
        for feat_xprn in self.fb_data_entities.values():
            if feat_xprn.is_problematic:
                continue
            xprn_pattern = self.expression_patterns[feat_xprn.expression_id]
            if xprn_pattern.is_problematic:
                continue
            elif not xprn_pattern.xprn_pattern_combos:
                continue
            for xp_combo in xprn_pattern.xprn_pattern_combos:
                xprn_tsv_dict = {
                    'feature_id': self.feature_lookup[feat_xprn.public_feature_id]['uniquename'],
                    'feature_symbol': self.feature_lookup[feat_xprn.public_feature_id]['name'],
                    'reference_id': feat_xprn.pub_curie,
                    'expression_type': feat_xprn.xprn_type,
                    'expression_id': xprn_pattern.db_primary_id,
                    'assay_term': self.cvterm_lookup[xp_combo['assay_cvterm_id']]['name'],
                    'stage_start': self.cvterm_lookup[xp_combo['stage_start_cvterm_id']]['name_plus_curie'],
                    'stage_end': self.cvterm_lookup[xp_combo['stage_end_cvterm_id']]['name_plus_curie'],
                    'stage_qualifiers': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in xp_combo['stage_qualifier_cvterm_ids']]),
                    'stage_slim_terms': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in xp_combo['stage_slim_cvterm_ids']]),
                    'anatomical_structure_term': self.cvterm_lookup[xp_combo['anatomical_structure_cvterm_id']]['name_plus_curie'],
                    'anatomical_structure_qualifiers': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in
                                                                  xp_combo['anatomical_structure_qualifier_cvterm_ids']]),
                    'anatomical_structure_slim_terms': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in
                                                                  xp_combo['anatomical_structure_slim_cvterm_ids']]),
                    'anatomical_substructure_term': self.cvterm_lookup[xp_combo['anatomical_substructure_cvterm_id']]['name_plus_curie'],
                    'anatomical_substructure_qualifiers': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in
                                                                     xp_combo['anatomical_substructure_qualifier_cvterm_ids']]),
                    'anatomical_substructure_slim_terms': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in
                                                                     xp_combo['anatomical_substructure_slim_cvterm_ids']]),
                    'cellular_component_term': self.cvterm_lookup[xp_combo['cellular_component_cvterm_id']]['name_plus_curie'],
                    'cellular_component_qualifiers': ' | '.join([self.cvterm_lookup[i]['name_plus_curie'] for i in
                                                                xp_combo['cellular_component_qualifier_cvterm_ids']]),
                    'notes': ' | '.join(feat_xprn.tap_stmt_notes),
                }
                self.export_data_for_tsv.append(xprn_tsv_dict)
                counter += 1
        self.log.info(f'Generated {counter} expression pattern TSV rows.')
        return

    # Elaborate on synthesize_info() for the ExpressionHandler.
    def synthesize_info(self, session):
        """Extend the method for the ExpressionHandler."""
        super().synthesize_info()
        self.assign_qualifiers()
        self.identify_stage_ranges()
        self.identify_tissue_ranges(session)
        self.identify_tissue_sub_parts()
        self.determine_public_feature_to_report()
        self.split_out_expression_patterns()
        return

    # Additional methods to be run by map_fb_data_to_alliance() below.
    # Placeholder
