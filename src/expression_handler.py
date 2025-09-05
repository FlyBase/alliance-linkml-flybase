"""Module:: expression_handler.

Synopsis:
    Data handler for FlyBase expression annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import re
from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.reporting import (
    Cv, Cvterm, Db, Dbxref, Expression, ExpressionCvterm, ExpressionCvtermprop,
    Feature, FeatureExpression    # FeatureExpressionprop
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
        self.expression_patterns = {}    # expression_id-keyed FBExpressionAnnotation objects.
        self.feat_xprn_annos = {}        # feature_expression_id-keyed FBFeatureExpressionAnnotation objects, for export.

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
    def regex_for_anatomical_terms_in_numerical_series(self, start_term, end_term):
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
            self.log.error(f'No number or letter+number pattern found for these terms: {start_term}, {end_term}')
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
        self.log.debug(f'For {start_term} and {end_term}, found this range: {start_position}--{end_position}')
        self.log.debug(f'For {start_term} and {end_term}, found this regex: {regex}')
        return regex, start_position, end_position

    def get_anatomical_terms_by_regex(self, session, term_regex):
        """Get anatomical terms matching a given regex."""
        self.log.debug(f'Find terms matching this regex: {term_regex}')
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
            self.log.debug(f'For {term_regex}, found these matching FBbt terms: {term_names}')
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
        self.log.debug(f'From {start}--{end}, look for numbers between {num_start} and {num_end}.')
        for term in terms:
            self.log.debug(f'Is "{term.name}" in range?')
            position = re.search(rgx, term.name).group(1)
            self.log.debug(f'Found this position "{position}".')
            num_position = int(re.search(num_rgx, position).group(1))
            self.log.debug(f'Term {term.name} is at position {position}, with number={num_position}')
            if num_position > num_start and num_position < num_end:
                filtered_terms.append(term)
        self.log.debug(f'Found these tissue range terms: {[i.name for i in filtered_terms]}')
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
            'curie': ''
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
                    self.log.debug(f'Slim term: cv_name={cv_name}, cvterm_name="{slim_term_name}", cvterm_id={slim_term.cvterm_id}')
                # Second, get the child terms under each slim term.
                child_cvterm_ids = self.get_child_cvterms(session, slim_term_name, cv_name)
                child_term_names = [self.cvterm_lookup[i]['name'] for i in child_cvterm_ids if i in self.cvterm_lookup.keys()]
                child_term_names.sort()
                child_term_name_str = '\n'.join(child_term_names)
                self.log.debug(f'Found {len(child_term_names)} child terms for the {cv_name} slim term "{slim_term_name}":\n{child_term_name_str}')
                for child_cvterm_id in child_cvterm_ids:
                    self.cvterm_lookup[child_cvterm_id]['slim_term_cvterm_ids'].append(slim_term.cvterm_id)
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

    # Elaborate on get_datatype_data() for the ExpressionHandler.
    def get_datatype_data(self, session):
        """Extend the method for the ExpressionHandler."""
        super().get_datatype_data(session)
        self.get_expression_patterns(session)
        self.get_expression_pattern_cvterms(session)
        self.get_expression_pattern_operators(session)
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
                self.log.debug(f'Evaluate type={slot_type} for xprn_id={xprn_pattern.db_primary_id}')
                # First, some QC on rank values to make sure they match expectation.
                observed_rank_list = [i.chado_obj.rank for i in xprn_pattern_slot.values()]
                observed_rank_list.sort()
                expected_rank_list = list(range(0, len(observed_rank_list)))
                if observed_rank_list != expected_rank_list:
                    msg = f'Expression pattern {xprn_pattern.db_primary_id} has unexpected ranks. '
                    msg += f'For {slot_name}, has unexpected ranks: {observed_rank_list}. Expected: {expected_rank_list}.'
                    self.log.warning(msg)
                    xprn_pattern.is_problematic = True
                    xprn_pattern.notes.append('Found problematic term ranks.')
                    continue
                # Sort qualifiers to their primary term, then delete them from the primary set.
                current_primary_cvt_id = None
                current_primary_cvterm_name = None
                qualifier_xprn_cvt_ids = []
                rank_sorted_xprn_cvts = {}
                for xprn_cvt in xprn_pattern_slot.values():
                    rank_sorted_xprn_cvts[xprn_cvt.chado_obj.rank] = xprn_cvt
                for rank in expected_rank_list:
                    this_xprn_cvt = rank_sorted_xprn_cvts[rank]
                    self.log.debug(f'Evaluate type={slot_type}, rank={this_xprn_cvt.chado_obj.rank}, term={this_xprn_cvt.cvterm_name}')
                    if this_xprn_cvt.cv_name == 'FlyBase miscellaneous CV':
                        xprn_pattern_slot[current_primary_cvt_id].qualifier_cvterm_ids.append(this_xprn_cvt.cvterm_id)
                        qualifier_xprn_cvt_ids.append(this_xprn_cvt.db_primary_id)
                        self.log.debug(f'Found qualifier {this_xprn_cvt.cvterm_name} for primary term="{current_primary_cvterm_name}"')
                    else:
                        current_primary_cvt_id = this_xprn_cvt.db_primary_id
                        current_primary_cvterm_name = this_xprn_cvt.cvterm_name
                        self.log.debug(f'Found primary term="{current_primary_cvterm_name}", xprn_cvterm_id={this_xprn_cvt.db_primary_id}')
                for qualifier_xprn_cvt_id in qualifier_xprn_cvt_ids:
                    del xprn_pattern_slot[qualifier_xprn_cvt_id]
                # Final review of terms.
                for xprn_cvt in xprn_pattern_slot.values():
                    qualifiers = [self.cvterm_lookup[i]['name'] for i in xprn_cvt.qualifier_cvterm_ids]
                    self.log.debug(f'Final term: type={slot_type}, rank={xprn_cvt.chado_obj.rank}, term={xprn_cvt.cvterm_name}, qualifiers={qualifiers}')
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
                stage_range_string = f'{start_terms[0].cvterm_name}--{end_terms[0].cvterm_name}'
                self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, found this stage range: {stage_range_string}')
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
                if 'FROM' in anatomy_term.operators:
                    anatomy_term.is_anat_start = True
                    start_terms.append(anatomy_term)
                if 'TO' in anatomy_term.operators:
                    anatomy_term.is_anat_end = True
                    end_terms.append(anatomy_term)
                    if anatomy_term.qualifier_cvterm_ids:
                        qualifiers = [self.cvterm_lookup[i]['name'] for i in anatomy_term.qualifier_cvterm_ids]
                        self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, "{anatomy_term.cvterm_name}" has qualifiers: {qualifiers}.')
            if not start_terms and not end_terms:
                continue
            elif len(start_terms) == 1 and len(end_terms) == 1:
                # Add start term to the range (under the end term object).
                end_terms[0].has_anat_terms.append(start_terms[0].cvterm_id)
                end_terms[0].operators.extend(start_terms[0].operators)    # Propagate operators from start of tissue range to end.
                # Get intervening tissue range terms, then add them to the list of terms in the tissue range.
                tissue_range_string = f'{start_terms[0].cvterm_name}--{end_terms[0].cvterm_name}'
                self.log.debug(f'For xprn_id={xprn_pattern.db_primary_id}, found this tissue range: {tissue_range_string}')
                rgx, start, end = self.regex_for_anatomical_terms_in_numerical_series(start_terms[0].cvterm_name, end_terms[0].cvterm_name)
                self.log.debug(f'Look for terms between positions {start} and {end} matching this regex: {rgx}')
                anatomical_series_terms = self.get_anatomical_terms_by_regex(session, rgx)
                filtered_terms = self.select_in_range_anatomical_terms(anatomical_series_terms, rgx, start, end)
                if filtered_terms:
                    anat_cvterm_ids = [i.cvterm_id for i in filtered_terms]
                else:
                    xprn_pattern.is_problematic = True
                    xprn_pattern.notes.append('Could not process tissue range.')
                end_terms[0].has_anat_terms.extend(anat_cvterm_ids)
                end_terms[0].has_anat_terms.sort()
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
                    self.log.error(f'For xprn_id={xprn_pattern.db_primary_id}, have "main" parts but no sub-parts.')
                elif len(potential_sub_parts) > 1:
                    self.log.error(f'For xprn_id={xprn_pattern.db_primary_id}, have MANY sub-parts.')
                else:
                    for sub_part in potential_sub_parts:
                        sub_part.is_sub_part = True
                        for main_part in main_parts:
                            main_part.has_sub_part = sub_part
                            msg = f'For xprn_id={xprn_pattern.db_primary_id}, found this part-sub_part pair: '
                            msg += f'part="{main_part.cvterm_name}", sub_part="{sub_part.cvterm_name}"'
                            self.log.debug(msg)
                            counter += 1
        self.log.info(f'Found {counter} part-sub_part annotations.')
        # BOB: Check xprn_id=42175, <a> cell | subset &&of mesoderm | dorsal &&of parasegment 2--12
        return

    def generate_xprn_pattern_dict(self, assay_term, stage_term, anatomy_term, cellular_term, anatomy_range_term_id):
        """Convert a specific combination of terms from an expression pattern into a simpler dict."""
        xprn_pattern_dict = {
            'assay': self.cvterm_lookup[assay_term.cvterm_id]['name'],
            'stage_start_name': self.cvterm_lookup[stage_term.cvterm_id]['name'],
            'stage_start_id': self.cvterm_lookup[stage_term.cvterm_id]['curie'],
            'stage_end_name': None,
            'stage_end_id': None,
            'stage_qualifier_names': None,
            'stage_qualifier_ids': None,
            'anatomical_structure_name': self.cvterm_lookup[anatomy_term.cvterm_id]['name'],
            'anatomical_structure_id': self.cvterm_lookup[anatomy_term.cvterm_id]['curie'],
            'anatomical_structure_qualifier_names': '|'.join([self.cvterm_lookup[i]['name'] for i in anatomy_term.qualifier_cvterm_ids]),
            'anatomical_structure_qualifier_ids': '|'.join([self.cvterm_lookup[i]['curie'] for i in anatomy_term.qualifier_cvterm_ids]),
            'anatomical_substructure_name': None,
            'anatomical_substructure_id': None,
            'anatomical_substructure_qualifier_names': None,
            'anatomical_substructure_qualifier_ids': None,
            'cellular_component_name': self.cvterm_lookup[cellular_term.cvterm_id]['curie'],
            'cellular_component_id': self.cvterm_lookup[cellular_term.cvterm_id]['curie'],
            'cellular_component_qualifier_names': None,
            'cellular_component_qualifier_ids': None,
        }
        # Stage end.
        if stage_term.has_stage_end:
            xprn_pattern_dict['stage_end_name'] = self.cvterm_lookup[stage_term.has_stage_end.cvterm_id]['name']
            xprn_pattern_dict['stage_end_id'] = self.cvterm_lookup[stage_term.has_stage_end.cvterm_id]['curie']
        # Stage qualifiers.
        stage_qualifier_cvterm_ids = []
        stage_qualifier_cvterm_ids.extend(stage_term.qualifier_cvterm_ids)
        if stage_term.has_stage_end:
            stage_qualifier_cvterm_ids.extend(stage_term.has_stage_end.qualifier_cvterm_ids)
        xprn_pattern_dict['stage_qualifier_names'] = '|'.join([self.cvterm_lookup[i]['name'] for i in stage_qualifier_cvterm_ids])
        xprn_pattern_dict['stage_qualifier_ids'] = '|'.join([self.cvterm_lookup[i]['curie'] for i in stage_qualifier_cvterm_ids])
        # If processing a term in a tissue range, replace the main anatomy term with the cvterm_id given for the term within the range.
        if anatomy_range_term_id:
            xprn_pattern_dict['anatomical_structure_name'] = self.cvterm_lookup[anatomy_range_term_id]['name']
            xprn_pattern_dict['anatomical_structure_id'] = self.cvterm_lookup[anatomy_range_term_id]['curie']
        # Anatomy sub_parts.
        if anatomy_term.has_sub_part:
            xprn_pattern_dict['anatomical_substructure_name'] = self.cvterm_lookup[anatomy_term.has_sub_part.cvterm_id]['name']
            xprn_pattern_dict['anatomical_substructure_id'] = self.cvterm_lookup[anatomy_term.has_sub_part.cvterm_id]['curie']
            qual_name_str = '|'.join([self.cvterm_lookup[i]['name'] for i in anatomy_term.has_sub_part.qualifier_cvterm_ids])
            qual_id_str = '|'.join([self.cvterm_lookup[i]['curie'] for i in anatomy_term.has_sub_part.qualifier_cvterm_ids])
            xprn_pattern_dict['anatomical_substructure_qualifier_names'] = qual_name_str
            xprn_pattern_dict['anatomical_substructure_qualifier_ids'] = qual_id_str
        # Cellular component qualifiers.
        xprn_pattern_dict['cellular_component_qualifier_names'] = '|'.join([self.cvterm_lookup[i]['name'] for i in cellular_term.qualifier_cvterm_ids])
        xprn_pattern_dict['cellular_component_qualifier_ids'] = '|'.join([self.cvterm_lookup[i]['curie'] for i in cellular_term.qualifier_cvterm_ids])
        return xprn_pattern_dict

    def split_out_expression_patterns(self):
        """Generate all combinations of anatomy/assay/cellular/stage terms for an expression pattern."""
        self.log.info('Generate all combinations of anatomy/assay/cellular/stage terms for an expression pattern.')
        for xprn_pattern in self.expression_patterns.values():
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
                    for cellular_term in xprn_pattern.cellular_terms.values():
                        for anatomy_term in xprn_pattern.anatomy_terms.values():
                            # Skip anatomy sub_part terms, as these are folded into reporting of the main_part term.
                            if anatomy_term.is_sub_part:
                                continue
                            # Skip anatomy range start terms, as these are folded into reporting of the stage range end term.
                            elif anatomy_term.is_anat_start:
                                continue
                            # If there is an anatomy range, first report the start and intermediate anatomy range terms.
                            for anatomy_range_term_id in anatomy_term.has_anat_terms:
                                xprn_pattern_dict = self.generate_xprn_pattern_dict(assay_term, stage_term, anatomy_term, cellular_term, anatomy_range_term_id)
                                xprn_pattern.xprn_pattern_combos.append(xprn_pattern_dict)
                            # Report for the main anatomy term, including end terms for anatomy ranges.
                            xprn_pattern_dict = self.generate_xprn_pattern_dict(assay_term, stage_term, anatomy_term, cellular_term, None)
                            xprn_pattern.xprn_pattern_combos.append(xprn_pattern_dict)
        return

    # Elaborate on synthesize_info() for the ExpressionHandler.
    def synthesize_info(self, session):
        """Extend the method for the ExpressionHandler."""
        super().synthesize_info()
        self.assign_qualifiers()
        self.identify_stage_ranges()
        self.identify_tissue_ranges(session)
        self.identify_tissue_sub_parts()
        self.split_out_expression_patterns()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    # Placeholder.

    # Elaborate on map_fb_data_to_alliance() for the ExpressionHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the ExpressionHandler."""
        super().map_fb_data_to_alliance()
        return
