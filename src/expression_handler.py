"""Module:: expression_handler.

Synopsis:
    Data handler for FlyBase expression annotations.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import csv
import re
from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.char_conversions import sgml_to_plain_text
from harvdev_utils.reporting import (
    Cv, Cvterm, Db, Dbxref, Expression, ExpressionCvterm, ExpressionCvtermprop,
    Feature, FeatureExpression, FeatureExpressionprop, Pub
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
        self.expression_patterns = {}    # expression_id-keyed FBExpressionAnnotation objects.
        self.feat_xprn_annos = {}        # feature_expression_id-keyed FBFeatureExpressionAnnotation objects, for export.

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
            self.log.debug(f'For {term_regex}, found these matching FBbt terms: {term_names}\n')
        else:
            self.log.error(f'For {term_regex}, found NO matching FBbt terms.\n')
        return terms

    def select_in_range_anatomical_terms(self, terms, rgx, start, end):
        """Filter for anatomical terms within a certain range."""
        filtered_terms = []
        # For intersegmental terms, this approach just uses the first number in the \d-\d pattern for range comparison.
        num_rgx = r'^[A-ZAa-z]{0,2}(\d+)'
        num_start = int(re.match(num_rgx, start).group(1))
        num_end = int(re.match(num_rgx, end).group(1))
        self.log.debug(f'From {start}--{end}, look for numbers between {num_start} and {num_end}.')
        for term in terms:
            self.log.debug(f'Is "{term.name}" in range?')
            position = re.search(rgx, term.name).group(1)
            self.log.debug(f'FOund this position "{position}".')
            num_position = int(re.match(num_rgx, position).group(1))
            self.log.debug(f'Term {term} is at position {position}, with number={num_position}')
            if num_position > num_start and num_position < num_end:
                filtered_terms.append(term)
        self.log.debug(f'Found these terms: {filtered_terms}')
        return filtered_terms

    # Add methods to be run by get_general_data() below.
    # Placeholder.

    # Elaborate on get_general_data() for the ExpressionHandler.
    def get_general_data(self, session):
        """Extend the method for the ExpressionHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
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
        term_types = ['assay', 'stage', 'anatomy', 'cellular']
        for xprn_pattern in self.expression_patterns.values():
            for term_type in term_types:
                slot_name = f'{term_type}_terms'
                xprn_pattern_slot = getattr(xprn_pattern, slot_name)
                if not xprn_pattern_slot:
                    continue
                self.log.debug(f'Evaluate type={term_type} for xprn_id={xprn_pattern.db_primary_id}')
                # First, some QC on rank values to make sure they match expectation.
                observed_rank_list = [i.chado_obj.rank for i in xprn_pattern_slot.values()]
                observed_rank_list.sort()
                expected_rank_list = list(range(0, len(observed_rank_list)))
                if observed_rank_list != expected_rank_list:
                    msg = f'Expression pattern {xprn_pattern.db_primary_id} has unexpected ranks. '
                    msg += f'For {slot_name}, has unexpected ranks: {observed_rank_list}. Expected: {expected_rank_list}.'
                    self.log.warning(msg)
                    xprn_pattern.is_problematic = True
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
                    self.log.debug(f'Evaluate type={term_type}, rank={this_xprn_cvt.chado_obj.rank}, term={this_xprn_cvt.cvterm_name}')
                    if this_xprn_cvt.cv_name == 'FlyBase miscellaneous CV':
                        xprn_pattern_slot[current_primary_cvt_id].qualifier_cvterm_ids.append(this_xprn_cvt.cvterm_id)
                        qualifier_xprn_cvt_ids.append(this_xprn_cvt.db_primary_id)
                        self.log.debug(f'Found qualifier {this_xprn_cvt.cvterm_name} for primary term="{current_primary_cvterm_name}"')
                    else:
                        current_primary_cvt_id = this_xprn_cvt.db_primary_id
                        current_primary_cvterm_name = this_xprn_cvt.cvterm_name
                        self.log.debug(f'Found primary term="{current_primary_cvterm_name}"')
                for qualifier_xprn_cvt_id in qualifier_xprn_cvt_ids:
                    del xprn_pattern_slot[qualifier_xprn_cvt_id]
                # Final review of terms.
                for xprn_cvt in xprn_pattern_slot.values():
                    qualifiers = [self.cvterm_lookup[i]['name'] for i in xprn_cvt.qualifier_cvterm_ids]
                    self.log.debug(f'Final term: type={term_type}, rank={xprn_cvt.chado_obj.rank}, term={xprn_cvt.cvterm_name}, qualifiers={qualifiers}')
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
                tissue_range_string = f'{start_terms[0].cvterm_name}--{end_terms[0].cvterm_name}'
                self.log.debug(f'\nFor {xprn_pattern.db_primary_id}, found this tissue range: {tissue_range_string}')
                rgx, start, end = self.regex_for_anatomical_terms_in_numerical_series(start_terms[0].cvterm_name, end_terms[0].cvterm_name)
                self.log.debug(f'Look for terms between positions {start} and {end} matching this regex: {rgx}')
                anatomical_series_terms = self.get_anatomical_terms_by_regex(session, rgx)
                # BILLY BOB: CONTINUE HERE, make sure range filtering works.
                filtered_terms = self.select_in_range_anatomical_terms(anatomical_series_terms, rgx, start, end)
                # BOB - add dummy FBExpressionCvterm objects to self.anatomy_terms for these interpolated terms.
                # BOB - propagate qualifiers from range end to start and intervening terms.
                counter += 1
            else:
                self.log.error(f'Many/partial tissue ranges found for xprn_id={xprn_pattern.db_primary_id}')
                xprn_pattern.is_problematic = True
                prob_counter += 1
        self.log.info(f'Found {counter} tissue ranges in expression patterns.')
        self.log.info(f'Found {prob_counter} expression patterns with many/partial tissue ranges that could not be processed.')
        return

    def identify_tissue_subparts(self):
        # BOB: Identify part/subpart tissues.
        # Check xprn_id=42175, cell | subset &&of mesoderm | dorsal &&of parasegment 2--12
        return

    # Elaborate on synthesize_info() for the ExpressionHandler.
    def synthesize_info(self, session):
        """Extend the method for the ExpressionHandler."""
        super().synthesize_info()
        self.assign_qualifiers()
        self.identify_stage_ranges()
        self.identify_tissue_ranges(session)
        self.identify_tissue_subparts()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    # Placeholder.

    # Elaborate on map_fb_data_to_alliance() for the ExpressionHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the ExpressionHandler."""
        super().map_fb_data_to_alliance()
        return
