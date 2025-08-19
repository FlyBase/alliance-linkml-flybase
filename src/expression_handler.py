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
        self.fb_export_type = fb_datatypes.FBExpressionAnnotation
        self.agr_export_type = None
        self.expression_patterns = {}    # expression_id-keyed FBExpressionAnnotation objects.
        self.feat_xprn_annos = {}        # feature_expression_id-keyed FBFeatureExpressionAnnotation objects, for export.


    # Utility functions.
    def regex_for_anatomical_terms_in_numerical_series(self, term):
        """For an anatomical term representing some numerical series, return a regex for all terms in that series, and the series position of the term."""
        # Use 1: For a term, find the regex to get all terms in the series.
        # Use 2: For a term, identify the numerical position in the series.
        # Finally, compare the number of each term in the series to the terms having the "FROM"/"TO" operators to determine which terms are in range.

        # First get the positional part of the term.
        # Case 1. Exceptional intersegmental region terms: e.g., "embryonic/larval abdominal 1-2 intersegmental apodeme".
        intersegmental_positions = ['1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9']
        if 'intersegmental' in term:
            rgx = r' ([0-9]-[0-9]) '
            match = re.search(rgx, term)
            # Filter out non-sensical intersegmental positions.
            if match and match.group(1) not in intersegmental_positions:
                match = None
        # Case 2. Standard number at end, which trumps any interal numbers: e.g., "1" in "medulla non-directional M6 local neuron 1".
        else:
            rgx = r' ([A-Z]{0,1}[0-9]{1,2})$'
            match = re.search(rgx, term)
            # Case 3. If no match, look for a number internally: e.g., "7" in "abdominal 7 neuroblast NB3-5", "A3" in "larval serotonergic A3 neuron".
            if not match:
                rgx = r' ([A-Z]{0,1}[0-9]{1,2})'
                match = re.search(rgx, term)
        # Build the regex.
        if not match:
            self.log.error(f'No number or letter+number pattern found in term: {term}')
            return None, None
        num_group = match.group(1)
        self.log.debug(f'Found series position "{match.group(1)}" for term: {term}')
        before = re.escape(term[:match.start()])
        after = re.escape(term[match.end():])
        
        # if starts with letters, keep them; replace numbers
        letter_part = re.match(r'([A-Za-z]*)', num_group).group(1)
        regex_number_part = r'\d+'
        group_regex = letter_part + regex_number_part
        
        regex = f'{before}{group_regex}{after}'
        return regex, num_group

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
        filters = (
            Feature.is_obsolete.is_(False),            
        )
        expression_patterns = session.query(Expression).\
            select_from(Expression).\
            join(FeatureExpression, (FeatureExpression.expression_id == Expression.expression_id)).\
            join(Feature, (Feature.feature_id == FeatureExpression.feature_id)).\
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
        xprn_cvterm = aliased(ExpressionCvterm, name='xprn_cvterm')
        type_cvterm = aliased(Cvterm, name='type_cvterm')
        operator_cvterm = aliased(Cvterm, name='operator_cvterm')
        operator_values = ['FROM', 'TO', 'OF']
        filters = (
            Cv.cv_name != 'FlyBase miscellaneous CV',
            xprn_cvterm.is_obsolete.is_(False),
            operator_cvterm.is_obsolete.is_(False),
            ExpressionCvtermprop.value.in_((operator_values)),
        )
        expression_cvterm_operators = session.query(type_cvterm, ExpressionCvterm, ExpressionCvtermprop).\
            select_from(ExpressionCvterm).\
            join(xprn_cvterm, (xprn_cvterm.cvterm_id == ExpressionCvterm.cvterm_id)).\
            join(Cv, (Cv.cv_id == xprn_cvterm.cv_id)).\
            join(type_cvterm, (type_cvterm.cvterm_id == ExpressionCvterm.cvterm_type_id)).\
            join(ExpressionCvtermprop, (ExpressionCvtermprop.expression_cvterm_id == ExpressionCvterm.expression_cvterm_id)).\
            join(operator_cvterm, (operator_cvterm.cvterm_id == ExpressionCvterm.type_id)).\
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
    # BOB - CONTINUE HERE.
    def assign_qualifiers(self):
        """Assign qualifiers to the appropriate primary CV terms."""
        self.log.info('Assign qualifiers to the appropriate primary CV terms.')
        term_types = ['assay', 'stage', 'anatomy', 'cellular']
        for xprn_pattern in self.expression_patterns.values():
            for term_type in term_types:
                slot_name = f'{term_type}_terms'
                if not xprn_pattern[slot_name]:
                    continue
                observed_rank_list = [i.rank for i in xprn_pattern[slot_name].values()]
                observed_rank_list.sort()
                expected_rank_list = list(range(0, len(observed_rank_list)))
                if observed_rank_list != expected_rank_list:
                    self.log.warning(f'Expression pattern {xprn_pattern.expression_id} has unexpected ranks: {observed_rank_list}. Expected: {expected_rank_list}.')
                    continue
        return

    # BOB: Group end stage with start stage.
    # BOB: Interpolate tissue ranges (and propagate qualifiers from end term to all others).

    # Elaborate on synthesize_info() for the ExpressionHandler.
    def synthesize_info(self):
        """Extend the method for the ExpressionHandler."""
        super().synthesize_info()
        self.assign_qualifiers()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    # Placeholder.

    # Elaborate on map_fb_data_to_alliance() for the ExpressionHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the ExpressionHandler."""
        super().map_fb_data_to_alliance()
        return
