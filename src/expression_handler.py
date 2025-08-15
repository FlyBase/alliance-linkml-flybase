"""Module:: expression_handler.

Synopsis:
    Data handler for FlyBase expression annotations.

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
    Cv, Cvterm, Db, Dbxref, Expression, ExpressionCvterm, ExpressionCvtermprop,
    Feature, FeatureExpression, FeatureExpressionprop, Pub
)
import agr_datatypes
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
        self.primary_export_set = 'temp_expression_ingest_set'
        self.expression_patterns = {}    # expression_id-keyed FBExpressionAnnotation objects.
        self.fb_data_entities = {}       # feature_expression_id-keyed FBFeatureExpressionAnnotation objects, for export.


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
            self.expression_patterns[expression_pattern] = fb_datatypes.FBExpressionAnnotation(expression_pattern)
            counter += 1
        self.log.info(f'Found {counter} distinct expression patterns in chado.')
        return

    # def get cvterms
    # def get qualifiers
    # other?


    # Elaborate on get_general_data() for the AGMDiseaseHandler.
    def get_general_data(self, session):
        """Extend the method for the ExpressionHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        # BOB: Query for FBco split system combnations is not WORKING!!!!
        self.build_feature_lookup(session, feature_types=['transcript', 'polypeptide', 'allele', 'gene', 'insertion', 'split system combination'])
        self.get_expression_patterns(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_expression_statements(self, session):
        return

    # Elaborate on get_datatype_data() for the AGMDiseaseHandler.
    def get_datatype_data(self, session):
        """Extend the method for the ExpressionHandler."""
        super().get_datatype_data(session)
        return

    # Add methods to be run by synthesize_info() below.
    # Placeholder.

    # Elaborate on synthesize_info() for the AGMDiseaseHandler.
    def synthesize_info(self):
        """Extend the method for the AGMDiseaseHandler."""
        super().synthesize_info()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    # Placeholder.

    # Elaborate on map_fb_data_to_alliance() for the AGMDiseaseHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AGMDiseaseHandler."""
        super().map_fb_data_to_alliance()
        return
