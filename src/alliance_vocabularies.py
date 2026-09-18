"""Module:: alliance_vocabularies.

Synopsis:
    Checks values that the Alliance validates against a controlled vocabulary, so that a term
    FlyBase has but the Alliance does not is caught here rather than by the load.

    The 2026_03 allele load stopped at its 1,000-failure cap having completed 157,680 of
    617,375 records. Every exception read:

        in_collection_name - Not a valid entry

    and every one carried the same value, "GV_FFL" - FBlc0008106, a current FlyBase reagent
    collection that is simply not a term in the Alliance's "allele_collection" vocabulary.
    3,008 records in that file carry it, so one missing term aborted the whole submission
    (SCRUM-6568 asks the Alliance to add it; FTA-260 has the submission context).

    The vocabulary is read from the A-Team curation API, which has no endpoint that filters
    terms by vocabulary - POST /api/vocabularyterm/find ignores a vocabulary filter and returns
    every term in the database - so this module takes the two-step route that does work:

        POST /api/vocabulary/find     -> all 48 vocabularies, matched on vocabularyLabel
        GET  /api/vocabulary/{id}     -> that vocabulary with its memberTerms

    Authentication is shared with resource_descriptors: a Cognito access token from the
    client_credentials exchange. See that module for the variables involved.

    When the vocabulary cannot be fetched, values pass through unchecked. That is deliberate:
    an unverifiable value is not evidence of a bad value, and silently dropping curated data
    on the strength of a failed API call would be worse than submitting it.

Author(s):
    Ian Longden ianlongden@morgan.harvard.edu

"""

from collections import Counter
from logging import Logger
from os import getenv

import requests

from resource_descriptors import (
    DEFAULT_BASE_URL, REQUEST_TIMEOUT, get_cognito_access_token
)

# Vocabulary labels this module knows about, keyed by the LinkML slot they validate.
ALLELE_COLLECTION_VOCABULARY = 'allele_collection'

VOCABULARY_FIND_PATH = 'api/vocabulary/find'
VOCABULARY_PATH = 'api/vocabulary'
# Comfortably above the 48 vocabularies that exist, since the endpoint cannot filter.
VOCABULARY_LIMIT = 200


def fetch_vocabulary_terms(log: Logger, vocabulary_label: str, token: str = None, base_url: str = None):
    """Fetch the non-obsolete term names of one Alliance vocabulary.

    Args:
        log (Logger): The handler's logger.
        vocabulary_label (str): The Alliance vocabularyLabel: e.g. "allele_collection".
        token (str): Bearer token to use as-is; omit to get one by Cognito exchange.
        base_url (str): API host; defaults to AGR_BASE_URL, then the public curation site.

    Returns:
        A set of term names, or None if the vocabulary could not be fetched. None means
        "unknown", and callers must let values through rather than treat them as invalid.

    """
    base = (base_url or getenv('AGR_BASE_URL', DEFAULT_BASE_URL)).rstrip('/')
    if not token:
        token = get_cognito_access_token(log) or getenv('ALLIANCETOKEN', None)
    if not token:
        log.warning(f'No curation API token available; cannot fetch the "{vocabulary_label}" vocabulary.')
        return None
    headers = {
        'Content-Type': 'application/json',
        'accept': 'application/json',
        'Authorization': f'Bearer {token}',
    }
    try:
        response = requests.post(f'{base}/{VOCABULARY_FIND_PATH}?limit={VOCABULARY_LIMIT}&page=0',
                                 headers=headers, json={}, timeout=REQUEST_TIMEOUT)
        if response.status_code != 200:
            log.error(f'Alliance vocabulary list returned {response.status_code}.')
            return None
        vocabulary_id = None
        for vocabulary in response.json().get('results', None) or []:
            if vocabulary.get('vocabularyLabel', None) == vocabulary_label:
                vocabulary_id = vocabulary.get('id', None)
                break
        if vocabulary_id is None:
            log.error(f'The Alliance has no vocabulary labelled "{vocabulary_label}".')
            return None
        response = requests.get(f'{base}/{VOCABULARY_PATH}/{vocabulary_id}',
                                headers=headers, timeout=REQUEST_TIMEOUT)
        if response.status_code != 200:
            log.error(f'Alliance vocabulary {vocabulary_id} returned {response.status_code}.')
            return None
        entity = response.json().get('entity', None) or {}
    except requests.RequestException as error:
        log.error(f'Could not reach the Alliance vocabulary API: {error}')
        return None
    except ValueError as error:
        log.error(f'Alliance vocabulary API returned unparseable JSON: {error}')
        return None
    terms = {term.get('name') for term in (entity.get('memberTerms', None) or [])
             if term.get('name') and term.get('obsolete', False) is not True}
    if not terms:
        log.error(f'The Alliance "{vocabulary_label}" vocabulary came back with no terms.')
        return None
    log.info(f'Fetched {len(terms)} terms of the Alliance "{vocabulary_label}" vocabulary.')
    return terms


class VocabularyGuard(object):
    """Checks slot values against Alliance vocabularies, fetching each one once.

    A guard with no API access checks nothing and says so in its report, so an export run
    without credentials submits exactly what it did before rather than dropping values it
    cannot verify.

    """
    def __init__(self, log: Logger):
        """Create a VocabularyGuard.

        Args:
            log (Logger): The handler's logger.

        """
        self.log = log
        self.terms_by_vocabulary = {}   # {vocabulary_label: set of term names, or None if unavailable}
        self.rejected = Counter()       # {(vocabulary_label, value): count} values dropped
        self.accepted = Counter()       # {vocabulary_label: count} values kept
        self.unchecked = Counter()      # {vocabulary_label: count} values passed through unverified

    def terms(self, vocabulary_label: str):
        """Return the term set for a vocabulary, fetching it on first use."""
        if vocabulary_label not in self.terms_by_vocabulary:
            self.terms_by_vocabulary[vocabulary_label] = fetch_vocabulary_terms(self.log, vocabulary_label)
        return self.terms_by_vocabulary[vocabulary_label]

    def check(self, vocabulary_label: str, value: str):
        """Return the value if the Alliance has it as a term, otherwise None.

        Args:
            vocabulary_label (str): The Alliance vocabularyLabel to check against.
            value (str): The value FlyBase would submit.

        Returns:
            The value unchanged when it is a term in that vocabulary, or when the vocabulary
            could not be fetched; None when the vocabulary is known and the value is not in
            it, meaning the caller must leave the slot empty rather than fail the record.

        """
        if value is None:
            return None
        terms = self.terms(vocabulary_label)
        if terms is None:
            self.unchecked[vocabulary_label] += 1
            return value
        if value in terms:
            self.accepted[vocabulary_label] += 1
            return value
        self.rejected[(vocabulary_label, value)] += 1
        return None

    def report(self):
        """Log what was dropped, so a curator can see it before the file is submitted."""
        for vocabulary_label, count in self.unchecked.items():
            self.log.warning(f'Submitted {count} "{vocabulary_label}" values unchecked: the Alliance '
                             'vocabulary could not be fetched. A value the Alliance does not have '
                             'will fail the load.')
        if self.rejected:
            total = sum(self.rejected.values())
            self.log.warning(f'Dropped {total} values absent from their Alliance vocabulary. Ask the '
                             'Alliance to add the term, or correct it in FlyBase:')
            for (vocabulary_label, value), count in self.rejected.most_common():
                self.log.warning(f'   {vocabulary_label}: "{value}" on {count} entities.')
        for vocabulary_label, count in self.accepted.items():
            self.log.info(f'Checked {count} "{vocabulary_label}" values against the Alliance vocabulary.')
        return
