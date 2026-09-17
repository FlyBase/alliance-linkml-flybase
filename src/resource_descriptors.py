"""Module:: resource_descriptors.

Synopsis:
    Resolves the "page_area" of an Alliance CrossReferenceDTO against the Alliance's own
    resource descriptors, fetched from the A-Team curation API.

    A page_area is only valid for a given prefix if that prefix's resource descriptor
    declares a resource page of that name. FlyBase has MOD pages per data type ("gene",
    "gene/expression", "allele", ...), but external resources generally do not: they have
    one landing page, addressed by the synthetic page area "default", which exists for
    every prefix carrying a default URL template. Deriving the page area from the FlyBase
    data type therefore produced page areas that no external prefix declares, and the
    Alliance rejected every affected record with:

        cross_reference_dtos - page_area - Not a valid entry (gene)

    That is FTA-263: it stopped the fb_2026_03 gene load after ~1,116 of 256,037 records,
    with 95% of what it processed rejected.

    The old resourceDescriptors.yaml is obsolete. Current data comes from the A-Team API
    endpoint that the shared agr_curation_api_client uses
    (src/agr_curation_api/api_methods.py, find_resource_descriptors_for_public):

        POST /api/resourcedescriptor/findForPublic?limit=5000&page=0&view=ResourceDescriptorView

    The ResourceDescriptorView is needed because the default ForPublic view omits
    resourcePages and synonyms.

    Authentication is NOT the ALLIANCETOKEN used for uploads. Verified 2026-09-17: the
    "Curation API Token" UUID from the curation site profile returns 401 with
    "www-authenticate: Bearer" on every endpoint, as does the "Cognito Id Token". The API
    validates AWS Cognito access tokens, so this module gets one the same way the shared
    agr_curation_api_client does, via agr_cognito_py: an OAuth client_credentials exchange
    against COGNITO_TOKEN_URL using COGNITO_ADMIN_CLIENT_ID/COGNITO_ADMIN_CLIENT_SECRET and
    COGNITO_ADMIN_SCOPE. The exchange is reimplemented here rather than adding the
    dependency, which pulls in linkml, sqlalchemy 2.x, elasticsearch and fastapi - sqlalchemy
    2.x conflicts with this repo's pin of <2.0. ALLIANCETOKEN is still honoured as a
    fallback, for a hand-exported Cognito access token during development.

Author(s):
    Ian Longden ianlongden@morgan.harvard.edu

"""

from collections import Counter
from logging import Logger
from os import getenv

import requests

# The synthetic page area every prefix with a default URL template resolves to.
DEFAULT_PAGE_AREA = 'default'

# The A-Team curation API endpoint and the view that includes resourcePages/synonyms.
RESOURCE_DESCRIPTOR_PATH = 'api/resourcedescriptor/findForPublic'
RESOURCE_DESCRIPTOR_VIEW = 'ResourceDescriptorView'
DEFAULT_BASE_URL = 'https://curation.alliancegenome.org'
# High enough to return every descriptor in one page: the endpoint only populates
# totalResults in count-only mode, so paging cannot be driven off the response.
RESOURCE_DESCRIPTOR_LIMIT = 5000
REQUEST_TIMEOUT = 60

# Environment variables for the Cognito client_credentials exchange, named as agr_cognito_py
# names them so that one set of GoCD variables serves both.
COGNITO_ENV_VARS = ('COGNITO_ADMIN_CLIENT_ID', 'COGNITO_ADMIN_CLIENT_SECRET', 'COGNITO_TOKEN_URL')


def get_cognito_access_token(log: Logger):
    """Get a Cognito access token for the curation API by client_credentials exchange.

    Mirrors agr_cognito_py.get_admin_token(): HTTP Basic with the admin client id and
    secret, form-encoded grant_type=client_credentials plus the scope, against
    COGNITO_TOKEN_URL.

    Args:
        log (Logger): The handler's logger.

    Returns:
        The access token string, or None if the credentials are absent or the exchange
        fails. Never logs the credentials or the token.

    """
    missing = [name for name in COGNITO_ENV_VARS if not getenv(name, None)]
    if missing:
        log.warning(f'Cannot get a Cognito token; these variables are not set: {", ".join(missing)}.')
        return None
    data = {'grant_type': 'client_credentials'}
    scope = getenv('COGNITO_ADMIN_SCOPE', None)
    if scope:
        data['scope'] = scope
    log.info('Requesting a Cognito access token for the curation API.')
    try:
        response = requests.post(
            getenv('COGNITO_TOKEN_URL'),
            auth=(getenv('COGNITO_ADMIN_CLIENT_ID'), getenv('COGNITO_ADMIN_CLIENT_SECRET')),
            headers={'Content-Type': 'application/x-www-form-urlencoded'},
            data=data,
            timeout=REQUEST_TIMEOUT,
        )
    except requests.RequestException as error:
        log.error(f'Cognito token request failed: {error}')
        return None
    if response.status_code != 200:
        # Deliberately not logging the body: a failed token response can echo the request.
        log.error(f'Cognito token request returned {response.status_code}.')
        return None
    try:
        access_token = response.json()['access_token']
    except (ValueError, KeyError):
        log.error('Cognito token response contained no access_token.')
        return None
    log.info('Obtained a Cognito access token.')
    return access_token


def fetch_resource_descriptors(log: Logger, token: str = None, base_url: str = None):
    """Fetch Alliance resource descriptors from the A-Team curation API.

    Args:
        log (Logger): The handler's logger.
        token (str): Bearer token to use as-is. Omit to get one by Cognito
            client_credentials exchange, falling back to the ALLIANCETOKEN variable.
        base_url (str): API host; defaults to the AGR_BASE_URL environment variable,
            then to the public curation site.

    Returns:
        A list of raw resource descriptor dicts, or None if they could not be fetched.
        None is distinct from an empty list: it means "unknown", and the caller falls
        back to conservative page areas rather than treating every prefix as unknown.

    """
    base = (base_url or getenv('AGR_BASE_URL', DEFAULT_BASE_URL)).rstrip('/')
    if not token:
        token = get_cognito_access_token(log)
    if not token:
        token = getenv('ALLIANCETOKEN', None)
        if token:
            log.warning('No Cognito credentials; trying ALLIANCETOKEN. Note the upload token is '
                        'NOT accepted by this API - only a Cognito access token is.')
    if not token:
        log.warning('No curation API token available; cannot fetch Alliance resource descriptors.')
        return None
    url = f'{base}/{RESOURCE_DESCRIPTOR_PATH}?limit={RESOURCE_DESCRIPTOR_LIMIT}&page=0&view={RESOURCE_DESCRIPTOR_VIEW}'
    headers = {
        'Content-Type': 'application/json',
        'accept': 'application/json',
        'Authorization': f'Bearer {token}',
    }
    log.info(f'Fetch Alliance resource descriptors from {base}/{RESOURCE_DESCRIPTOR_PATH}.')
    try:
        response = requests.post(url, headers=headers, json={}, timeout=REQUEST_TIMEOUT)
    except requests.RequestException as error:
        log.error(f'Could not reach the Alliance resource descriptor API: {error}')
        return None
    if response.status_code != 200:
        log.error(f'Alliance resource descriptor API returned {response.status_code}: {response.text[:200]}')
        return None
    try:
        results = response.json().get('results', None)
    except ValueError as error:
        log.error(f'Alliance resource descriptor API returned unparseable JSON: {error}')
        return None
    if not results:
        log.error('Alliance resource descriptor API returned no descriptors.')
        return None
    log.info(f'Fetched {len(results)} Alliance resource descriptors.')
    return results


class PageAreaResolver(object):
    """Resolves cross-reference page areas against Alliance resource descriptors.

    Built from API data, the resolver knows each prefix's declared resource pages and
    whether it has a default URL template, so it can answer both "is this page area valid
    for this prefix?" and "is this prefix known at all?".

    Built without API data (no token, or the API is unreachable), it falls back to the
    rule the FTA-263 evidence supports: FlyBase keeps its data type pages, every other
    prefix gets "default". It then reports no prefix as unknown, because in that state it
    cannot tell an unknown prefix from one it simply has no data for.

    """
    def __init__(self, log: Logger, descriptors: list = None):
        """Create a PageAreaResolver.

        Args:
            log (Logger): The handler's logger.
            descriptors (list): Raw descriptor dicts from fetch_resource_descriptors(),
                or None to build a fallback resolver.

        """
        self.log = log
        self.api_backed = descriptors is not None
        self.pages_by_prefix = {}        # {prefix: set of non-obsolete resource page names}
        self.has_default_by_prefix = {}  # {prefix: bool, whether a default URL template exists}
        self.synonym_to_prefix = {}      # {lowercased prefix or synonym: canonical prefix}
        self.substitutions = Counter()   # {(prefix, requested page area): count} swapped for "default"
        self.dropped_xrefs = Counter()   # {prefix: count} xrefs dropped for an unknown prefix
        self.prefix_corrections = Counter()  # {(submitted, canonical): count} prefixes respelled
        self.unresolved = Counter()      # {(prefix, requested page area): count} with no valid page area
        if descriptors is not None:
            self._index_descriptors(descriptors)

    def _index_descriptors(self, descriptors: list):
        """Index raw descriptor dicts by prefix."""
        for descriptor in descriptors:
            prefix = descriptor.get('prefix', None)
            if not prefix:
                continue
            pages = set()
            for page in descriptor.get('resourcePages', None) or []:
                if page.get('obsolete', False) is True:
                    continue
                page_name = page.get('name', None)
                if page_name:
                    pages.add(page_name)
            self.pages_by_prefix[prefix] = pages
            self.has_default_by_prefix[prefix] = bool(descriptor.get('defaultUrlTemplate', None))
            self.synonym_to_prefix[prefix.lower()] = prefix
            for synonym in descriptor.get('synonyms', None) or []:
                self.synonym_to_prefix.setdefault(synonym.lower(), prefix)
        self.log.info(f'Indexed resource pages for {len(self.pages_by_prefix)} Alliance prefixes.')

    def canonical_prefix(self, prefix: str):
        """Return the prefix spelled as the Alliance spells it, or None if unrecognized.

        The Alliance validates the literal prefix string, so a prefix that differs only in
        case or that is a declared synonym is rejected as submitted even though the resource
        exists: "dgrc" is rejected where "DGRC" loads. Returning the canonical spelling lets
        the caller correct the cross-reference instead of dropping it.

        A fallback resolver echoes the prefix back unchanged: with no descriptor data it can
        neither correct a spelling nor judge one unrecognized.

        """
        if not self.api_backed:
            return prefix
        if prefix in self.pages_by_prefix:
            return prefix
        canonical = self.synonym_to_prefix.get(prefix.lower(), None)
        if canonical is not None and canonical != prefix:
            self.prefix_corrections[(prefix, canonical)] += 1
        return canonical

    def knows_prefix(self, prefix: str):
        """Return True if the Alliance recognizes the prefix, by name, case or synonym.

        Always True for a fallback resolver, which has no basis on which to call a prefix
        unknown; callers must not drop cross-references on the strength of that answer.

        """
        return self.canonical_prefix(prefix) is not None

    def resolve(self, prefix: str, page_area: str):
        """Return a page area valid for the prefix, or None if none is.

        Args:
            prefix (str): The Alliance prefix of the cross-reference: e.g., "FB", "interpro".
            page_area (str): The page area the caller would like: usually a FlyBase data
                type page such as "gene" or "gene/expression".

        Returns:
            The requested page area if the prefix declares it; otherwise DEFAULT_PAGE_AREA
            if the prefix has a default URL template; otherwise None, meaning the
            cross-reference cannot be given a valid page area and should be dropped.

        """
        if not self.api_backed:
            # Evidence-based fallback: only FlyBase has per-data-type pages.
            if prefix == 'FB':
                return page_area
            return DEFAULT_PAGE_AREA
        if page_area in self.pages_by_prefix.get(prefix, set()):
            return page_area
        canonical = self.synonym_to_prefix.get(prefix.lower(), prefix)
        if page_area in self.pages_by_prefix.get(canonical, set()):
            return page_area
        if self.has_default_by_prefix.get(canonical, False):
            if page_area != DEFAULT_PAGE_AREA:
                self.substitutions[(prefix, page_area)] += 1
            return DEFAULT_PAGE_AREA
        self.unresolved[(prefix, page_area)] += 1
        return None

    def note_dropped_xref(self, prefix: str):
        """Record that a cross-reference was dropped because its prefix is unknown."""
        self.dropped_xrefs[prefix] += 1

    def report(self):
        """Log what the resolver changed, so a run can be checked before submission."""
        if self.api_backed is False:
            self.log.warning('Page areas resolved WITHOUT Alliance resource descriptors: '
                             'FlyBase page areas kept, all other prefixes set to '
                             f'"{DEFAULT_PAGE_AREA}". Unknown prefixes were NOT detected.')
        if self.substitutions:
            total = sum(self.substitutions.values())
            self.log.info(f'Substituted "{DEFAULT_PAGE_AREA}" for {total} cross-reference page areas '
                          f'not declared by their prefix ({len(self.substitutions)} prefix/page_area combinations).')
            for (prefix, page_area), count in self.substitutions.most_common():
                self.log.debug(f'PAGE_AREA: {prefix} does not declare "{page_area}"; used '
                               f'"{DEFAULT_PAGE_AREA}" for {count} cross-references.')
        if self.prefix_corrections:
            total = sum(self.prefix_corrections.values())
            self.log.info(f'Corrected the spelling of {total} cross-reference prefixes to match the Alliance: '
                          f'{ {f"{sub} -> {canon}": n for (sub, canon), n in self.prefix_corrections.items()} }.')
        if self.dropped_xrefs:
            total = sum(self.dropped_xrefs.values())
            self.log.warning(f'Dropped {total} cross-references whose prefix the Alliance does not recognize: '
                             f'{dict(self.dropped_xrefs)}. Fix the prefix in fb_agr_db_dict or remove the mapping.')
        if self.unresolved:
            total = sum(self.unresolved.values())
            self.log.warning(f'Dropped {total} cross-references with no valid page area, because the prefix has '
                             f'neither the requested page nor a default URL: {dict(self.unresolved)}.')
        return
