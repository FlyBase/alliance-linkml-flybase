"""Offline harness for AberrationHandler.add_fbba_to_fbab() and the carries injection (FTA-236).

Extracts method source with ast and drives it with duck-typed stand-ins, since sqlalchemy and
harvdev_utils cannot be imported here (see docs/FTA-236/plan.md).

Run: python3 docs/FTA-236/verify_graft.py
"""

import ast
import re
from os import environ

SRC = '/Users/ilongden/harvard/alliance-linkml-flybase/src/allele_handlers.py'

tree = ast.parse(open(SRC).read())
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'AberrationHandler')
ns = {'re': re, 'getenv': lambda k, d=None: environ.get(k, d)}
wanted_attrs = ('balancer_graft_prop_types', 'balancer_carries_exclusions')
for node in cls.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') in wanted_attrs:
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)
# Later tasks add more methods; the guard lets the harness run before they exist.
for name in ('add_fbba_to_fbab', 'merge_balancers_into_aberrations',
             'synthesize_balancer_carries_associations'):
    node = next((n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name), None)
    if node is not None:
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)


class FakeLog:
    """Collects log calls by level, so the checks can assert on what was reported."""

    def __init__(self):
        """Start with empty lists for each level the handler methods use."""
        self.lines, self.warnings = [], []

    def info(self, msg):
        """Record an info line."""
        self.lines.append(msg)

    def debug(self, msg):
        """Discard debug output; the checks never assert on it."""
        pass

    def warning(self, msg):
        """Record a warning."""
        self.warnings.append(msg)

    def error(self, msg):
        """Record an error alongside the warnings; the checks treat them together."""
        self.warnings.append(msg)


class FakeProp:
    """Stand-in for an FBProp: only .chado_obj.value and .pubs are read."""

    def __init__(self, value):
        """Wrap a featureprop value the way FBProp exposes it."""
        self.chado_obj = type('O', (), {'value': value})()
        self.pubs = []


class FakeSynonym:
    """Stand-in for a FeatureSynonym association, which is grafted by synonym_id."""

    def __init__(self, synonym_id, name):
        """Record the synonym_id the demotion logic keys on, plus its text."""
        self.synonym_id = synonym_id
        self.name = name


class FakeRel:
    """Stand-in for an FBRelationship: subject/object ids plus supporting pubs."""

    def __init__(self, rel_id, subject_id, object_id, pubs=()):
        """Build a relationship between two feature_ids."""
        self.chado_obj = type('O', (), {'subject_id': subject_id, 'object_id': object_id})()
        self.rel_id = rel_id
        self.pubs = list(pubs)


class FakeEntity:
    """Stand-in for an FBAberration or FBBalancer entity."""

    def __init__(self, feature_id, uniquename, name='', rels=None):
        """Build an entity with the source-data lists the graft extends."""
        self.db_primary_id = feature_id
        self.uniquename = uniquename
        self.name = name
        self.synonyms = []
        self.merged_synonym_ids = set()
        self.demoted_synonym_ids = set()
        self.alt_fb_ids = []
        self.fb_sec_dbxrefs = []
        self.pub_associations = []
        self.props_by_type = {}
        self.linkmldto = object()
        self._rels = rels or {}

    def recall_relationships(self, log, **kwargs):
        """Return relationships keyed by (entity_role, rel_types).

        Deliberately ignores rel_entity_types, and records whether it was passed: the nested
        BalancerHandler has no feature_lookup, so its by-related-feature-type buckets are empty and
        that filter would silently match nothing (the 2026-08-14 run moved 0 associations because of
        it). The caller must filter partners against the AberrationHandler's own lookup instead.
        """
        if 'rel_entity_types' in kwargs:
            RECALL_MISUSE.append(kwargs['rel_entity_types'])
        return self._rels.get((kwargs.get('entity_role'), kwargs.get('rel_types')), [])

    def __str__(self):
        """Match FBDataEntity's "name (uniquename)" form, which log messages rely on."""
        return f'{self.name} ({self.uniquename})'


class FakeHandler:
    """Minimal AberrationHandler: the real methods, bound to duck-typed state."""

    balancer_graft_prop_types = ns['balancer_graft_prop_types']
    balancer_carries_exclusions = ns['balancer_carries_exclusions']
    add_fbba_to_fbab = ns.get('add_fbba_to_fbab')
    merge_balancers_into_aberrations = ns.get('merge_balancers_into_aberrations')
    synthesize_balancer_carries_associations = ns.get('synthesize_balancer_carries_associations')

    def __init__(self, aberrations, balancers, merge_map, testing=False):
        """Key the given entities by feature_id and start with empty report/rel collections."""
        self.log = FakeLog()
        self.testing = testing
        self.fb_data_entities = {a.db_primary_id: a for a in aberrations}
        self.fbba_entities = {b.db_primary_id: b for b in balancers}
        self.balancer_merge_map = merge_map
        self.balancer_merge_report = []
        self.aberration_allele_rels = {}
        self.feature_subtypes = {'allele': ['allele'], 'insertion': ['insertion']}
        # The nested BalancerHandler builds no feature_lookup, so partner types must be resolved
        # against THIS handler's lookup (the one FTA-218 builds under the same gate).
        self.feature_lookup = {
            5001: {'type': 'allele', 'name': 'Bar[1]', 'uniquename': 'FBal0000817'},
            5002: {'type': 'allele', 'name': 'sc[8]', 'uniquename': 'FBal0015197'},
            5003: {'type': 'allele', 'name': 'held-back', 'uniquename': 'FBal0000003'},
            5004: {'type': 'allele', 'name': 'cassette', 'uniquename': 'FBal0009999'},
            5005: {'type': 'gene', 'name': 'wrong-type', 'uniquename': 'FBgn0000001'},
            6001: {'type': 'insertion', 'name': 'P{x}1', 'uniquename': 'FBti0000001'},
        }
        self.cassette_ignore_list = {5004}


results = []
RECALL_MISUSE = []


def check(name, condition):
    """Record and print one assertion."""
    results.append((name, condition))
    print(f'  {"PASS" if condition else "FAIL"}: {name}')


print('--- graft: names, identifiers, whitelisted notes, references ---')
parent = FakeEntity(1, 'FBab0004219', 'In(1)Basc')
parent.synonyms = [FakeSynonym(900, 'In(1)Basc')]
parent.props_by_type = {'misc': [FakeProp('Parent comment.')]}
parent.pub_associations = ['parent_pub']

balancer = FakeEntity(101, 'FBba0000014', 'Basc')
balancer.synonyms = [FakeSynonym(910, 'Basc'), FakeSynonym(911, 'M5')]
balancer.fb_sec_dbxrefs = ['fbba_sec_dbxref']
balancer.pub_associations = ['fbba_pub_1', 'fbba_pub_2']
balancer.props_by_type = {
    'misc': [FakeProp('N[opa33b] is a natural variant ...')],
    'internal_notes': [FakeProp('FTA: Balancer - merge with parent In(1)Basc (FBab0004219).')],
    'availability': [FakeProp('Stated to be lost.')],      # must NOT move: would set is_extinct
    'derived_stock_count': [FakeProp('12')],               # must NOT move
}

h = FakeHandler([parent], [balancer], {101: 1})
h.merge_balancers_into_aberrations()

check('balancer synonyms appended', {s.synonym_id for s in parent.synonyms} == {900, 910, 911})
check('balancer synonym ids marked as merged', parent.merged_synonym_ids == {910, 911})
check('parent own synonym not marked merged', 900 not in parent.merged_synonym_ids)
check('balancer synonym ids also demoted (protects the full name slot)', parent.demoted_synonym_ids == {910, 911})
check('balancer primary ID added as secondary', 'FB:FBba0000014' in parent.alt_fb_ids)
check('balancer secondary dbxrefs moved', 'fbba_sec_dbxref' in parent.fb_sec_dbxrefs)
check('misc props merged, parent prop kept', len(parent.props_by_type['misc']) == 2)
check('internal_notes props merged', len(parent.props_by_type['internal_notes']) == 1)
check('availability NOT merged', 'availability' not in parent.props_by_type)
check('derived_stock_count NOT merged', 'derived_stock_count' not in parent.props_by_type)
check('balancer pubs merged', parent.pub_associations == ['parent_pub', 'fbba_pub_1', 'fbba_pub_2'])

row = h.balancer_merge_report[0]
check('report identifies both entities', row['fbba_id'] == 'FBba0000014' and row['fbab_id'] == 'FBab0004219')
check('report counts synonyms', row['synonyms'] == 2)
check('report counts comments', row['comments'] == 1)
check('report counts references', row['references'] == 2)
check('report counts secondary ids (primary + own)', row['secondary_ids'] == 2)

print('--- two balancers into one parent (the SM6 case) ---')
parent2 = FakeEntity(2, 'FBab0004818', 'In(2LR)SM6')
sm6a = FakeEntity(102, 'FBba0000039', 'SM6a')
sm6a.synonyms = [FakeSynonym(920, 'SM6a')]
sm6b = FakeEntity(103, 'FBba0000040', 'SM6b')
sm6b.synonyms = [FakeSynonym(930, 'SM6b')]
h2 = FakeHandler([parent2], [sm6a, sm6b], {102: 2, 103: 2})
h2.merge_balancers_into_aberrations()
check('both balancers merged into one parent', {s.synonym_id for s in parent2.synonyms} == {920, 930})
check('two report rows written', len(h2.balancer_merge_report) == 2)

print('--- carries associations: injected under the parent, exclusions honoured ---')
basc = FakeEntity(101, 'FBba0000014', 'Basc', rels={
    ('object', 'carried_on'): [FakeRel(1, 5001, 101), FakeRel(2, 5002, 101),
                               FakeRel(5, 5004, 101),    # cassette FBal: must be skipped with a warning
                               FakeRel(6, 5005, 101),    # wrong partner type: must be skipped
                               FakeRel(7, 5999, 101)],   # not in feature_lookup at all: must be skipped
    ('subject', 'associated_with'): [FakeRel(3, 101, 6001)],
})
fm1 = FakeEntity(111, 'FBba0000011', 'FM1', rels={
    ('object', 'carried_on'): [FakeRel(4, 5003, 111)],
})
parent_basc = FakeEntity(1, 'FBab0004219', 'In(1)Basc')
parent_fm1 = FakeEntity(10, 'FBab0010486', 'In(1)FM1')

h3 = FakeHandler([parent_basc, parent_fm1], [basc, fm1], {101: 1, 111: 10})
h3.merge_balancers_into_aberrations()
h3.synthesize_balancer_carries_associations()

keys = set(h3.aberration_allele_rels.keys())
check('carried_on allele keyed under the parent as "carries"', (1, 5001, 'carries') in keys)
check('second carried_on allele present', (1, 5002, 'carries') in keys)
check('associated_with insertion keyed as "carries"', (1, 6001, 'carries') in keys)
check('excluded balancer contributes nothing', not any(k[0] == 10 for k in keys))
check('nothing keyed under the balancer itself', not any(k[0] in (101, 111) for k in keys))
basc_row = next(r for r in h3.balancer_merge_report if r['fbba_id'] == 'FBba0000014')
fm1_row = next(r for r in h3.balancer_merge_report if r['fbba_id'] == 'FBba0000011')
check('report counts moved associations', basc_row['carries_alleles'] == 3)
check('cassette partner skipped with a warning', any('cassette' in w.lower() for w in h3.log.warnings))
check('wrong-type partner not moved', (1, 5005, 'carries') not in keys)
check('partner missing from feature_lookup not moved', (1, 5999, 'carries') not in keys)
check('report counts excluded associations', fm1_row['carries_excluded'] == 1)
check('excluded balancer moved nothing', fm1_row['carries_alleles'] == 0)

print('--- the relationship objects themselves are kept, for evidence lookup ---')
rel_list = h3.aberration_allele_rels[(1, 5001, 'carries')]
check('rel_key maps to the FBRelationship list', rel_list[0].rel_id == 1)

check('rel_entity_types not passed to recall_relationships (empty buckets in the nested handler)',
      RECALL_MISUSE == [])

print(f'\n{sum(1 for _, ok in results if ok)}/{len(results)} checks PASSED')
