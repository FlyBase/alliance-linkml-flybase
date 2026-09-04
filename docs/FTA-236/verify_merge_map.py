"""Offline harness for AberrationHandler.synthesize_balancer_merge_map() (FTA-236).

Extracts the method source with ast and drives it with duck-typed stand-ins, since
sqlalchemy/harvdev_utils cannot be imported here (see docs/FTA-236/plan.md).

Run: python3 docs/FTA-236/verify_merge_map.py
"""

import ast
import re

SRC = '/Users/ilongden/harvard/alliance-linkml-flybase/src/allele_handlers.py'

tree = ast.parse(open(SRC).read())
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'AberrationHandler')
ns = {'re': re}
for node in cls.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') == 'balancer_merge_regex':
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)
for name in ('synthesize_balancer_merge_map', '_resolve_parent_feature_id'):
    node = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)
for node in tree.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') == 'BALANCER_MERGE_EXPECTED':
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)
print(f'Regex in use: {ns["balancer_merge_regex"].pattern}')
print(f'Expected flag count: {ns["BALANCER_MERGE_EXPECTED"]}\n')


class FakeLog:
    """Collects log calls by level, so the checks can assert on what was reported."""

    def __init__(self):
        """Start with empty lists for each level the handler methods use."""
        self.info_lines, self.warnings, self.errors = [], [], []

    def info(self, msg):
        """Record an info line."""
        self.info_lines.append(msg)

    def debug(self, msg):
        """Discard debug output; the checks never assert on it."""
        pass

    def warning(self, msg):
        """Record a warning."""
        self.warnings.append(msg)

    def error(self, msg):
        """Record an error."""
        self.errors.append(msg)


class FakeProp:
    """Stand-in for an FBProp: only .chado_obj.value is read by the code under test."""

    def __init__(self, value):
        """Wrap a featureprop value the way FBProp exposes it."""
        self.chado_obj = type('O', (), {'value': value})()


class FakeEntity:
    """Stand-in for an FBAberration or FBBalancer entity."""

    def __init__(self, feature_id, uniquename, name='', internal_notes=(), is_obsolete=False):
        """Build an entity with optional internal_notes props; obsolete entities are exported too."""
        self.db_primary_id = feature_id
        self.uniquename = uniquename
        self.name = name
        self.is_obsolete = is_obsolete
        self.props_by_type = {'internal_notes': [FakeProp(v) for v in internal_notes]} if internal_notes else {}

    def __str__(self):
        """Match FBDataEntity's "name (uniquename)" form, which log messages rely on."""
        return f'{self.name} ({self.uniquename})'


class FakeHandler:
    """Minimal AberrationHandler: the real methods, bound to duck-typed state."""

    balancer_merge_regex = ns['balancer_merge_regex']
    synthesize_balancer_merge_map = ns['synthesize_balancer_merge_map']
    _resolve_parent_feature_id = ns['_resolve_parent_feature_id']

    def __init__(self, aberrations, balancers, testing=False):
        """Key the given aberrations and balancers by feature_id, as the handler does."""
        self.log = FakeLog()
        self.testing = testing
        self.fb_data_entities = {a.db_primary_id: a for a in aberrations}
        self.fbba_entities = {b.db_primary_id: b for b in balancers}
        self.balancer_merge_map = {}


# The real chado values, verified 2026-08-14.
FLAG = 'FTA: Balancer - merge with parent In(1)Basc (FBab0004219).'
FLAG_SM6A = 'FTA: Balancer - merge with parent In(2LR)SM6 (FBab0004818).'
FLAG_SM6B = 'FTA: Balancer - merge with parent In(2LR)SM6 (FBab0004818).'
# AM1's note before the week-of-2026-08-17 correction: names the obsolete FBab0007127.
FLAG_OBSOLETE = 'FTA: Balancer - merge with parent T(2;3)A1-W (FBab0007127).'
FLAG_UNKNOWN = 'FTA: Balancer - merge with parent Nope (FBab9999999).'
# Sibling flag families that must not be picked up here.
FLAG_RENAME = 'FTA: Balancer - use balancer symbol and fullname for parent In(1)Basc (FBab0004219).'
MARK = "FTA: Balancer - mark this aberration as 'balancer'."
PROSE = 'Used as a balancer; merge with parent was discussed but never agreed.'

results = []


def check(name, condition):
    """Record and print one assertion."""
    results.append((name, condition))
    print(f'  {"PASS" if condition else "FAIL"}: {name}')


print('--- happy path: 3 flagged balancers, 2 parents (SM6 takes two) ---')
aberrations = [
    FakeEntity(1, 'FBab0004219', 'In(1)Basc'),
    FakeEntity(2, 'FBab0004818', 'In(2LR)SM6'),
    FakeEntity(3, 'FBab0049550', 'T(2;3)A1-W'),
]
balancers = [
    FakeEntity(101, 'FBba0000014', 'Basc', [FLAG]),
    FakeEntity(102, 'FBba0000039', 'SM6a', [FLAG_SM6A]),
    FakeEntity(103, 'FBba0000040', 'SM6b', [FLAG_SM6B]),
    FakeEntity(104, 'FBba0000999', 'X', [PROSE, FLAG_RENAME, MARK]),
]
h = FakeHandler(aberrations, balancers)
h.synthesize_balancer_merge_map()
check('3 flagged balancers mapped', h.balancer_merge_map == {101: 1, 102: 2, 103: 2})
check('rename/mark/prose balancer ignored', 104 not in h.balancer_merge_map)
check('no errors logged', h.log.errors == [])

print('--- obsolete parent is an error, not a silent redirect (Steven: option 1) ---')
h = FakeHandler([FakeEntity(3, 'FBab0049550', 'T(2;3)A1-W')],
                [FakeEntity(105, 'FBba0000688', 'AM1', [FLAG_OBSOLETE])])
h.synthesize_balancer_merge_map()
check('obsolete/absent parent not merged', h.balancer_merge_map == {})
check('error names the balancer', any('FBba0000688' in e for e in h.log.errors))
check('error names the stale parent ID', any('FBab0007127' in e for e in h.log.errors))
check('current T(2;3)A1-W is NOT silently substituted', 3 not in h.balancer_merge_map.values())

print('--- obsolete parent PRESENT in fb_data_entities is still an error (2026-08-14 run bug) ---')
# fb_data_entities holds obsolete aberrations too - they export as internal/obsolete - so membership
# alone is not enough. Before the fix, AM1's data merged into the obsolete FBab0007127.
h = FakeHandler([FakeEntity(7, 'FBab0007127', 'T(2;3)A1-W', is_obsolete=True),
                 FakeEntity(3, 'FBab0049550', 'T(2;3)A1-W')],
                [FakeEntity(105, 'FBba0000688', 'AM1', [FLAG_OBSOLETE])])
h.synthesize_balancer_merge_map()
check('obsolete parent in fb_data_entities is rejected', h.balancer_merge_map == {})
check('error raised for the obsolete parent', any('FBab0007127' in e for e in h.log.errors))
check('current same-named aberration not substituted', 3 not in h.balancer_merge_map.values())

print('--- unknown parent ID is an error ---')
h = FakeHandler([FakeEntity(1, 'FBab0004219', 'In(1)Basc')],
                [FakeEntity(106, 'FBba0000777', 'Y', [FLAG_UNKNOWN])])
h.synthesize_balancer_merge_map()
check('unknown parent not merged', h.balancer_merge_map == {})
check('error names the unknown ID', any('FBab9999999' in e for e in h.log.errors))

print('--- count warning off the expected 38, suppressed in testing mode ---')
check('count warning present', any('38' in w for w in h.log.warnings))
h = FakeHandler([FakeEntity(1, 'FBab0004219', 'In(1)Basc')],
                [FakeEntity(101, 'FBba0000014', 'Basc', [FLAG])], testing=True)
h.synthesize_balancer_merge_map()
check('no count warning in testing mode', not any('38' in w for w in h.log.warnings))

print(f'\n{sum(1 for _, ok in results if ok)}/{len(results)} checks PASSED')
