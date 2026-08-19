"""Offline harness for the FTA-237 balancer rename methods.

Extracts method source with ast and drives it with duck-typed stand-ins, since sqlalchemy and
harvdev_utils cannot be imported on this machine (see docs/FTA-237/plan.md).

Run: python3 docs/FTA-237/verify_renames.py
"""

import ast
import re

SRC = '/Users/ilongden/harvard/alliance-linkml-flybase/src/allele_handlers.py'

tree = ast.parse(open(SRC).read())
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'AberrationHandler')


class FakeNameSlotAnnotationDTO:
    """Stand-in for agr_datatypes.NameSlotAnnotationDTO (Task 3's hard-coded path builds these)."""

    def __init__(self, name_type_name, format_text, display_text, evidence_curies):
        """Build the stand-in with the same four positional args as the real DTO."""
        self.dto = {'name_type_name': name_type_name, 'format_text': format_text,
                    'display_text': display_text, 'synonym_scope_name': 'exact',
                    'evidence_curies': list(evidence_curies), 'internal': False, 'obsolete': False}

    def dict_export(self):
        """Return the DTO as a plain dict, as the real class does."""
        return dict(self.dto)


ns = {
    're': re,
    'sub_sup_sgml_to_html': lambda x: x,
    'sub_sup_to_sgml': lambda x: x,
    'agr_datatypes': type('agr_datatypes', (), {'NameSlotAnnotationDTO': FakeNameSlotAnnotationDTO}),
}
wanted_attrs = ('balancer_rename_regex', 'balancer_hardcoded_renames')
for node in cls.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') in wanted_attrs:
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)
# Later tasks add more methods; the guard below lets the harness run before they exist.
for name in ('synthesize_balancer_rename_map', 'map_balancer_renames', 'drop_duplicate_synonyms',
             'apply_rename', 'map_hardcoded_balancer_renames'):
    node = next((n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name), None)
    if node is not None:
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)

# The module-level expected count is read by synthesize_balancer_rename_map().
for node in tree.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') == 'BALANCER_RENAME_EXPECTED':
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)


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


def name_dto(name_type, text):
    return {'name_type_name': name_type, 'format_text': text, 'display_text': text,
            'synonym_scope_name': 'exact', 'evidence_curies': [], 'internal': False, 'obsolete': False}


class FakeDTO:
    """Stand-in for AlleleDTO, holding just the three name slots the rename touches."""

    def __init__(self, symbol=None, full_name=None, synonyms=None):
        """Set the symbol and full-name slots, and copy the synonym list."""
        self.allele_symbol_dto = symbol
        self.allele_full_name_dto = full_name
        self.allele_synonym_dtos = list(synonyms or [])


class FakeEntity:
    """Stand-in for an FBAberration or FBBalancer entity."""

    def __init__(self, feature_id, uniquename, name='', internal_notes=(), dto=None):
        """Build an entity with optional internal_notes props and an optional LinkML DTO."""
        self.db_primary_id = feature_id
        self.uniquename = uniquename
        self.name = name
        self.is_obsolete = False
        self.props_by_type = {'internal_notes': [FakeProp(v) for v in internal_notes]} if internal_notes else {}
        self.linkmldto = dto

    def __str__(self):
        """Match FBDataEntity's "name (uniquename)" form, which log messages rely on."""
        return f'{self.name} ({self.uniquename})'


class FakeHandler:
    """Minimal AberrationHandler: the real methods, bound to duck-typed state."""

    balancer_rename_regex = ns['balancer_rename_regex']
    balancer_hardcoded_renames = ns.get('balancer_hardcoded_renames', {})
    synthesize_balancer_rename_map = ns['synthesize_balancer_rename_map']
    map_balancer_renames = ns.get('map_balancer_renames')
    drop_duplicate_synonyms = ns.get('drop_duplicate_synonyms')
    apply_rename = ns.get('apply_rename')
    map_hardcoded_balancer_renames = ns.get('map_hardcoded_balancer_renames')

    def __init__(self, aberrations, balancers, testing=False):
        """Key the given aberrations and balancers by feature_id, as the handler does."""
        self.log = FakeLog()
        self.testing = testing
        self.fb_data_entities = {a.db_primary_id: a for a in aberrations}
        self.fbba_entities = {b.db_primary_id: b for b in balancers}
        self.balancer_rename_map = {}
        self.balancer_rename_report = []

    def _resolve_parent_feature_id(self, fbab_uniquename):
        """Mirror the real helper: find an exportable parent by uniquename, else None."""
        for aberration in self.fb_data_entities.values():
            if aberration.uniquename == fbab_uniquename:
                return aberration.db_primary_id
        return None


RENAME = 'FTA: Balancer - use balancer symbol and fullname for parent In(1)Basc (FBab0004219).'
MERGE = 'FTA: Balancer - merge with parent In(1)Basc (FBab0004219).'
MARK = "FTA: Balancer - mark this aberration as 'balancer'."
RENAME_UNKNOWN = 'FTA: Balancer - use balancer symbol and fullname for parent Nope (FBab9999999).'

results = []


def check(name, condition):
    results.append((name, condition))
    print(f'  {"PASS" if condition else "FAIL"}: {name}')


print('--- flag family discrimination ---')
basc_parent = FakeEntity(1, 'FBab0004219', 'In(1)Basc')
basc = FakeEntity(101, 'FBba0000014', 'Basc', [MERGE, RENAME])
merge_only = FakeEntity(102, 'FBba0000039', 'SM6a', [MERGE])
marker_only = FakeEntity(103, 'FBba0000999', 'X', [MARK])
h = FakeHandler([basc_parent], [basc, merge_only, marker_only])
h.synthesize_balancer_rename_map()
check('rename flag mapped', h.balancer_rename_map == {101: 1})
check('merge-only balancer not renamed', 102 not in h.balancer_rename_map)
check('mark flag not treated as a rename', 103 not in h.balancer_rename_map)
check('no errors', h.log.errors == [])

print('--- unknown parent is an error ---')
h = FakeHandler([basc_parent], [FakeEntity(104, 'FBba0000777', 'Y', [RENAME_UNKNOWN])])
h.synthesize_balancer_rename_map()
check('unknown parent not renamed', h.balancer_rename_map == {})
check('error names the unknown parent', any('FBab9999999' in e for e in h.log.errors))

print('--- count warning away from 24, suppressed in testing mode ---')
check('count warning present', any('24' in w for w in h.log.warnings))
h = FakeHandler([basc_parent], [FakeEntity(101, 'FBba0000014', 'Basc', [RENAME])], testing=True)
h.synthesize_balancer_rename_map()
check('no count warning in testing mode', not any('24' in w for w in h.log.warnings))

print('--- rename: balancer names take the slots, aberration names become synonyms ---')
parent_dto = FakeDTO(
    symbol=name_dto('nomenclature_symbol', 'In(1)FM3'),
    full_name=name_dto('full_name', 'Inversion (1) First Multiple'),
    synonyms=[name_dto('nomenclature_symbol', 'FM3'),        # grafted by FTA-236, must not duplicate
              name_dto('nomenclature_symbol', 'dm-FM3')],
)
parent = FakeEntity(1, 'FBab0003926', 'In(1)FM3', dto=parent_dto)
balancer_dto = FakeDTO(symbol=name_dto('nomenclature_symbol', 'FM3'),
                       full_name=name_dto('full_name', 'First Multiple 3'))
balancer = FakeEntity(101, 'FBba0000002', 'FM3', dto=balancer_dto)
h = FakeHandler([parent], [balancer])
h.balancer_rename_map = {101: 1}
h.map_balancer_renames()

check('parent symbol is now the balancer symbol', parent_dto.allele_symbol_dto['format_text'] == 'FM3')
check('parent full name is now the balancer full name', parent_dto.allele_full_name_dto['format_text'] == 'First Multiple 3')
syn_texts = [i['format_text'] for i in parent_dto.allele_synonym_dtos]
check('old symbol kept as a synonym', 'In(1)FM3' in syn_texts)
check('old full name kept as a synonym', 'Inversion (1) First Multiple' in syn_texts)
check('promoted symbol removed from synonyms (no duplicate)', syn_texts.count('FM3') == 0)
check('unrelated synonym untouched', 'dm-FM3' in syn_texts)
check('old symbol synonym keeps its name type', next(i['name_type_name'] for i in parent_dto.allele_synonym_dtos
                                                     if i['format_text'] == 'In(1)FM3') == 'nomenclature_symbol')
check('slot dicts are copies, not aliases', parent_dto.allele_symbol_dto is not balancer_dto.allele_symbol_dto)
row = h.balancer_rename_report[0]
check('report captures both names', row['new_symbol'] == 'FM3' and row['old_symbol'] == 'In(1)FM3')
check('report marks the source', row['source'] == 'flag')

print('--- parent with no full name of its own (18 of the 24) ---')
parent_dto2 = FakeDTO(symbol=name_dto('nomenclature_symbol', 'In(1)Basc'), full_name=None, synonyms=[])
parent2 = FakeEntity(2, 'FBab0004219', 'In(1)Basc', dto=parent_dto2)
balancer2 = FakeEntity(102, 'FBba0000014', 'Basc',
                       dto=FakeDTO(symbol=name_dto('nomenclature_symbol', 'Basc'),
                                   full_name=name_dto('full_name', 'Bar-apricot-scute')))
h2 = FakeHandler([parent2], [balancer2])
h2.balancer_rename_map = {102: 2}
h2.map_balancer_renames()
check('empty full name slot just gets filled', parent_dto2.allele_full_name_dto['format_text'] == 'Bar-apricot-scute')
check('no phantom full-name synonym added', [i['format_text'] for i in parent_dto2.allele_synonym_dtos] == ['In(1)Basc'])
check('report records an empty old full name', h2.balancer_rename_report[0]['old_full_name'] == '')

print('--- a parent whose DTO was never built is skipped, not crashed on ---')
parent3 = FakeEntity(3, 'FBab0005463', 'In(3LR)TM3', dto=None)
h3 = FakeHandler([parent3], [FakeEntity(103, 'FBba0000047', 'TM3', dto=FakeDTO(symbol=name_dto('nomenclature_symbol', 'TM3')))])
h3.balancer_rename_map = {103: 3}
h3.map_balancer_renames()
check('no DTO means no rename and no exception', h3.balancer_rename_report == [])

print('--- hard-coded SM6: names come from the ticket, not from chado ---')
sm6_dto = FakeDTO(symbol=name_dto('nomenclature_symbol', 'In(2LR)SM6'), full_name=None,
                  synonyms=[name_dto('nomenclature_symbol', 'SM6'),        # non-current synonym in chado
                            name_dto('nomenclature_symbol', 'SM6a'),       # grafted by FTA-236
                            name_dto('full_name', 'Second Multiple 6a')])  # grafted by FTA-236
sm6_parent = FakeEntity(4, 'FBab0004818', 'In(2LR)SM6', dto=sm6_dto)
h4 = FakeHandler([sm6_parent], [])
h4.balancer_rename_map = {}
h4.map_balancer_renames()

check('SM6 symbol applied with no flag present', sm6_dto.allele_symbol_dto['format_text'] == 'SM6')
check('SM6 full name fabricated', sm6_dto.allele_full_name_dto['format_text'] == 'Second Multiple 6')
sm6_syns = [i['format_text'] for i in sm6_dto.allele_synonym_dtos]
check('own symbol demoted to synonym', 'In(2LR)SM6' in sm6_syns)
check('pre-existing SM6 synonym de-duplicated', sm6_syns.count('SM6') == 0)
check('balancer names left as synonyms', 'SM6a' in sm6_syns and 'Second Multiple 6a' in sm6_syns)
sm6_row = next(r for r in h4.balancer_rename_report if r['fbab_id'] == 'FBab0004818')
check('report marks it hard-coded', sm6_row['source'] == 'hard-coded' and sm6_row['fbba_id'] == 'n/a')
check('hard-coded rename has no evidence curies', sm6_dto.allele_symbol_dto['evidence_curies'] == [])

print('--- a hard-coded target that is not in the export is reported, not applied ---')
h5 = FakeHandler([], [])
h5.balancer_rename_map = {}
h5.map_balancer_renames()
check('absent hard-coded parent logs an error', any('FBab0004818' in e for e in h5.log.errors))

print(f'\n{sum(1 for _, ok in results if ok)}/{len(results)} checks PASSED')
