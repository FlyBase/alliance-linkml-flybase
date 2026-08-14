# FTA-237 — Rename flagged aberrations to the balancer symbol and full name

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For the 24 FBba balancers carrying `FTA: Balancer - use balancer symbol and fullname for parent <SYMBOL> (<FBab_ID>).`, plus the one hard-coded SM6 case, export the parent FBab aberration under the balancer's current symbol and full name, keeping the aberration's own symbol and full name as synonyms.

**Architecture:** The nested `BalancerHandler` run (shared with FTA-236) already produces a complete `AlleleDTO` per FBba, including a correct `allele_symbol_dto` and `allele_full_name_dto` — right `display_text` SGML→HTML conversion, right evidence curies. So the rename is a DTO-level swap performed after `map_synonyms()`: copy the balancer's two name DTOs into the parent's slots, push the parent's previous ones into `allele_synonym_dtos`, and de-duplicate. No new name text is constructed for the 24, and no part of the synonym-synthesis machinery is modified.

**Tech Stack:** Python 3.11, SQLAlchemy ORM over FlyBase chado (`harvdev_utils.reporting` models), Alliance LinkML DTOs in `src/agr_datatypes.py`.

## Global Constraints

- **PEP8, `max-line-length = 160`** from `.flake8`. `ignore = D100,D103,D104` **replaces** pycodestyle's defaults, so **W503 and W504 are both active**: never break a line before *or* after a boolean operator. Docstrings are required on new public methods (D102 is active).
- **No local script runs.** The repo venv is x86_64 and unusable on this arm64 Mac. Verification is `python3 -m py_compile`, the offline `ast`-extraction harness in this plan, SQL against chado, and a real run on flysql26.
- **chado access:** `ssh flysql26` then `PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado`. Set `standard_conforming_strings = on` before using `\y` in a regex.
- **Run artifacts:** `flysql26:/data/alliance/PERSISTENT_*`. A full allele run takes ~85 minutes.
- **No gate.** This ticket only touches `allele_symbol_dto`, `allele_full_name_dto` and `allele_synonym_dtos`, all long-released LinkML slots. Nothing here needs `ADD_IS_ABERRATION` or `ADD_ALLELE_ALLELE_ASSOC`.
- **Never commit without explicit permission from Ian.** Commit steps are prepared but only run on request.
- **Conventional Commits**, scope `FTA-237`, `git add <specific paths>` only.

## Relationship to FTA-236 — read this first

FTA-236 (plan: `docs/FTA-236/plan.md`, branch `FTA-236`, commit `ec4ce3b`) is **parked** pending a curator answer, but its **Task 1 is the shared foundation**: the nested `BalancerHandler` run (`get_balancer_entities()`) and the flag-parsing/parent-resolution helpers (`_resolve_parent_feature_id()`).

- If FTA-236's Task 1 has already landed, **Task 1 here is just the second regex plus one map method** — skip its Steps 3-5 and say so in the commit.
- If it has not, Task 1 here includes the nested retrieval verbatim from `docs/FTA-236/plan.md` Task 1 Steps 3-6. FTA-236's parked question does not affect this ticket: **all 24 rename parents are current and non-obsolete**, so both candidate parent-resolution policies behave identically here.

The two tickets must stay in this order of operations, and the interaction is the main design risk:

1. FTA-236 grafts the balancer's names onto the parent and marks them all non-current (`demoted_synonym_ids`), so they land in `allele_synonym_dtos`.
2. FTA-237 then promotes two of them — the current symbol and current full name — into the parent's name slots, and removes the now-duplicated synonym entries.

Implemented this way, **either ticket works alone**: with FTA-236 absent the balancer names simply are not in the synonym list, and the de-duplication finds nothing to remove.

---

## Verified input data (production_chado, 2026-08-14)

24 rename flags → 24 distinct parents, **all current and non-obsolete**; all 24 balancers have both a current symbol and a current full name; the 24 are a strict subset of FTA-236's 38 merge-flagged balancers (14 merge without renaming, none rename without merging).

| FBba | new symbol | new full name | parent FBab | parent symbol → synonym | parent full name → synonym |
|---|---|---|---|---|---|
| FBba0000002 | FM3 | First Multiple 3 | FBab0003926 | In(1)FM3 | Inversion (1) First Multiple |
| FBba0000070 | FM4 | First Multiple 4 | FBab0003927 | In(1)FM4 | *(none)* |
| FBba0000003 | FM6 | First Multiple 6 | FBab0003928 | In(1)FM6 | *(none)* |
| FBba0000007 | FM7a | First Multiple 7a | FBab0003929 | In(1)FM7a | *(none)* |
| FBba0000014 | Basc | Bar-apricot-scute | FBab0004219 | In(1)Basc | *(none)* |
| FBba0000025 | CyO | Curly of Oster | FBab0004786 | In(2LR)CyO | Inversion (2LR) Curly of Oster |
| FBba0000037 | SM1 | Second Multiple 1 | FBab0004815 | In(2LR)SM1 | Inversion (2LR) Second Multiple |
| FBba0000038 | SM5 | Second Multiple 5 | FBab0004817 | In(2LR)SM5 | *(none)* |
| FBba0000104 | DcxF | Dichaete crossover suppressor of Federova | FBab0005386 | In(3LR)DcxF | Inversion (3LR) Dichaete crossover suppressor of Federova |
| FBba0000042 | TM1 | Third Multiple 1 | FBab0005462 | In(3LR)TM1 | Inversion (3LR) Third Multiple |
| FBba0000047 | TM3 | Third Multiple 3 | FBab0005463 | In(3LR)TM3 | *(none)* |
| FBba0000056 | TM6 | Third Multiple 6 | FBab0005464 | In(3LR)TM6 | *(none)* |
| FBba0000057 | TM6B | Third Multiple 6B | FBab0005465 | In(3LR)TM6B | *(none)* |
| FBba0000071 | TM6C | Third Multiple 6C | FBab0005466 | In(3LR)TM6C | *(none)* |
| FBba0000060 | TM8 | Third Multiple 8 | FBab0005467 | In(3LR)TM8 | *(none)* |
| FBba0000061 | TM9 | Third Multiple 9 | FBab0005468 | In(3LR)TM9 | *(none)* |
| FBba0000072 | TMS | Third Multiple of Singson | FBab0005469 | In(3LR)TMS | *(none)* |
| FBba0000062 | TM2 | Third Multiple 2 | FBab0005473 | In(3LR)TM2 | Inversion (3LR) Ultrabithorax |
| FBba0000068 | MRS | Minute-rosy-Stubble | FBab0010017 | Tp(3;3)MRS | *(none)* |
| FBba0000011 | FM1 | First Multiple 1 | FBab0010486 | In(1)FM1 | *(none)* |
| FBba0000021 | Insc | Inversion scute | FBab0010488 | In(1)Insc | *(none)* |
| FBba0000001 | CyO-Df(2R)B80 | Curly of Oster-Df(2R)B80 | FBab0022240 | In(2LR)CyO-Df(2R)B80 | *(none)* |
| FBba0000237 | CyO-Tp(2;2)pk-sple[26] | Curly of Oster-Tp(2;2)pk-sple[26] | FBab0024785 | In(2LR)CyO-Tp(2;2)pk-sple[26] | *(none)* |
| FBba0000659 | FM8 | First Multiple 8 | FBab0049294 | In(1)FM8 | *(none)* |

Only 6 of the 24 parents have a full name of their own to demote; the other 18 have an empty slot that the balancer's full name fills.

## The SM6 case is not flag-driven (correction to the ticket text)

The ticket presents SM6 as an exception "because SM6a and SM6b share a parent". The real reason it must be hard-coded is different, and it changes what Task 3 has to do: **neither `FBba0000039` (SM6a) nor `FBba0000040` (SM6b) carries a rename flag at all.** They are merge-flagged only, so `FBab0004818` is not among the 24 and no flag can produce its new names.

In chado:

- `FBab0004818` current symbol is `In(2LR)SM6`, and it has **no** current full name.
- `SM6` exists as a **non-current** symbol synonym on `FBab0004818` itself (and on both balancers).
- **`Second Multiple 6` does not exist anywhere in chado** — not as a synonym of the aberration or of either balancer.
- The balancers' own current names are `SM6a` / `Second Multiple 6a` and `SM6b` / `Second Multiple 6b`, none of which is wanted.

So both target strings come from the ticket, and Task 3 constructs the two DTOs itself. Total: **24 flag-driven renames + 1 hard-coded = 25 renamed aberrations.**

## File Structure

- **Modify `src/allele_handlers.py`** — everything lands in `class AberrationHandler`: two class attributes, one synthesis method, one mapping method, and two call sites in the existing `synthesize_info()` / `map_fb_data_to_alliance()` hooks.
- **Modify `src/AGR_data_retrieval_curation_allele.py`** — a curator TSV of the renames, plus a line in the `--help` epilog.
- **Create `docs/FTA-237/verify_renames.py`** — offline harness for Tasks 1-3.
- **Create `docs/FTA-237/expected_renames.sql`** — the SQL producing the table above (written and verified 2026-08-14), so a run's TSV can be diffed against chado. It is the oracle for `source=flag` rows, since the curator spreadsheet needs auth.

---

### Task 1: Parse the rename flags

**Files:**
- Modify: `src/allele_handlers.py` (class `AberrationHandler`: new `balancer_rename_regex` attribute beside the FTA-235/236 regexes; new `synthesize_balancer_rename_map()`; one call in `synthesize_info()`)
- Create: `docs/FTA-237/verify_renames.py`
- Create: `docs/FTA-237/expected_renames.sql`

**Interfaces:**
- Consumes: `self.fbba_entities` and `self._resolve_parent_feature_id()` from FTA-236 Task 1 (see "Relationship to FTA-236" — implement that task's Steps 3-6 first if they are not already present), `self.fb_data_entities`, `self.testing`, `self.log`.
- Produces:
  - `self.balancer_rename_map` — `{FBba feature_id: parent FBab feature_id}`, 24 entries.
  - `self.balancer_rename_regex` — compiled pattern with named groups `symbol` and `fbab`.
  - `BALANCER_RENAME_EXPECTED = 24` — module-level constant.

- [ ] **Step 1: Write the failing test**

Create `docs/FTA-237/verify_renames.py`:

```python
"""Offline harness for the FTA-237 balancer rename methods.

Extracts method source with ast and drives it with duck-typed stand-ins, since sqlalchemy and
harvdev_utils cannot be imported on this machine (see docs/FTA-237/plan.md).
"""

import ast
import re

SRC = '/Users/ilongden/harvard/alliance-linkml-flybase/src/allele_handlers.py'

tree = ast.parse(open(SRC).read())
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'AberrationHandler')



class FakeNameSlotAnnotationDTO:
    """Stand-in for agr_datatypes.NameSlotAnnotationDTO (Task 3's hard-coded path builds these)."""

    def __init__(self, name_type_name, format_text, display_text, evidence_curies):
        self.dto = {'name_type_name': name_type_name, 'format_text': format_text,
                    'display_text': display_text, 'synonym_scope_name': 'exact',
                    'evidence_curies': list(evidence_curies), 'internal': False, 'obsolete': False}

    def dict_export(self):
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


class FakeLog:
    def __init__(self):
        self.info_lines, self.warnings, self.errors = [], [], []

    def info(self, msg):
        self.info_lines.append(msg)

    def debug(self, msg):
        pass

    def warning(self, msg):
        self.warnings.append(msg)

    def error(self, msg):
        self.errors.append(msg)


class FakeProp:
    def __init__(self, value):
        self.chado_obj = type('O', (), {'value': value})()


def name_dto(name_type, text):
    return {'name_type_name': name_type, 'format_text': text, 'display_text': text,
            'synonym_scope_name': 'exact', 'internal': False, 'obsolete': False}


class FakeDTO:
    def __init__(self, symbol=None, full_name=None, synonyms=None):
        self.allele_symbol_dto = symbol
        self.allele_full_name_dto = full_name
        self.allele_synonym_dtos = list(synonyms or [])


class FakeEntity:
    def __init__(self, feature_id, uniquename, name='', internal_notes=(), dto=None):
        self.db_primary_id = feature_id
        self.uniquename = uniquename
        self.name = name
        self.is_obsolete = False
        self.props_by_type = {'internal_notes': [FakeProp(v) for v in internal_notes]} if internal_notes else {}
        self.linkmldto = dto

    def __str__(self):
        return f'{self.name} ({self.uniquename})'


class FakeHandler:
    balancer_rename_regex = ns['balancer_rename_regex']
    balancer_hardcoded_renames = ns.get('balancer_hardcoded_renames', {})
    synthesize_balancer_rename_map = ns['synthesize_balancer_rename_map']
    map_balancer_renames = ns.get('map_balancer_renames')
    drop_duplicate_synonyms = ns.get('drop_duplicate_synonyms')
    apply_rename = ns.get('apply_rename')
    map_hardcoded_balancer_renames = ns.get('map_hardcoded_balancer_renames')

    def __init__(self, aberrations, balancers, testing=False):
        self.log = FakeLog()
        self.testing = testing
        self.fb_data_entities = {a.db_primary_id: a for a in aberrations}
        self.fbba_entities = {b.db_primary_id: b for b in balancers}
        self.balancer_rename_map = {}
        self.balancer_rename_report = []

    def _resolve_parent_feature_id(self, fbab_uniquename):
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

print(f'\n{sum(1 for _, ok in results if ok)}/{len(results)} checks PASSED')
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python3 docs/FTA-237/verify_renames.py`
Expected: `KeyError: 'balancer_rename_regex'` — the attribute does not exist yet.

- [ ] **Step 3: Add the constant and the regex**

At module level in `src/allele_handlers.py`, beside `BALANCER_MERGE_EXPECTED` (add that too if FTA-236 Task 1 has not landed):

```python
BALANCER_RENAME_EXPECTED = 24    # FTA-237: curated rename flags in chado as of 2026-08-14.
```

In `class AberrationHandler`, below the FTA-236 attributes:

```python
    # FTA-237: this flag says the parent FBab should be exported under the balancer's current symbol
    # and full name, e.g.
    #   FTA: Balancer - use balancer symbol and fullname for parent In(1)Basc (FBab0004219).
    # 24 of these exist, one per parent. As in FTA-235/236 the match runs through the instruction
    # clause, since all three flag families open with "FTA: Balancer -".
    balancer_rename_regex = re.compile(r'^\s*:?\s*FTA:\s*balancer\s*-\s*use balancer symbol and fullname for parent\s+(?P<symbol>.+?)\s*\((?P<fbab>FBab[0-9]+)\)\s*\.?\s*$',
                                       re.IGNORECASE)
```

- [ ] **Step 4: Add the rename-map method**

In the "Additional sub-methods to be run by synthesize_info()" area of `class AberrationHandler`:

```python
    def synthesize_balancer_rename_map(self):
        """Map each rename-flagged FBba balancer to the parent FBab it renames (FTA-237)."""
        self.log.info('Map rename-flagged FBba balancers to their parent FBab aberrations.')
        flag_counter = 0
        for balancer in self.fbba_entities.values():
            for fb_prop in balancer.props_by_type.get('internal_notes', []):
                prop_value = fb_prop.chado_obj.value
                if not prop_value:
                    continue
                match = self.balancer_rename_regex.match(prop_value)
                if match is None:
                    continue
                flag_counter += 1
                fbab_uniquename = match.group('fbab')
                parent_feature_id = self._resolve_parent_feature_id(fbab_uniquename)
                if parent_feature_id is None:
                    self.log.error(f'FTA-237: balancer {balancer} would rename parent {fbab_uniquename} '
                                   f'("{match.group("symbol")}"), which is not an exportable aberration '
                                   f'(obsolete or unknown). No rename applied; please fix the internal note.')
                    break
                self.balancer_rename_map[balancer.db_primary_id] = parent_feature_id
                break
        self.log.info(f'Found {flag_counter} balancer rename flags; resolved {len(self.balancer_rename_map)} renames.')
        if self.testing is False and flag_counter != BALANCER_RENAME_EXPECTED:
            self.log.warning(f'Expected {BALANCER_RENAME_EXPECTED} balancer rename flags per FTA-237, but found '
                             f'{flag_counter}. Check the "internal_notes" flag text in chado.')
        return
```

- [ ] **Step 5: Add the class attribute for the report and wire the call**

Beside the other FTA-23x attributes:

```python
    balancer_rename_map = {}      # FTA-237: {FBba feature_id: parent FBab feature_id} for resolved rename flags.
    balancer_rename_report = []   # FTA-237: dicts of old/new names per renamed aberration, for the curator TSV.
```

In `AberrationHandler.synthesize_info()`, immediately after `self.synthesize_balancer_merge_map()` if that exists, otherwise immediately after `super().synthesize_info()`:

```python
        self.synthesize_balancer_rename_map()
```

- [ ] **Step 6: Run the harness to verify it passes**

Run: `python3 docs/FTA-237/verify_renames.py`
Expected: `8/8 checks PASSED`

- [ ] **Step 7: Compile check and line length**

Run:
```bash
python3 -m py_compile src/allele_handlers.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py
grep -n ' $' src/allele_handlers.py
```
Expected: `COMPILE OK` and no other output.

- [ ] **Step 8: Record the chado expectations**

`docs/FTA-237/expected_renames.sql` is already written and was verified against production_chado on 2026-08-14 — every number in this plan came from it. Re-run it to confirm chado has not moved:

```bash
scp docs/FTA-237/expected_renames.sql flysql26:/tmp/
ssh flysql26 "PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado -f /tmp/expected_renames.sql"
```
Expected: 24 rows matching the table; `SM6` present only as a non-current synonym; `Second Multiple 6` absent.

- [ ] **Step 9: Commit (only when Ian asks)**

```bash
git add src/allele_handlers.py docs/FTA-237/plan.md docs/FTA-237/verify_renames.py docs/FTA-237/expected_renames.sql
git commit -m "feat(FTA-237): map rename-flagged balancers to their parent aberrations"
```

---

### Task 2: Apply the rename at the DTO level

**Files:**
- Modify: `src/allele_handlers.py` (class `AberrationHandler`: new `map_balancer_renames()`; one call in `map_fb_data_to_alliance()`)
- Modify: `docs/FTA-237/verify_renames.py` (add the rename checks)

**Interfaces:**
- Consumes: `self.balancer_rename_map` from Task 1; `balancer.linkmldto.allele_symbol_dto` / `.allele_full_name_dto` — the nested `BalancerHandler` builds a full `AlleleDTO` per FBba (its `agr_export_type` is `AlleleDTO`, and `map_synonyms()` derives the same `allele_*` slot names), so these dicts already hold correct `format_text`, `display_text` and `evidence_curies`.
- Produces: `map_balancer_renames()` — mutates each renamed parent's `linkmldto`; appends one dict per rename to `self.balancer_rename_report` with keys `fbab_id`, `fbba_id`, `new_symbol`, `old_symbol`, `new_full_name`, `old_full_name`, `source` (`'flag'` or `'hard-coded'`).

- [ ] **Step 1: Write the failing test**

Append to `docs/FTA-237/verify_renames.py`:

```python
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
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python3 docs/FTA-237/verify_renames.py`
Expected: `TypeError: 'NoneType' object is not callable` from `h.map_balancer_renames()` — the harness sets `map_balancer_renames = ns.get(...)`, which is still `None`.

- [ ] **Step 3: Write the method**

```python
    def map_balancer_renames(self):
        """Export a rename-flagged aberration under its balancer's symbol and full name (FTA-237).

        The nested BalancerHandler has already built a full AlleleDTO for every FBba, so the balancer's
        allele_symbol_dto and allele_full_name_dto are correct - right display_text conversion, right
        evidence curies - and are simply copied into the parent's slots. The aberration's own names
        move to allele_synonym_dtos, as the ticket asks.

        Must run after map_synonyms(), which is what populates the parent's name slots in the first
        place. Any synonym matching a newly promoted name is dropped, so a name grafted by FTA-236
        does not appear twice.
        """
        self.log.info('Rename flagged aberrations to use their balancer symbol and full name.')
        NAME_SLOTS = (('allele_symbol_dto', 'symbol'), ('allele_full_name_dto', 'full_name'))
        renamed = 0
        for fbba_feature_id, fbab_feature_id in self.balancer_rename_map.items():
            balancer = self.fbba_entities[fbba_feature_id]
            aberration = self.fb_data_entities[fbab_feature_id]
            if aberration.linkmldto is None:
                self.log.warning(f'FTA-237: {aberration} has no LinkML object to rename; skipping.')
                continue
            report = {
                'fbab_id': aberration.uniquename,
                'fbba_id': balancer.uniquename,
                'source': 'flag',
                'new_symbol': '',
                'old_symbol': '',
                'new_full_name': '',
                'old_full_name': '',
            }
            for slot_name, label in NAME_SLOTS:
                new_dto = getattr(balancer.linkmldto, slot_name, None)
                if new_dto is None:
                    continue
                old_dto = getattr(aberration.linkmldto, slot_name)
                report[f'new_{label}'] = new_dto['format_text']
                if old_dto is not None:
                    report[f'old_{label}'] = old_dto['format_text']
                    aberration.linkmldto.allele_synonym_dtos.append(old_dto)
                # Copy, so later edits to either DTO cannot bleed into the other entity's export.
                setattr(aberration.linkmldto, slot_name, dict(new_dto))
                self.drop_duplicate_synonyms(aberration, new_dto)
            self.balancer_rename_report.append(report)
            renamed += 1
            self.log.debug(f'FTA-237: renamed {aberration} to "{report["new_symbol"]}" '
                           f'("{report["new_full_name"]}") after balancer {balancer}.')
        self.log.info(f'Renamed {renamed} aberrations to their balancer symbol and full name.')
        return

    def drop_duplicate_synonyms(self, aberration, promoted_dto):
        """Remove synonyms that duplicate a name just promoted into a name slot (FTA-237).

        FTA-236 grafts the balancer's names onto the parent as synonyms; once one of them becomes the
        parent's symbol or full name, the synonym entry is redundant. Matching is on name type plus
        format_text, so a same-text synonym of a different type is left alone.
        """
        keep = []
        for synonym_dto in aberration.linkmldto.allele_synonym_dtos:
            is_duplicate = synonym_dto['format_text'] == promoted_dto['format_text']
            if is_duplicate and synonym_dto['name_type_name'] == promoted_dto['name_type_name']:
                continue
            keep.append(synonym_dto)
        aberration.linkmldto.allele_synonym_dtos = keep
        return
```

- [ ] **Step 4: Wire it into `map_fb_data_to_alliance()`**

In `AberrationHandler.map_fb_data_to_alliance()`, immediately **after** `self.map_synonyms()`:

```python
        self.map_balancer_renames()
```

- [ ] **Step 5: Run the harness to verify it passes**

Run: `python3 docs/FTA-237/verify_renames.py`
Expected: `22/22 checks PASSED` (8 from Task 1 plus 14 new)

- [ ] **Step 6: Compile check and line length**

Run:
```bash
python3 -m py_compile src/allele_handlers.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py
```
Expected: `COMPILE OK` and no `awk` output.

- [ ] **Step 7: Commit (only when Ian asks)**

```bash
git add src/allele_handlers.py docs/FTA-237/verify_renames.py
git commit -m "feat(FTA-237): rename flagged aberrations to the balancer symbol and full name"
```

---

### Task 3: The hard-coded SM6 rename

**Files:**
- Modify: `src/allele_handlers.py` (class `AberrationHandler`: new `balancer_hardcoded_renames` attribute; extend `map_balancer_renames()`)
- Modify: `docs/FTA-237/verify_renames.py` (add the SM6 checks)

**Interfaces:**
- Consumes: `sub_sup_sgml_to_html` and `sub_sup_to_sgml` — **not currently imported by `src/allele_handlers.py`**, so add them to the `harvdev_utils.char_conversions` import list at the top of the module (`src/entity_handler.py:16` shows the same import, and `entity_handler.py:1280` the same fabrication idiom).
- Produces: renames driven by `balancer_hardcoded_renames`, reported with `source = 'hard-coded'` and `fbba_id = 'n/a'`.

- [ ] **Step 1: Write the failing test**

Append to `docs/FTA-237/verify_renames.py`:

```python
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
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python3 docs/FTA-237/verify_renames.py`
Expected: FAIL on `SM6 symbol applied with no flag present` — with no rename flag and no hard-coded table, `map_balancer_renames()` does nothing for `FBab0004818`.

- [ ] **Step 3: Add the import**

At the top of `src/allele_handlers.py`, add the char-conversion helpers (the module currently imports none):

```python
from harvdev_utils.char_conversions import sub_sup_sgml_to_html, sub_sup_to_sgml
```

- [ ] **Step 4: Add the hard-coded table**

In `class AberrationHandler`, beside the other FTA-237 attributes:

```python
    # FTA-237: one aberration is renamed without a curated flag. Per the ticket, In(2LR)SM6 takes the
    # symbol "SM6" and full name "Second Multiple 6", because its two balancers SM6a (FBba0000039) and
    # SM6b (FBba0000040) share it as a parent and neither of their names should win. Neither balancer
    # carries a rename flag, so nothing in chado can drive this: "SM6" exists only as a non-current
    # synonym of the aberration, and "Second Multiple 6" does not exist in chado at all. Keyed by FBab
    # uniquename; the values are (symbol, full_name).
    balancer_hardcoded_renames = {
        'FBab0004818': ('SM6', 'Second Multiple 6'),
    }
```

- [ ] **Step 5: Extend `map_balancer_renames()`**

Refactor the per-slot work into a helper both paths share. `map_balancer_renames()` from Task 2 becomes exactly this — the slot loop is gone, replaced by an `apply_rename()` call, and the hard-coded table runs at the end so one call site in `map_fb_data_to_alliance()` still covers everything:

```python
    def map_balancer_renames(self):
        """Export a rename-flagged aberration under its balancer's symbol and full name (FTA-237).

        The nested BalancerHandler has already built a full AlleleDTO for every FBba, so the balancer's
        allele_symbol_dto and allele_full_name_dto are correct - right display_text conversion, right
        evidence curies - and are simply copied into the parent's slots. The aberration's own names
        move to allele_synonym_dtos, as the ticket asks.

        Must run after map_synonyms(), which is what populates the parent's name slots in the first
        place. Any synonym matching a newly promoted name is dropped, so a name grafted by FTA-236
        does not appear twice.
        """
        self.log.info('Rename flagged aberrations to use their balancer symbol and full name.')
        renamed = 0
        for fbba_feature_id, fbab_feature_id in self.balancer_rename_map.items():
            balancer = self.fbba_entities[fbba_feature_id]
            aberration = self.fb_data_entities[fbab_feature_id]
            if aberration.linkmldto is None:
                self.log.warning(f'FTA-237: {aberration} has no LinkML object to rename; skipping.')
                continue
            report = {
                'fbab_id': aberration.uniquename,
                'fbba_id': balancer.uniquename,
                'source': 'flag',
                'new_symbol': '',
                'old_symbol': '',
                'new_full_name': '',
                'old_full_name': '',
            }
            new_dtos = {
                'allele_symbol_dto': getattr(balancer.linkmldto, 'allele_symbol_dto', None),
                'allele_full_name_dto': getattr(balancer.linkmldto, 'allele_full_name_dto', None),
            }
            self.apply_rename(aberration, new_dtos, report)
            self.balancer_rename_report.append(report)
            renamed += 1
            self.log.debug(f'FTA-237: renamed {aberration} to "{report["new_symbol"]}" '
                           f'("{report["new_full_name"]}") after balancer {balancer}.')
        self.map_hardcoded_balancer_renames()
        self.log.info(f'Renamed {renamed} aberrations from curated flags, plus '
                      f'{len(self.balancer_rename_report) - renamed} hard-coded; '
                      f'{len(self.balancer_rename_report)} in total.')
        return
```

`drop_duplicate_synonyms()` from Task 2 is unchanged. Add the two new methods:

```python
    def apply_rename(self, aberration, new_dtos, report):
        """Move an aberration's own names to synonyms and install the new ones (FTA-237).

        Args:
            aberration (FBAberration): the aberration being renamed.
            new_dtos (dict): {slot_name: NameSlotAnnotationDTO dict or None} for the new names.
            report (dict): the report row to fill in, mutated in place.

        """
        NAME_SLOT_LABELS = {'allele_symbol_dto': 'symbol', 'allele_full_name_dto': 'full_name'}
        for slot_name, label in NAME_SLOT_LABELS.items():
            new_dto = new_dtos.get(slot_name)
            if new_dto is None:
                continue
            old_dto = getattr(aberration.linkmldto, slot_name)
            report[f'new_{label}'] = new_dto['format_text']
            if old_dto is not None:
                report[f'old_{label}'] = old_dto['format_text']
                aberration.linkmldto.allele_synonym_dtos.append(old_dto)
            setattr(aberration.linkmldto, slot_name, dict(new_dto))
            self.drop_duplicate_synonyms(aberration, new_dto)
        return

    def map_hardcoded_balancer_renames(self):
        """Apply the FTA-237 renames that no curated flag can drive (currently only In(2LR)SM6)."""
        for fbab_uniquename, names in self.balancer_hardcoded_renames.items():
            new_symbol, new_full_name = names
            parent_feature_id = self._resolve_parent_feature_id(fbab_uniquename)
            if parent_feature_id is None:
                self.log.error(f'FTA-237: hard-coded rename target {fbab_uniquename} ("{new_symbol}") is not an '
                               f'exportable aberration; rename not applied.')
                continue
            aberration = self.fb_data_entities[parent_feature_id]
            if aberration.linkmldto is None:
                self.log.warning(f'FTA-237: {aberration} has no LinkML object to rename; skipping.')
                continue
            new_dtos = {
                'allele_symbol_dto': agr_datatypes.NameSlotAnnotationDTO(
                    'nomenclature_symbol', new_symbol, sub_sup_sgml_to_html(sub_sup_to_sgml(new_symbol)), []).dict_export(),
                'allele_full_name_dto': agr_datatypes.NameSlotAnnotationDTO(
                    'full_name', new_full_name, sub_sup_sgml_to_html(sub_sup_to_sgml(new_full_name)), []).dict_export(),
            }
            report = {
                'fbab_id': aberration.uniquename,
                'fbba_id': 'n/a',
                'source': 'hard-coded',
                'new_symbol': '',
                'old_symbol': '',
                'new_full_name': '',
                'old_full_name': '',
            }
            self.apply_rename(aberration, new_dtos, report)
            self.balancer_rename_report.append(report)
            self.log.info(f'FTA-237: applied the hard-coded rename of {aberration} to '
                          f'"{new_symbol}" ("{new_full_name}").')
        return
```

No change is needed in `map_fb_data_to_alliance()` — the Task 2 call site now covers both paths, and the closing log line counts 24 + 1 = 25.

- [ ] **Step 6: Run the harness to verify it passes**

Run: `python3 docs/FTA-237/verify_renames.py`
Expected: `30/30 checks PASSED` (22 from Tasks 1-2 plus 8 new)

- [ ] **Step 7: Compile check and line length**

Run:
```bash
python3 -m py_compile src/allele_handlers.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py
```
Expected: `COMPILE OK` and no `awk` output.

- [ ] **Step 8: Commit (only when Ian asks)**

```bash
git add src/allele_handlers.py docs/FTA-237/verify_renames.py
git commit -m "feat(FTA-237): hard-code the In(2LR)SM6 rename to SM6 / Second Multiple 6"
```

---

### Task 4: Curator TSV, docs, and full-run verification

**Files:**
- Modify: `src/AGR_data_retrieval_curation_allele.py` (TSV emission in `main()`; epilog note)
- Create: `docs/FTA-237/verify_run.sh`

**Interfaces:**
- Consumes: `aberration_handler.balancer_rename_report`.
- Produces: `<output>_balancer_renames.tsv` beside the other allele TSVs.

- [ ] **Step 1: Emit the TSV**

In `main()`, after the `write_notes_tsv(...)` call (and after FTA-236's merge TSV if that has landed):

```python
        # FTA-237: one row per renamed aberration, so curators can check the 24 flag-driven renames and
        # the hard-coded In(2LR)SM6 case against their spreadsheet.
        renames_filename = tsv_filename.replace('.tsv', '_balancer_renames.tsv')
        curation_tsv.write_association_tsv(
            filename=renames_filename,
            rows=aberration_handler.balancer_rename_report,
            first_field='fbab_id',
            second_field='fbba_id',
            extra_fields=[
                ('source', 'source', ''),
                ('new_symbol', 'new_symbol', ''),
                ('old_symbol', 'old_symbol', ''),
                ('new_full_name', 'new_full_name', ''),
                ('old_full_name', 'old_full_name', ''),
            ],
        )
        log.info(f'Generated TSV: {renames_filename} ({len(aberration_handler.balancer_rename_report)} renames)')
```

- [ ] **Step 2: Document it in the epilog**

Add to the `Environment variables:` epilog, after the existing entries (it is behaviour rather than a variable, so keep it in its own short paragraph):

```
Notes:
  FTA-237: 24 aberrations are exported under their balancer's symbol and full name, driven by curated
  "FTA: Balancer - use balancer symbol and fullname for parent ..." internal notes, plus In(2LR)SM6
  (FBab0004818), which is renamed to "SM6" / "Second Multiple 6" from a hard-coded table because no
  flag exists for it. The aberration's own names are kept as synonyms. Not gated: the slots involved
  are long-released. See the *_balancer_renames.tsv output for the full list.
```

- [ ] **Step 3: Compile check both files**

Run:
```bash
python3 -m py_compile src/allele_handlers.py src/AGR_data_retrieval_curation_allele.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py src/AGR_data_retrieval_curation_allele.py
```
Expected: `COMPILE OK` and no `awk` output.

- [ ] **Step 4: Write the run-verification script**

Create `docs/FTA-237/verify_run.sh`:

```bash
#!/usr/bin/env bash
# Usage: verify_run.sh /data/alliance/PERSISTENT_main_FTA-237_balancer_rename   (run on flysql26)
set -euo pipefail
DIR="$1"
LOG=$(ls "$DIR"/allele_curation_*.log)
TSV=$(ls "$DIR"/allele_curation_*.tsv | grep -v -e _notes -e _renames -e _merges -e _failures)
RENAMES=$(ls "$DIR"/*_balancer_renames.tsv)
echo '--- rename tallies from the log (expect 24 flags, 25 renamed) ---'
grep -E 'balancer rename flags|Renamed .* aberrations to their balancer|hard-coded rename of' "$LOG"
echo '--- errors and rename warnings (expect none) ---'
grep -c ' ERROR ' "$LOG" || true
grep -E 'many current symbols|many current full_names' "$LOG" | wc -l
echo '--- rename TSV row count (expect 25) ---'
tail -n +2 "$RENAMES" | wc -l
echo '--- the worked example: In(1)Basc must export as Basc / Bar-apricot-scute ---'
grep -m1 'FBab0004219' "$TSV"
echo '--- the SM6 case: symbol SM6, full name Second Multiple 6, source hard-coded ---'
grep 'FBab0004818' "$RENAMES"
grep -m1 'FBab0004818' "$TSV"
echo '--- every renamed aberration must keep its old symbol as a synonym ---'
awk -F'\t' 'NR>1 {print $1"\t"$4"\t"$5}' "$RENAMES" | while IFS=$'\t' read -r fbab new old; do
  grep -m1 "^$fbab" "$TSV" | grep -q "$old" || echo "MISSING old symbol synonym: $fbab ($old)"
done
echo '--- done ---'
```

- [ ] **Step 5: Ask Ian to run the export, then verify**

The run cannot happen on this machine. Ask Ian for an allele run against `production_chado` into `PERSISTENT_main_FTA-237_balancer_rename`, then:

```bash
ssh flysql26 "bash -s" < docs/FTA-237/verify_run.sh /data/alliance/PERSISTENT_main_FTA-237_balancer_rename
```

Expected: 24 rename flags, 25 renamed aberrations, 25 TSV rows, 0 errors, **0** `many current symbols` / `many current full_names` lines, `FB:FBab0004219` exporting as `Basc`, `FB:FBab0004818` as `SM6`, and no `MISSING old symbol synonym` lines.

- [ ] **Step 6: Diff the whole rename list against chado**

```bash
ssh flysql26 "PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado -f /tmp/expected_renames.sql" > /tmp/chado_renames.txt
# compare the 24 flag-driven rows (new_symbol, new_full_name per FBab) with the TSV's source=flag rows
```
Expected: the 24 flag rows match `docs/FTA-237/expected_renames.sql` exactly; the only `source=hard-coded` row is `FBab0004818`.

- [ ] **Step 7: Commit (only when Ian asks)**

```bash
git add src/AGR_data_retrieval_curation_allele.py docs/FTA-237/verify_run.sh
git commit -m "feat(FTA-237): curator TSV of balancer renames, plus docs"
```

---

## Interaction checks to run once both tickets are implemented

If FTA-236 lands too, re-run its harness and this one, then check on a combined run:

- No aberration logs `many current symbols` or `many current full_names`. FTA-236 demotes all grafted balancer names; FTA-237 then promotes exactly two per renamed parent.
- For the 14 merge-only balancers (merge-flagged but not rename-flagged), the parent keeps its own symbol and full name, with the balancer's as synonyms.
- For the 24 renamed parents, `allele_synonym_dtos` contains the aberration's old symbol **and** the balancer's other synonyms, with no entry duplicating the promoted names.
- `FBab0004818` shows one FTA-237 row (`source=hard-coded`) and two FTA-236 merge rows (SM6a and SM6b), with `carries_excluded` non-zero for both and `carries_alleles` zero.

## Deferred, with reasons

- **Symbol collisions at the Alliance** — renaming `In(2LR)CyO` to `CyO` could collide with another MOD's or another FlyBase entity's symbol. Out of scope: the ticket specifies the names, and FlyBase does not police Alliance-side uniqueness here.
- **The curator spreadsheet as an automated oracle** — the "expected final submission data" Google Sheet needs auth to fetch, so `expected_renames.sql` is the machine-checkable oracle instead. If a discrepancy shows up, compare against the sheet by hand.
