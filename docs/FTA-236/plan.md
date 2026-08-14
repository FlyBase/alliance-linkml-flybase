# FTA-236 — Merge flagged FBba balancer data into parent FBab aberration entries

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For each of the 38 FBba balancers carrying the curated `FTA: Balancer - merge with parent <SYMBOL> (<FBab_ID>).` internal note, fold that balancer's synonyms, secondary IDs, comments, internal notes, references and "carries alleles"/"transposon insertions" relationships into its parent FBab aberration's Alliance allele record.

**Architecture:** The `AberrationHandler` gains a nested `BalancerHandler` run (exactly as `AlleleHandler.get_insertion_entities()` runs a nested `InsertionHandler`), parses the merge flags into a `{FBba feature_id: parent FBab feature_id}` map, and then grafts whitelisted FBba data onto the parent `FBAberration` object **before** `synthesize_*` runs — the same technique `add_fbal_to_fbti()` uses for FBti absorbing FBal. Grafting means every existing mapping method (`map_synonyms`, `map_secondary_ids`, `map_entity_props_to_notes`, `map_pubs`) produces the merged output with no changes. The "carries alleles" relationships are the one exception: they are injected directly into `self.aberration_allele_rels` rather than grafted, so they cannot leak into aberration-**gene** association synthesis.

**Tech Stack:** Python 3.11, SQLAlchemy ORM over FlyBase chado (`harvdev_utils.reporting` models), Alliance LinkML DTOs in `src/agr_datatypes.py`.

## Global Constraints

- **PEP8, `max-line-length = 160`** from `.flake8`. `ignore = D100,D103,D104` **replaces** pycodestyle's defaults, so **W503 and W504 are both active**: never break a line before *or* after a boolean operator — wrap the whole condition in a helper or a single line. Docstrings are required on new public methods (D102 is active).
- **No local script runs.** The repo venv is x86_64 and unusable on this arm64 Mac; there is no `flake8`, `sqlalchemy` or `harvdev_utils` on any system Python. Verification is: `python3 -m py_compile`, the offline `ast`-extraction harnesses defined in this plan, SQL against chado, and a real run on flysql26.
- **chado access:** `ssh flysql26` then `PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado`. Run `set standard_conforming_strings = on;` first if a regex uses `\y`, or the escape is silently eaten and the pattern matches nothing.
- **Run artifacts** live in `flysql26:/data/alliance/PERSISTENT_*`. A full allele run takes ~85 minutes.
- **Gating.** Everything that emits an `allele_allele_association_ingest_set` entry stays behind `ADD_ALLELE_ALLELE_ASSOC=YES` (FTA-218: the Alliance has no such ingest set and lacks the `carries`/`breakpoint_allele` CV terms). Everything else in this ticket uses slots that already exist in released LinkML and needs **no** gate.
- **Never commit without explicit permission from Ian.** Each task's commit step is prepared, but only run when he asks.
- **Conventional Commits**, scope `FTA-236`, and `git add <specific paths>` only.

## Open dependency — resolve before Task 2 lands

`FBba0000688`'s note names `FBab0007127`, which is **obsolete** (a non-current secondary ID of `FBab0049550`, the current T(2;3)A1-W and the aberration that carries the FTA-235 `is_balancer` flag). Asked on FTA-236 (comment 43812, 2026-08-14):

- **Option 1 (preferred, assumed by this plan):** Steven fixes the note; the code treats an obsolete/unknown parent as a reported error and skips the merge. Task 1 implements this.
- **Option 2:** the code resolves parent IDs through `feature_dbxref` to the current feature. If Steven picks this, Task 1's `_resolve_parent_feature_id()` gains a secondary-ID fallback — the hook is called out in that task, so the change is one method, not a redesign.

Until Steven answers, expect **37 of 38** merges to succeed and one logged error naming `FBba0000688`.

---

## Verified input data (production_chado, 2026-08-14)

| fact | value |
|---|---|
| FBba with merge flag | 38 (one `internal_notes` prop each) |
| distinct parent FBab named | 37 (`FBab0004818` named twice: `FBba0000039` SM6a + `FBba0000040` SM6b) |
| parents that also carry the FTA-235 flag | 36 (the T(2;3)A1-W pair is the mismatch above) |
| symbol synonyms to move | 96 rows over 38 FBba |
| fullname synonyms to move | 35 rows over 28 FBba |
| FBba secondary IDs (`feature_dbxref`) | 10 rows over 7 FBba (plus the 38 current FBba IDs) |
| `carried_on` alleles (FBal subject → FBba object) | 326 rows, 54 distinct FBal, over 37 FBba |
| `associated_with` insertions (FBba subject → FBti object) | 10 rows over 6 FBba |
| `misc` props (AB6) | 27 rows over 22 FBba |
| `internal_notes` props | 68 rows over 38 FBba (includes the merge flag itself) |
| `feature_pub` references | 440 rows over 38 FBba |
| rows removed by the 3 hard-coded exclusions | 53 of the 326 (`FBba0000011` FM1 = 15, `FBba0000039` SM6a = 20, `FBba0000040` SM6b = 18) |

## The full-name trap (drives a shared-file change in Task 2)

`map_synonyms()` treats a second *current symbol* and a second *current full name* very differently:

- **symbols** (`entity_handler.py:1282-1298`): a chooser prefers the entity's own name, demotes the rest to synonyms, and logs `Found many current symbols`. Recoverable.
- **full names** (`entity_handler.py:1300-1304`): with more than one current full name it logs `Found many current full_names` and **sets no `allele_full_name_dto` at all** — the aberration silently loses its own full name from the export.

28 of the 38 balancers have current `fullname` synonyms, so a naive graft would blank the full name on most merged aberrations. `merged_synonym_ids` does not help: `entity_handler.py:1071` only demotes types in `symbol_type_names`, and that narrowness is deliberate (FTA-234: FBti *should* inherit a merged FBal's full name).

FTA-236 wants neither slot overwritten — the ticket's target for both AB1b and AB2b is `allele_synonyms`. So Task 2 adds one new entity attribute, `demoted_synonym_ids`, honoured for **all** name types, and populates it from the grafted balancer synonyms. Existing behaviour is unchanged because the set is empty for every other entity.

## File Structure

- **Modify `src/allele_handlers.py`** — all handler work lands in `class AberrationHandler`. New class attributes, one retrieval method, one synthesis-map method, one graft method, one association-injection method, and four one-line calls into the existing `get_datatype_data()` / `synthesize_info()` hooks. This file is already ~1800 lines and organised by handler class; follow that structure rather than splitting.
- **Modify `src/fb_datatypes.py`** — one new `FBDataEntity` attribute, `demoted_synonym_ids`, beside `merged_synonym_ids` (line 95).
- **Modify `src/entity_handler.py`** — three lines in `synthesize_synonyms()` honouring that set for every name type.
- **Modify `src/AGR_data_retrieval_curation_allele.py`** — document the behaviour in the `--help` epilog and add the balancer-merge curator TSV.
- **Modify `src/curation_tsv.py`** — nothing. The new TSV is written by a small local writer in the retrieval script only if `write_association_tsv()` cannot be reused; Task 4 checks that first and reuses it if it fits.
- **Create `docs/FTA-236/verify_merge_map.py`** — offline harness for Task 1 (flag parsing + parent resolution).
- **Create `docs/FTA-236/verify_graft.py`** — offline harness for Tasks 2 and 3 (grafting and association injection).
- **Create `docs/FTA-236/expected_counts.sql`** — the SQL that produces the expectation table above (written and verified against production_chado on 2026-08-14; every number in the table came from it), so any run's log/TSV can be re-checked against chado.

---

### Task 1: Parse the merge flags and resolve each parent aberration

**Files:**
- Modify: `src/allele_handlers.py` (class `AberrationHandler`: new class attributes next to `balancer_flag_regex` at ~line 1225; new methods `get_balancer_entities()` and `synthesize_balancer_merge_map()`; calls added in `get_datatype_data()` ~line 1345 and `synthesize_info()` ~line 1548)
- Create: `docs/FTA-236/verify_merge_map.py`
- Create: `docs/FTA-236/expected_counts.sql`

**Interfaces:**
- Consumes: `self.fb_data_entities` (feature_id-keyed `FBAberration`), `self.testing`, `self.log`, `export_chado_data` (already imported at the top of the module), `BalancerHandler` (defined later in the same module — reference it inside the method body, not at import time).
- Produces:
  - `self.fbba_entities` — `{FBba feature_id: FBBalancer}`, every FBba in chado, from the nested handler.
  - `self.balancer_merge_map` — `{FBba feature_id: parent FBab feature_id}`, only successfully resolved merges.
  - `self.balancer_merge_regex` — compiled pattern with named groups `symbol` and `fbab`.
  - `BALANCER_MERGE_EXPECTED = 38` — module-level constant, used by later log warnings.

- [ ] **Step 1: Write the failing test**

Create `docs/FTA-236/verify_merge_map.py`. It extracts the real method source with `ast` and drives it with duck-typed stand-ins, because no package imports work on this machine.

```python
"""Offline harness for AberrationHandler.synthesize_balancer_merge_map() (FTA-236).

Extracts the method source with ast and drives it with duck-typed stand-ins, since
sqlalchemy/harvdev_utils cannot be imported here (see docs/FTA-236/plan.md).
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


class FakeEntity:
    def __init__(self, feature_id, uniquename, internal_notes=(), is_obsolete=False):
        self.db_primary_id = feature_id
        self.uniquename = uniquename
        self.is_obsolete = is_obsolete
        self.props_by_type = {'internal_notes': [FakeProp(v) for v in internal_notes]} if internal_notes else {}

    def __str__(self):
        return self.uniquename


class FakeHandler:
    balancer_merge_regex = ns['balancer_merge_regex']
    synthesize_balancer_merge_map = ns['synthesize_balancer_merge_map']
    _resolve_parent_feature_id = ns['_resolve_parent_feature_id']

    def __init__(self, aberrations, balancers, testing=False):
        self.log = FakeLog()
        self.testing = testing
        self.fb_data_entities = {a.db_primary_id: a for a in aberrations}
        self.fbba_entities = {b.db_primary_id: b for b in balancers}
        self.balancer_merge_map = {}


FLAG = "FTA: Balancer - merge with parent In(1)Basc (FBab0004219)."
FLAG_SM6A = "FTA: Balancer - merge with parent In(2LR)SM6 (FBab0004818)."
FLAG_SM6B = "FTA: Balancer - merge with parent In(2LR)SM6 (FBab0004818)."
FLAG_OBSOLETE = "FTA: Balancer - merge with parent T(2;3)A1-W (FBab0007127)."
FLAG_UNKNOWN = "FTA: Balancer - merge with parent Nope (FBab9999999)."
FLAG_OTHER = "FTA: Balancer - use balancer symbol and fullname for parent In(1)Basc (FBab0004219)."
PROSE = 'Used as a balancer; merge with parent was discussed but never agreed.'

results = []


def check(name, condition):
    results.append((name, condition))
    print(f'  {"PASS" if condition else "FAIL"}: {name}')


print('--- happy path: 4 balancers, 3 parents (SM6 takes two) ---')
aberrations = [
    FakeEntity(1, 'FBab0004219'),
    FakeEntity(2, 'FBab0004818'),
    FakeEntity(3, 'FBab0049550'),
]
balancers = [
    FakeEntity(101, 'FBba0000014', [FLAG]),
    FakeEntity(102, 'FBba0000039', [FLAG_SM6A]),
    FakeEntity(103, 'FBba0000040', [FLAG_SM6B]),
    FakeEntity(104, 'FBba0000999', [PROSE, FLAG_OTHER]),
]
h = FakeHandler(aberrations, balancers)
h.synthesize_balancer_merge_map()
check('3 flagged balancers mapped', h.balancer_merge_map == {101: 1, 102: 2, 103: 2})
check('unflagged balancer ignored', 104 not in h.balancer_merge_map)
check('no errors logged', h.log.errors == [])

print('--- obsolete parent is an error, not a silent skip ---')
aberrations = [FakeEntity(3, 'FBab0049550')]
balancers = [FakeEntity(105, 'FBba0000688', [FLAG_OBSOLETE])]
h = FakeHandler(aberrations, balancers)
h.synthesize_balancer_merge_map()
check('obsolete/absent parent not merged', h.balancer_merge_map == {})
check('error names the balancer', any('FBba0000688' in e for e in h.log.errors))
check('error names the parent ID', any('FBab0007127' in e for e in h.log.errors))

print('--- unknown parent ID is an error ---')
h = FakeHandler([FakeEntity(1, 'FBab0004219')], [FakeEntity(106, 'FBba0000777', [FLAG_UNKNOWN])])
h.synthesize_balancer_merge_map()
check('unknown parent not merged', h.balancer_merge_map == {})
check('error names the unknown ID', any('FBab9999999' in e for e in h.log.errors))

print('--- count warning fires away from the expected 38 ---')
check('count warning present', any('38' in w for w in h.log.warnings))

print('--- testing mode suppresses the count warning ---')
h = FakeHandler([FakeEntity(1, 'FBab0004219')], [FakeEntity(101, 'FBba0000014', [FLAG])], testing=True)
h.synthesize_balancer_merge_map()
check('no count warning in testing mode', not any('38' in w for w in h.log.warnings))

print(f'\n{sum(1 for _, ok in results if ok)}/{len(results)} checks PASSED')
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python3 docs/FTA-236/verify_merge_map.py`
Expected: `StopIteration` from the `next(...)` that looks for `synthesize_balancer_merge_map` — the method does not exist yet.

- [ ] **Step 3: Add the class attributes**

In `src/allele_handlers.py`, at module level near the top (after the imports), add:

```python
BALANCER_MERGE_EXPECTED = 38    # FTA-236: curated merge flags in chado as of 2026-08-14.
```

In `class AberrationHandler`, immediately below the existing `balancer_flag_expected_count` line, add:

```python
    # FTA-236: the merge flag on an FBba balancer names its parent FBab aberration, e.g.
    #   FTA: Balancer - merge with parent In(1)Basc (FBab0004219).
    # 38 of these exist, naming 37 parents (In(2LR)SM6 takes both SM6a and SM6b). Note that FBba
    # balancers also carry a "use balancer symbol and fullname for parent ..." flag (FTA-237, 24 of
    # them) with the same opening, so the match runs through the instruction clause, as in FTA-235.
    balancer_merge_regex = re.compile(r'^\s*:?\s*FTA:\s*balancer\s*-\s*merge with parent\s+(?P<symbol>.+?)\s*\((?P<fbab>FBab[0-9]+)\)\s*\.?\s*$',
                                      re.IGNORECASE)
    # FTA-236: 'carries alleles' must NOT move for these balancers (curator decision, reasons not given).
    balancer_carries_exclusions = {
        'FBba0000011',    # FM1 -> In(1)FM1 (FBab0010486); 15 alleles held back.
        'FBba0000039',    # SM6a -> In(2LR)SM6 (FBab0004818); 20 alleles held back.
        'FBba0000040',    # SM6b -> In(2LR)SM6 (FBab0004818); 18 alleles held back.
    }
```

And with the other "additional export sets" attributes:

```python
    fbba_entities = {}        # FTA-236: feature_id-keyed FBBalancer objects from the nested BalancerHandler.
    balancer_merge_map = {}   # FTA-236: {FBba feature_id: parent FBab feature_id} for resolved merge flags.
    balancer_merge_report = []  # FTA-236: dicts of what moved, for the curator TSV.
```

- [ ] **Step 4: Add the nested retrieval method**

Add to `class AberrationHandler`, in its "Additional sub-methods for get_datatype_data()" area:

```python
    def get_balancer_entities(self, session):
        """Have the AberrationHandler run the BalancerHandler to get FBba data (FTA-236).

        Mirrors AlleleHandler.get_insertion_entities(): the nested handler does its own full
        retrieval, and only its fb_data_entities are kept. Nothing it maps or exports is used, so
        FBba balancers still do not appear in the allele_ingest_set - their data is grafted onto the
        parent FBab aberrations instead.
        """
        separator = 80 * '#'
        self.log.info(f'Have the AberrationHandler run the BalancerHandler.\n{separator}')
        balancer_handler = BalancerHandler(self.log, self.testing)
        export_chado_data(session, self.log, balancer_handler)
        self.fbba_entities = balancer_handler.fb_data_entities
        self.log.info(f'The AberrationHandler obtained {len(self.fbba_entities)} FBba balancers from chado.\n{separator}')
        return
```

- [ ] **Step 5: Add the merge-map method**

Add to `class AberrationHandler`, in its "Additional sub-methods to be run by synthesize_info()" area:

```python
    def _resolve_parent_feature_id(self, fbab_uniquename):
        """Return the feature_id of an exportable parent FBab, or None (FTA-236).

        Only aberrations in self.fb_data_entities can receive merged data: that dict holds the
        current, non-obsolete FBab features this handler exports. An obsolete or unrecognised ID
        gets no fallback on purpose - see the "Open dependency" section of docs/FTA-236/plan.md. If
        the curators decide stale IDs should be followed instead, resolve them through
        feature_dbxref here and nowhere else.
        """
        for aberration in self.fb_data_entities.values():
            if aberration.uniquename == fbab_uniquename:
                return aberration.db_primary_id
        return None

    def synthesize_balancer_merge_map(self):
        """Map each flagged FBba balancer to its parent FBab aberration (FTA-236)."""
        self.log.info('Map flagged FBba balancers to their parent FBab aberrations.')
        flag_counter = 0
        for balancer in self.fbba_entities.values():
            for fb_prop in balancer.props_by_type.get('internal_notes', []):
                prop_value = fb_prop.chado_obj.value
                if not prop_value:
                    continue
                match = self.balancer_merge_regex.match(prop_value)
                if match is None:
                    continue
                flag_counter += 1
                fbab_uniquename = match.group('fbab')
                parent_feature_id = self._resolve_parent_feature_id(fbab_uniquename)
                if parent_feature_id is None:
                    self.log.error(f'FTA-236: balancer {balancer} names parent {fbab_uniquename} '
                                   f'("{match.group("symbol")}"), which is not an exportable aberration '
                                   f'(obsolete or unknown). No data merged; please fix the internal note.')
                    break
                self.balancer_merge_map[balancer.db_primary_id] = parent_feature_id
                self.log.debug(f'FTA-236: {balancer} merges into {self.fb_data_entities[parent_feature_id]}.')
                break
        parents = set(self.balancer_merge_map.values())
        self.log.info(f'Found {flag_counter} balancer merge flags; resolved {len(self.balancer_merge_map)} '
                      f'balancers onto {len(parents)} parent aberrations.')
        if self.testing is False and flag_counter != BALANCER_MERGE_EXPECTED:
            self.log.warning(f'Expected {BALANCER_MERGE_EXPECTED} balancer merge flags per FTA-236, but found '
                             f'{flag_counter}. Check the "internal_notes" flag text in chado.')
        return
```

- [ ] **Step 6: Wire both methods into the handler hooks**

In `AberrationHandler.get_datatype_data()`, add as the **last** call before `return`:

```python
        self.get_balancer_entities(session)
```

In `AberrationHandler.synthesize_info()`, add immediately after `super().synthesize_info()` and **before** `self.flag_new_additions_and_obsoletes()`:

```python
        self.synthesize_balancer_merge_map()
```

- [ ] **Step 7: Run the harness to verify it passes**

Run: `python3 docs/FTA-236/verify_merge_map.py`
Expected: `10/10 checks PASSED`

- [ ] **Step 8: Compile check and line length**

Run:
```bash
python3 -m py_compile src/allele_handlers.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py
grep -n ' $' src/allele_handlers.py
```
Expected: `COMPILE OK`, and no output from `awk` or `grep`.

- [ ] **Step 9: Record the chado expectations**

Create `docs/FTA-236/expected_counts.sql` holding the queries that produced the "Verified input data" table (flag counts by family, parents per balancer, per-field row counts, exclusion counts). Run it and confirm the numbers still match:

```bash
scp docs/FTA-236/expected_counts.sql flysql26:/tmp/
ssh flysql26 "PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado -f /tmp/expected_counts.sql"
```
Expected: 38 merge flags, 37 distinct parents, `FBab0004818` with 2 balancers.

- [ ] **Step 10: Commit (only when Ian asks)**

```bash
git add src/allele_handlers.py docs/FTA-236/plan.md docs/FTA-236/verify_merge_map.py docs/FTA-236/expected_counts.sql
git commit -m "feat(FTA-236): map flagged FBba balancers to their parent FBab aberrations"
```

---

### Task 2: Graft balancer names, secondary IDs, notes and references onto the parent

**Files:**
- Modify: `src/fb_datatypes.py:95` (new `demoted_synonym_ids` attribute)
- Modify: `src/entity_handler.py:1071` (honour it in `synthesize_synonyms()`)
- Modify: `src/allele_handlers.py` (class `AberrationHandler`: new methods `add_fbba_to_fbab()` and `merge_balancers_into_aberrations()`; one call in `synthesize_info()`)
- Create: `docs/FTA-236/verify_graft.py`

**Interfaces:**
- Consumes: `self.balancer_merge_map` and `self.fbba_entities` from Task 1.
- Produces: `FBDataEntity.demoted_synonym_ids` — a `set()` of `synonym_id`s that `synthesize_synonyms()` forces to `is_current = False` regardless of name type.
- Produces: `add_fbba_to_fbab(balancer, aberration)` — mutates the aberration in place; `merge_balancers_into_aberrations()` — walks the map and calls it; appends one dict per balancer to `self.balancer_merge_report` with keys `fbba_id`, `fbba_symbol`, `fbab_id`, `synonyms`, `secondary_ids`, `comments`, `internal_notes`, `references`, `carries_alleles`, `carries_excluded` (the last two are filled by Task 3, defaulting to 0).

- [ ] **Step 1: Write the failing test**

Create `docs/FTA-236/verify_graft.py`:

```python
"""Offline harness for AberrationHandler.add_fbba_to_fbab() (FTA-236)."""

import ast
import re

SRC = '/Users/ilongden/harvard/alliance-linkml-flybase/src/allele_handlers.py'

tree = ast.parse(open(SRC).read())
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'AberrationHandler')
ns = {'re': re}
for node in cls.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') in ('balancer_graft_prop_types',):
        exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)
for name in ('add_fbba_to_fbab', 'merge_balancers_into_aberrations'):
    node = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[node], type_ignores=[]), SRC, 'exec'), ns)


class FakeLog:
    def __init__(self):
        self.lines, self.warnings = [], []

    def info(self, msg):
        self.lines.append(msg)

    def debug(self, msg):
        pass

    def warning(self, msg):
        self.warnings.append(msg)

    def error(self, msg):
        self.warnings.append(msg)


class FakeProp:
    def __init__(self, value):
        self.chado_obj = type('O', (), {'value': value})()
        self.pubs = []


class FakeSynonym:
    def __init__(self, synonym_id, name):
        self.synonym_id = synonym_id
        self.name = name


class FakeEntity:
    def __init__(self, feature_id, uniquename, name=''):
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

    def __str__(self):
        return f'{self.name} ({self.uniquename})'


class FakeHandler:
    balancer_graft_prop_types = ns['balancer_graft_prop_types']
    add_fbba_to_fbab = ns['add_fbba_to_fbab']
    merge_balancers_into_aberrations = ns['merge_balancers_into_aberrations']

    def __init__(self, aberrations, balancers, merge_map, testing=False):
        self.log = FakeLog()
        self.testing = testing
        self.fb_data_entities = {a.db_primary_id: a for a in aberrations}
        self.fbba_entities = {b.db_primary_id: b for b in balancers}
        self.balancer_merge_map = merge_map
        self.balancer_merge_report = []


results = []


def check(name, condition):
    results.append((name, condition))
    print(f'  {"PASS" if condition else "FAIL"}: {name}')


parent = FakeEntity(1, 'FBab0004219', 'In(1)Basc')
parent.synonyms = [FakeSynonym(900, 'In(1)Basc')]
parent.props_by_type = {'misc': [FakeProp('Parent comment.')]}
parent.pub_associations = ['parent_pub']

balancer = FakeEntity(101, 'FBba0000014', 'Basc')
balancer.synonyms = [FakeSynonym(910, 'Basc'), FakeSynonym(911, 'M5')]
balancer.alt_fb_ids = []
balancer.fb_sec_dbxrefs = ['fbba_sec_dbxref']
balancer.pub_associations = ['fbba_pub_1', 'fbba_pub_2']
balancer.props_by_type = {
    'misc': [FakeProp('N[opa33b] is a natural variant ...')],
    'internal_notes': [FakeProp("FTA: Balancer - merge with parent In(1)Basc (FBab0004219).")],
    'availability': [FakeProp('Stated to be lost.')],      # must NOT move: would set is_extinct
    'derived_stock_count': [FakeProp('12')],               # must NOT move
}

h = FakeHandler([parent], [balancer], {101: 1})
h.merge_balancers_into_aberrations()

print('--- names ---')
check('balancer synonyms appended', {s.synonym_id for s in parent.synonyms} == {900, 910, 911})
check('balancer synonym ids marked as merged', parent.merged_synonym_ids == {910, 911})
check('parent own synonym not marked merged', 900 not in parent.merged_synonym_ids)
check('balancer synonym ids also demoted (protects the full name slot)', parent.demoted_synonym_ids == {910, 911})

print('--- identifiers ---')
check('balancer primary ID added as secondary', 'FB:FBba0000014' in parent.alt_fb_ids)
check('balancer secondary dbxrefs moved', 'fbba_sec_dbxref' in parent.fb_sec_dbxrefs)

print('--- notes: whitelist only ---')
check('misc props merged, parent prop kept', len(parent.props_by_type['misc']) == 2)
check('internal_notes props merged', len(parent.props_by_type['internal_notes']) == 1)
check('availability NOT merged', 'availability' not in parent.props_by_type)
check('derived_stock_count NOT merged', 'derived_stock_count' not in parent.props_by_type)

print('--- references ---')
check('balancer pubs merged', parent.pub_associations == ['parent_pub', 'fbba_pub_1', 'fbba_pub_2'])

print('--- report row ---')
row = h.balancer_merge_report[0]
check('report identifies both entities', row['fbba_id'] == 'FBba0000014' and row['fbab_id'] == 'FBab0004219')
check('report counts synonyms', row['synonyms'] == 2)
check('report counts comments', row['comments'] == 1)
check('report counts references', row['references'] == 2)

print('--- two balancers into one parent (SM6 case) ---')
parent2 = FakeEntity(2, 'FBab0004818', 'In(2LR)SM6')
sm6a = FakeEntity(102, 'FBba0000039', 'SM6a')
sm6a.synonyms = [FakeSynonym(920, 'SM6a')]
sm6b = FakeEntity(103, 'FBba0000040', 'SM6b')
sm6b.synonyms = [FakeSynonym(930, 'SM6b')]
h2 = FakeHandler([parent2], [sm6a, sm6b], {102: 2, 103: 2})
h2.merge_balancers_into_aberrations()
check('both balancers merged into one parent', {s.synonym_id for s in parent2.synonyms} == {920, 930})
check('two report rows written', len(h2.balancer_merge_report) == 2)

print(f'\n{sum(1 for _, ok in results if ok)}/{len(results)} checks PASSED')
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python3 docs/FTA-236/verify_graft.py`
Expected: `StopIteration` — `add_fbba_to_fbab` does not exist yet.

- [ ] **Step 3: Add the `demoted_synonym_ids` attribute and honour it**

In `src/fb_datatypes.py`, immediately after the `merged_synonym_ids` line (95):

```python
        self.demoted_synonym_ids = set()  # synonym_ids to force non-current for ALL name types (FTA-236 balancer merges).
```

In `src/entity_handler.py`, in `synthesize_synonyms()`, immediately after the existing `is_merged_symbol` demotion block (the two lines at 1071-1073), add:

```python
                # FTA-236: a balancer's names must never rival the parent aberration's own. Unlike
                # merged_synonym_ids this covers full names too: map_synonyms() sets no
                # allele_full_name_dto at all when it finds two current full names, so an
                # undemoted balancer full name would delete the aberration's own from the export.
                if syno_id in fb_data_entity.demoted_synonym_ids:
                    syno_dict['is_current'] = False
```

- [ ] **Step 4: Add the prop whitelist attribute**

In `class AberrationHandler`, next to the other FTA-236 attributes:

```python
    # FTA-236: only these FBba prop types move to the parent aberration. A blanket props_by_type merge
    # (as add_fbal_to_fbti does) would also carry 'availability' and 'derived_stock*' props, which
    # map_extinction_info() reads to set is_extinct - a balancer's stock status is not the
    # aberration's. The Alliance note type for each comes from aberration_prop_to_note_mapping.
    balancer_graft_prop_types = ('misc', 'internal_notes')
```

- [ ] **Step 5: Write the graft methods**

```python
    def add_fbba_to_fbab(self, balancer, aberration):
        """Graft one FBba balancer's data onto its parent FBab aberration (FTA-236).

        Follows add_fbal_to_fbti(): extend the parent's source-data lists before synthesize_* runs, so
        every existing mapping method produces the merged output unchanged. Only whitelisted prop
        types move (see balancer_graft_prop_types), and relationships are handled separately by
        synthesize_balancer_carries_associations() so they cannot reach aberration-gene synthesis.
        """
        report = {
            'fbba_id': balancer.uniquename,
            'fbba_symbol': balancer.name,
            'fbab_id': aberration.uniquename,
            'synonyms': len(balancer.synonyms),
            'secondary_ids': 1 + len(balancer.fb_sec_dbxrefs),
            'comments': len(balancer.props_by_type.get('misc', [])),
            'internal_notes': len(balancer.props_by_type.get('internal_notes', [])),
            'references': len(balancer.pub_associations),
            'carries_alleles': 0,
            'carries_excluded': 0,
        }
        # Names. Per the ticket, the balancer's symbols and full names become synonyms of the parent;
        # neither may rival the parent's own. merged_synonym_ids covers symbols (a rival current
        # symbol makes map_synonyms() pick one arbitrarily and log "Found many current symbols");
        # demoted_synonym_ids also covers full names, where a rival is worse - map_synonyms() sets no
        # allele_full_name_dto at all and the aberration's own full name vanishes from the export.
        # FTA-237 deliberately renames 24 of these aberrations to the balancer symbol/full name; it
        # must do that by setting the slots, not by skipping this demotion.
        aberration.synonyms.extend(balancer.synonyms)
        balancer_synonym_ids = [i.synonym_id for i in balancer.synonyms]
        aberration.merged_synonym_ids.update(balancer_synonym_ids)
        aberration.demoted_synonym_ids.update(balancer_synonym_ids)
        # Identifiers: the balancer's current ID plus any of its own secondary IDs.
        aberration.alt_fb_ids.append(f'FB:{balancer.uniquename}')
        aberration.fb_sec_dbxrefs.extend(balancer.fb_sec_dbxrefs)
        # Notes (whitelisted prop types only).
        for prop_type in self.balancer_graft_prop_types:
            balancer_props = balancer.props_by_type.get(prop_type, [])
            if not balancer_props:
                continue
            try:
                aberration.props_by_type[prop_type].extend(balancer_props)
            except KeyError:
                aberration.props_by_type[prop_type] = list(balancer_props)
        # References: synthesize_pubs() reads pub_associations, so graft before it runs.
        aberration.pub_associations.extend(balancer.pub_associations)
        self.balancer_merge_report.append(report)
        return report

    def merge_balancers_into_aberrations(self):
        """Graft every flagged FBba balancer onto its parent FBab aberration (FTA-236)."""
        self.log.info('Merge flagged FBba balancer data into parent FBab aberrations.')
        totals = {'synonyms': 0, 'secondary_ids': 0, 'comments': 0, 'internal_notes': 0, 'references': 0}
        for fbba_feature_id, fbab_feature_id in self.balancer_merge_map.items():
            balancer = self.fbba_entities[fbba_feature_id]
            aberration = self.fb_data_entities[fbab_feature_id]
            report = self.add_fbba_to_fbab(balancer, aberration)
            for key in totals.keys():
                totals[key] += report[key]
        parents = set(self.balancer_merge_map.values())
        self.log.info(f'Merged {len(self.balancer_merge_map)} balancers into {len(parents)} aberrations.')
        for key in sorted(totals.keys()):
            self.log.info(f'Moved {totals[key]} {key} from balancers to aberrations.')
        return
```

- [ ] **Step 6: Wire it into `synthesize_info()`**

Add immediately after `self.synthesize_balancer_merge_map()` (from Task 1) and **before** `self.synthesize_secondary_ids()`:

```python
        self.merge_balancers_into_aberrations()
```

Order matters and is the whole reason grafting works: `synthesize_secondary_ids()`, `synthesize_synonyms()` and `synthesize_pubs()` all run later in the same method and consume the lists just extended.

- [ ] **Step 7: Run the harness to verify it passes**

Run: `python3 docs/FTA-236/verify_graft.py`
Expected: `17/17 checks PASSED`

- [ ] **Step 8: Compile check and line length**

Run:
```bash
python3 -m py_compile src/allele_handlers.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py
```
Expected: `COMPILE OK` and no `awk` output.

- [ ] **Step 9: Commit (only when Ian asks)**

```bash
git add src/allele_handlers.py docs/FTA-236/verify_graft.py
git commit -m "feat(FTA-236): graft balancer names, IDs, notes and refs onto parent aberrations"
```

---

### Task 3: Move the balancer's "carries alleles" and transposon insertions (gated)

**Files:**
- Modify: `src/allele_handlers.py` (class `AberrationHandler`: new method `synthesize_balancer_carries_associations()`; one gated call in `synthesize_info()`)
- Modify: `docs/FTA-236/verify_graft.py` (add the association checks)

**Interfaces:**
- Consumes: `self.balancer_merge_map`, `self.fbba_entities`, `self.balancer_carries_exclusions`, `self.aberration_allele_rels` and `self.feature_subtypes` (all existing), plus `recall_relationships()` on the FBba entity.
- Produces: extra keys in `self.aberration_allele_rels`, shaped exactly as `synthesize_aberration_allele_associations()` builds them — `(parent FBab feature_id, partner feature_id, 'carries')` keying a list of `FBRelationship` objects — so `map_aberration_allele_associations()` needs no change. Updates `carries_alleles` / `carries_excluded` in each `balancer_merge_report` row.

- [ ] **Step 1: Write the failing test**

Append to `docs/FTA-236/verify_graft.py` (extract the new method the same way, then):

```python
print('--- carries associations: injected under the parent, exclusions honoured ---')


class FakeRel:
    def __init__(self, rel_id, subject_id, object_id, pubs=()):
        self.chado_obj = type('O', (), {'subject_id': subject_id, 'object_id': object_id})()
        self.rel_id = rel_id
        self.pubs = list(pubs)


class RelEntity(FakeEntity):
    def __init__(self, feature_id, uniquename, name='', rels=None):
        super().__init__(feature_id, uniquename, name)
        self._rels = rels or {}

    def recall_relationships(self, log, **kwargs):
        key = (kwargs.get('entity_role'), kwargs.get('rel_types'))
        return self._rels.get(key, [])


basc = RelEntity(101, 'FBba0000014', 'Basc', rels={
    ('object', 'carried_on'): [FakeRel(1, 5001, 101), FakeRel(2, 5002, 101)],
    ('subject', 'associated_with'): [FakeRel(3, 101, 6001)],
})
fm1 = RelEntity(111, 'FBba0000011', 'FM1', rels={
    ('object', 'carried_on'): [FakeRel(4, 5003, 111)],
})
parent_basc = FakeEntity(1, 'FBab0004219', 'In(1)Basc')
parent_fm1 = FakeEntity(10, 'FBab0010486', 'In(1)FM1')

h3 = FakeHandler([parent_basc, parent_fm1], [basc, fm1], {101: 1, 111: 10})
h3.feature_subtypes = {'allele': ['allele'], 'insertion': ['insertion']}
h3.aberration_allele_rels = {}
h3.balancer_carries_exclusions = {'FBba0000011', 'FBba0000039', 'FBba0000040'}
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
check('report counts excluded associations', fm1_row['carries_excluded'] == 1)
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python3 docs/FTA-236/verify_graft.py`
Expected: `StopIteration` — `synthesize_balancer_carries_associations` does not exist yet.

- [ ] **Step 3: Write the method**

```python
    def synthesize_balancer_carries_associations(self):
        """Move a flagged balancer's carried alleles and insertions onto its parent (FTA-236).

        In chado the two sources are:
            1. type="carried_on",      subject=FBal, object=FBba (proforma AB5b) -> "carries"
            2. type="associated_with", subject=FBba, object=FBti (proforma AB5a) -> "carries"
        Both become Alliance "carries" associations from the parent FBab to the partner, so the keys
        built here match synthesize_aberration_allele_associations() exactly and the existing
        map_aberration_allele_associations() emits them with no changes.

        These relationships are deliberately NOT grafted onto the parent entity: the parent's own
        relationships also drive aberration-GENE association synthesis, and a balancer's links would
        show up there as phantom gene associations.

        Gated with the rest of the allele-allele associations (FTA-218): the Alliance has no
        "allele_allele_association_ingest_set" yet.
        """
        self.log.info('Move balancer carried alleles and insertions onto parent aberrations.')
        rel_specs = [
            ('object', 'carried_on', self.feature_subtypes['allele'], 'AB5b'),
            ('subject', 'associated_with', self.feature_subtypes['insertion'], 'AB5a'),
        ]
        partner_attr = {'subject': 'object_id', 'object': 'subject_id'}
        moved = 0
        excluded = 0
        for fbba_feature_id, fbab_feature_id in self.balancer_merge_map.items():
            balancer = self.fbba_entities[fbba_feature_id]
            report_rows = [i for i in self.balancer_merge_report if i['fbba_id'] == balancer.uniquename]
            is_excluded = balancer.uniquename in self.balancer_carries_exclusions
            for entity_role, fb_rel_type, rel_entity_types, proforma_field in rel_specs:
                relevant_rels = balancer.recall_relationships(self.log, entity_role=entity_role, rel_types=fb_rel_type,
                                                              rel_entity_types=rel_entity_types)
                for feat_rel in relevant_rels:
                    partner_id = getattr(feat_rel.chado_obj, partner_attr[entity_role])
                    if is_excluded:
                        excluded += 1
                        for row in report_rows:
                            row['carries_excluded'] += 1
                        continue
                    rel_key = (fbab_feature_id, partner_id, 'carries')
                    try:
                        self.aberration_allele_rels[rel_key].append(feat_rel)
                    except KeyError:
                        self.aberration_allele_rels[rel_key] = [feat_rel]
                        moved += 1
                        for row in report_rows:
                            row['carries_alleles'] += 1
        self.log.info(f'Moved {moved} balancer "carries" associations onto parent aberrations.')
        self.log.info(f'Held back {excluded} associations for the {len(self.balancer_carries_exclusions)} '
                      f'balancers excluded by FTA-236.')
        return
```

- [ ] **Step 4: Wire it into `synthesize_info()`**

Inside the existing `ADD_ALLELE_ALLELE_ASSOC` block, so it runs **after** `synthesize_aberration_allele_associations()`:

```python
        if getenv('ADD_ALLELE_ALLELE_ASSOC', None) == 'YES':
            self.synthesize_aberration_allele_associations()
            self.synthesize_balancer_carries_associations()
        else:
            self.log.info('ADD_ALLELE_ALLELE_ASSOC not set to "YES"; skipping aberration-allele association synthesis.')
```

Also extend the `get_general_data()` comment above the gated `build_feature_lookup(...)` call to note that FTA-236's balancer partners rely on the same lookup — the FBal/FBti partners of a balancer must be in `feature_lookup` for `map_aberration_allele_associations()` to resolve their uniquenames.

- [ ] **Step 5: Run the harness to verify it passes**

Run: `python3 docs/FTA-236/verify_graft.py`
Expected: `24/24 checks PASSED` (17 from Task 2 plus 7 new)

- [ ] **Step 6: Compile check and line length**

Run:
```bash
python3 -m py_compile src/allele_handlers.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py
```
Expected: `COMPILE OK` and no `awk` output.

- [ ] **Step 7: Commit (only when Ian asks)**

```bash
git add src/allele_handlers.py docs/FTA-236/verify_graft.py
git commit -m "feat(FTA-236): move balancer carried alleles and insertions to parent aberrations"
```

---

### Task 4: Curator TSV, docs, and full-run verification

**Files:**
- Modify: `src/AGR_data_retrieval_curation_allele.py` (epilog text; new TSV emission in `main()` after `write_notes_tsv`)
- Modify: `src/curation_tsv.py` (only if a new writer is needed — check `write_association_tsv()` first)
- Create: `docs/FTA-236/verify_run.sh`

**Interfaces:**
- Consumes: `aberration_handler.balancer_merge_report` (list of dicts from Tasks 2-3).
- Produces: `<output>_balancer_merges.tsv` beside the other allele TSVs.

- [ ] **Step 1: Check whether an existing writer fits**

Run: `grep -n "def write_" src/curation_tsv.py`
`write_association_tsv()` takes `first_field`/`second_field`/`extra_fields` over a list of dicts, so it fits: `first_field='fbab_id'`, `second_field='fbba_id'`, and the counts as `extra_fields`. Reuse it — do not add a new writer.

**But it reads `entity_dict['relation_name']` unconditionally**, so a report row without that key raises `KeyError` (found while implementing FTA-237, which hit exactly this). Pass a constant relation and blank the evidence sentinel:

```python
            rows=[dict(row, relation_name='merged_from_balancer') for row in aberration_handler.balancer_merge_report],
            no_pubs_sentinel='',
```

- [ ] **Step 2: Emit the TSV**

In `main()`, inside the `else` branch that writes the other TSVs, after the `write_notes_tsv(...)` call:

```python
        # FTA-236: one row per FBba balancer merged into a parent FBab aberration, so curators can
        # check what moved. Written unconditionally; the 'carries' counts are 0 unless
        # ADD_ALLELE_ALLELE_ASSOC=YES, since that gate controls the association export.
        merges_filename = tsv_filename.replace('.tsv', '_balancer_merges.tsv')
        curation_tsv.write_association_tsv(
            filename=merges_filename,
            rows=[dict(row, relation_name='merged_from_balancer') for row in aberration_handler.balancer_merge_report],
            first_field='fbab_id',
            second_field='fbba_id',
            extra_fields=[
                ('fbba_symbol', 'fbba_symbol', ''),
                ('synonyms', 'synonyms', 0),
                ('secondary_ids', 'secondary_ids', 0),
                ('comments', 'comments', 0),
                ('internal_notes', 'internal_notes', 0),
                ('references', 'references', 0),
                ('carries_alleles', 'carries_alleles', 0),
                ('carries_excluded', 'carries_excluded', 0),
            ],
            no_pubs_sentinel='',
        )
        log.info(f'Generated TSV: {merges_filename} ({len(aberration_handler.balancer_merge_report)} balancer merges)')
```

- [ ] **Step 3: Document the behaviour in the epilog**

Extend the `ADD_ALLELE_ALLELE_ASSOC` entry to mention FTA-236, and add a plain-English note that balancer merges themselves are ungated:

```
  ADD_ALLELE_ALLELE_ASSOC   Set to 'YES' to emit the FTA-218 'allele_allele_association_ingest_set' (aberration
                            'carries'/'breakpoint_allele' relations, plus the FTA-236 alleles and insertions
                            carried by a balancer and moved to its parent aberration). Off by default: the
                            Alliance schema has no such ingest set and lacks the two CV terms, so the data
                            cannot yet be loaded. FTA-236's names, IDs, notes and references are NOT gated -
                            they use slots that already exist in released LinkML.
```

- [ ] **Step 4: Compile check both files**

Run:
```bash
python3 -m py_compile src/allele_handlers.py src/AGR_data_retrieval_curation_allele.py && echo COMPILE OK
awk 'length > 160 {print FILENAME": "FNR": "length}' src/allele_handlers.py src/AGR_data_retrieval_curation_allele.py
```
Expected: `COMPILE OK` and no `awk` output.

- [ ] **Step 5: Write the run-verification script**

Create `docs/FTA-236/verify_run.sh`, taking the run directory as `$1`, to be pointed at a `flysql26:/data/alliance/PERSISTENT_*` copy. It must check:

```bash
#!/usr/bin/env bash
# Usage: verify_run.sh /data/alliance/PERSISTENT_main_FTA-236_balancer_merge   (run on flysql26)
set -euo pipefail
DIR="$1"
LOG="$DIR"/allele_curation_*.log
echo '--- merge tallies from the log ---'
grep -E 'balancer merge flags|Merged .* balancers into|Moved .* from balancers|Moved .* balancer "carries"|Held back' $LOG
echo '--- errors (expect one only if Steven has not fixed FBba0000688) ---'
grep -c ' ERROR ' $LOG || true
grep ' ERROR ' $LOG || true
echo '--- merge TSV row count (expect 38, or 37 pending the FBba0000688 fix) ---'
tail -n +2 "$DIR"/*_balancer_merges.tsv | wc -l
echo '--- In(1)Basc must gain M5 as a synonym and FBba0000014 as a secondary ID ---'
grep -m1 'FBab0004219' "$DIR"/allele_curation_*.tsv
echo '--- In(2LR)SM6 must show two merge rows ---'
grep -c 'FBab0004818' "$DIR"/*_balancer_merges.tsv
echo '--- excluded balancers must move 0 carries and hold back 15/20/18 ---'
grep -E 'FBba0000011|FBba0000039|FBba0000040' "$DIR"/*_balancer_merges.tsv
echo '--- no aberration should report many current symbols (the merged-synonym trap) ---'
grep -c 'many current symbols' $LOG || true
```

- [ ] **Step 6: Ask Ian to run the export, then verify**

The run cannot happen on this machine. Ask Ian to run the allele export against `production_chado` twice, into `PERSISTENT_main_FTA-236_balancer_merge` (default gate) and `PERSISTENT_main_FTA-236_balancer_merge_assoc` (`ADD_ALLELE_ALLELE_ASSOC=YES ADD_IS_ABERRATION=YES`). Then:

```bash
ssh flysql26 "bash -s" < docs/FTA-236/verify_run.sh /data/alliance/PERSISTENT_main_FTA-236_balancer_merge
```

Expected, per the verified-input table: 38 merge flags (37 merges + 1 error until the note is fixed), 96 synonyms, 45 secondary IDs (38 primary + 7 balancers' own), 27 comments, 68 internal notes, 440 references, `many current symbols` count of 0, and in the gated run 273 moved `carries` associations with 53 held back.

- [ ] **Step 7: Cross-check the JSON for one worked example**

The ticket's worked example is Basc → In(1)Basc. Confirm in the gated run's JSON that `FB:FBab0004219` has: `M5` among `allele_synonym_dtos`, `FB:FBba0000014` among `allele_secondary_id_dtos`, the `N[opa33b]` comment as a `comment` note, and `carries` associations to `Bar[1]`, `sc[8]`, `sc[S1]` and `w[a]`:

```bash
RUN=/data/alliance/PERSISTENT_main_FTA-236_balancer_merge_assoc
ssh flysql26 "grep -n 'FB:FBab0004219' $RUN/allele_curation_production_chado.json | head -3"
# take the first line number N from that output, then read the record around it:
ssh flysql26 "sed -n '53583740,53583860p' $RUN/allele_curation_production_chado.json"   # substitute N-400,N+50
```

`docs/FTA-236/expected_counts.sql` query 7 prints the same facts straight from chado, so the two can be
compared line by line: Basc carries exactly `Bar[1]`, `sc[8]`, `sc[S1]` and `w[a]`, and its AB6 comment is
the `N[opa33b]` sentence (stored with `@...@` FlyBase markup, which `clean_note_free_text()` strips).

- [ ] **Step 8: Commit (only when Ian asks)**

```bash
git add src/AGR_data_retrieval_curation_allele.py docs/FTA-236/verify_run.sh
git commit -m "feat(FTA-236): curator TSV of balancer merges, plus env var docs"
```

---

## Deferred, with reasons

- **`is_balancer` on merged parents** — already true for 36 of the 37 parents via FTA-235 (merged, `d8155b6`). The 37th is the T(2;3)A1-W mismatch; fixing `FBba0000688`'s note resolves both tickets at once, so no code here.
- **FTA-237 (rename 24 aberrations to the balancer symbol/full name)** — a separate ticket, but it collides with Task 2: this plan demotes every merged balancer synonym so the parent keeps its own symbol. FTA-237 must set the symbol slot explicitly for its 24 aberrations rather than removing the demotion, or the export will log "many current symbols" and pick a symbol arbitrarily. Note this in FTA-237 before starting it.
- **Uploading the result** — blocked with the rest of FTA-218 until the Alliance has an `allele_allele_association_ingest_set` and the `carries` CV term.
