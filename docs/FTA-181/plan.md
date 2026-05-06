# FTA-181 — `ConstructCassetteAssociationDTO` linking anon `<FBti>_con` → anon `<FBti>_cas`

## Context

FTA-181 (sub-task of FTA-136, priority Highest, "Ready for work") asks for two
things, per Gillian's writeup (ticket description edited 2026-04-17, with
attachment `FTA-181_anon_con_anon_cas_associations.tsv`, 5409 rows):

**A. Anonymous Cassette per anonymous Construct** — for each `<FBti>_con`
Construct created by FTA-180, emit an anonymous Cassette with:

- `placeholder = true`
- `cassette_symbol_dto` symbol = `<FBti>_cas`
- `primary_external_id = <FBti>_cas`
- `obsolete = false`
- timestamps "same as FTA-141"

**B. `ConstructCassetteAssociationDTO`** linking the two anons:

```
relation_name   = 'has_component'
cons_curie      = <FBti>_con
cassette_curie  = <FBti>_cas
pub_curies      = pubs on the `FBti producedby FBtp` feature_relationship
                  (same set used in the FTA-180 AlleleConstructAssociationDTO
                   for FBti contains FBti_con)
```

**Status of each part in the current code:**

- **Part A is already done.** The cassette retrieval script routes generic-TI
  anon data into the existing anon-cassette pipeline
  (`AGR_data_retrieval_curation_cassette.py:228-233`), which calls
  `CassetteHandler.map_anon_cassette_basic()` at `cassette_handler.py:619`.
  That method already produces a `CassetteDTO` with
  `primary_external_id = FB:{insertion_uniquename}_cas`, `placeholder=True`,
  `obsolete=False`, plus `cassette_symbol_dto` and `data_provider_dto`, and
  pushes it into `anon_cassettes` → `cassette_ingest_set`. No changes
  required for A.

- **Part B is missing.** No code currently emits
  `ConstructCassetteAssociationDTO(FB:<FBti>_con, 'has_component', FB:<FBti>_cas, pubs)`.
  - `map_anon_cassette_to_construct_association` (`construct_handler.py:1191`) only
    runs for FBtp-derived anon cassettes (`construct.needs_anon_cassette`), not
    generic-TI. Its relation_name branch is `has_functional_unit` / `has_component`
    keyed off `tool_uses_data`, and it points at the FBtp's own curie, not the
    `<FBti>_con` anon construct.
  - `map_generic_ti_anon_marker_associations` (`construct_handler.py:1238`) is
    FTA-183 (markers on the anon construct) — different DTOs,
    `has_selectable_marker` / `has_transcriptional_unit`, no cassette.

Part B is the only new work.

## Recommended approach

Mirror the existing `map_generic_ti_anon_marker_associations` pattern in
`ConstructHandler`: one new emitter method, reusing
`construct_cassette_associations` (the list that already drains into
`construct_cassette_association_ingest_set` at `construct_handler.py:1313-1314`).

Collect the `producedby` pubs inside `ConstructHandler` so the cassette
retrieval script remains the single owner of generic-TI cassette plumbing —
do **not** reach across to `AlleleHandler.map_generic_ti_anon_construct_associations`
even though it already walks the same `FBti producedby FBtp` chain (different
script, different handler instance).

### Files modified

| File | Change |
|---|---|
| `src/construct_handler.py` | Extend `get_generic_ti_insertions` to carry `feature_relationship_id`; new method `get_generic_ti_producedby_pubs(session)`; new `__init__` field `ti_producedby_pubs: {ti_fid: [pub_id, ...]}`; new emitter method `map_generic_ti_anon_con_cas_associations(anon_data)`; three wiring lines in `query_chado_and_export`. |

No change to `cassette_handler.py`, `AGR_data_retrieval_curation_cassette.py`,
`agr_datatypes.py`, or `allele_handlers.py`. `ConstructCassetteAssociationDTO`
already exists (`agr_datatypes.py:345`).

### 1. Extend `get_generic_ti_insertions` with `feature_relationship_id`

At `construct_handler.py:644-687`, add one column to the query and one field to
the dict:

```python
results = session.query(
    Feature.uniquename,
    Feature.feature_id,
    FeatureRelationship.object_id,
    FeatureRelationship.feature_relationship_id,   # NEW
).\
    ...
for ins_uname, ins_fid, construct_fid, fr_id in results:
    ...
    insertions_by_construct[construct_fid].append({
        'uniquename': ins_uname,
        'feature_id': ins_fid,
        'producedby_fr_id': fr_id,                  # NEW
    })
```

### 2. New `__init__` attribute

In `ConstructHandler.__init__` (alongside existing instance caches), add:

```python
self.ti_producedby_pubs = {}   # FTA-181: {ti_feature_id: [pub_id, ...]}
```

### 3. New method: collect pubs for the producedby FRs

Add next to `attach_generic_ti_insertions` (after `construct_handler.py:697`):

```python
def get_generic_ti_producedby_pubs(self, session):
    """Load pubs on each FBti producedby FBtp feature_relationship (FTA-181)."""
    self.log.info('Get producedby pubs for generic TI insertions.')
    ti_by_fr = {}
    for construct in self.fb_data_entities.values():
        if not construct.is_generic_ti:
            continue
        for insertion in getattr(construct, 'generic_ti_insertions', []):
            ti_by_fr[insertion['producedby_fr_id']] = insertion['feature_id']
    if not ti_by_fr:
        return
    results = session.query(
        FeatureRelationshipPub.feature_relationship_id,
        FeatureRelationshipPub.pub_id,
    ).filter(
        FeatureRelationshipPub.feature_relationship_id.in_(ti_by_fr.keys()),
    ).distinct()
    for fr_id, pub_id in results:
        ti_fid = ti_by_fr[fr_id]
        self.ti_producedby_pubs.setdefault(ti_fid, []).append(pub_id)
    self.log.info(
        f'Cached producedby pubs for {len(self.ti_producedby_pubs)} FBti.')
```

`FeatureRelationshipPub` must be added to the `harvdev_utils.reporting` import
block at the top of `construct_handler.py`.

### 4. Enrich `get_generic_ti_anon_construct_data`

At `construct_handler.py:911-924`, add one field to the anon_data entry so the
emitter does not need a second lookup:

```python
anon_data.append({
    ...
    'producedby_pub_ids': list(self.ti_producedby_pubs.get(ti_fid, [])),
})
```

### 5. New emitter method

Insert after `map_generic_ti_anon_marker_associations`
(`construct_handler.py:1262`), before `map_fb_data_to_alliance`:

```python
def map_generic_ti_anon_con_cas_associations(self, anon_data):
    """Emit ConstructCassetteAssociationDTOs linking <FBti>_con → <FBti>_cas (FTA-181)."""
    self.log.info('Map anon construct to anon cassette associations.')
    counter = 0
    for data in anon_data:
        ins_uname = data['insertion_uniquename']
        cons_curie = f'FB:{ins_uname}_con'
        cassette_curie = f'FB:{ins_uname}_cas'
        pub_curies = self.lookup_pub_curies(data['producedby_pub_ids'])
        fb_rel = fb_datatypes.FBExportEntity()
        rel_dto = agr_datatypes.ConstructCassetteAssociationDTO(
            cons_curie, 'has_component', cassette_curie, pub_curies)
        # Parent FBtp is effectively obsolete for generic-TI (FTA-179).
        if data['parent_is_obsolete']:
            rel_dto.obsolete = True
            rel_dto.internal = True
        fb_rel.linkmldto = rel_dto
        self.construct_cassette_associations.append(fb_rel)
        counter += 1
    self.log.info(f'Created {counter} anon con-cas ConstructCassetteAssociationDTOs.')
```

### 6. Wire into `query_chado_and_export`

In the generic-TI block (`construct_handler.py:1300-1310`):

```python
# FTA-136: Create anonymous constructs for generic TI insertions.
insertions_by_construct = self.get_generic_ti_insertions(session)
self.attach_generic_ti_insertions(insertions_by_construct)
self.get_generic_ti_producedby_pubs(session)           # NEW - before anon_data built
# FTA-182: load associated-allele tool data before building anon_data.
self.get_associated_allele_tool_data(session)
generic_ti_data = self.get_generic_ti_anon_construct_data()
if generic_ti_data:
    self.map_generic_ti_anon_constructs(generic_ti_data)
    self.map_generic_ti_anon_marker_associations(generic_ti_data)   # FTA-183
    self.map_generic_ti_anon_con_cas_associations(generic_ti_data)  # NEW - FTA-181
    self.export_generic_ti_anon_constructs()
```

The existing `flag_unexportable_entities` / `generate_export_dict` calls at
`construct_handler.py:1313-1314` already drain `construct_cassette_associations`
into `construct_cassette_association_ingest_set`, so the new DTOs ride out on
that pipeline — no extra wiring.

### 7. Obsolete/internal inheritance — decide

Gillian's TSV spec does not ask for obsolete/internal flags, but the established
pattern for every other generic-TI anon DTO (FTA-179 obsolete propagation —
see construct `map_generic_ti_anon_constructs` not doing it, but
`map_generic_ti_anon_marker_associations:1255-1257` does, and the FTA-180 Part B
AlleleConstructAssociation does too) sets `obsolete=true, internal=true` when
`parent_is_obsolete`. For a generic-TI FBtp, `parent_is_obsolete` is always true
because `is_generic_ti` feeds into it (`construct_handler.py:919`).

Decision: follow the pattern — set obsolete+internal when parent is
effectively obsolete. The `<FBti>_con` construct and `<FBti>_cas` cassette are
themselves `placeholder=True` so they're invisible on the web regardless; this
just keeps the association metadata consistent with its endpoints.

## Reused utilities

- `ConstructCassetteAssociationDTO` — `agr_datatypes.py:345-372`
- `self.lookup_pub_curies(pub_ids)` — handler base
- `FBExportEntity()` wrapper — `fb_datatypes.py`
- `FeatureRelationshipPub` model — `harvdev_utils.reporting`
- `self.construct_cassette_associations` drain into
  `construct_cassette_association_ingest_set` — already wired at
  `construct_handler.py:1313-1314`

## Verification

1. Run the cassette retrieval script (same invocation used for FTA-136/180/182):
   ```
   python src/AGR_data_retrieval_curation_cassette.py ...
   ```

2. In the output `construct_cassette_association_ingest_set` file, extract the
   subset with `relation_name == 'has_component'` where both ends start `FB:FBti`
   and end `_con` / `_cas`, and diff vs
   `docs/FTA-181/FTA-181_anon_con_anon_cas_associations.tsv`:
   - Row count: **5409**.
   - Same `(<FBti>_con, has_component, <FBti>_cas)` triples.
   - `pub_curies` match (some rows expected to have multi-pub via `|`-list if
     the producedby FR has multiple pubs — mirror of FTA-180 Part B which had 4
     such rows).

3. Cross-check the producedby pubs against the FTA-180 Part B file
   `docs/FTA-180/FTA-180_anon_con_associations.tsv`: for each FBti, the
   `pub_curies` column of the FTA-181 row must equal the `pub_curies` column of
   the FTA-180 row. (Both ultimately come from `FBti producedby FBtp`
   `feature_relationship_pub`.)

4. Every new DTO has `obsolete=true, internal=true` (parent generic-TI is
   effectively obsolete).

5. Schema validation:
   ```
   python src/validate_agr_schema.py <new-file.json>
   ```

6. flake8 on the modified file:
   ```
   flake8 src/construct_handler.py
   ```

7. Spot-check: `FBti0167947_con has_component FBti0167947_cas FBrf0225747`
   appears in the output (first row of Gillian's TSV).

## Out of scope / noted

- Timestamps ("do the same as FTA-141"): the existing
  `map_anon_cassette_basic` doesn't set any timestamps on the anon Cassette
  entities — `ConstructCassetteAssociationDTO` inherits from
  `EvidenceAssociationDTO` and likewise has no timestamp slots the ticket would
  apply to. Nothing to do unless Gillian clarifies otherwise.
- `primary_external_id` question in the description ("is that
  primary_external_id ??"): yes — the Cassette already uses
  `primary_external_id = <FBti>_cas` in `map_anon_cassette_basic`, so this
  question is resolved in the existing Part A code.
- TSV attachment `FTA-181_anon_con_anon_cas_associations.tsv` is copied into
  `docs/FTA-181/` so the verification diff is in-repo.
