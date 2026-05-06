# FTA-180 Part B — `AlleleConstructAssociationDTO` for anon `<FBti>_con` constructs

## Context

Part A of FTA-180 (create anonymous `<FBti>_con` ConstructDTOs per generic-TI
insertion) shipped in commits `ccea8af`, `574d3bc`. Part B — emitting an
`AlleleConstructAssociationDTO` linking each anon construct back to its FBti —
was never implemented.

Per the ticket:

```
allele_curie    = FB:<FBti>
rel_type        = 'contains'
construct_curie = FB:<FBti>_con
pub_curies      = pubs on the producedby feature_relationship
                  (FBti subject -> generic-TI FBtp object)
```

Expected shape (from `docs/FTA-180/FTA-180_anon_con_associations.tsv`, 5409 rows):

```
#allele_curie    rel_type   construct_curie       pub_curies
FBti0167947      contains   FBti0167947_con       FBrf0225747
...
FBti0168382      contains   FBti0168382_con       FBrf0223447|FBrf0223769
```

Four rows have multiple producedby pubs; most have one.

## Current state

- `AlleleConstructAssociationDTO` class exists (`agr_datatypes.py:345-362`).
- Ingest set `allele_construct_association_ingest_set` is already used (emitted
  by the allele retrieval script at `AGR_data_retrieval_curation_allele.py:137-141`
  from `allele_handlers.py:980-981`).
- The anon-construct pipeline lives in `construct_handler.py` and is driven by the
  cassette retrieval script (`AGR_data_retrieval_curation_cassette.py:225-235`).
- `get_generic_ti_insertions()` (`construct_handler.py:632-675`) already
  identifies the FBti↔FBtp producedby relationships but **doesn't collect the
  pubs** on that rel — only `uniquename`, `feature_id`, `object_id`.

## Approach

Own the new DTOs inside `ConstructHandler` (adjacent to the anon pipeline),
export them into `allele_construct_association_ingest_set` via the existing
pipeline, and have the cassette retrieval script write a new output JSON file
for the set. No cross-handler data handoff needed.

### Files modified

- `src/construct_handler.py` — new chado query for producedby pubs, new emitter
  method, new instance list + attribute, wire into export.
- `src/AGR_data_retrieval_curation_cassette.py` — write a new JSON file
  containing the anon `allele_construct_association_ingest_set` entries.

No change to `allele_handlers.py` (the FBal-construct associations it already
emits are independent).

### 1. Collect producedby pubs

Extend `get_generic_ti_insertions()` (`construct_handler.py:632`) to also return
`FeatureRelationship.feature_relationship_id` per row. The query shape is simple
— add the column and pass it through to the insertion dict:

```python
session.query(
    Feature.uniquename,
    Feature.feature_id,
    FeatureRelationship.object_id,
    FeatureRelationship.feature_relationship_id,    # NEW
)
...
insertions_by_construct[construct_fid].append({
    'uniquename': ins_uname,
    'feature_id': ins_fid,
    'producedby_fr_id': fr_id,                      # NEW
})
```

Add a new method `get_generic_ti_producedby_pubs(session)` that loads pubs for
the collected `fr_id`s (drive off the FRs collected by the previous step). Store
as `self.ti_producedby_pubs: {ti_feature_id: [pub_id, ...]}`.

```python
def get_generic_ti_producedby_pubs(self, session):
    fr_ids = []
    ti_by_fr = {}
    for construct in self.fb_data_entities.values():
        for ins in getattr(construct, 'generic_ti_insertions', []):
            fr_ids.append(ins['producedby_fr_id'])
            ti_by_fr[ins['producedby_fr_id']] = ins['feature_id']
    if not fr_ids:
        return
    results = session.query(
        FeatureRelationshipPub.feature_relationship_id,
        FeatureRelationshipPub.pub_id,
    ).filter(
        FeatureRelationshipPub.feature_relationship_id.in_(fr_ids),
    ).distinct()
    for fr_id, pub_id in results:
        ti_fid = ti_by_fr[fr_id]
        self.ti_producedby_pubs.setdefault(ti_fid, []).append(pub_id)
```

Add `self.ti_producedby_pubs = {}` to `__init__`.

### 2. Enrich `get_generic_ti_anon_construct_data`

Add one field per anon_data entry:

```python
'producedby_pub_ids': list(self.ti_producedby_pubs.get(ti_fid, [])),
```

### 3. Add emitter method

New instance list attribute: `self.anon_construct_allele_associations = []`.

New method in construct_handler.py (after the FTA-183 marker associations
method):

```python
def map_generic_ti_anon_construct_allele_associations(self, anon_data):
    """Emit AlleleConstructAssociationDTOs linking each <FBti>_con to its FBti (FTA-180 Part B)."""
    self.log.info('Map anon construct -> FBti AlleleConstructAssociationDTOs.')
    counter = 0
    for data in anon_data:
        ins_uname = data['insertion_uniquename']
        allele_curie = f'FB:{ins_uname}'
        construct_curie = f'FB:{ins_uname}_con'
        pub_curies = self.lookup_pub_curies(data['producedby_pub_ids'])
        fb_rel = fb_datatypes.FBExportEntity()
        rel_dto = agr_datatypes.AlleleConstructAssociationDTO(
            allele_curie, 'contains', construct_curie, pub_curies)
        # FTA-179: if parent is effectively obsolete, flag the association too.
        if data['parent_is_obsolete']:
            rel_dto.obsolete = True
            rel_dto.internal = True
        fb_rel.linkmldto = rel_dto
        self.anon_construct_allele_associations.append(fb_rel)
        counter += 1
    self.log.info(f'Created {counter} anon AlleleConstructAssociationDTOs.')
```

Note: `parent_is_obsolete` already factors in `is_generic_ti` (from FTA-179
commit `c611fa3`), so this inherits the "hidden parent → hidden association"
rule automatically. Since the parent generic-TI FBtp is always obsolete+internal
in our model, every new DTO will come out obsolete+internal — but the logic is
correct regardless of future changes to `parent_is_obsolete`.

### 4. Wire into `query_chado_and_export`

Two additions inside the generic-TI block (after `map_generic_ti_anon_constructs`,
adjacent to the FTA-183 marker-association call):

```python
self.get_generic_ti_producedby_pubs(session)
...
self.map_generic_ti_anon_construct_allele_associations(generic_ti_data)
...
self.flag_unexportable_entities(
    self.anon_construct_allele_associations,
    'allele_construct_association_ingest_set')
self.generate_export_dict(
    self.anon_construct_allele_associations,
    'allele_construct_association_ingest_set')
```

Call `get_generic_ti_producedby_pubs(session)` right after
`attach_generic_ti_insertions` (so the pubs are cached before
`get_generic_ti_anon_construct_data` reads them).

### 5. Write output file from cassette retrieval script

`AGR_data_retrieval_curation_cassette.py` already builds an
`association_export_dict` for cassette-* ingest sets (lines 260-288). Extend
that block to also carry the anon allele-construct set and write a dedicated
file:

```python
# FTA-180 Part B: anon construct -> FBti associations written by this script
ingest_name = 'allele_construct_association_ingest_set'
if cons_handler.export_data.get(ingest_name):
    aca_export_dict = {
        'linkml_version': linkml_release,
        'alliance_member_release_version': database_release,
        ingest_name: list(cons_handler.export_data[ingest_name]),
    }
    aca_output_filename = output_filename.replace('cassette', 'anon_allele_construct_association')
    generate_export_file(aca_export_dict, log, aca_output_filename)
```

Guard on `get(...)` because this only runs when `ADD_CASS_TO_CONSTRUCT=YES`.

## Files modified

| File | Change |
|---|---|
| `src/construct_handler.py` | `__init__` adds `ti_producedby_pubs` + `anon_construct_allele_associations`; extend `get_generic_ti_insertions` with `feature_relationship_id`; new `get_generic_ti_producedby_pubs`; enrich `get_generic_ti_anon_construct_data` with `producedby_pub_ids`; new emitter method; export wiring in `query_chado_and_export`. |
| `src/AGR_data_retrieval_curation_cassette.py` | After existing anon-cassette export, write `anon_allele_construct_association_*.json` with the new ingest set. |

No new DTOs. No changes to `allele_handlers.py`.

## Reused utilities

- `AlleleConstructAssociationDTO` — `agr_datatypes.py:345`
- `self.lookup_pub_curies(pub_ids)` — handler base
- `FBExportEntity()` wrapper
- `flag_unexportable_entities`, `generate_export_dict` — handler base
- `FeatureRelationshipPub` model — already imported in construct_handler

## Verification

1. Run the cassette retrieval script with `ADD_CASS_TO_CONSTRUCT=YES`:
   ```
   ADD_CASS_TO_CONSTRUCT=YES python src/AGR_data_retrieval_curation_cassette.py ...
   ```
2. Confirm a new file `anon_allele_construct_association_*.json` is produced
   with the `allele_construct_association_ingest_set` populated.
3. Diff the extracted rows against `docs/FTA-180/FTA-180_anon_con_associations.tsv`:
   - Same count (5409 rows for fb_2026_02_07).
   - Same `allele_curie`, `relation_name='contains'`, `construct_curie` for
     every row.
   - Pub unions match (4 rows with `|`-separated FBrfs).
4. Confirm every new DTO has `obsolete=true, internal=true` (parent is
   effectively obsolete via FTA-179 logic).
5. Schema validation:
   ```
   python src/validate_agr_schema.py <new-file.json>
   ```
6. Spot-check that the existing allele-side `allele_construct_association_ingest_set`
   emitted by the allele script is unchanged (no overlap in keys — anon
   constructs never appear in the allele side).
