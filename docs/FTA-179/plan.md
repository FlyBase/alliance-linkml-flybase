# FTA-179 — Submit generic-TI FBtps as obsolete constructs

## Context

Sub-task of FTA-136. Each "generic TI" FBtp identified by FTA-178 should still be
submitted to the Alliance as a `ConstructDTO`, but flagged `obsolete=True` AND
`internal=True`. This keeps them visible to curators while preventing them from
leaking onto the Alliance website — the real user-facing constructs for each
generic-TI insertion are the anonymous `<FBti>_con` constructs from FTA-180.

Ticket text explicitly says "no special code" for the payload — the normal
construct mapping (symbols, synonyms, data provider, pubs, timestamps) should
still run; only `obsolete` and `internal` need to be forced true.

## Current state

**The core requirement is already implemented** on this branch, via the existing
generic-TI plumbing:

1. `construct_handler.py:1020` in `map_construct_basic()`:
   ```python
   agr_construct.obsolete = construct.chado_obj.is_obsolete or construct.is_generic_ti
   ```
   → forces `obsolete=True` on the `ConstructDTO` for every generic-TI FBtp.

2. `handler.py:1009-1035` `flag_internal_fb_entities()` runs after mapping
   (invoked at `construct_handler.py:1272` and `1278`):
   ```python
   if fb_data_entity.linkmldto.obsolete is True:
       fb_data_entity.linkmldto.internal = True
       fb_data_entity.internal_reasons.append('Obsolete')
   ```
   → auto-sets `internal=True` for any obsolete construct.

Result: generic-TI FBtps are exported with their normal synonyms, data provider,
pubs, and timestamps — but with `obsolete=True, internal=True`. That matches the
ticket's "follow what normally happens for obsolete things" note.

## Small consistency gap to close

The downstream construct-scoped association DTOs check
`construct.is_obsolete` (the raw chado `feature.is_obsolete` — False for a
generic-TI FBtp whose underlying chado row is non-obsolete) rather than the
effective-obsolete status. Affected call sites:

| File:line | DTO | Current check |
|---|---|---|
| `construct_handler.py:1045` | `ConstructComponentSlotAnnotationDTO` | `self.feature_lookup[feature_id]['is_obsolete']` only |
| `construct_handler.py:1102` | `ConstructGenomicEntityAssociationDTO` | `construct.is_obsolete or feature_lookup[...]['is_obsolete']` |
| `construct_handler.py:1133` | `ConstructCassetteAssociationDTO` (marked_with) | same pattern |
| `construct_handler.py:1173` | `ConstructCassetteAssociationDTO` (associated_with) | same pattern |
| `construct_handler.py:1246` | `ConstructCassetteAssociationDTO` (FTA-183 anon marker) | uses `data['parent_is_obsolete']` set at line 913 from `construct.is_obsolete` |

Because these use the chado flag directly, a generic-TI FBtp's associations are
emitted with `obsolete=False, internal=False` even though the construct they
reference is now hidden. Alliance would then surface associations pointing at an
internal/obsolete construct, which is at minimum confusing.

The ticket doesn't mandate flipping these, but the safe reading of "submit them
as obsolete" is that anything downstream of the obsolete construct also gets
flagged. Cost is small (5 one-line changes) and keeps Alliance output internally
consistent.

## Approach

### 1. Add a computed helper on `FBConstruct` (or inline per callsite)

Simplest: inline the same expression already used at `construct_handler.py:1020`
— `construct.chado_obj.is_obsolete or construct.is_generic_ti` — at each
callsite. Five one-line edits:

```python
# construct_handler.py:1045
if (construct.chado_obj.is_obsolete or construct.is_generic_ti
        or self.feature_lookup[feature_id]['is_obsolete']):
    ...

# construct_handler.py:1102, 1133, 1173  (same shape, same fix)

# construct_handler.py:1246 (FTA-183 marker-on-anon path)
# data['parent_is_obsolete'] already covers chado obsolete; need to also include is_generic_ti.
```

For line 1246, the simplest fix is to update line 913 in
`get_generic_ti_anon_construct_data`:

```python
'parent_is_obsolete': construct.is_obsolete or construct.is_generic_ti,
```

That single change propagates the effective-obsolete status through `anon_data`
without touching the downstream callsite.

Note the anon `<FBti>_con` construct itself (emitted by
`map_generic_ti_anon_constructs`) must remain non-obsolete — FTA-180 emits it as
the *replacement* public face of the insertion. Don't touch
`map_generic_ti_anon_constructs`.

### 2. (Optional) Add `FBConstruct.is_effectively_obsolete` property

Alternative to inline: add a computed property on `FBConstruct` in
`fb_datatypes.py` that returns `self.is_obsolete or self.is_generic_ti`. Use it
at the 5 callsites. Marginally cleaner but adds a new concept for a one-branch
situation. Inline is fine.

## Files modified

- `src/construct_handler.py` — 5 one-line edits (or 4 + 1 via anon_data line 913).

No changes to `cassette_handler.py`, schema, or new DTOs.

## Verification

1. Run the construct retrieval script on a dataset containing `FBtp0143312`
   (the current generic-TI example):
   ```
   python src/AGR_data_retrieval_curation_construct.py
   ```
2. In the JSON output, find the `construct_ingest_set` entry for
   `FB:FBtp0143312` and confirm:
   ```json
   {"primary_external_id": "FB:FBtp0143312",
    "obsolete": true,
    "internal": true,
    "construct_symbol_dto": {...},   # still populated
    "data_provider_dto": {...}}      # still populated
   ```
3. Find the cassette-association DTO(s) for FBtp0143312 in
   `construct_cassette_association_ingest_set` (the marked_with → FBal0361621
   link from FTA-139) and confirm it too has `obsolete: true, internal: true`
   (after consistency fix).
4. Confirm the anonymous `<FBti>_con` constructs for the same parent (e.g.
   `FB:FBti0212460_con`) are emitted with `obsolete: false, internal: false`
   (they are the public face).
5. Confirm existing non-generic FBtp constructs are unchanged.
6. Schema validation:
   ```
   python src/validate_agr_schema.py <output.json>
   ```

## Notes

- No new chado queries.
- No risk to the FTA-180/181/182/183 anon pipelines — those remain the public
  surface. This only flips a flag on the already-obsolete parent.
- This plan does NOT suppress the parent FBtp from any ingest_set, because
  Alliance's obsolete+internal behavior already hides it from the website. The
  ticket explicitly wants curators to still see the data.
