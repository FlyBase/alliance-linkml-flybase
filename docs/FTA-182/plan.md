# FTA-182 — Add tool-related data to anonymous `<FBti>_cas` cassettes

## Context

Sub-task of FTA-136. FTA-181 already creates anonymous `<FBti>_cas` cassettes for
each generic-TI insertion, driven by `construct_handler.generic_ti_data_for_cassette_handler`
→ `cassette_handler.receive_anon_cassette_data` → `map_anon_cassettes`. Today the
anon cassette inherits `direct_rels` (`encodes_tool`, `has_reg_region`,
`tagged_with`, `carries_tool`) and `tool_uses_data` **from the parent FBtp**, using
the parent's pub_ids as pub_curies — which is wrong per the ticket.

FTA-182 replaces the single-source logic with a gated two-source model:

1. Tool data from the **parent generic-TI FBtp** is emitted with **blank** `pub_curies`,
   AND only if `propagate_transgenic_uses` allows.
2. Tool data from each **`associated_with` FBal** is emitted with the FBal's own
   `pub_curies`, gated by an `is_represented_at_alliance_as` cross-check between that
   FBal and the FBti.
3. If the cross-check fails for an FBal that *has* tool data, the cassette gets a
   NoteDTO (`internal_note`) explaining why the tool data was withheld.

Authoritative algorithm: `docs/FTA-182/text_for_FTA-182.v2.txt`.
Expected output shape: `docs/FTA-182/FTA-182_output.v2.txt` (9856 rows).

## Algorithm (from ticket, restated)

```
for each <FBti>_cas:
    associated_fbals = FeatureRelationship(type='associated_with',
                                           object_id=FBti, subject_id~FBal)

    if len(associated_fbals) == 0:
        # Use parent FBtp data, pub_curies BLANK.
        emit_parent_tool_data()

    else:
        propagate_transgenic_uses = True
        for fbal in associated_fbals:
            if fbal has Featureprop(type='propagate_transgenic_uses'):
                propagate_transgenic_uses = False
                break

        if propagate_transgenic_uses:
            emit_parent_tool_data()  # pub_curies BLANK

        for fbal in associated_fbals:
            if fbal has tool data:  # direct_rels OR tool_uses cvtermprops
                if FeatureRelationship(type='is_represented_at_alliance_as',
                                       subject=fbal, object=FBti) exists:
                    emit_fbal_tool_data(pub_curies=<FBal pubs>)
                else:
                    cassette.note_dtos.append(NoteDTO(
                        'internal_note',
                        f"FTA: unable to add tool information from {fbal.uniquename} "
                        f"(failed 'is_represented_at_alliance_as' cross-check)",
                        []))
                    log_warn_and_counter()
```

"Tool data" = the 4 feature_relationships (`encodes_tool`, `has_reg_region`,
`tagged_with`, `carries_tool`) + `tool_uses` FeatureCvtermprops.

## Approach

Data loading stays in `ConstructHandler` (consistent with how the generic-TI
pipeline already works). `CassetteHandler` remains the DTO builder but learns to
distinguish the generic-TI path from the FTA-181 non-generic path via an explicit
flag on each `anon_cassette_data` entry.

### Files modified

- `src/construct_handler.py` — new chado queries + enriched anon_data.
- `src/cassette_handler.py` — refactor `map_anon_cassette_*` to honor the new
  gating, emit FBal-sourced rels, blank parent pub_curies, and emit NoteDTOs.

No schema changes; no new DTO classes.

### 1. New chado data loads in `ConstructHandler`

Add a single orchestrator method `get_associated_allele_tool_data(session, ti_feature_ids)`
that issues the following queries once per run, keyed for cheap downstream lookups.
Called from `query_chado_and_export` after `get_generic_ti_insertions` and before
`get_generic_ti_anon_construct_data`.

| Data | Query (chado) | Output shape |
|---|---|---|
| Associated FBals per FBti | `FeatureRelationship` where `type='associated_with'`, `object_id ∈ ti_feature_ids`, `subject_id` matches allele regex `^FBal[0-9]{7}$` | `self.ti_associated_alleles: {ti_feature_id: [fbal_feature_id, ...]}` |
| Block flag per FBal | `Featureprop` where `type.name='propagate_transgenic_uses'`, `feature_id ∈ all_fbal_ids` | `self.fbal_block_propagation: set[fbal_feature_id]` (presence = block) |
| `is_represented_at_alliance_as` pairs | `FeatureRelationship` where `type='is_represented_at_alliance_as'`, `subject_id ∈ fbal_ids`, `object_id ∈ ti_feature_ids` | `self.fbal_ti_represented: set[(fbal_fid, ti_fid)]` |
| FBal tool-rels (4 types) | `FeatureRelationship` where `type.name ∈ {encodes_tool, has_reg_region, tagged_with, carries_tool}`, `subject_id ∈ fbal_ids`; join `FeatureRelationshipPub` | `self.fbal_tool_rels: {fbal_fid: [{rel_type, object_feature_id, pub_ids}, ...]}` |
| FBal tool_uses | Mirror existing `get_construct_tool_uses` pattern: `FeatureCvterm` + `FeatureCvtermprop` where `cvtermprop.type.name='tool_uses'`, `feature_id ∈ fbal_ids` | `self.fbal_tool_uses: {fbal_fid: [{cvterm_name, accession, pub_id}, ...]}` |

Rationale for batch loading: there can be thousands of generic-TI FBtis; a single
bulk query per datum (filtered by `IN (...)` on the cached fbal_ids) is much cheaper
than per-insertion queries, and matches the pattern used in `get_construct_tool_uses`
(`construct_handler.py:278-310`).

Reuse existing helpers:
- `self.regex['allele']` — already defined in `handler.py` (uniquename regex pattern)
- `self.feature_lookup[feature_id]` — uniquename + is_obsolete for any FB feature
- `self.lookup_pub_curies(pub_ids)` — pub_id → FBrf curie conversion

### 2. Extend `get_generic_ti_anon_construct_data` (`construct_handler.py:682`)

Per insertion, attach the new data so downstream code doesn't need the handler's
own caches. Add these keys to each anon_data dict:

```python
'is_generic_ti': True,             # FTA-182: differentiates from FTA-181 anon path
'associated_alleles': [            # from self.ti_associated_alleles[ti_fid]
    {
        'feature_id': int,
        'uniquename': str,
        'blocks_propagation': bool, # fbal_fid in self.fbal_block_propagation
        'is_represented_at_alliance_as': bool,  # (fbal_fid, ti_fid) in self.fbal_ti_represented
        'direct_rels': [{rel_type, object_feature_id, pub_ids}, ...],
        'tool_uses_data': [{cvterm_name, accession, pub_id}, ...],
    },
    ...
],
'propagate_transgenic_uses': bool, # computed: True unless any associated FBal blocks
```

`propagate_transgenic_uses` is computed here rather than in cassette_handler so
the logic stays adjacent to the loaded chado data.

Simultaneously, tag existing FTA-181 anon_data entries (produced by
`get_anon_cassette_data` at `construct_handler.py:599`) with `'is_generic_ti': False`
so the cassette handler can branch cleanly.

### 3. Refactor `map_anon_cassette_*` methods (`cassette_handler.py:642-746`)

The three mapping methods (`map_anon_cassette_simple_components`,
`map_anon_cassette_encodes_tool`, `map_anon_cassette_tool_uses`) currently iterate
`data['direct_rels']` / `data['tool_uses_data']` per entry. Change each to:

```python
for data in self.anon_cassette_data:
    if data.get('is_generic_ti'):
        self._emit_generic_ti_tool_data(data)   # NEW: FTA-182 logic
    else:
        self._emit_parent_tool_data_with_pubs(data)  # existing FTA-181 behavior, extracted
```

The extraction into two helpers preserves FTA-181 semantics verbatim and keeps the
new gating isolated. Shared DTO-building code (pick FBto / FBgn / FBsf branch,
instantiate the association DTO, append to the right list) factors into a small
shared helper `_append_tool_component_dto(cassette_id, rel_type, object_feature_id,
pub_curies, cassette_dto)`.

### 4. New helper: `_emit_generic_ti_tool_data(data)` in `CassetteHandler`

Implements the algorithm:

```python
def _emit_generic_ti_tool_data(self, data):
    cassette_id = f'FB:{data["construct_uniquename"]}_cas'
    cassette_dto = data['anon_cassette_dto']
    propagate = data['propagate_transgenic_uses']

    # 1. Parent FBtp data — only if propagate, and pub_curies BLANK.
    if propagate:
        for rel in data['direct_rels']:
            self._append_tool_component_dto(cassette_id, rel, [], cassette_dto)
        for tu in data['tool_uses_data']:
            self._append_tool_use_dto(cassette_id, tu, [], cassette_dto)

    # 2. Per associated FBal — emit iff represented, else NoteDTO.
    for fbal in data['associated_alleles']:
        has_tool_data = bool(fbal['direct_rels']) or bool(fbal['tool_uses_data'])
        if not has_tool_data:
            continue
        if fbal['is_represented_at_alliance_as']:
            for rel in fbal['direct_rels']:
                pub_curies = self.lookup_pub_curies(rel['pub_ids'])
                self._append_tool_component_dto(cassette_id, rel, pub_curies, cassette_dto)
            # Group tool_uses by accession so each DTO gets its own pubs.
            by_acc = {}
            for tu in fbal['tool_uses_data']:
                by_acc.setdefault(tu['accession'], []).append(tu['pub_id'])
            for accession, pub_ids in by_acc.items():
                pub_curies = self.lookup_pub_curies(pub_ids)
                cassette_dto.cassette_use_dtos.append(
                    agr_datatypes.CassetteUseSlotAnnotationDTO(
                        pub_curies, [f'FBcv:{accession}']).dict_export())
        else:
            self.log.warning(
                f"Skipping tool data from {fbal['uniquename']} for "
                f"{data['construct_uniquename']}_cas: failed "
                f"'is_represented_at_alliance_as' cross-check")
            note_dto = agr_datatypes.NoteDTO(
                'internal_note',
                f"FTA: unable to add tool information from {fbal['uniquename']} "
                f"(failed 'is_represented_at_alliance_as' cross-check)",
                []).dict_export()
            cassette_dto.note_dtos.append(note_dto)
            self._ir_failure_counter += 1
```

`_ir_failure_counter` is summed + logged at the end of `map_anon_cassettes`.

### 5. Helper: `_append_tool_component_dto`

Consolidates the FBto/FBgn/FBsf branches from
`map_anon_cassette_simple_components` (`cassette_handler.py:660-681`) and
`map_anon_cassette_encodes_tool` (`684-718`). Handles the rel_type → alliance_rel
mapping (`has_reg_region` → `is_regulated_by`, `tagged_with` → `tagged_with`,
`carries_tool` → `contains`, `encodes_tool` → `expresses`). A single function used
by both FTA-181 and FTA-182 paths eliminates the current near-duplication.

### 6. Deduplication

The ticket notes that the sample output has duplicate rows when the same tool is
carried by multiple FBals (`FBti0168384_cas carries_tool FBto0000349` appears twice,
once per source FBal). LinkML expects a single cassette-component entry. Use a
`set` keyed by `(cassette_id, rel_type, object_feature_id)` inside
`_emit_generic_ti_tool_data` to suppress repeats from multiple FBals, merging
their `pub_curies`.

For `tool_uses`, the grouping-by-accession pattern already does this naturally.

## Touched files summary

| File | Changes |
|---|---|
| `src/construct_handler.py` | New method `get_associated_allele_tool_data`; extend `get_generic_ti_anon_construct_data` to emit `is_generic_ti`, `associated_alleles`, `propagate_transgenic_uses`; tag FTA-181 path entries with `is_generic_ti: False` in `get_anon_cassette_data`; wire call into `query_chado_and_export`. |
| `src/cassette_handler.py` | Split `map_anon_cassette_simple_components` / `_encodes_tool` / `_tool_uses` into FTA-181 + FTA-182 branches; add `_emit_generic_ti_tool_data`, `_emit_parent_tool_data_with_pubs`, `_append_tool_component_dto`, `_append_tool_use_dto`; add per-run IR-failure counter. |

## Reused utilities / patterns

- `FeatureRelationship` / `FeatureRelationshipPub` query pattern —
  `construct_handler.py:627-670` (generic TI insertions)
- `FeatureCvterm` + `FeatureCvtermprop` `tool_uses` pattern —
  `construct_handler.py:278-310`
- `Featureprop` query pattern — `construct_handler.py:314-336` (marked_with pub lookup)
- `self.feature_lookup`, `self.organism_lookup`, `self.lookup_pub_curies` — handler base
- Existing DTOs in `agr_datatypes.py`: `CassetteTransgenicToolAssociationDTO`,
  `CassetteGenomicEntityAssociationDTO`, `CassetteComponentSlotAnnotationDTO`,
  `CassetteUseSlotAnnotationDTO`, `NoteDTO`

## Verification

1. Run the cassette retrieval script with `ADD_CASS_TO_CONSTRUCT=YES`:
   ```bash
   ADD_CASS_TO_CONSTRUCT=YES python src/AGR_data_retrieval_curation_cassette.py
   ```
2. Grep the output JSON for `_cas` entries and compare against
   `docs/FTA-182/FTA-182_output.v2.txt`. For each insertion in the sample:
   - Confirm rel_type + tool ID + pub_curies match (after dedup).
   - Confirm `pub_curies` is `[]` for every parent-FBtp-sourced row.
   - Confirm `pub_curies` is populated for FBal-sourced rows.
3. Targeted test cases (all currently in `construct_handler.test_set`):
   - `FBti0212460` (parent FBtp0143312) — 0 associated FBals → parent-only path,
     blank pub_curies.
   - A generic-TI FBti with associated FBals where all pass `is_represented_at_alliance_as`
     — no NoteDTOs, FBal pubs populated.
   - A generic-TI FBti with a `propagate_transgenic_uses` featureprop on one of
     its associated FBals — no parent tool data emitted.
   - A generic-TI FBti with an associated FBal that has tool data but FAILS the
     cross-check — one `NoteDTO` on the cassette with the exact text
     `FTA: unable to add tool information from <FBal> (failed 'is_represented_at_alliance_as' cross-check)`.
4. Check logs for the new IR-failure counter summary line.
5. FTA-181 regression: for constructs in test_set that have `needs_anon_cassette`
   (not `is_generic_ti`), confirm the emitted cassette DTOs are identical before
   and after this change (use `git stash` + diff the output JSON).
6. Schema validation:
   ```bash
   python src/validate_agr_schema.py <output.json>
   ```

## Open questions to resolve before coding

1. **Dedup policy for `pub_curies` when merging multiple FBals.** Ticket shows the
   raw (pre-dedup) rows; when two FBals both supply `carries_tool FBto0000349`
   with different FBrfs, should the merged DTO carry the union of FBrfs, or just
   the first? Proposal: union (set-deduped). Cheap and matches what a curator
   would expect.
2. **Whether FTA-181 path also needs pub_curies blanked for parent data.** The
   ticket text restricts the "no parent pubs" rule to generic-TI cassettes, so
   FTA-181 stays as-is. Confirm this with Gillian if unsure.
3. **Allele regex for `associated_with` subject filter.** Use the strict
   `^FBal[0-9]{7}$` from `self.regex['allele']`? Or accept anything subject_id
   returns (chado integrity should already constrain)? Proposal: use the regex for
   defensive correctness, matching how `get_generic_ti_insertions` filters FBti.
