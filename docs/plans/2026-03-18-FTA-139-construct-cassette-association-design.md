# FTA-139: Construct-Cassette Association Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Export `marked_with` relationships between constructs (FBtp) and cassettes (FBal) as `ConstructCassetteAssociationDTO` for Alliance LinkML submission.

**Architecture:** The `marked_with` relationships (FBtp subject → FBal object) are already captured by `get_entity_relationships(session, 'subject')` in `construct_handler.py:207`. We add a DTO class, a mapping method that uses `recall_relationships()` to filter for `marked_with`, and wire up the export pipeline.

**Tech Stack:** Python, SQLAlchemy, chado PostgreSQL database

---

### Task 1: Create branch FTA-139

**Step 1: Create and switch to new branch**

Run: `cd /Users/ilongden/harvard/alliance-linkml-flybase && git checkout -b FTA-139`

**Step 2: Commit** — N/A (no files changed yet)

---

### Task 2: Add ConstructCassetteAssociationDTO to agr_datatypes.py

**Files:**
- Modify: `src/agr_datatypes.py:354` (insert before `ConstructGenomicEntityAssociationDTO`)

**Step 1: Add the new DTO class**

Insert after line 372 (end of `ConstructGenomicEntityAssociationDTO`), before `class AnnotationDTO`:

```python
class ConstructCassetteAssociationDTO(EvidenceAssociationDTO):
    """ConstructCassetteAssociationDTO class."""
    def __init__(self, construct_id: str, rel_type: str, cassette_id: str, evidence_curies: list):
        """Create ConstructCassetteAssociationDTO for FlyBase object.

        Args:
            construct_id (str): The FB:FBtp curie for the construct.
            rel_type (str): A relation name: has_selectable_marker or has_transcriptional_unit.
            cassette_id (str): The FB:FBal curie for the cassette.
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.construct_identifier = construct_id
        self.relation_name = rel_type
        self.cassette_identifier = cassette_id
        self.required_fields.extend(['construct_identifier', 'relation_name', 'cassette_identifier'])
```

**Step 2: Run flake8 check**

Run: `cd /Users/ilongden/harvard/alliance-linkml-flybase && flake8 src/agr_datatypes.py`
Expected: No errors

**Step 3: Commit**

```bash
git add src/agr_datatypes.py
git commit -m "feat(FTA-139): add ConstructCassetteAssociationDTO"
```

---

### Task 3: Add cassette association list and mapping method to construct_handler.py

**Files:**
- Modify: `src/construct_handler.py`

**Step 1: Add class attribute for cassette associations**

After line 62 (`construct_associations = []`), add:

```python
    construct_cassette_associations = []    # Will be a list of FBExportEntity objects, map to ConstructCassetteAssociationDTO.
```

**Step 2: Replace the TODO comment with a pass-through comment**

At line 216, replace:
```python
        # BOB: get_allele_marked_with()
```
with:
```python
        # marked_with rels already captured by get_entity_relationships(session, 'subject') above.
```

**Step 3: Add `map_construct_cassette_associations()` method**

Insert after `map_construct_genomic_associations()` (after line 490), before `map_fb_data_to_alliance()`:

```python
    def map_construct_cassette_associations(self):
        """Map construct marked_with cassette relations to ConstructCassetteAssociationDTO."""
        self.log.info('Map construct marked_with cassette associations.')
        counter = 0
        for construct in self.fb_data_entities.values():
            marked_with_rels = construct.recall_relationships(
                self.log, entity_role='subject', rel_types='marked_with', rel_entity_types='allele')
            for rel in marked_with_rels:
                cassette_feature_id = rel.chado_obj.object_id
                cassette_uniquename = self.feature_lookup[cassette_feature_id]['uniquename']
                cons_curie = f'FB:{construct.uniquename}'
                cassette_curie = f'FB:{cassette_uniquename}'
                pub_curies = self.lookup_pub_curies(rel.pubs)
                # FBal0345196 is a special case per FTA-139.
                if cassette_uniquename == 'FBal0345196':
                    relation_name = 'has_transcriptional_unit'
                else:
                    relation_name = 'has_selectable_marker'
                fb_rel = fb_datatypes.FBExportEntity()
                rel_dto = agr_datatypes.ConstructCassetteAssociationDTO(
                    cons_curie, relation_name, cassette_curie, pub_curies)
                if construct.is_obsolete is True or self.feature_lookup[cassette_feature_id]['is_obsolete'] is True:
                    rel_dto.obsolete = True
                    rel_dto.internal = True
                fb_rel.linkmldto = rel_dto
                self.construct_cassette_associations.append(fb_rel)
                counter += 1
        self.log.info(f'Mapped {counter} ConstructCassetteAssociationDTOs.')
        return
```

**Step 4: Run flake8 check**

Run: `flake8 src/construct_handler.py`
Expected: No errors

**Step 5: Commit**

```bash
git add src/construct_handler.py
git commit -m "feat(FTA-139): add map_construct_cassette_associations method"
```

---

### Task 4: Wire up mapping and export in construct_handler.py

**Files:**
- Modify: `src/construct_handler.py`

**Step 1: Update `map_fb_data_to_alliance()`**

In `map_fb_data_to_alliance()`, after line 509 (`self.flag_internal_fb_entities('construct_associations')`), add:

```python
        self.map_construct_cassette_associations()
        self.flag_internal_fb_entities('construct_cassette_associations')
```

**Step 2: Update `query_chado_and_export()`**

In `query_chado_and_export()`, after line 517 (`self.generate_export_dict(self.construct_associations, ...)`), add:

```python
        self.flag_unexportable_entities(self.construct_cassette_associations, 'construct_cassette_association_ingest_set')
        self.generate_export_dict(self.construct_cassette_associations, 'construct_cassette_association_ingest_set')
```

**Step 3: Run flake8 check**

Run: `flake8 src/construct_handler.py`
Expected: No errors

**Step 4: Commit**

```bash
git add src/construct_handler.py
git commit -m "feat(FTA-139): wire cassette associations into export pipeline"
```

---

### Task 5: Update entry point script to output cassette associations

**Files:**
- Modify: `src/AGR_data_retrieval_curation_construct.py`

**Step 1: Add cassette association export**

In `main()`, after line 116 (`generate_export_file(association_export_dict, ...)`), and still inside the `if not reference_session:` block, add:

```python
        # Export the construct-cassette associations to a separate file.
        cassette_assoc_output_filename = output_filename.replace('construct', 'construct_cassette_association')
        cassette_assoc_export_dict = {
            'linkml_version': linkml_release,
            'alliance_member_release_version': database_release,
        }
        cassette_assoc_export_dict['construct_cassette_association_ingest_set'] = cons_handler.export_data['construct_cassette_association_ingest_set']
        if len(cassette_assoc_export_dict['construct_cassette_association_ingest_set']) == 0:
            log.warning('The "construct_cassette_association_ingest_set" is empty.')
        else:
            generate_export_file(cassette_assoc_export_dict, log, cassette_assoc_output_filename)
```

Note: Using `log.warning` instead of raising ValueError since it's plausible that in test mode there may be no marked_with relationships in the test set.

**Step 2: Run flake8 check**

Run: `flake8 src/AGR_data_retrieval_curation_construct.py`
Expected: No errors

**Step 3: Commit**

```bash
git add src/AGR_data_retrieval_curation_construct.py
git commit -m "feat(FTA-139): export construct-cassette associations to JSON"
```

---

### Task 6: Final verification

**Step 1: Run flake8 on all modified files**

Run: `flake8 src/agr_datatypes.py src/construct_handler.py src/AGR_data_retrieval_curation_construct.py`
Expected: No errors

**Step 2: Review diff**

Run: `git diff main --stat`
Expected: 3 files changed

---

## Key Decisions

- **No new query method needed**: `get_entity_relationships(session, 'subject')` already captures `marked_with` relationships. We use `recall_relationships()` to filter them in the mapping step.
- **Separate output file**: Cassette associations get their own JSON file (`construct_cassette_association`), following the existing pattern where genomic entity associations have a separate file.
- **Warning not error for empty set**: Test mode may not have `marked_with` data in the test_set, so we warn instead of raising ValueError.
- **FBal0345196 hardcoded exception**: Per the Jira ticket, this specific allele uses `has_transcriptional_unit` instead of `has_selectable_marker`.