# Fix Hardcoded `feature_cvterm` in `get_entity_cvterms()` Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace hardcoded `feature_cvterm` references in `entity_handler.py:get_entity_cvterms()` with generic `chado_type`-based attribute access so the method works for strains, genotypes, and all other entity types.

**Architecture:** The method already uses a `chado_type` variable (set to `'feature'`, `'strain'`, `'genotype'`, etc.) and `getattr()` for dynamic table/column access throughout most of the method. Lines 579-591 break this pattern by hardcoding `feature_cvterm`. The fix makes these lines follow the same `getattr()` pattern, using the fact that all chado `*Cvtermprop` models have a `{chado_type}_cvterm` relationship, and all `*Cvterm` models have a `{chado_type}` relationship back to the parent entity plus `cvterm` and `pub` relationships.

**Tech Stack:** Python, SQLAlchemy (harvdev_utils chado models)

---

## Root Cause

In `src/entity_handler.py`, the `get_entity_cvterms()` method processes cvterm property results (lines 574-598). While the query building (lines 491-573) correctly uses `chado_type` and `getattr()` to be generic, the result-processing loop hardcodes `feature_cvterm`:

```python
# Line 579 - hardcoded
entity_id = cvtermprop_result.feature_cvterm.feature.feature_id
# Lines 588-591 - hardcoded
prop_data = {'name': cvtermprop_result.feature_cvterm.cvterm.name, ...}
```

When `chado_type = 'strain'`, the `StrainCvtermprop` object has no `feature_cvterm` attribute. It has `strain_cvterm`, which points to `StrainCvterm` (which has `.strain`, `.cvterm`, `.pub`).

### Model relationship pattern (all types follow this):

| Cvtermprop model | Relationship to Cvterm | Cvterm → Entity | Cvterm → Cvterm | Cvterm → Pub |
|---|---|---|---|---|
| `FeatureCvtermprop` | `.feature_cvterm` | `.feature` | `.cvterm` | `.pub` |
| `StrainCvtermprop` | `.strain_cvterm` | `.strain` | `.cvterm` | `.pub` |
| `GenotypeCvtermprop` | `.genotype_cvterm` | `.genotype` | `.cvterm` | `.pub` |

---

### Task 1: Replace hardcoded `feature_cvterm` with generic `chado_type`-based access

**Files:**
- Modify: `src/entity_handler.py:574-598`

**Step 1: Create the reusable cvterm accessor variable**

At line 574 (just before the loop begins), add a line to get the cvterm relationship name:

```python
        cvterm_rel_name = f'{chado_type}_cvterm'
```

**Step 2: Replace line 579 (entity_id lookup)**

Replace:
```python
            # Strain does not have feature_cvterm !!
            entity_id = cvtermprop_result.feature_cvterm.feature.feature_id
```

With:
```python
            cvterm_obj = getattr(cvtermprop_result, cvterm_rel_name)
            entity_obj = getattr(cvterm_obj, chado_type)
            entity_id = getattr(entity_obj, f'{chado_type}_id')
```

**Step 3: Replace line 584 (error logging)**

Replace:
```python
                    self.log.error(f"Entity_id:{entity_id} not in list of data_entities feature {cvtermprop_result.feature_cvterm.feature}")
```

With:
```python
                    self.log.error(f"Entity_id:{entity_id} not in list of data_entities {chado_type} {entity_obj}")
```

**Step 4: Replace lines 588-591 (prop_data dict)**

Replace:
```python
                prop_data = {'name': cvtermprop_result.feature_cvterm.cvterm.name,
                             'type': cvtermprop_result.feature_cvterm.cvterm.cv.name,
                             'pub': cvtermprop_result.feature_cvterm.pub.uniquename,
                             'accession': cvtermprop_result.feature_cvterm.cvterm.dbxref.accession}
```

With:
```python
                prop_data = {'name': cvterm_obj.cvterm.name,
                             'type': cvterm_obj.cvterm.cv.name,
                             'pub': cvterm_obj.pub.uniquename,
                             'accession': cvterm_obj.cvterm.dbxref.accession}
```

Note: `cvterm_obj` is already set earlier in the loop iteration (Step 2), so we reuse it here. The `.cvterm` and `.pub` relationships exist on all `*Cvterm` models (FeatureCvterm, StrainCvterm, GenotypeCvterm, etc.).

**Step 5: Verify the final code block looks correct**

The full loop (lines 574-598) should read:

```python
        cvterm_prop_counter = 0
        cvterm_rel_name = f'{chado_type}_cvterm'
        for cvtermprop_result in cvtermprop_results:
            entity_cvterm_id = getattr(cvtermprop_result, f'{chado_type}_cvterm_id')
            entity_prop_type_name = self.cvterm_lookup[cvtermprop_result.type_id]['name']
            cvterm_obj = getattr(cvtermprop_result, cvterm_rel_name)
            entity_obj = getattr(cvterm_obj, chado_type)
            entity_id = getattr(entity_obj, f'{chado_type}_id')
            if entity_id in self.ignore_list:
                continue
            elif entity_id not in self.fb_data_entities:
                if not self.testing:
                    self.log.error(f"Entity_id:{entity_id} not in list of data_entities {chado_type} {entity_obj}")
                    self.log.error(f"Ignore_list is {self.ignore_list}")
                continue
            if entity_prop_type_name in self.fb_data_entities[entity_id].prop_data:
                prop_data = {'name': cvterm_obj.cvterm.name,
                             'type': cvterm_obj.cvterm.cv.name,
                             'pub': cvterm_obj.pub.uniquename,
                             'accession': cvterm_obj.cvterm.dbxref.accession}
                self.fb_data_entities[entity_id].prop_data[entity_prop_type_name].append(prop_data)
            if entity_prop_type_name in cvterm_annotation_dict[entity_cvterm_id].props_by_type.keys():
                cvterm_annotation_dict[entity_cvterm_id].props_by_type[entity_prop_type_name].append(fb_datatypes.FBProp(cvtermprop_result))
                cvterm_prop_counter += 1
            else:
                cvterm_annotation_dict[entity_cvterm_id].props_by_type[entity_prop_type_name] = [fb_datatypes.FBProp(cvtermprop_result)]
                cvterm_prop_counter += 1
```

**Step 6: Commit**

```bash
git add src/entity_handler.py
git commit -m "fix: replace hardcoded feature_cvterm with generic chado_type access in get_entity_cvterms

The cvtermprop result processing loop hardcoded feature_cvterm relationships,
causing AttributeError for StrainCvtermprop (and any non-feature entity type).
Now uses getattr() with chado_type, matching the pattern used elsewhere in the method."
```