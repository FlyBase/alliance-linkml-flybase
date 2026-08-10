# INTRODUCTION

FlyBase expression curation has three main branches:

- 9.5K genes with 109K annotations.
- 18.1K transgenic alleles with 120K annotations.
- 1.8K split system combinations with 4K annotations.

The current Alliance model for expression requires a gene to be the subject of the annotation.
However, at FlyBase, only about half of annotations are for endogenous gene products.
We need to figure out how to export split system combination data (later discussion) and transgenic allele expression (this discussion).
For transgenic alleles, curators are able to indicate when a transgenic gene product reflects the expression of an endogenous gene.
However, this covers only ~ 1/4 of transgenic product annotations.
In about 5% of cases, the transgenic product expression is aberrant relative to the wildtype expression of the associated gene.

An overview of non-melanogaster transgenic product expression annotation by species.
We have to be careful how these are handled at the Alliance (prevent leaking into Alliance gene pages of non-Dros species).

| Species | n_alleles | n_xprn_annotations | Alliance species |
| --- | ---: | ---: | --- |
| Scer (yeast) | 12988 | 100956 | yes |
| Ecol (E. coli) | 3032 | 12826 | |
| Hsap (human) | 831 | 2712 | yes |
| Avic (jellyfish) | 726 | 1938 | |
| Mmus (mouse) | 48 | 132 | yes |
| Ncra (Neurospora) | 35 | 121 | |
| Disc (soft coral) | 32 | 91 | |
| Hsim (herpes simplex virus) | 12 | 58 | |
| Btau (cow) | 2 | 24 | |
| Rnor (rat) | 4 | 6 | yes |
| Tn10 (bact. TE) | 2 | 6 | |
| Uuuu (unknown) | 1 | 2 | |

We need to figure out how to map transgenic expression annotations to the Alliance model.
We should submit in a way that matches how other MODs handle similar data (ensure that all MODs are harmonized).

# STRATEGY FOR HANDLING TRANSGENIC EXPRESSION ANNOTATIONS

A gene product that has expression data can be related, directly or indirectly, to a Dmel gene.
For gene products that are "associated_with" an allele (`uniquename ~ '^FBal[0-9]{7}$'`), there can be many different pathways from allele to gene(s), even for a single allele: a single path usually leads to one gene, but sometimes many; a single allele may have many paths to a gene.
In some cases, the path to a gene dead ends at an FBti insertion site.

Review the categories below, then write a python script (following conventions in this repo) to query the chado database and figure out how many gene products (and related annotations in the `feature_expression` table) fall into each category.
Also, assess and report all pairwise category overlaps (# gene products, # annotations).
Create two output files, a log file and a report file.
The log file should document the script run, and provide bulk counts requested.
The report file should list all gene products (`feature_id`, `uniquename`, `name` in columns 1-3), and list categories (column 4 as a comma-separated list of category labels).
For the bulk counts, indicate:

1. how many gene products map to MANY genes;
2. how many gene products represent aberrant gene expression (see CATEGORY 7 below).

In the descriptions below, `feature_relationship` chains are indicated as follows:

1. features are described in square brackets with info on uniquename patterns and organism_id.
2. `feature_relationship` type names are indicated in round brackets flanked by "-" on the left and "->" on the right, with the ">" (subject and object indicated left and right, respectively, of the type).

## CATEGORY 0

- **Label:** UNKNOWN
- **Description:** A gene product fits none of the following categories.

## CATEGORY 1

- **Label:** SPLIT_COMBO
- **Description:** The annotation subject is the result of a split system combination. No gene tracing is required.
- **Scope:** 1.8K FBco products with 4.0K annotations.
- **Examples:**
    - `GAL4[DBD.R10G01]+RELA[AD.R49F02]` (FBco0001980) has 5 expression annotations.
    - `GAL4[DBD.R11A03]+RELA[AD.R58E02]` (FBco0000101) has 13 expression annotations.
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotations
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
WHERE gp.is_obsolete IS FALSE
  AND gp.uniquename ~ '^FBco[0-9]{7}$';
```

## CATEGORY 2A

- **Label:** ENDO
- **Description:** A gene product directly "associated_with" a Dmel gene.
- **Chain:** `[gene_product (XR/XP name suffix, FBtr/FBpp uniquename)] -(associated_with)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 13.9K gene products, 108.5K annotations.
- **Examples:**
    - `en-XP` is a generic polypeptide of `en` (FBgn0000577).
    - `wg-XR` is a generic transcript of `wg` (FBgn0284084).
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotations
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature g ON g.feature_id = gpg.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.name ~ '-X(R|P)$'
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$'
  AND g.organism_id = 1;
```

## CATEGORY 2B

- **Label:** ENDO_ISO
- **Description:** A gene product "isa" gene product that is directly "associated_with" a Dmel gene.
- **Chain:** `[gene_product (FBtr/FBpp uniquename)] -(isa)-> [gene_product (XR/XP name suffix, FBtr/FBpp uniquename)] -(associated_with)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 0.3K gene products, 0.7K annotations
- **Examples:**
    - `Abd-B[+]P493` is a polypeptide isoform of `Abd-B` (FBgn0000015).
    - `Antp[+]R5.1` is an RNA isoform of `Antp` (FBgn0260642).
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotations
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'isa')
JOIN feature gp2 ON gp2.feature_id = gpg.object_id
JOIN feature_relationship gp2g ON gp2g.subject_id = gp2.feature_id AND gp2g.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature g ON g.feature_id = gp2g.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp2.is_obsolete IS FALSE
  AND gp2.name ~ '-X(R|P)$'
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$'
  AND g.organism_id = 1;
```

## CATEGORY 3

- **Label:** TAGGED
- **Description:** A gene product is directly "associated_with" a Dmel allele, which is in turn directly "alleleof" a Dmel gene.
- **Chain:** `[gene_product (FBtr/FBpp uniquename)] -(associated_with)-> [allele (FBal uniquename, organism_id=1)] -(alleleof)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 0.4K gene products, 1.3K annotations
- **Examples:**
    - `Bsg[G00413]` (FBal0243207) is an allele of `Bsg` (FBgn0261822).
    - `esg[G66]` (FBal0039323) is an allele of `esg` (FBgn0287768).
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotations
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature a ON a.feature_id = gpa.object_id
JOIN feature_relationship ag ON ag.subject_id = a.feature_id AND ag.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'alleleof')
JOIN feature g ON g.feature_id = ag.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
  AND a.is_obsolete IS FALSE
  AND a.uniquename ~ '^FBal[0-9]{7}$'
  AND a.organism_id = 1
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$';
```

## CATEGORY 4A

- **Label:** PROMOTER_CHAR_GENE
- **Description:** A transgenic non-Dmel allele has a Dmel gene's regulatory region.
- **Chain:** `[gene_product (FBtr/FBpp uniquename)] -(associated_with)-> [allele (FBal uniquename, organism_id!=1)] -(has_reg_region)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 1.1K gene products, 7.9K annotations
- **Examples:**
    - `Avic\GFP[EGFP.lbe.K]` (FBal0296331) has a regulatory region from `lbe` (FBgn0011278).
    - `Avic\GFP[mb247.Tag:TS(brpS)]` (FBal0397200) has a regulatory region from `Mef2` (FBgn0011656).
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotation
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature a ON a.feature_id = gpa.object_id
JOIN feature_relationship ag ON ag.subject_id = a.feature_id AND ag.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'has_reg_region')
JOIN feature g ON g.feature_id = ag.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
  AND a.is_obsolete IS FALSE
  AND a.uniquename ~ '^FBal[0-9]{7}$'
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$'
  AND a.organism_id != 1
  AND g.organism_id = 1;
```

## CATEGORY 4B

- **Label:** PROMOTER_CHAR_FBSF
- **Description:** A transgenic non-Dmel allele has a Dmel regulatory region that is associated with a gene.
- **Chain:** `[gene_product (FBtr/FBpp uniquename)] -(associated_with)-> [allele (FBal uniquename, organism_id!=1)] -(has_reg_region)-> [seqfeat (FBsf uniquename, organism_id=1] -(associated_with)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 10.8K gene products, 82.1K annotations
- **Examples:**
    - `Avic\GFP[Hand.1054-1400]` (FBal0263675) has `Hand_HV-element` (FBsf0000435371), related to `Hand` (FBgn0032209).
    - `Avic\GFP[Obp57d.Y]` (FBal0248806) has `Obp57d_Dmel-obp57d` (FBsf0000435801), related to `Obp57d` (FBgn0043536).
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotation
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature a ON a.feature_id = gpa.object_id
JOIN feature_relationship asf ON asf.subject_id = a.feature_id AND asf.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'has_reg_region')
JOIN feature sf ON sf.feature_id = asf.object_id
JOIN feature_relationship sfg ON sfg.subject_id = sf.feature_id AND sfg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature g ON g.feature_id = sfg.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
  AND a.is_obsolete IS FALSE
  AND a.uniquename ~ '^FBal[0-9]{7}$'
  AND sf.is_obsolete IS FALSE
  AND sf.uniquename ~ '^FBsf[0-9]{10}$'
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$'
  AND a.organism_id != 1
  AND g.organism_id = 1;
```

## CATEGORY 5A

- **Label:** INS_TRAP
- **Description:** An insertion traps nearby regulatory elements. Expression annotated under the reporter (e.g., GFP gene), but there is an indirectly associated Dmel allele.
- **Chain:** `[gene_product (FBtr/FBpp uniquename)] -(associated_with)-> [allele (FBal uniquename, organism_id!=1)] -(associated_with)-> [insertion (FBti uniquename, organism_id=1] <-(associated_with)- [allele, FBal uniquename, organism_id=1) -(alleleof)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 2.9K gene products, 15.7K annotations
- **Examples:**
    - `Avic\GFP::Ppyr\LUC[lola-226]` (FBal0387102), related to the same FBti insertion as `lola[226]` (FBal0179030).
    - `Scer\GAL4[zfh2-MS209]` (FBal0040477), related to the same FBti insertion as `zfh2[MS209]` (FBal0143099).
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotation
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature a ON a.feature_id = gpa.object_id
JOIN feature_relationship ai ON ai.subject_id = a.feature_id AND ai.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature i ON i.feature_id = ai.object_id
JOIN feature_relationship ia ON ia.object_id = i.feature_id AND ia.subject_id != ai.subject_id AND ia.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature a2 ON a2.feature_id = ia.subject_id AND a2.feature_id != a.feature_id
JOIN feature_relationship a2g ON a2g.subject_id = a2.feature_id AND a2g.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'alleleof')
JOIN feature g ON g.feature_id = a2g.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
  AND a.is_obsolete IS FALSE
  AND a.uniquename ~ '^FBal[0-9]{7}$'
  AND i.is_obsolete IS FALSE
  AND i.uniquename ~ '^FBti[0-9]{7}$'
  AND a.organism_id != 1
  AND i.organism_id = 1
  AND a2.is_obsolete IS FALSE
  AND a2.uniquename ~ '^FBal[0-9]{7}$'
  AND a2.organism_id = 1
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$'
  AND g.organism_id = 1;
```

## CATEGORY 5B

- **Label:** INS_TRAP_UNK
- **Description:** An insertion traps nearby regulatory elements. Expression annotated under the reporter (e.g., GFP gene), but there is no known Dmel allele/gene.
- **Chain:** `[gene_product (FBtr/FBpp uniquename)] -(associated_with)-> [allele (FBal uniquename, organism_id!=1)] -(associated_with)-> [insertion (FBti uniquename, organism_id=1]`
- **Scope:** 2.3K gene products, 11.3K annotations
- **Examples:**
    - `Avic\GFP[Bacc-YC0022]` (FBal0230215); FBti is not localized to the genome
    - `Avic\GFP[Nrv2-PTT]` (FBal0241865); FBti is not localized to the genome; name suggests a relationship to `nrv2` (FBgn0015777) but not formally traceable in database.
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotation
--SELECT DISTINCT gp.name, a.name, i.uniquename
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature a ON a.feature_id = gpa.object_id
JOIN feature_relationship ai ON ai.subject_id = a.feature_id AND ai.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
JOIN feature i ON i.feature_id = ai.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
  AND a.is_obsolete IS FALSE
  AND a.uniquename ~ '^FBal[0-9]{7}$'
  AND i.is_obsolete IS FALSE
  AND i.uniquename ~ '^FBti[0-9]{7}$'
  AND a.organism_id != 1
  AND i.organism_id = 1
  AND i.feature_id NOT IN
  (
      SELECT DISTINCT i.feature_id
      FROM feature a
      JOIN feature_relationship fr ON fr.subject_id = a.feature_id AND fr.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
      JOIN feature i ON i.feature_id = fr.object_id
      WHERE a.is_obsolete IS FALSE
        AND a.uniquename ~ '^FBal[0-9]{7}$'
        AND a.organism_id = 1
        AND i.is_obsolete IS FALSE
        AND i.uniquename ~ '^FBti[0-9]{7}$'
        AND i.organism_id = 1
   );
```

## CATEGORY 6

- **Label:** CURATED
- **Description:** Curators specify, in some cases, that an allelic gene product reflects expression of a gene. Only applies to `](R|P)A$` products of alleles. Excludes any that are marked as "aberrant" by curator.
- **Chain:** `[gene_product (RA/PA name suffix, FBtr/FBpp uniquename)] -(attributed_as_expression_of)-> [gene (FBgn uniquename, organism_id=1)]`
- **Scope:** 4.1K gene products, 23.6K annotations
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotations
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'attributed_as_expression_of')
JOIN feature g ON g.feature_id = gpg.object_id
WHERE gp.is_obsolete IS FALSE
  AND gp.name ~ '](R|P)A$'
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$'
  AND gpg.feature_relationship_id NOT IN
  (
      SELECT DISTINCT feature_relationship_id FROM feature_relationshipprop WHERE type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'is_relative_wildtype')
  );
```

## CATEGORY 7

- **Label:** ABERRANT
- **Description:** When curators relate a gene product to a gene using the "attributed_as_expression_of" relationship, they can add an "is_relative_wildtype" `feature_relationshipprop` (value='y' always) that indicates that the transgenic product is aberrant in some way when compared to wild type expression.
- **Chain:** `[gene_product (RA/PA name suffix, FBtr/FBpp uniquename)] -(attributed_as_expression_of)-> [gene (FBgn uniquename, organism_id=1)]` WITH "is_relative_wildtype" `feature_relationshipprop.type`
- **Scope:** 0.1K gene products, 1.2K annotations
- **Query:**

```sql
SELECT COUNT(DISTINCT gp.feature_id) AS n_gene_product, COUNT(DISTINCT fe.feature_expression_id) AS n_xprn_annotations
FROM feature gp
JOIN feature_expression fe ON fe.feature_id = gp.feature_id
JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'attributed_as_expression_of')
JOIN feature g ON g.feature_id = gpg.object_id
JOIN feature_relationshipprop frp ON frp.feature_relationship_id = gpg.feature_relationship_id AND frp.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'is_relative_wildtype')
WHERE gp.is_obsolete IS FALSE
  AND gp.name ~ '](R|P)A$'
  AND g.is_obsolete IS FALSE
  AND g.uniquename ~ '^FBgn[0-9]{7}$';
```
