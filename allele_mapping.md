TO DO
correct reason for no mapping
mention genotype_cvterm marking compliance
- similar for alleles?



# Reporting Build SOP
- [Overview](#Overview)
- [AlleleConversion](#AlleleConversion)
- [GenotypeConversion](#GenotypeConversion)
- [AssociatedAnnotations](#AssociatedAnnotations)
- [FutureWork](#FutureWork)

## Overview
FlyBase has had a complicated way of representing transgenic alleles, with a single insertional event being represented by many FBal and FBti features.  
To align with the Alliance, FlyBase is switching to representing these events by a single FBti insertion entity at the Alliance (submitted as an allele).  
FlyBase FBal alleles that are superseded in this way are still submitted to the Alliance, but marked as internal and obsolete. They should not be used in curation at the Alliance.  
- Aspects of these FBal alleles are propagated to the superseding FBti insertion entity.  
CRITICAL - Given all of this, an FBal allele may be current in chado, but obsolete at the Alliance (deprecated in favor of an FBti).

This decision on how to represent FlyBase transgenic FBal alleles at the Alliance impacts genotypes.  
Genotypes composed of FBal alleles that are deprecated at the Alliance must be updated to a form that is composed of the favored FBti insertions (as applicable).  

This decision also impacts annotations that refer to alleles or genotypes, as these alleles and genotypes must be reported as the "Alliance-compliant" form upon export to the Alliance.  
This was the case for disease annotations, which have been exported to the Alliance, and the Alliance is now the source of truth for these annotations.  
Expression, phenotype and genetic interaction data will also require some conversion of the alleles and/or genotypes involved upon Alliance export.  

## AlleleConversion
The logic for these allele conversions is as follows.

Case 1. Non-converted alleles.  
FBal alleles are exported "as themselves", and not converted, in the following cases:  
The allele is not the result of an FBtp construct genomic insertion: e.g., point mutation.  
The allele represents two or more closeby FBti insertions.  
The allele represents a mix of non-insertion and insertion events.  
The FBal allele is directly "associated_with" more than one FBtp construct (the FBal-FBti mapping is not one-to-one).  

Case 2. Converted at-locus insertions.  
Here we consider cases where an FBtp construct is inserted into the genome, disrupting a given locus at the insertion site ("at-locus", insertion is inside the allele).  
In chado, this is modeled by an FBti insertion to represent the insertion, and an FBal allele for every gene disrupted by that insertion (can be many where genes overlap the insertion).  
The FBal allele is "associated_with" the FBti insertion, and the FBti feature is "producedby" the FBtp construct (all in the feature_relationship table).  
For these, there is a script "map_alleles_to_insertions_for_alliance_genotypes.py" from the "alliance-linkml-flybase" repo.  
This script creates an "is_represented_at_alliance_as" feature_relationship between the FBal allele (subject) and the superseding FBti insertion (object).  
This script is run every epicycle (data_validation_p2 > updateAlleleInsertionRelationships).  
Upon export to the Alliance, the AlleleHandler looks for this relationship and propagates select information from the FBal allele to the FBti insertion, submitting the FBal allele as internal/obsolete.  

Case 3. Constructs carrying engineered genes.  
Here we consider cases where an FBtp construct carries engineered genes. The location of the construct insertion is generally not relevant to the function of the transgenic allele.  
In chado, this is modeled by an FBal allele for each gene carried in the FBtp construct (FBal "associated_with" FBtp in the feature_relationship).  
In chado, there may or may not be an insertion site curated for the FBtp construct insertion site (because it's not so important in this case).  
For these, there is a script "make_unspecified_construct_insertions.py" from the "harvdev-production-bulk-update" repo.  
This script creates a new "unspecified" FBti feature representing a generic insertion for every FBtp. The name of the FBti is the name of the FBtp with the "unspecified" suffix.  
This new "unspecified" FBti insertion is in addition to any specific FBti insertions (at curated locations) that may already exist.  
This script is run every epicycle (data_validation_p2 > updateAlleleInsertionRelationships).  
Upon export to the Alliance, information from the FBal is propagated to all related FBti insertions (the "unspecified" one, plus any curated ones).  

Every epicycle, the "alliance_allele_conversion_production_chado" SVN file is generated, which reports to curators the FB-to-Alliance FBal-FBti conversions.  
This is file is generated every epicycle by the "report_alliance_allele_mappings.py" script from the "harvdev-epicycle" repo.  
This script is run in the "Curation_Data > runCurationDataPython" pipeline step.  

## GenotypeConversion
In addition to script that set up allele-to-insertion conversions, there is a similar script for genotypes.  
For these, there is a script "make_alliance_genotypes.py" from the "harvdev-production-bulk-update" repo.  
This script is run every epicycle (data_validation_p2 > makeAllianceGenotypes).  
If a genotype is fine as is (it contains no FBal alleles that require conversion to an FBti insertion), then the genotype is marked as `alliance_compliant` in the `genotype_cvterm` table.  
If any features in a genotype require FBal-to-FBti conversion:  
- If all features that require conversion can be converted, a new genotype is made using the converted features. The original genotype (subject) is associated with the new `alliance_compliant` genotype (object) in the `genotype_relationship` table by a `is_represented_at_alliance_as` relation type.  
- If there are any transgenic FBal directly `associated_with` many FBtp constructs, those FBal features cannot be converted, and an `alliance_compliant` genotype cannot be made. So genotypes that cannot be converted lack both the `alliance_compliant` CV term annotation and a `is_represented_at_alliance_as` genotype association.  
Every release, the "alliance_genotype_conversion_production_chado" SVN file is generated, which reports to curators the FB-to-Alliance FBal-FBti conversions.  
This is file is generated every epicycle by the "report_alliance_genotype_mappings.py" script from the "harvdev-epicycle" repo.  
This script is run in the "Curation_Data > runCurationDataPython" pipeline step.  

## AssociatedAnnotations
In the "alliance-linkml-flybase" repo, the AGMDiseaseHandler provides an example of how to convert genotypes cited in annotations.
- Use this as an example for performing similar conversions in the export of phenotype data.
Basically, the "derive_genotypes()" method takes the set of alleles mentioned in the feature_cvterm disease annotations (including related feature_cvtermprops), and derives the name of the genotype.
Then, "get_genotypes()" method takes the name of that genotype, and passes it into the GenotypeAnnotation() object (from "harvdev-utils").
This GenotypeAnnotation() object takes the genotype name, identifies the features, then performs all necessary allele conversions to determine the final genotype required. It then either gets that genotype from chado, or it creates that genotype in chado.
CRITICAL - the `GenotypeAnnotation` object has the logic required to convert each genotype component, and create an `alliance_compliant` genotype if necessary.  
CRITICAL - The AGMDiseaseHandler should be run on a production version of chado (not a release), such that any new genotypes are created in production (and not an offshoot reporting instance).

## FutureWork
Just looking at production_chado, it still takes many steps to figure out if an allele needs to be converted, because FBal-is_represented_at_alliance_as-FBti relationships only cover case 1.
- So, at the moment, if an allele lacks such an "is_represented_at_alliance_as" relationship to an FBti, it could represent case 1 (should not be converted) or a subset of case 3.
- In retrospect, it might be good to have some sort of relationship created when an FBal is superceded unambiguously by a single "unspecified" FBti insertion.
- In retrospect, it might be good to have some way to flag alleles that should not be converted (case 1), or could not be converted (case 3 subset).
- These changes would make it easier to determine if FBal-FBti conversion is required, not required, or not possible.
