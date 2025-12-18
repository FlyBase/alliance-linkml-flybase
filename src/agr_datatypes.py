"""Module:: Alliance datatypes.

Synopsis:
    Objects representing Alliance LinkML objects, tailored to FlyBase.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""


# Base class for all Alliance DTO (data transfer object) classes.
class AuditedObjectDTO(object):
    """Base Alliance DTO class."""
    def __init__(self):
        """Create AuditedObjectDTO for FlyBase objects."""
        self.internal = False
        self.obsolete = False
        self.date_created = None
        self.date_updated = None
        self.created_by_curie = 'FB:FB_curator'
        self.updated_by_curie = 'FB:FB_curator'
        self.required_fields = ['internal']
        self.internal_fields = ['internal_fields', 'required_fields']

    def dict_export(self):
        """Return a JSON-friendly dict for cases where inlined object needed."""
        export_dict = {}
        for k, v in self.__dict__.items():
            if k not in self.internal_fields and v is not None and v != []:
                export_dict[k] = v
        return export_dict


# Primary Alliance DTO Classes for FlyBase primary entities (having curies, web pages).
class SubmittedObjectDTO(AuditedObjectDTO):
    """SubmittedObjectDTO class."""
    def __init__(self):
        """Create SubmittedObjectDTO for FlyBase objects."""
        super().__init__()
        self.primary_external_id = None
        self.mod_internal_id = None
        self.data_provider_dto = None
        self.required_fields.extend(['data_provider_dto'])


class BiologicalEntityDTO(SubmittedObjectDTO):
    """BiologicalEntityDTO class."""
    def __init__(self):
        """Create BiologicalEntityDTO for FlyBase object."""
        super().__init__()
        self.taxon_curie = None
        self.required_fields.extend(['taxon_curie'])


class GenomicEntityDTO(BiologicalEntityDTO):
    """GenomicEntityDTO class."""
    def __init__(self):
        """Create GenomicEntityDTO for FlyBase object."""
        super().__init__()
        self.cross_reference_dtos = []
        self.required_fields.extend([])


class AffectedGenomicModelDTO(GenomicEntityDTO):
    """AffectedGenomicModelDTO class."""
    def __init__(self):
        """Create AffectedGenomicModelDTO for FlyBase object."""
        super().__init__()
        # Legacy "name" field retired in LinkML v2.12.0 - replaced with structured name DTOs
        self.agm_full_name_dto = None      # One NameSlotAnnotationDTO for the current full name
        self.agm_synonym_dtos = []         # List of NameSlotAnnotationDTO objects for synonyms and aliases
        self.subtype_name = None           # "strain" or "genotype".
        self.agm_secondary_id_dtos = []    # Secondary IDs.
        self.reference_curies = []         # Publication curies (PMID or FBrf).
        # self.component_dtos = []         # Retired in LinkML v2.9.0.
        self.required_fields.extend(['subtype_name'])


class AlleleDTO(GenomicEntityDTO):
    """AlleleDTO class."""
    def __init__(self):
        """Create AlleleDTO for FlyBase object."""
        super().__init__()
        self.allele_symbol_dto = None                          # One NameSlotAnnotationDTO.
        self.allele_full_name_dto = None                       # One NameSlotAnnotationDTO.
        self.allele_synonym_dtos = []                          # Many NameSlotAnnotationDTO objects.
        self.in_collection_name = None                         # Will be the name of a FlyBase library/collection.
        self.is_extinct = None                                 # Make True if extinction reported; make False is stock exists; leave as None otherwise.
        self.allele_mutation_type_dtos = []                    # AlleleMutationTypeSlotAnnotationDTOs.
        self.allele_inheritance_mode_dtos = []                 # AlleleInheritanceModeSlotAnnotationDTOs.
        self.allele_database_status_dto = None                 # AlleleDatabaseStatusSlotAnnotationDTOs.
        self.allele_secondary_id_dtos = []                     # SecondaryIdSlotAnnotationDTO (for 2o FB IDs).
        self.reference_curies = []                             # Will be a list of reference curies (directly or indirectly related).
        # Additional fields that may be used in the future.
        self.allele_functional_impact_dtos = []                # ToDo - Waiting on "Functional Impact" CV. Get feature_cvterm, child of "allele class" term.
        self.transgene_chromosome_location_curie = None        # ToDo - get chr via FBtp from FBti floc, derived_chromosome_location featureprop, or dock site.
        self.note_dtos = []                                    # ToDo - Waiting on "Allele Note Type" CV. Get from featureprop.
        # Additional unused fields.
        self.allele_germline_transmission_status_dto = None    # N/A (MGI).
        self.allele_nomenclature_event_dtos = []               # N/A.
        self.is_extrachromosomal = None                        # N/A (WB).
        self.is_integrated = None                              # N/A (WB).
        self.laboratory_of_origin_curie = None                 # N/A (WB).
        self.required_fields.extend(['allele_symbol_dto', 'primary_external_id'])


class CassetteDTO(GenomicEntityDTO):
    """CassetteDTO class."""
    def __init__(self):
        """Create CassetteDTO for FlyBase object."""
        super().__init__()
        self.cassette_symbol_dto = None      # One NameSlotAnnotationDTO.
        self.cassette_full_name_dto = None   # One NameSlotAnnotationDTO.
        self.cassette_synonym_dtos = []      # Many NameSlotAnnotationDTO objects.
        self.gene_type_curie = None                 # SO term ID for gene's promoted_gene_type.
        self.secondary_identifiers = []             # Annotation IDs and 2o FlyBase IDs.
        self.note_dtos = []                         # Will be NoteDTO objects.
        # self.cross_reference_dtos = []
        self.required_fields.extend(['cassette_symbol_dto'])
        self.required_fields.remove('taxon_curie')  # Does not have it!


class CassetteTransgenicToolAssociationDTO(AuditedObjectDTO):
    """CassetteTransgenicToolAssociationDTO class."""
    def __init__(self, cassette_association_subject, cassette_association_object,
                 pub_curies, obsolete, relation):
        """Create CassetteAssociationDTO for FlyBase object."""
        super().__init__()
        self.cassette_identifier = cassette_association_subject
        self.transgenic_tool_identifier = cassette_association_object
        # self.evidence = pub_curies
        self.obsolete = obsolete
        self.relation_name = relation


class CassetteGenomicEntityAssociationDTO(AuditedObjectDTO):
    """CassetteGenomicEntityAssociationDTO class."""
    def __init__(self, cassette_association_subject, cassette_association_object,
                 pub_curies, obsolete, relation):
        """Create CassetteAssociationDTO for FlyBase object."""
        super().__init__()
        self.cassette_identifier = cassette_association_subject
        self.genomic_entity_identifier = cassette_association_object
        # self.evidence = pub_curies
        self.obsolete = obsolete
        self.relation_name = relation


class CassetteStrAssociationDTO(AuditedObjectDTO):
    """CassetteStrAssociationDTO class."""
    def __init__(self, cassette_association_subject, cassette_association_object,
                 pub_curies, obsolete, relation):
        """Create CassetteStrAssociationDTO for FlyBase object."""
        super().__init__()
        self.cassette_identifier = cassette_association_subject
        self.sequence_targeting_reagent_identifier = cassette_association_object
        # self.evidence = pub_curies
        self.obsolete = obsolete
        self.relation_name = relation


class GeneDTO(GenomicEntityDTO):
    """GeneDTO class."""
    def __init__(self):
        """Create GeneDTO for FlyBase object."""
        super().__init__()
        self.gene_symbol_dto = None             # One NameSlotAnnotationDTO.
        self.gene_full_name_dto = None          # One NameSlotAnnotationDTO.
        self.gene_systematic_name_dto = None    # One NameSlotAnnotationDTO.
        self.gene_synonym_dtos = []             # Many NameSlotAnnotationDTO objects.
        self.gene_type_curie = None             # SO term ID for gene's promoted_gene_type.
        self.gene_secondary_id_dtos = []        # Annotation IDs and 2o FlyBase IDs.
        # self.reference_curies = []              # Not yet part of LinkML, so not exported - should be added to LinkML model?
        self.note_dtos = []                     # Will be NoteDTO objects.
        self.gcrp_cross_reference_dto = None    # Will be a single CrossReferenceDTO object for UniProt/GCRP xref, if any.
        self.required_fields.extend(['gene_symbol_dto'])


class TransgenicToolDTO(GenomicEntityDTO):
    """TransgenicToolDTO class."""
    def __init__(self):
        """Create TransgenicToolDTO for FlyBase object."""
        super().__init__()
        self.transgenic_tool_symbol_dto = None      # One NameSlotAnnotationDTO.
        self.transgenic_tool_full_name_dto = None   # One NameSlotAnnotationDTO.
        self.transgenic_tool_synonym_dtos = []      # Many NameSlotAnnotationDTO objects.
        self.gene_type_curie = None                 # SO term ID for gene's promoted_gene_type.
        self.secondary_identifiers = []             # Annotation IDs and 2o FlyBase IDs.
        self.note_dtos = []                         # Will be NoteDTO objects.
        self.cross_reference_dtos = []
        self.required_fields.extend(['transgenic_tool_symbol_dto'])
        self.required_fields.remove('taxon_curie')  # Does not have it!


class TransgenicToolAssociationDTO(AuditedObjectDTO):
    """TransgenicToolAssociationDTO class."""
    def __init__(self, transgenic_tool_subject_identifier, transgenic_tool_object_identifier,
                 pub_curies, obsolete, relation):
        """Create TransgenicToolAssociationDTO for FlyBase object."""
        super().__init__()
        self.transgenic_tool_subject_identifier = transgenic_tool_subject_identifier
        self.transgenic_tool_object_identifier = transgenic_tool_object_identifier
        # self.evidence_curies = pub_curies  # ??
        # self.obsolete = obsolete
        self.relation_name = relation
        self.required_fields.extend(['transgenic_tool_subject_identifier',
                                     'transgenic_tool_object_identifier'])


class ReagentDTO(SubmittedObjectDTO):
    """ReagentDTO class."""
    def __init__(self):
        """Create ReagentDTO for FlyBase object."""
        super().__init__()
        self.secondary_identifiers = []    # Will be list of 2o FB curies (strings).
        self.required_fields.extend([])


class ConstructDTO(ReagentDTO):
    """ConstructDTO class."""
    def __init__(self):
        """Create ConstructDTO for FlyBase object."""
        super().__init__()
        self.construct_symbol_dto = None       # One NameSlotAnnotationDTO.
        self.construct_full_name_dto = None    # One NameSlotAnnotationDTO.
        self.construct_synonym_dtos = None     # NameSlotAnnotationDTOs.
        self.construct_component_dtos = []     # ConstructComponentSlotAnnotationDTOs.
        self.reference_curies = []
        self.required_fields.extend(['construct_symbol_dto'])


# Primary Alliance DTO Classes for FlyBase relationships/annotations.
class EvidenceAssociationDTO(AuditedObjectDTO):
    """EvidenceAssociationDTO class."""
    def __init__(self, evidence_curies):
        """Create EvidenceAssociationDTO for FlyBase object.

        Args:
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__()
        self.evidence_curies = evidence_curies
        self.required_fields.extend([])


class SingleReferenceAssociationDTO(AuditedObjectDTO):
    """SingleReferenceAssociationDTO class."""
    def __init__(self):
        """Create SingleReferenceAssociationDTO for FlyBase object."""
        super().__init__()
        self.evidence_curie = None
        self.required_fields.extend([])


class AgmAlleleAssociationDTO(AuditedObjectDTO):
    """AgmAlleleAssociationDTO class."""
    def __init__(self, genotype_curie, component_curie, zygosity):
        """Create AgmAlleleAssociationDTO for FlyBase object."""
        super().__init__()
        self.agm_subject_identifier = genotype_curie
        self.allele_identifier = component_curie
        self.zygosity_curie = self.zygosity_id[zygosity]
        self.relation_name = 'contains'
        self.required_fields.extend(['agm_subject_identifier', 'allele_identifier', 'zygosity_curie', 'relation_name'])
    # Zygosity mapping to GENO IDs.
    # https://github.com/monarch-initiative/GENO-ontology/blob/develop/geno-base.obo
    zygosity_id = {
        # 'heterozygous': 'GENO:0000135',             # Retired in favor of more specific terms.
        # 'hemizygous': 'GENO:0000134_hemizygous',    # Not yet implemented in FB code, may never be.
        'simple heterozygous': 'GENO:0000458',
        'compound heterozygous': 'GENO:0000402',
        'homozygous': 'GENO:0000136',
        'unspecified zygosity': 'GENO:0000137',
        # 'homoplasmic': 'GENO:0000602',
        # 'heteroplasmic': 'GENO:0000603',
        # 'hemizygous X-linked': 'GENO:0000604',
        # 'hemizygous Y-linked': 'GENO:0000605',
        # 'hemizygous insertion-linked': 'GENO:0000606'
    }


class AlleleGenomicEntityAssociationDTO(EvidenceAssociationDTO):
    """AlleleGenomicEntityAssociationDTO class."""
    def __init__(self, evidence_curies):
        """Create AlleleGenomicEntityAssociationDTO for FlyBase object.

        Args:
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.allele_identifier = None
        self.relation_name = None
        self.evidence_code_curie = None
        self.note_dto = None
        self.required_fields.extend(['allele_identifier', 'relation_name'])


class AlleleConstructAssociationDTO(AlleleGenomicEntityAssociationDTO):
    """AlleleConstructAssociationDTO class."""
    def __init__(self, allele_id: str, rel_type: str, construct_id: str, evidence_curies: list):
        """Create AlleleConstructAssociationDTO for FlyBase object.

        Args:
            allele_id (str): The FB:FBal curie for the allele subject.
            rel_type (str): A CV term: TBD.
            construct_id (str): The FB:FBtp curie for the construct object.
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.allele_identifier = allele_id
        self.relation_name = rel_type
        self.construct_identifier = construct_id
        self.evidence_curies = evidence_curies
        self.required_fields.extend(['construct_identifier'])


class AlleleGeneAssociationDTO(AlleleGenomicEntityAssociationDTO):
    """AlleleGeneAssociationDTO class."""
    def __init__(self, allele_id: str, rel_type: str, gene_id: str, evidence_curies: list):
        """Create AlleleGeneAssociationDTO for FlyBase object.

        Args:
            allele_id (str): The FB:FBal curie for the allele subject.
            rel_type (str): A CV term: TBD.
            gene_id (str): The FB:FBgn curie for the gene object.
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.allele_identifier = allele_id
        self.relation_name = rel_type
        self.gene_identifier = gene_id
        self.evidence_curies = evidence_curies
        self.required_fields.extend(['gene_identifier'])


class ConstructGenomicEntityAssociationDTO(EvidenceAssociationDTO):
    """ConstructGenomicEntityAssociationDTO class."""
    def __init__(self, construct_id: str, rel_type: str, genomic_id: str, evidence_curies: list):
        """Create ConstructGenomicEntityAssociationDTO for FlyBase object.

        Args:
            construct_id (str): The FB:FBtp curie for the construct subject.
            rel_type (str): A CV term from "Construct Genomic Entity Association Relation": expressed, targets, or, is_regulated_by.
            genomic_id (str): The FB:FB curie for the construct object, limited to LinkML exported entities.
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.construct_identifier = construct_id
        self.genomic_entity_relation_name = rel_type
        self.genomic_entity_identifier = genomic_id
        # self.evidence_curies = evidence_curies
        self.note_dtos = []
        self.required_fields.extend(['construct_identifier', 'genomic_entity_relation_name', 'genomic_entity_identifier'])


class AnnotationDTO(SingleReferenceAssociationDTO):
    """AnnotationDTO class."""
    def __init__(self, evidence_curie):
        """Create AnnotationDTO for FlyBase object."""
        super().__init__()
        self.evidence_curie = evidence_curie
        self.primary_external_id = None
        self.mod_internal_id = None
        self.data_provider_dto = None
        self.note_dtos = []
        self.condition_relation_dtos = []
        self.required_fields.extend(['data_provider_dto', 'evidence_curie'])


class DiseaseAnnotationDTO(AnnotationDTO):
    """DiseaseAnnotationDTO class."""
    def __init__(self, do_term_curie, evidence_curie):
        """Create DiseaseAnnotationDTO for FlyBase object."""
        super().__init__(evidence_curie)
        self.do_term_curie = do_term_curie
        self.annotation_type_name = 'manually_curated'
        self.negated = False
        self.evidence_code_curies = []
        self.disease_genetic_modifier_identifiers = None
        self.disease_genetic_modifier_relation_name = None
        self.required_fields.extend(['disease_relation_name', 'do_term_curie', 'evidence_code_curies'])


class AlleleDiseaseAnnotationDTO(DiseaseAnnotationDTO):
    """AlleleDiseaseAnnotationDTO class."""
    def __init__(self, allele_identifier, do_term_curie, evidence_curie):
        """Create AlleleDiseaseAnnotationDTO for FlyBase object."""
        super().__init__(do_term_curie, evidence_curie)
        self.disease_relation_name = 'is_implicated_in'
        self.allele_identifier = allele_identifier
        self.inferred_gene_identifier = None
        self.required_fields.extend(['allele_identifier'])


class AGMDiseaseAnnotationDTO(DiseaseAnnotationDTO):
    """AGMDiseaseAnnotationDTO class."""
    def __init__(self, agm_identifier, do_term_curie, evidence_curie):
        """Create AGMDiseaseAnnotationDTO for FlyBase object."""
        super().__init__(do_term_curie, evidence_curie)
        self.disease_relation_name = None
        self.agm_identifier = agm_identifier
        self.inferred_gene_identifier = None
        self.asserted_gene_identifiers = []
        self.inferred_allele_identifier = None
        self.asserted_allele_identifiers = []
        self.required_fields.extend(['agm_identifier'])


# Secondary Alliance DTO Classes for data associated with first class objects.
# Note - when attaching these annotations to first class entities above, make
#        sure to first convert them to dict objects using the built-in
#        dict_export() method. Otherwise, upon attempting to print out the JSON
#        file, there will be a TypeError: Object ... is not JSON serializable.
class AffectedGenomicModelComponentDTO(AuditedObjectDTO):
    """AffectedGenomicModelComponentDTO class."""
    def __init__(self, component_curie, zygosity):
        """Create AffectedGenomicModelComponentDTO for FlyBase object."""
        super().__init__()
        self.allele_identifier = component_curie
        self.zygosity_curie = self.zygosity_id[zygosity]
        self.required_fields.extend(['allele_identifier', 'zygosity_curie'])
    # Zygosity mapping to GENO IDs.
    # https://github.com/monarch-initiative/GENO-ontology/blob/develop/geno-base.obo
    zygosity_id = {
        # 'hemizygous': 'GENO:0000134_hemizygous',                          # Not yet implemented in FB code.
        # 'heterozygous': 'GENO:0000135',                                   # Retire use of generic parent term.
        'simple heterozygous': 'GENO:0000135',                              # Map specific term to ID of generic parent term for now.
        'compound heterozygous': 'GENO:0000135',                            # Map specific term to ID of generic parent term for now.
        # 'simple heterozygous': 'GENO:0000458_simple_heterozygous',        # Implement this only when this GENOTerm is allowed at the Alliance.
        # 'compound heterozygous': 'GENO:0000402_compound_heterozygous',    # Implement this only when this GENOTerm is allowed at the Alliance.
        'homozygous': 'GENO:0000136_homozygous',
        'unspecified zygosity': 'GENO:0000137_unspecified_zygosity',
        # 'homoplasmic': 'GENO:0000602',
        # 'heteroplasmic': 'GENO:0000603',
        # 'hemizygous X-linked': 'GENO:0000604',
        # 'hemizygous Y-linked': 'GENO:0000605',
        # 'hemizygous insertion-linked': 'GENO:0000606'
    }


class CrossReferenceDTO(AuditedObjectDTO):
    """CrossReferenceDTO class."""
    def __init__(self, prefix, referenced_curie, page_area, display_name):
        """Create a simple CrossReferenceDTO for FlyBase object."""
        super().__init__()
        self.prefix = prefix
        self.referenced_curie = referenced_curie
        self.page_area = page_area
        self.display_name = display_name
        self.required_fields.extend(['prefix', 'referenced_curie', 'page_area', 'display_name'])


class DataProviderDTO(AuditedObjectDTO):
    """DataProviderDTO class."""
    def __init__(self, data_provider_xref):
        """Create a simple DataProviderDTO for FlyBase objects.

        Args:
            data_provider_xref (dict): A dict generated by CrossReferenceDTO.dict_export().

        """
        super().__init__()
        self.source_organization_abbreviation = 'FB'
        self.cross_reference_dto = data_provider_xref
        self.required_fields.extend(['source_organization_abbreviation'])


class NoteDTO(AuditedObjectDTO):
    """NoteDTO class."""
    def __init__(self, note_type_name: str, free_text: str, evidence_curies: list):
        """Create a NoteDTO for a FlyBase statement.

        Args:
            note_type_name (str): The type of note.
            free_text (str): The note itself.
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__()
        self.note_type_name = note_type_name
        self.free_text = free_text
        self.evidence_curies = evidence_curies
        self.required_fields.extend(['note_type_name', 'free_text'])


class SlotAnnotationDTO(AuditedObjectDTO):
    """SlotAnnotationDTO class."""
    def __init__(self, evidence_curies: list):
        """Create a SlotAnnotationDTO for FlyBase object.

        Args:
            evidence_curies (list): A list of FB:FBrf or PMID curies.

        """
        super().__init__()
        self.evidence_curies = evidence_curies


class AlleleDatabaseStatusSlotAnnotationDTO(SlotAnnotationDTO):
    """AlleleDatabaseStatusSlotAnnotationDTO class."""
    def __init__(self, database_status_name: str, evidence_curies: list):
        """Create an AlleleDatabaseStatusSlotAnnotationDTO for FlyBase allele.

        Args:
            database_status_name (str): The name of the database status term.
            evidence_curies (list): A list of FB:FBrf or PMID curies, or, an empty list.

        """
        super().__init__(evidence_curies)
        self.database_status_name = database_status_name
        self.required_fields.extend(['database_status_name'])


class AlleleInheritanceModeSlotAnnotationDTO(SlotAnnotationDTO):
    """AlleleInheritanceModeSlotAnnotationDTO class."""
    def __init__(self, inheritance_mode_name: str, phenotype_term_curie: str, phenotype_statement: str, evidence_curies: list):
        """Create an AlleleInheritanceModeSlotAnnotationDTO for FlyBase allele.

        Args:
            inheritance_mode_name (str): The name of the inheritance mode term.
            phenotype_term_curie (str): The curie of the FB term representing the phenotype.
            phenotype_statement (str): The phenotype context for the inheritance mode.
            evidence_curies (list): A list of FB:FBrf or PMID curies, or, an empty list.

        """
        super().__init__(evidence_curies)
        self.inheritance_mode_name = inheritance_mode_name
        self.phenotype_term_curie = phenotype_term_curie    # Cannot be used until FB terms are in the persistent store.
        self.phenotype_statement = phenotype_statement
        self.required_fields.extend(['inheritance_mode_name'])
        self.internal_fields.extend(['phenotype_term_curie'])


class AlleleMutationTypeSlotAnnotationDTO(SlotAnnotationDTO):
    """AlleleMutationTypeSlotAnnotationDTO class."""
    def __init__(self, mutation_type_curie: str, evidence_curies: list):
        """Create an AlleleMutationTypeSlotAnnotationDTO for FlyBase allele.

        Args:
            mutation_type_curie (str): The curie of the SO term representing the mutation type.
            evidence_curies (list): A list of FB:FBrf or PMID curies, or, an empty list.

        """
        super().__init__(evidence_curies)
        self.mutation_type_curies = [mutation_type_curie]
        self.required_fields.extend(['mutation_type_curies'])


class ConstructComponentSlotAnnotationDTO(SlotAnnotationDTO):
    """ConstructComponentSlotAnnotationDTO class."""
    def __init__(self, rel_type: str, component_symbol: str, taxon_curie: str, taxon_text: str, evidence_curies: list):
        """Create a ConstructComponentSlotAnnotation for a FlyBase construct component.

        Args:
            rel_type (str): A CV term from "Construct Genomic Entity Association Relation": expressed, targets, or, is_regulated_by.
            component_symbol (str): The symbol for the component.
            taxon_curie (str): The NCBITaxon ID of the component.
            taxon_text (str): The species name of the component.
            evidence_curies (list): A list of FB:FBrf or PMID curies.

        """
        super().__init__(evidence_curies)
        self.relation_name = rel_type
        self.component_symbol = component_symbol
        self.taxon_curie = taxon_curie
        self.taxon_text = taxon_text
        self.note_dtos = []
        self.required_fields.extend([])


class NameSlotAnnotationDTO(SlotAnnotationDTO):
    """NameSlotAnnotationDTO class."""
    def __init__(self, name_type_name: str, format_text: str, display_text: str, evidence_curies: list):
        """Create a NameSlotAnnotationDTO for a FlyBase name.

        Args:
            name_type_name (str): The type of name.
            format_text (str): The ASCII-only version of the name.
            display_text (str): The UTF-8 version of the name (using HTML superscript/subscript tags).
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.name_type_name = name_type_name
        self.format_text = format_text
        self.display_text = display_text
        self.synonym_scope_name = 'exact'
        self.required_fields.extend(['name_type_name', 'format_text', 'display_text', 'synonym_scope'])


class SecondaryIdSlotAnnotationDTO(SlotAnnotationDTO):
    """SecondaryIdSlotAnnotationDTO class."""
    def __init__(self, secondary_id: str, evidence_curies: list):
        """Create a SecondaryIdSlotAnnotationDTO for FlyBase object.

        Args:
            evidence_curies (list): A list of FB:FBrf or PMID curies.
            secondary_id (str): The secondary ID string (should include the db prefix).

        """
        super().__init__(evidence_curies)
        self.secondary_id = secondary_id
        self.required_fields.extend(['secondary_id'])
