"""Module:: datatypes.

Synopsis:
    Objects representing Alliance LinkML objects, tailored to FlyBase.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""


# Organized hierarchically, then alphabetically.
class AuditedObjectDTO(object):
    """Base Alliance class."""
    def __init__(self):
        """Create base AuditedObjectDTO for FlyBase objects."""
        self.internal = False
        self.obsolete = False
        self.date_created = None
        self.date_updated = None
        self.created_by_curie = 'FB:FB_curator'
        self.updated_by_curie = 'FB:FB_curator'
        self.required_fields = ['internal']

    def dict_export(self):
        """Return a JSON-friendly dict for cases where inlined object needed."""
        export_dict = {}
        for k, v in self.__dict__.items():
            if k != 'required_fields' and v is not None and v != []:
                export_dict[k] = v
        return export_dict


class SubmittedObjectDTO(AuditedObjectDTO):
    """SubmittedObjectDTO Alliance class."""
    def __init__(self):
        """Create base AuditedObjectDTO for FlyBase objects."""
        super().__init__()
        self.mod_entity_id = None
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


class ReagentDTO(AuditedObjectDTO):
    """ReagentDTO class."""
    def __init__(self):
        """Create ReagentDTO for FlyBase object."""
        super().__init__()
        self.mod_entity_id = None          # Will be the MOD curie.
        self.mod_internal_id = None        # Will be the MOD internal db id.
        self.secondary_identifiers = []    # Will be list of 2o FB IDs (strings).
        self.required_fields.extend(['mod_entity_id'])


class GenomicEntityDTO(BiologicalEntityDTO):
    """GenomicEntityDTO class."""
    def __init__(self):
        """Create GenomicEntityDTO for FlyBase object."""
        super().__init__()
        self.cross_reference_dtos = []


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


class AffectedGenomicModelDTO(GenomicEntityDTO):
    """AffectedGenomicModelDTO class."""
    def __init__(self):
        """Create AffectedGenomicModelDTO for FlyBase object."""
        super().__init__()
        self.name = None                   # Current fullname synonym (ASCII).
        self.subtype_name = None           # "strain" or "genotype".
        self.agm_secondary_id_dtos = []    # Secondary IDs.
        self.reference_curies = []         # Publication curies (PMID or FBrf).
        self.component_dtos = []           # AffectedGenomicModelComponentDTOs.
        self.required_fields.extend(['subtype_name'])


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
        self.reference_curies = []              # Not yet part of LinkML, so not exported - should be added to LinkML model?
        self.related_notes = []                 # Will be NoteDTO objects.
        self.required_fields.extend(['gene_symbol_dto'])


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


class AlleleGenomicEntityAssociationDTO(EvidenceAssociationDTO):
    """AlleleGenomicEntityAssociationDTO class."""
    def __init__(self, evidence_curies):
        """Create EvidenceAssociationDTO for FlyBase object.

        Args:
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__(evidence_curies)
        self.allele_identifier = None
        self.relation_name = None
        self.evidence_code_curie = None
        self.note_dto = None
        self.required_fields.extend(['allele_identifier', 'relation_name'])


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


# Secondary Alliance DTO Classes for FlyBase data.
# Use the export() method to get a dict for these (easier to print out).
class AffectedGenomicModelComponentDTO(AuditedObjectDTO):
    """AffectedGenomicModelComponentDTO class."""
    def __init__(self, component_curie, zygosity):
        """Create AffectedGenomicModelComponentDTO for FlyBase object."""
        super().__init__()
        self.allele_curie = component_curie
        self.zygosity_curie = self.zygosity_id[zygosity]
    # Zygosity mapping.
    zygosity_id = {
        'hemizygous': 'GENO:0000134',
        'heterozygous': 'GENO:0000135',
        'homozygous': 'GENO:0000136',
        'unspecified zygosity': 'GENO:0000137',
        'homoplasmic': 'GENO:0000602',
        'heteroplasmic': 'GENO:0000603',
        'hemizygous X-linked': 'GENO:0000604',
        'hemizygous Y-linked': 'GENO:0000605',
        'hemizygous insertion-linked': 'GENO:0000606'
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
