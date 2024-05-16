"""Module:: datatypes.

Synopsis:
    Ojects representing FlyBase chado data types, for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""


# FlyBase Classes
# Attributes are sorted into "primary chado data" (i.e., query results) and
# processed data (synthesis of sql results).
class FBEntity(object):
    """A generic FlyBase entity."""
    def __init__(self):
        """Create a FBEntity object with bins for Alliance mapping."""
        self.db_primary_id = None     # Table primary key (or concatenation).
        self.uniq_key = None          # Uniquely identifying string.
        self.org_abbr = None          # Organism.abbreviation, if applicable.
        self.org_genus = None         # Organism.genus, if applicable.
        self.org_species = None       # Organism.species, if applicable.
        self.linkmldto = None         # Alliance LinkML object for mapped data.
        self.for_export = True        # Made False to prevent Alliance export.
        self.internal_reasons = []    # Reasons an object was marked internal.
        self.export_warnings = []     # Reasons an object was not exported.
        self.entity_desc = None       # A succinct description for this entity.

    def __str__(self):
        """Basic descriptive info for the object."""
        return self.entity_desc


# A placeholder for development of association/annotation export.
class FBAssociation(FBEntity):
    """An abstract, generic FlyBase association/annotation."""
    def __init__(self):
        """Create a FBAssociation object."""
        super().__init__()


# Objects below exclude associations/annotations.
class FBDataEntity(FBEntity):
    """An abstract, generic FlyBase data entity with all it related data."""
    def __init__(self, chado_obj):
        """Create a FBDataEntity object from the main db entry.

        Args:
            chado_obj (SQLAlchemy object): The object representing the table entry.

        Returns:
            A FBDataEntity object.
        """
        super().__init__()
        self.chado_obj = chado_obj
        table_columns = [
            'uniquename',
            'name',
            'type_id',
            'organism_id',
            'is_obsolete',
            'is_analysis',
            'timeaccessioned',
            'timelastmodified']
        for column_name in table_columns:
            try:
                setattr(self, column_name, getattr(chado_obj, column_name))
            except AttributeError:
                pass
        try:
            self.org_abbr = chado_obj.organism.abbreviation
            self.org_genus = chado_obj.organism.genus
            self.org_species = chado_obj.organism.species
        except AttributeError:
            pass
        self.entity_desc = f'{self.name} ({self.uniquename})'
        # Primary FB chado data - direct db query results, no processing.
        self.pubs = []                  # Pub associations: e.g., FeaturePub.
        self.synonyms = []              # Synonym associations: e.g., FeatureSynonym.
        self.fb_sec_dbxrefs = []        # Dbxref non-current associations for "FlyBase" db: e.g., FeatureDbxref.
        self.dbxrefs = []               # Dbxref associations: e.g., FeatureDbxref.
        self.props = []                 # entity props: e.g., Featureprop.
        self.prop_pubs = []             # entity prop_pubs: e.g., Featureprop_pub.
        self.cvterms = []               # Cvterm associations: e.g., FeatureCvterm.
        self.cvtermprops = []           # Cvtermprop associations: e.g., FeatureCvtermprop.
        self.timestamps = []            # FB timestamps.
        # Processed FB data - processed from primary FB chado data above.
        self.ncbi_taxon_id = None       # The NCBITaxon dbxref.accession (str).
        self.synonym_dict = {}          # Will be synonym_id-keyed dicts of processed synonym info, similar to NameDTO.
        self.curr_fb_symbol = None      # The current symbol for the entity (UTF-8).
        self.alt_fb_ids = []            # Secondary FB IDs (including the "FB:" prefix).
        self.all_pub_ids = []           # Pub.pub_ids for pubs associated in any way with the entity.
        self.prop_dict = {}             # cvterm name-keyed lists of prop objects: e.g., Featureprop.
        self.prop_pub_dict = {}         # prop_id-keyed lists of pub_ids.


class FBFeature(FBDataEntity):
    """An abstract, generic FlyBase Feature entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.db_primary_id = chado_obj.feature_id
        self.chr_flocs = []          # Will be chromosomal Featureloc objects for the entity.
        self.fb_anno_dbxrefs = []    # Will be "FlyBase Annotation IDs" FeatureDbxref objects.
        # Processed FB data.
        self.curr_anno_id = None     # Will be current annotation ID for the gene (str).
        self.alt_anno_ids = []       # Will be list of non-current annotation IDs for the gene (str).


class FBGene(FBFeature):
    """A FlyBase Gene entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.gene_type_names = []     # Will be "promoted_gene_type" Featureprops.
        self.gene_snapshots = []      # Will be "gene_summary_text" Featureprops.
        # Processed FB data.
        self.gene_type_name = None    # Will be SO term name from "promoted_gene_type" Featureprop.
        self.gene_type_id = None      # Will be SO term ID from "promoted_gene_type" Featureprop.


class FBConstruct(FBFeature):
    """A FlyBase Construct entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBConstruct object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        # For constructs, relationships to components are key; may be direct or indirect via allele.
        # Direct relationships (FeatureRelationship objects); ignoring carries_tool and tagged_with.
        self.parent_allele_rels = []      # Direct "associated_with" relationships to FBal allele(s).
        self.encodes_tool_rels = []       # Direct "encodes_tool" relationships to FBto tools.
        self.reg_region_rels = []         # Direct "has_reg_region" relationships to FBto/FBgn/FBsf tools/genes/seqfeats.
        self.seqfeat_rels = []            # Direct FBsf-FBtp "associated_with" relationships (old style).
        # Indirect relationships (lists of (FeatureRelationship, FeatureRelationshipPub) results) via the allele.
        self.al_encodes_tool_rels = []    # Indirect "encodes" relationships: a list of allele-to-FBto/FBsf FeatureRelationship objects.
        self.al_reg_region_rels = []      # Indirect "has_reg_region" relationships: a list of allele-to-FBto/FBsf/FBgn FeatureRelationship objects.
        self.al_genes = []                # Indirect gene relationships: a list of allele-to-FBgn FeatureRelationship objects.
        # Processed FB data.
        # Final relationship assessments for ConstructComponentSlotAnnotationDTO mapping.
        self.expressed_features = {}      # Will be list of feature_id-keyed pub_id list for expressed things: FBgn and FBto.
        self.targeted_features = {}       # Will be list of feature_id-keyed pub_id list for targeted things: FBgn.
        self.regulating_features = {}     # Will be list of feature_id-keyed pub_id list for things that regulate the construct: FBgn, FBto and FBsf.
        self.expressed_tool_genes = []    # Will be list of feature_ids for genes for which related tools are also associated in construct.expressed_features.
        self.regulating_tool_genes = []   # Will be list of feature_ids for genes for which related tools are also associated in construct.targeted_features.


class FBStrain(FBDataEntity):
    """A FlyBase Strain entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBStrain object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.db_primary_id = chado_obj.strain_id
        # Processed FB data.


class FBGenotype(FBDataEntity):
    """A FlyBase Genotype entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBStrain object."""
        super().__init__(chado_obj)
        self.db_primary_id = chado_obj.genotype_id


# Primary Alliance DTO Classes for FB entities.
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


class BiologicalEntityDTO(AuditedObjectDTO):
    """BiologicalEntityDTO class."""
    def __init__(self):
        """Create BiologicalEntityDTO for FlyBase object."""
        super().__init__()
        self.curie = None
        self.taxon_curie = None
        self.data_provider_dto = None
        self.required_fields.extend(['curie', 'taxon_curie', 'data_provider_dto'])


class ReagentDTO(AuditedObjectDTO):
    """ReagentDTO class."""
    def __init__(self):
        """Create ReagentDTO for FlyBase object."""
        super().__init__()
        self.mod_entity_id = None          # Will be the MOD curie.
        self.mod_internal_id = None        # Will be the MOD internal db id.
        self.secondary_identifiers = []    # Will be list of 2o FB IDs (strings)
        self.data_provider_dto = None
        self.required_fields.extend(['mod_entity_id', 'data_provider_dto'])


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
    def __init__(self, evidence_curies: list):
        """Create EvidenceAssociationDTO for FlyBase object.

        Args:
            evidence_curies (list): A list of FB:FBrf or PMID:### curies.

        """
        super().__init__()
        self.evidence_curies = evidence_curies
        self.required_fields.extend([])


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
        self.genomic_entity_curie = genomic_id
        self.note_dtos = []
        self.required_fields.extend(['construct_identifier', 'genomic_entity_relation_name', 'genomic_entity_curie'])


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
