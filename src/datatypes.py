"""Module:: datatypes.

Synopsis:
    Datatype objects representing FlyBaseData retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""


# FlyBase Classes
class FBEntity(object):
    """An abstract, generic FlyBase entity."""
    def __init__(self):
        """Create the generic FlyBase entity with bins for Alliance mapping."""
        self.db_primary_id = None     # The chado table primary key (or concatenation of primary keys).
        self.uniquename = None        # The FlyBase uniquename, if applicable.
        self.uniq_key = None          # A string derived from the uniquely defining properties of the entity.
        self.linkmldto = None         # The Alliance LinkML object containing mapped data.
        self.for_export = True        # Change to False if object should be excluded from export.
        self.internal_reasons = []    # Reasons for marking an object as internal in the export file.
        self.export_warnings = []     # Reasons for suppressing an object from the export file.


class FBAssociation(FBEntity):
    """An abstract, generic FlyBase association/annotation."""
    def __init__(self, chado_objs):
        """Create the generic FlyBase association/annotation object from the main db entry/entries."""
        super().__init__()
        self.chado_objs = chado_objs    # The defining SQLAlchemy chado object(s).


class FBDataEntity(FBEntity):
    """An abstract, generic FlyBase data entity with all it related data, excluding associations/annotations."""
    def __init__(self, chado_obj):
        """Create the generic FlyBase data entity object from the main db entry."""
        super().__init__()
        self.uniquename = chado_obj.uniquename
        try:
            self.organism_abbr = chado_obj.organism.abbreviation
        except AttributeError:
            self.organism_abbr = 'Dmel'
        # Primary FB chado data - direct db query results, no processing.
        self.chado_obj = chado_obj      # The primary SQLAlchemy chado object.
        self.pubs = []                  # Pub associations: e.g., FeaturePub, StrainPub.
        self.synonyms = []              # Synonym associations: e.g., FeatureSynonym, StrainSynonym.
        self.dbxrefs = []               # Dbxref associations: e.g., FeatureDbxref, StrainDbxref.
        self.props = []                 # entity props: e.g., Featureprop, Strainprop.
        self.prop_pubs = []             # entity prop_pubs: e.g., Featureprop_pub, Strainprop_pub.
        self.cvterms = []               # Cvterm associations: e.g., FeatureCvterm, StrainCvterm.
        self.cvtermprops = []           # Cvtermprop associations: e.g., FeatureCvtermprop, StrainDbxref.
        self.timestamps = []            # FB timestamps.
        # Processed FB data - processed from primary FB chado data above.
        self.organism_abbr = None       # The organism.abbreviation (str).
        self.ncbi_taxon_id = None       # The NCBITaxon dbxref.accession (str).
        self.curr_fb_symbol = None      # The current symbol Synonym chado for the entity.
        self.curr_fb_fullname = None    # The current fullname Synonym chado object for the entity.
        self.alt_fb_ids = []            # Secondary FB IDs.
        self.all_pub_ids = []           # Pub.pub_id db IDs for pubs associated in any way with the entity.
        self.prop_dict = {}             # cvterm name-keyed lists of prop objects: e.g., Featureprop, Strainprop.
        self.prop_pub_dict = {}         # prop_id-keyed lists of pub_ids.

    def __str__(self):
        """Basic descriptive info for the object."""
        entity_desc = f'{self.chado_obj.name} ({self.chado_obj.uniquename})'
        return entity_desc


class FBFeature(FBDataEntity):
    """An abstract, generic FlyBase Feature entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.db_primary_id = chado_obj.feature_id
        self.chr_flocs = []        # Will be chromosomal Featureloc objects for the entity.
        self.fb_anno_xrefs = []    # Will be "FlyBase Annotation IDs" FeatureDbxref objects.
        # Processed FB data.


class FBGene(FBFeature):
    """A FlyBase Gene entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.gene_type_names = []    # Will be "promoted_gene_type" Featureprops.
        self.gene_snapshots = []     # Will be "gene_summary_text" Featureprops.
        # Processed FB data.
        self.curr_anno_id = None     # Will be current annotation ID for the gene (str).
        self.alt_anno_ids = []       # Will be list of non-current annotation IDs for the gene (str).


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


# Primary Alliance DTO Classes for FlyBase data.
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
        """Return a JSON-friendly dict for cases where inlined object is required."""
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


class GenomicEntityDTO(BiologicalEntityDTO):
    """GenomicEntityDTO class."""
    def __init__(self):
        """Create GenomicEntityDTO for FlyBase object."""
        super().__init__()
        self.cross_reference_dtos = []


class GeneDTO(GenomicEntityDTO):
    """GeneDTO class."""
    def __init__(self):
        """Create GeneDTO for FlyBase object."""
        super().__init__()
        self.gene_symbol_dto = None             # Will be a single SymbolSlotAnnotationDTO.
        self.gene_full_name_dto = None          # Will be a single FullNameSlotAnnotation.
        self.gene_systematic_name_dto = None    # Will be a single SystematicNameSlotAnnotation.
        self.gene_synonym_dtos = []             # Will be list of NameSlotAnnotationDTO objects.
        self.gene_type_curie = None             # Will be the SO term ID corresponding to the gene's promoted_gene_type.
        self.gene_secondary_id_dtos = []        # Annotation IDs and 2o FlyBase IDs.


class AffectedGenomicModelDTO(GenomicEntityDTO):
    """AffectedGenomicModelDTO class."""
    def __init__(self):
        """Create AffectedGenomicModelDTO for FlyBase object."""
        super().__init__()
        self.name = None                   # The current fullname synonym, ASCII format.
        self.subtype_name = None           # "strain" or "genotype".
        self.agm_secondary_id_dtos = []    # Secondary IDs.
        self.reference_curies = []         # Publication curies (PMID or FBrf).
        self.component_dtos = []           # AffectedGenomicModelComponentDTO objects.
        self.required_fields.extend(['subtype_name'])


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


class SlotAnnotationDTO(AuditedObjectDTO):
    """SlotAnnotationDTO class."""
    def __init__(self):
        """Create a SlotAnnotationDTO for FlyBase object."""
        self.evidence_curies = []


class SecondaryIdSlotAnnotationDTO(SlotAnnotationDTO):
    """SecondaryIdSlotAnnotationDTO class."""
    def __init__(self, secondary_id):
        """Create a SecondaryIdSlotAnnotationDTO for FlyBase object."""
        self.secondary_id = secondary_id
        self.required_fields.extend(['secondary_id'])
