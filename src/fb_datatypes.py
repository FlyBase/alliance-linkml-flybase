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
        self.pubs = []                  # Pub associations: e.g., FeaturePub, StrainPub.
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
        self.allele_rels = []         # Direct FBal "alleleof" FBgn FeatureRelationships.
        # Processed FB data.
        self.gene_type_name = 'gene'        # Will be SO term name from "promoted_gene_type" Featureprop, if available.
        self.gene_type_id = 'SO:0000704'    # Will be SO term ID from "promoted_gene_type" Featureprop, if available.
        self.alleles = {}                   # Will be allele_id-keyed list of pub curies.


class FBConstruct(FBFeature):
    """A FlyBase Construct entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBConstruct object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        # For constructs, relationships to components are key; may be direct or indirect via allele.
        # Direct relationships (FeatureRelationship objects); ignoring carries_tool and tagged_with.
        self.parent_allele_rels = []      # Direct FBal "associated_with" FBtp FeatureRelationships.
        self.encodes_tool_rels = []       # Direct FBtp "encodes_tool" FBto FeatureRelationships.
        self.reg_region_rels = []         # Direct FBtp "has_reg_region" FBto/FBgn/FBsf FeatureRelationships.
        self.seqfeat_rels = []            # Direct FBsf "associated_with" FBtp FeatureRelationships (old style).
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
