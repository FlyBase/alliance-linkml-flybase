"""Module:: FlyBase datatypes.

Synopsis:
    Ojects representing FlyBase chado data types, for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""


# FlyBase Classes
# Attributes are sorted into "primary chado data" (i.e., query results) and
# processed data (synthesis of sql results).
class FBExportEntity(object):
    """A base, generic FlyBase data export entity."""
    def __init__(self):
        """Create a FBExportEntity object with bins for Alliance mapping."""
        self.db_primary_id = None     # Table primary key (or concatenation of keys).
        self.uniq_key = None          # Uniquely identifying string.
        self.uniquename = None        # FB uniquename, if applicable.
        self.name = None              # FB name, if applicable.
        self.type_id = None           # FB type_id, if applicable.
        self.is_obsolete = None       # FB is_obsolete, if applicable.
        self.is_analysis = None       # FB is_analysis, if applicable.
        self.organism_id = None       # FB organism_id, if applicable.
        self.org_abbr = None          # Organism.abbreviation, if applicable.
        self.org_genus = None         # Organism.genus, if applicable.
        self.org_species = None       # Organism.species, if applicable.
        self.timeaccessioned = None   # FB timeaccessioned, if applicable.
        self.timelastmodified = None  # FB timelastmodified, if applicable.
        self.timestamps = []          # FB timestamps.
        self.linkmldto = None         # Alliance LinkML object for mapped data.
        self.for_export = True        # Made False to prevent Alliance export.
        self.internal_reasons = []    # Reasons an object was marked internal.
        self.export_warnings = []     # Reasons an object was not exported.
        self.entity_desc = None       # A succinct description for this entity.

    def __str__(self):
        """Basic descriptive info for the object."""
        return self.entity_desc


# Primary data entities (have curie and usually a FB web page).
class FBDataEntity(FBExportEntity):
    """An abstract, generic FlyBase data export entity for first class entities.

    First class entities refer to object that typically have a FlyBase curie
    and a dedicated web report: e.g., gene, strain, genotype, gene group. These
    entities are the subjects of relationships and annotations; they typically
    have related "_cvterm", "prop", "_pub" tables, etc.

    """
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
            'timelastmodified'
        ]
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
        # Primary FB chado data from direct db query results, no or minimal processing.
        # These attributes apply to various FlyBase entities: e.g., gene, strain, genotype, gene group, etc.
        self.pubs = []                  # Pub associations: e.g., FeaturePub, StrainPub.
        self.synonyms = []              # Synonym associations: e.g., FeatureSynonym.
        self.fb_sec_dbxrefs = []        # 2o/non-current FlyBase xref objects: e.g., FeatureDbxref.
        self.dbxrefs = []               # Current xref objects: e.g., FeatureDbxref.
        self.props_by_type = {}         # Lists of FBProp objects keyed by prop type name.
        # Processed FB data - processed from primary FB chado data above.
        self.ncbi_taxon_id = None       # The NCBITaxon dbxref.accession (str).
        self.synonym_dict = {}          # Will be synonym_id-keyed dicts of processed synonym info.
        self.curr_fb_symbol = None      # The current symbol for the entity (UTF-8).
        self.curr_fb_fullname = None    # The current fullname for the entity (UTF-8).
        self.alt_fb_ids = []            # Secondary FB IDs (including the "FB:" prefix).
        self.all_pubs = []           # Pub.pub_ids for pubs associated in any way with the entity.
        self.prop_dict = {}             # cvterm name-keyed lists of FBProp objects.


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
        self.curr_anno_id = None     # Will be current annotation ID for the gene, transcript or protein (str).
        self.alt_anno_ids = []       # Will be list of non-current annotation IDs for the gene, transcript or protein (str).


class FBGene(FBFeature):
    """A FlyBase Gene entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.allele_rels = []               # Direct FBal "alleleof" FBgn FeatureRelationships.
        # Processed FB data.
        self.gene_type_name = 'gene'        # Update default gene to SO term name from "promoted_gene_type" Featureprop, if available.
        self.gene_type_id = 'SO:0000704'    # Update default gene ID to SO term ID from "promoted_gene_type" Featureprop, if available.
        self.allele_pubs = {}               # Will be allele_id-keyed list of pub_ids that support the allele-gene relationship.


class FBAllele(FBFeature):
    """A FlyBase Allele entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBAllele object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.parent_gene_rels = []              # Direct FBal "alleleof" FBgn FeatureRelationships.
        self.constructs = []                    # List of FeatureRelationships to associated constructs (FBtp).
        self.dmel_insertions = []               # List of FeatureRelationships to associated Dmel insertions (FBti).
        self.non_dmel_insertions = []           # List of FeatureRelationships to associated non-Dmel insertions (FBti).
        self.args = []                          # List of FeatureRelationships to associated ARG (variation) features.
        self.direct_colls = []                  # List of collecionts (Library objects) directly associated with the allele.
        self.ins_colls = []                     # List of collections (Library objects) indirectly associated with the allele via an FBti insertion.
        self.cons_colls = []                    # List of collections (Library objects) indirectly associated with the allele via an FBtp construct.
        self.sf_colls = []                      # List of collections (Library objects) indirectly associated with the allele via an FBsf feature.
        self.phenstatements = []                # List of SQLAlchemy (Feature, Genotype, Phenotype, Cvterm, Pub) results from Phenstatements.
        # Processed FB data.
        self.parent_gene_id = None              # The FBgn ID fo the allele's parent gene.
        self.adj_org_abbr = 'Dmel'              # Assume allele is Dmel (classical/transgenic) unless it can be shown to be a non-Dmel classical allele.
        self.in_vitro = False                   # Change to True if the allele is associated with an "in vitro%" term.
        self.allele_of_internal_gene = False    # Change to True if the allele is related to an internal-type gene (e.g., origin of replication).


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
        """Create the FBGenotype object."""
        super().__init__(chado_obj)
        self.db_primary_id = chado_obj.genotype_id


# First class associations and annotations.
class FBRelationship(FBExportEntity):
    """FBRelationship class."""
    def __init__(self, chado_table, subject_id, object_id, rel_type):
        """Create a FBRelationship object.

        Args:
            chado_table (str): Name of chado table holding relationship: e.g., feature_relationship, strain_feature.
            subject_id (int): Internal primary_key id for subject (e.g., strain_id for strain_feature table).
            object_id (int): Internal primary_key id for object (e.g., feature_id for strain_feature table).
            rel_type (str): CV term name for the relationship type. Use "unspecified" if not available.

        """
        super().__init__()
        # Primary FB chado data.
        self.chado_table_name = chado_table
        self.subject_id = subject_id
        self.object_id = object_id
        self.rel_type = rel_type
        self.pub_ids = []                      # Will be list of pub_ids supporting the relationship.


class FBAlleleDiseaseAnnotation(FBExportEntity):
    """FBAlleleDiseaseAnnotation class."""
    def __init__(self, feature_cvterm, provenance_prop):
        """Create a FBAlleleDiseaseAnnotation object.

        Args:
            feature_cvterm (FeatureCvterm): FeatureCvterm object.
            provenance_prop (FeatureCvtermprop): FeatureCvtermprop object of type "provenance".

        """
        super().__init__()
        # Primary FB chado data.
        self.feature_cvterm = feature_cvterm
        self.provenance = provenance_prop
        self.db_primary_id = f'{feature_cvterm.feature_cvterm_id}_{provenance_prop.rank}'
        self.evidence_code = None               # Will be the "evidence_code" FeatureCvtermprop.
        self.qualifier = None                   # Will be the "qualifier" FeatureCvtermprop.
        # Processed FB data.
        self.preferred_gene_curie = None        # Get the most appropriate curie for the allele's parental gene.
        self.fb_modifier_type = None
        self.fb_modifier_id = None
        self.modifier_id_was_updated = False    # Change to true if modifier ID in evidence text was updated.
        self.modifier_problem = False           # Change to true if there's a problem finding the modifier allele.
        self.is_redundant = False

    def set_entity_desc(self):
        """Succinct text string describing the disease annotation."""
        desc = 'feature_cvterm_id={}\trank={}\tfeature={} ({})\tcvterm={} (DOID:{})\tpub={}\tqual={}\tevidence={}'.\
            format(self.feature_cvterm.feature_cvterm_id,
                   self.provenance.rank,
                   self.feature_cvterm.feature.name,
                   self.feature_cvterm.feature.uniquename,
                   self.feature_cvterm.cvterm.name,
                   self.feature_cvterm.cvterm.dbxref.accession,
                   self.feature_cvterm.pub.uniquename,
                   self.qualifier.value,
                   self.evidence_code.value)
        self.entity_desc = desc
        return


# Second class annotations (submitted as part of other objects).
class FBProp(object):
    """A generic FlyBase property annotation for a first class entity."""
    def __init__(self, chado_obj):
        """Create a FBProp object from the main db entry.

        Args:
            chado_obj (SQLAlchemy object): The object representing the prop table entry.

        Returns:
            A FBProp object.
        """
        self.chado_obj = chado_obj    # The prop object: e.g., Featureprop, Strainprop, FeatureRelationshipprop, etc.
        self.pubs = []                # Will be a list of pub_ids.
