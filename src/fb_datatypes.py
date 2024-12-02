"""Module:: FlyBase datatypes.

Synopsis:
    Objects representing FlyBase chado data types, for export to the Alliance in
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
        self.entity_desc = ''         # A succinct description for this entity.

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
        self.pub_associations = []        # Pub associations: e.g., FeaturePub, StrainPub.
        self.synonyms = []                # Synonym associations: e.g., FeatureSynonym.
        self.fb_sec_dbxrefs = []          # 2o/non-current FlyBase xref objects: e.g., FeatureDbxref.
        self.dbxrefs = []                 # Current xref objects: e.g., FeatureDbxref.
        self.props_by_type = {}           # Lists of FBProp objects keyed by prop type name.
        self.cvt_annos_by_id = {}         # entity_cvterm_id-keyed dict of FBCVtermAnnotation objects.
        self.cvt_anno_ids_by_cv = {}      # Cv.name-keyed lists of entity_cvterm_ids.
        self.cvt_anno_ids_by_term = {}    # Cvterm.name-keyed lists of entity_cvterm_ids.
        self.cvt_anno_ids_by_prop = {}    # Cvtermprop type (name) keyed lists of entity_cvterm_ids.
        self.rels_by_id = {}              # Internal relationship_id-keyed dict of FBRelationship entities.
        self.sbj_rel_ids_by_type = {}     # Lists of relationship_ids where this entity is the subject; keyed by relationship type name.
        self.obj_rel_ids_by_type = {}     # Lists of relationship_ids where this entity is the object; keyed by relationship type name.
        # Processed FB data - processed from primary FB chado data above.
        self.ncbi_taxon_id = None       # The NCBITaxon dbxref.accession (str).
        self.synonym_dict = {}          # Will be synonym_id-keyed dicts of processed synonym info.
        self.curr_fb_symbol = None      # The current symbol for the entity (UTF-8).
        self.curr_fb_fullname = None    # The current fullname for the entity (UTF-8).
        self.alt_fb_ids = []            # Secondary FB IDs (including the "FB:" prefix).
        self.all_pubs = []              # Pub.pub_ids for pubs associated in any way with the entity.

    def recall_relationships(self, log, **kwargs):
        """Recall FBRelationship objects with specified attributes from previous db queries for entity relationships.

        Args:
            log (Logger): The logger for any messages.

        Keyword Args:
            entity_role (str): Indicate if entity is the subject or object in the relationship.
            rel_types (str, list): The CV term name(s) for the relationship of interest.
            rel_entity_types (str, list): The CV term name(s) for the related entity type; only useful for features.
            At least one keyword argument must be specified.

        Raises:
            Raises an error if no keyword argument is specified.

        Returns:
            A list of FBRelationship objects that meet the specified criteria.

        """
        rel_ids_of_interest = []
        rels_of_interest = []
        rel_buckets = {
            'subject': 'sbj_rel_ids_by_type',
            'object': 'obj_rel_ids_by_type',
        }
        additional_rel_buckets = {
            'subject': 'sbj_rel_ids_by_obj_type',
            'object': 'obj_rel_ids_by_sbj_type',
        }
        # Check for allowable entity_role value.
        entity_role = None
        if 'entity_role' in kwargs.keys():
            entity_role = kwargs['entity_role']
        if entity_role not in rel_buckets.keys():
            log.error(f'For recall_relationships(), the entity role given is not allowed: "{entity_role}"; must be "subject" or "object" only.')
            raise ValueError
        # Convert rel_types to list if applicable.
        rel_types = None
        if 'rel_types' in kwargs.keys():
            rel_types = kwargs['rel_types']
            if type(rel_types) is str:
                rel_types = [rel_types]
        # Check that rel_entity_types is allowed for the relevant FB datatype.
        rel_entity_types = None
        if 'rel_entity_types' in kwargs.keys():
            if not isinstance(self, FBFeature):
                log.error(f'For recall_relationships(), cannot specify "rel_entity_types" for {type(self)} entities, only for feature types.')
                raise ValueError
            else:
                rel_entity_types = kwargs['rel_entity_types']
                if type(rel_entity_types) is str:
                    rel_entity_types = [rel_entity_types]
        # Check that rel_types and/or rel_entity_types are defined.
        if rel_types is None and rel_entity_types is None:
            log.error('For recall_relationships(), must provide rel_types and/or rel_entity_types keyword arguments: neither were given.')
            raise ValueError
        # Check for unexpected kwargs.
        good_kwargs = {'entity_role', 'rel_types', 'rel_entity_types'}
        bad_kwargs = kwargs.keys() - good_kwargs
        if bad_kwargs:
            log.error(f'For recall_relationships(), found unexpected kwargs: {bad_kwargs}. Only these are allowed: {good_kwargs}.')
            raise ValueError
        # Get the initial list of relevant relationship_ids.
        if entity_role:
            relevant_bucket_names = [rel_buckets[entity_role]]
        else:
            relevant_bucket_names = list(rel_buckets.values())
        for bucket_name in relevant_bucket_names:
            relevant_rel_dict = getattr(self, bucket_name)
            if rel_types:
                for rel_type in rel_types:
                    if rel_type in relevant_rel_dict.keys():
                        rel_ids_of_interest.extend(relevant_rel_dict[rel_type])
            else:
                for rel_list in relevant_rel_dict.values():
                    rel_ids_of_interest.extend(rel_list)
        # Further filter the list by related entity type if applicable.
        if rel_entity_types:
            additional_rel_ids_of_interest = []
            if entity_role:
                additional_relevant_bucket_names = [additional_rel_buckets[entity_role]]
            else:
                additional_relevant_bucket_names = list(additional_rel_buckets.values())
            for additional_bucket_name in additional_relevant_bucket_names:
                additional_relevant_rel_dict = getattr(self, additional_bucket_name)
                for rel_entity_type in rel_entity_types:
                    if rel_entity_type in additional_relevant_rel_dict.keys():
                        additional_rel_ids_of_interest.extend(additional_relevant_rel_dict[rel_entity_type])
            rel_ids_of_interest = list(set(rel_ids_of_interest).intersection(set(additional_rel_ids_of_interest)))
        for rel_id in rel_ids_of_interest:
            rels_of_interest.append(self.rels_by_id[rel_id])
        return rels_of_interest

    def recall_cvterm_annotations(self, log, **kwargs):
        """Recall FBCVTermAnnotation objects with specified attributes from previous db queries for CV term annotations.

        Args:
            log (Logger): The logger for any messages.

        Keyword Args:
            cv_names (str, list): The CV name(s) for annotations of interest.
            cvterm_names (str, list): The CV term name(s) for annotations of interest.
            prop_type_names (str, list): The name of prop type(s) for annotations of interest.
            One or both keyword arguments must be specified.

        Raises:
            Raises an error if no keyword argument is specified.

        Returns:
            A list of FBCVTermAnnotation objects that meet the specified criteria.

        """
        anno_bin_types = {
            'cv_names': 'cvt_anno_ids_by_cv',
            'cvterm_names': 'cvt_anno_ids_by_term',
            'prop_type_names': 'cvt_anno_ids_by_prop',
        }
        if not kwargs:
            log.error(f'Must give at least one keyword argument for recall_cvterm_annotations() method: {anno_bin_types.keys()}.')
            raise ValueError
        unrecognized_kwargs = set(kwargs.keys()) - set(anno_bin_types.keys())
        if unrecognized_kwargs:
            log.error(f'Unrecognized keyword args were given: {unrecognized_kwargs}. Only these are accepted: {anno_bin_types.keys()}.')
            raise ValueError
        anno_ids_of_interest = set(self.cvt_annos_by_id.keys())
        log.info(f'BOB: Start with {len(anno_ids_of_interest)} annotations.')
        annos_of_interest = []
        for kwarg_name in anno_bin_types.keys():
            if kwarg_name not in kwargs.keys():
                continue
            if type(kwargs[kwarg_name]) is str:
                kwargs[kwarg_name] = [kwargs[kwarg_name]]
            log.debug(f'BOB: Get annotations for {kwarg_name}={kwargs[kwarg_name]}')
            relevant_anno_dict = getattr(self, anno_bin_types[kwarg_name])
            log.debug(f'BOB: Look for annotations in this bucket: {anno_bin_types[kwarg_name]}')
            relevant_anno_ids = []
            # Take union of all specified values for a given attribute.
            for anno_attr_name in kwargs[kwarg_name]:
                relevant_anno_ids.extend(relevant_anno_dict[anno_attr_name])
            relevant_anno_ids = set(relevant_anno_ids)
            # Apply the filter.
            anno_ids_of_interest = anno_ids_of_interest.intersection(relevant_anno_ids)
            log.info(f'BOB: Down to {len(anno_ids_of_interest)} annotations after {kwarg_name} filter.')
        for anno_id in anno_ids_of_interest:
            annos_of_interest.append(self.cvt_annos_by_id[anno_id])
        return annos_of_interest


class FBFeature(FBDataEntity):
    """An abstract, generic FlyBase feature entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.db_primary_id = chado_obj.feature_id
        self.chr_flocs = []                  # Will be chromosomal Featureloc objects for the entity.
        self.fb_anno_dbxrefs = []            # Will be "FlyBase Annotation IDs" FeatureDbxref objects.
        self.sbj_rel_ids_by_obj_type = {}    # Lists of relationship_ids where this entity is the subject; keyed by object feature type.
        self.obj_rel_ids_by_sbj_type = {}    # Lists of relationship_ids where this entity is the object; keyed by subject feature type.
        self.reagent_colls = []              # List of reagent collections (Library objects) directly associated with the feature.
        self.al_reagent_colls = []           # List of reagent collections (Library objects) indirectly associated with the feature via an FBal allele.
        self.ti_reagent_colls = []           # List of reagent collections (Library objects) indirectly associated with the feature via an FBti insertion.
        self.tp_reagent_colls = []           # List of reagent collections (Library objects) indirectly associated with the feature via an FBtp construct.
        self.sf_reagent_colls = []           # List of reagent collections (Library objects) indirectly associated with the feature via an FBsf feature.
        # Processed FB data.
        self.curr_anno_id = None     # Will be current annotation ID for the gene, transcript or protein (str).
        self.alt_anno_ids = []       # Will be list of non-current annotation IDs for the gene, transcript or protein (str).


class FBAberration(FBFeature):
    """A FlyBase aberration entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBAberration object."""
        super().__init__(chado_obj)


class FBAllele(FBFeature):
    """A FlyBase allele entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBAllele object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.phenstatements = []                # List of SQLAlchemy (Feature, Genotype, Phenotype, Cvterm, Pub) results from Phenstatements.
        # Processed FB data.
        self.parent_gene_id = None              # The FBgn ID for the allele's parent gene.
        self.cons_rels = []                     # List of current cons FBRelationships.
        self.dmel_ins_rels = []                 # List of current Dmel FBti FBRelationships.
        self.non_dmel_ins_rels = []             # List of current non-Dmel FBti FBRelationships.
        self.arg_rels = []                      # List of current ARG FBRelationships.
        self.adj_org_abbr = 'Dmel'              # Assume allele is Dmel (classical/transgenic) unless it can be shown to be a non-Dmel classical allele.
        self.in_vitro = False                   # Change to True if the allele is associated with an "in vitro%" term.
        self.allele_of_internal_gene = False    # Change to True if the allele is related to an internal-type gene (e.g., origin of replication).


class FBBalancer(FBFeature):
    """A FlyBase balancer entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBBalancer object."""
        super().__init__(chado_obj)


class FBConstruct(FBFeature):
    """A FlyBase construct entity with all its related data."""
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


class FBGene(FBFeature):
    """A FlyBase gene entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGene object."""
        super().__init__(chado_obj)
        # Processed FB data.
        self.gene_type_name = 'gene'        # Update this default gene to SO term name from "promoted_gene_type" Featureprop, if available.
        self.gene_type_id = 'SO:0000704'    # Update this default gene ID to SO term ID from "promoted_gene_type" Featureprop, if available.


class FBInsertion(FBFeature):
    """A FlyBase insertion entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBInsertion object."""
        super().__init__(chado_obj)
        # Processed FB data.
        self.parent_gene_ids = []                  # List of what? BOB


class FBStrain(FBDataEntity):
    """A FlyBase strain entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBStrain object."""
        super().__init__(chado_obj)
        # Primary FB chado data.
        self.db_primary_id = chado_obj.strain_id
        # Processed FB data.


class FBGenotype(FBDataEntity):
    """A FlyBase genotype entity with all its related data."""
    def __init__(self, chado_obj):
        """Create the FBGenotype object."""
        super().__init__(chado_obj)
        self.db_primary_id = chado_obj.genotype_id


# First class associations and annotations.
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


class FBRelationship(FBExportEntity):
    """FBRelationship class."""
    def __init__(self, chado_obj, table_name):
        """Create a FBRelationship object, limited to associations of entities in the same table.

        Args:
            chado_obj (SQLAlchemy object): The Chado object representing the relationship: e.g., FeatureRelationship, StrainRelationship.
            table_name (str): The relationship table name for the relationship.

        """
        super().__init__()
        self.chado_obj = chado_obj
        self.db_primary_id = getattr(chado_obj, f'{table_name}_id')
        self.entity_desc = f'{table_name}_id={self.db_primary_id}'
        self.props_by_type = {}    # Lists of FBProp objects keyed by prop type name.
        self.pubs = []    # Will be list of Pub.pub_ids supporting the relationship.


class FBCVTermAnnotation(FBExportEntity):
    """FBCVTermAnnotation class."""
    def __init__(self, chado_obj, table_name):
        """Create a FBCVTermAnnotation object.

        Args:
            chado_obj (SQLAlchemy object): The Chado object representing the CV term annotation: e.g., FeatureCvterm, GrpCvterm.
            table_name (str): The relationship table name for the relationship: e.g., feature_cvterm, strain_cvterm.

        """
        super().__init__()
        self.chado_obj = chado_obj
        self.db_primary_id = getattr(chado_obj, f'{table_name}_id')
        self.entity_desc = f'{table_name}_id={self.db_primary_id}'
        self.props_by_type = {}    # Lists of FBProp objects keyed by prop type name.
        self.pub_id = self.chado_obj.pub_id


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
