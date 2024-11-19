"""Module:: insertion_handler.

Synopsis:
    A data handler that exports FlyBase data for insertions to Alliance Gene LinkML
    objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
import agr_datatypes
from fb_datatypes import FBInsertion
from feature_handler import FeatureHandler
from harvdev_utils.production import (
    Cvterm, Feature, FeatureRelationship, Library, LibraryFeature, LibraryFeatureprop
)


class InsertionHandler(FeatureHandler):
    """This object gets, synthesizes and filters insertion data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the InsertionHandler object."""
        super().__init__(log, testing)
        self.datatype = 'insertion'
        self.fb_export_type = FBInsertion
        self.agr_export_type = agr_datatypes.AlleleDTO
        self.primary_export_set = 'allele_ingest_set'

    # Types: 228747 transposable_element_insertion_site; 7726 insertion_site; 5753 transposable element; 3573 match (internal).
    # Relationships: 234754 FBti(producedby)FBtp; 64920 FBal(associated_with)FBti.
    test_set = {
        'FBti0000040': 'P{hsneo}Xrp1[142]',     # type=transposable_element_insertion_site. Location trap. FBal-(associated_with)->FBti-(producedby)->FBtp.
        'FBti0151770': 'P{UAS-stnB.M}vl',       # type=transposable_element_insertion_site. FBti-(producedby)->FBtp<-(associated_with)-FBal.
        'FBti0167947': 'TI{TI}wg[GFP]',         # type=insertion_site. FBti-(producedby)->FBtp<-(associated_with)-FBal.
        'FBti0018862': '17.6{}804',             # type=17.6{}804; this insertion shares its uniquename with two internal "match" features.
        'FBti0016979': 'P{PZ}Vha44[06072b]',    # type=transposable_element_insertion_site. Direct "associated_with" association to a gene.
    }

    # Additional export sets.
    gene_insertion_rels = {}            # Will be (insertion feature_id, gene feature_id) tuples keying lists of FBRelationships.
    gene_insertion_associations = []    # Will be the final list of gene-insertion FBRelationships to export.

    # Elaborate on get_general_data() for the InsertionHandler.
    def get_general_data(self, session):
        """Extend the method for the InsertionHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        self.build_feature_lookup(session, feature_types=['gene', 'allele', 'construct'])
        self.build_allele_gene_lookup(session)
        self.get_internal_genes(session)
        return

    # Additional sub-methods for get_datatype_data().
    def get_direct_collections(self, session):
        """Find collections directly related to insertions."""
        self.log.info('Find collections directly related to insertions.')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        filters = (
            Feature.uniquename.op('~')(self.regex['insertion']),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.regex['library']),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection'
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        collections = session.query(Feature, Library).\
            select_from(Feature).\
            join(LibraryFeature, (LibraryFeature.feature_id == Feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in collections:
            self.fb_data_entities[result.Feature.feature_id].direct_colls.append(result.Library)
            counter += 1
        self.log.info(f'Found {counter} direct insertion-collection associations.')
        return

    def get_indirect_collections(self, session):
        # Placeholder.
        return

    def get_sf_collections(self, session):
        """Find collections indirectly related to insertions via sequence features."""
        self.log.info('Find collections indirectly related to insertions via sequence features.')
        # Placeholder.
        # insertion = aliased(Feature, name='insertion')
        # construct = aliased(Feature, name='construct')
        # seqfeat = aliased(Feature, name='seqfeat')
        # libtype = aliased(Cvterm, name='libtype')
        # libfeattype = aliased(Cvterm, name='libfeattype')
        # featreltype = aliased(Cvterm, name='featreltype')
        # insertion_construct = aliased(FeatureRelationship, name='insertion_construct')
        # seqfeat_construct = aliased(FeatureRelationship, name='seqfeat_construct')
        # filters = (
        #     insertion.uniquename.op('~')(self.regex['insertion']),
        #     construct.uniquename.op('~')(self.regex['construct']),
        #     seqfeat.uniquename.op('~')(self.regex['seqfeat']),
        #     construct.is_obsolete.is_(False),
        #     seqfeat.is_obsolete.is_(False),
        #     Library.is_obsolete.is_(False),
        #     Library.uniquename.op('~')(self.regex['library']),
        #     libtype.name == 'reagent collection',
        #     libfeattype.name == 'member_of_reagent_collection',
        #     featreltype.name == 'associated_with'
        # )
        # if self.testing:
        #     self.log.info(f'TESTING: limit to these entities: {self.test_set}')
        #     filters += (insertion.uniquename.in_((self.test_set.keys())), )
        # sf_collections = session.query(insertion, Library).\
        #     select_from(insertion).\
        #     join(insertion_construct, (insertion_construct.subject_id == insertion.feature_id)).\
        #     join(construct, (construct.feature_id == insertion_construct.object_id)).\
        #     join(featreltype, (featreltype.cvterm_id == insertion_construct.type_id)).\
        #     join(seqfeat_construct, (seqfeat_construct.object_id == construct.feature_id)).\
        #     join(seqfeat, (seqfeat.feature_id == seqfeat_construct.subject_id)).\
        #     join(LibraryFeature, (LibraryFeature.feature_id == seqfeat.feature_id)).\
        #     join(Library, (Library.library_id == LibraryFeature.library_id)).\
        #     join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
        #     join(libtype, (libtype.cvterm_id == Library.type_id)).\
        #     join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
        #     filter(*filters).\
        #     distinct()
        # counter = 0
        # for result in sf_collections:
        #     self.fb_data_entities[result.insertion.feature_id].sf_colls.append(result.Library)
        #     counter += 1
        # self.log.info(f'Found {counter} sequence feature-mediated insertion-library associations.')
        return

    # Elaborate on get_datatype_data() for the InsertionHandler.
    def get_datatype_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject')
        self.get_entity_relationships(session, 'object')
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_direct_collections(session)
        # self.get_indirect_collections(session)    # BOB: Need to update this for insertions.
        # self.get_sf_collections(session)          # BOB: Need to update this for insertions.
        return

    # Add sub-methods to be run by synthesize_info() below.
    def synthesize_parent_genes(self):
        """Get insertion parent gene IDs."""
        self.log.info('Get insertion parent gene IDs.')
        # Placeholder.
        return

    def flag_insertions_of_internal_genes(self):
        """Flag insertions of internal genes."""
        self.log.info('Flag insertions of internal genes.')
        # Placeholder - insertion can have many genes associated - make internal only if all related genes are internal.
        return

    def synthesize_gene_insertions(self):
        """Synthesize gene-insertion relationships."""
        self.log.info('Synthesize gene-insertion relationships.')
        # Placeholder.
        return

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_info(self):
        """Extend the method for the GeneHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_parent_genes()
        self.flag_insertions_of_internal_genes()
        self.synthesize_gene_insertions()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_insertion_basic(self):
        """Map basic FlyBase insertion data to the Alliance LinkML object."""
        self.log.info('Map basic insertion info to Alliance object.')
        for insertion in self.fb_data_entities.values():
            agr_insertion = self.agr_export_type()
            agr_insertion.obsolete = insertion.chado_obj.is_obsolete
            agr_insertion.mod_entity_id = f'FB:{insertion.uniquename}'
            agr_insertion.mod_internal_id = str(insertion.chado_obj.feature_id)
            agr_insertion.taxon_curie = insertion.ncbi_taxon_id
            insertion.linkmldto = agr_insertion
        return

    def map_insertion_database_status(self):
        """Map insertion database status."""
        self.log.info('Map insertion database status.')
        evidence_curies = []
        for insertion in self.fb_data_entities.values():
            if insertion.is_obsolete is False:
                db_status = 'approved'
            else:
                db_status = 'deleted'
            db_status_annotation = agr_datatypes.AlleleDatabaseStatusSlotAnnotationDTO(db_status, evidence_curies)
            insertion.linkmldto.allele_database_status_dto = db_status_annotation.dict_export()
        return

    def map_collections(self):
        """Map insertion collections."""
        self.log.info('Map insertion collections.')
        for insertion in self.fb_data_entities.values():
            collections = []
            if insertion.direct_colls:
                collections.extend(insertion.direct_colls)
            elif insertion.ins_colls:
                collections.extend(insertion.ins_colls)
            elif insertion.cons_colls:
                collections.extend(insertion.cons_colls)
            elif insertion.sf_colls:
                collections.extend(insertion.sf_colls)
            if collections:
                collections = list(set(collections))
                insertion.linkmldto.in_collection_name = collections[0].name
                if len(collections) > 1:
                    self.log.warning(f'{insertion} has many relevant collections: {[i.name for i in collections]}')
        return

    def map_extinction_info(self):
        """Map extinction info."""
        self.log.info('Map extinction info.')
        counter = 0
        for insertion in self.fb_data_entities.values():
            if 'availability' in insertion.props_by_type.keys():
                for prop in insertion.props_by_type['availability']:
                    if prop.chado_obj.value == 'Stated to be lost.':
                        insertion.linkmldto.is_extinct = True
            for prop_type in insertion.props_by_type.keys():
                if prop_type.startswith('derived_stock'):
                    # Stock availability trumps curated extinction comment.
                    if insertion.linkmldto.is_extinct is True:
                        self.log.warning(f'{insertion} is stated to be lost but has stocks.')
                    insertion.linkmldto.is_extinct = False
            if insertion.linkmldto.is_extinct is True:
                counter += 1
        self.log.info(f'Flagged {counter} insertions as extinct.')
        return

    def map_internal_insertion_status(self):
        """Flag internal insertions using insertion-specific criteria."""
        self.log.info('Flag internal insertions using insertion-specific criteria.')
        internal_gene_counter = 0
        non_dmel_drosophilid_counter = 0
        for insertion in self.fb_data_entities.values():
            if insertion.insertion_of_internal_gene is True:
                insertion.linkmldto.internal = True
                insertion.internal_reasons.append('Insertion related to an internal type FB gene')
                internal_gene_counter += 1
            if insertion.org_abbr != 'Dmel':
                insertion.linkmldto.internal = True
                insertion.internal_reasons.append('A non-Dmel insertion')
                non_dmel_drosophilid_counter += 1
        self.log.info(f'Flagged {internal_gene_counter} insertions of internal-type genes as internal.')
        self.log.info(f'Flagged {non_dmel_drosophilid_counter} non-Dmel Drosophilid insertions as internal.')
        return

    def map_mutation_types(self):
        """Map mutation types."""
        self.log.info('Map mutation types.')
        counter = 0
        for insertion in self.fb_data_entities.values():
            mutation_type_curie = self.cvterm_lookup[insertion.type_id]['curie']
            pub_curies = []
            mutant_type_annotation = agr_datatypes.AlleleMutationTypeSlotAnnotationDTO(mutation_type_curie, pub_curies)
            insertion.linkmldto.allele_mutation_type_dtos.append(mutant_type_annotation.dict_export())
            counter += 1
        self.log.info(f'Mapped {counter} mutation type annotations.')
        return

    def map_gene_insertion_associations(self):
        """Map gene-insertion associations to Alliance object."""
        self.log.info('Map gene-insertion associations to Alliance object.')
        # Placeholder
        # INSERTION = 0
        # GENE = 1
        # counter = 0
        # for insertion_gene_key, insertion_gene_rels in self.gene_insertion_rels.items():
        #     insertion = self.fb_data_entities[insertion_gene_key[INSERTION]]
        #     insertion_curie = f'FB:{insertion.uniquename}'
        #     gene = self.feature_lookup[insertion_gene_key[GENE]]
        #     gene_curie = f'FB:{gene["uniquename"]}'
        #     first_feat_rel = insertion_gene_rels[0]
        #     all_pub_ids = []
        #     for insertion_gene_rel in insertion_gene_rels:
        #         all_pub_ids.extend(insertion_gene_rel.pubs)
        #     first_feat_rel.pubs = all_pub_ids
        #     pub_curies = self.lookup_pub_curies(all_pub_ids)
        #     rel_dto = agr_datatypes.AlleleGeneAssociationDTO(insertion_curie, 'is_allele_of', gene_curie, pub_curies)
        #     if insertion.is_obsolete is True or gene['is_obsolete'] is True:
        #         rel_dto.obsolete = True
        #         rel_dto.internal = True
        #     first_feat_rel.linkmldto = rel_dto
        #     self.gene_insertion_associations.append(first_feat_rel)
        #     counter += 1
        # self.log.info(f'Generated {counter} insertion-gene unique associations.')
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_insertion_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_internal_insertion_status()
        self.map_extinction_info()
        self.map_collections()
        self.map_insertion_database_status()
        self.map_mutation_types()
        # self.map_pubs()    # TEMPORARILY SUPPRESS UNTIL LOAD SPEED IMPROVES
        self.map_timestamps()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        # self.map_gene_insertion_associations()    # BOB
        # self.flag_internal_fb_entities('gene_insertion_associations')    # BOB
        return

    # Elaborate on query_chado_and_export() for the GeneHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the GeneHandler."""
        super().query_chado_and_export(session)
        # self.flag_unexportable_entities(self.gene_insertion_associations, 'allele_gene_association_ingest_set')
        # self.generate_export_dict(self.gene_insertion_associations, 'allele_gene_association_ingest_set')
        return
