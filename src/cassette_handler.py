"""Module:: cassette_handler.

Synopsis:
    A data handler that exports FlyBase data for cassette to Alliance
    TransgenicTool LinkML objects.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

# import csv
# import re
from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class CassetteHandler(FeatureHandler):
    """This object gets, synthesizes and filters cassette data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the CassetteHandler object."""
        super().__init__(log, testing)
        self.datatype = 'cassette'
        self.fb_export_type = fb_datatypes.FBCassette
        self.agr_export_type = agr_datatypes.CassetteDTO
        self.primary_export_set = 'cassette_ingest_set'
        self.incremental_update = False

    test_set = {
        'FBal0386858': 'SppL[CR70402-TG4.1]',   # Insertion allele superceded by FBti0226866 (superseded_by_at_locus_insertion).
        'FBal0386745': 'Scer_GAL4[SppL...]',    # Insertion allele superceded by FBti0226866 too (superseded_by_at_locus_insertion).
        'FBal0137236': 'gukh[142]',             # Insertion allele superceded by FBti0000040 (superseded_by_at_locus_insertion).
        'FBal0137618': 'Xrp1[142]',             # Insertion allele superceded by FBti0000040 too (superseded_by_at_locus_insertion).
        'FBal0036007': 'wg[en11]',              # Insertion allele superceded by FBti0002065 (superseded_by_at_locus_insertion).
        'FBal0018482': 'wg[1]',                 # X-ray mutation.
        'FBal0015148': 'Sb[Spi]',               # point mutation.
        'FBal0043981': 'Ecol_lacZ[en-14]',      # Has an allele full name. Relationship to ARG has no pub support, superceded by FBti0002067.
        'FBal0279489': 'Scer_GAL4[how-24B]',    # Has a 2o ID, superceded by FBti0150063.
        'FBal0000010': 'alphaTub67C[3]',        # Has semidominant annotation from single locus homozygous genotype.
        'FBal0403680': 'Atg8a[3-31]',           # Has recessive annotation from single locus homozygous genotype.
        'FBal0011189': 'sti[1]',                # Has recessive annotation from single locus hemizygous genotype over a deficiency.
        'FBal0007942': '18w[00053]',            # Has recessive annotation from single locus unspecified zygosity genotype.
        'FBal0015410': 'sei[2]',                # Has codominant annotation from single locus unspecified zygosity genotype.
        'FBal0198096': 'tal[1]',                # Allele of internal type gene tal (gene_with_polycistronic_transcript).
        'FBal0055793': 'wg[UAS.Tag:HA]',        # Allele is directly related to a construct, superceded by FBti0256568 (superseded_by_transgnc_insertions).
        'FBal0048226': 'Dmau_w[a23]',           # Non-Dmel allele related to non-Dmel insertion, superceded by FBti0014970 (superseded_by_at_locus_insertion).
        'FBal0011649': 'Dsim_Lhr[1]',           # Non-Dmel classical allele.
        'FBal0043132': 'Hsap_MAPT[UAS.cAa]',    # Transgenic, superceded by FBti0000969, FBti0249419 (superseded_by_transgnc_insertions).
        'FBal0062057': 'Scer_CDC42[V12.hs]',    # Transgenic, superceded by FBti0012506, FBti0249909 (superseded_by_transgnc_insertions).
        'FBal0198528': 'CG33269[HMJ22303]',     # Transegnic, superceded by four FBti (superseded_by_transgnc_insertions).
        'FBal0322755': 'Mcm3[+tBa]',
        'FBal0322754': 'flfl[DeltaRanBD.UAS.Venus]',
        'FBal0296109': 'sSemp1[R41G.UAS]',
        'FBal0193766': 'Gr63a[UAS.cJa]',
        'FBal0239883': 'sd[RNAi.N.UAS]',
        'FBal0000531': 'Amy-p[IX]',
        'FBal0028742': 'Act88F[E334K]',
        'FBal0212171': r'Avic\GFP[UAS.FRT1]',  # in vitro only
        'FBal0290956': 'Csas[21]',  # in vitro only
        'FBal0392043': r'Avic\GFP[EYFP.3xP3.cUa]',  # in vitro only
    }

    cassette_prop_to_note_mapping = {
        'description': ('summary', 'note_dtos'),
        'misc': ('comment', 'note_dtos'),
        'aminoacid_rep': ('comment', 'note_dtos'),
        'molecular_info': ('comment', 'note_dtos'),
        'nucleotide_sub': ('comment', 'note_dtos'),
        # At the moment, just for code development. (line below)
        # 'internal_notes': ('internal_note', 'note_dtos'),
    }
    cassette_associations = []
    cassette_cassette_rels = {}

    def get_general_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_feature_lookup(session, feature_types=['cassette', 'construct', 'allele', 'tool', 'gene', 'seqfeat'])

    def get_entities(self, session, **kwargs):
        """Extend the method for the CassetteHandler."""
        reference_set = False
        if 'reference' in kwargs.keys() and kwargs['reference'] is True:
            reference_set = True
            self.incremental_update = True

        # Get the main set of cassettes
        self.get_main_entities(session, reference_set)
        # Get in vitro set of cassettes
        self.add_in_vitro_allele_entries(session, reference_set)
        if self.testing:
            self.log.debug("BOB: print list")
            for bob in self.fb_data_entities:
                self.log.debug(f"BOB: {bob}")

    def add_in_vitro_allele_entries(self, session, reference_set):
        """Extend list of entities."""
        self.log.info('Add entities for alleles having "in vitro construct" annotations.')
        counter = 0
        results = self.get_cassettes_in_vitro_entries(session)
        for result in results:
            pkey_id = getattr(result, 'feature_id')
            if reference_set is True:
                self.fb_reference_entity_ids.append(pkey_id)
            else:
                if pkey_id not in self.fb_data_entities:
                    self.fb_data_entities[pkey_id] = self.fb_export_type(result)
                    counter += 1
        if reference_set is True:
            self.log.info(f'Found {counter} FlyBase {self.datatype} in vitro entities in reference chado instance not in main set.')
        else:
            self.log.info(f'Found {counter} FlyBase {self.datatype} in vitro entities in chado not in main set.')

    def get_main_entities(self, session, reference_set):
        """Get simple FlyBase cassette/allele data entities.

        So Subject FBal linked to
           Object FBtp  linked via
           Cvterm 'associated_with'
        Args:
            session (Session): SQLAlchemy session for the query.
            reference_set (bool): If True, retrieves only non-obsolete objects from
                              a previous reference database; for incremental
                              updates.

        """
        if self.datatype in self.feat_type_export.keys():
            chado_type = 'feature'
        else:
            chado_type = self.datatype

        pkey_name = f'{chado_type}_id'
        self.log.info(f'Have this primary_key name: {pkey_name}')
        counter = 0
        results = self.get_cassette_main_entities(session, reference_set)
        for result in results:
            pkey_id = getattr(result, pkey_name)
            if self.testing:
                self.log.debug(f"BOB: {pkey_id}: {result}")
            if reference_set is True:
                self.fb_reference_entity_ids.append(pkey_id)
            else:
                self.fb_data_entities[pkey_id] = self.fb_export_type(result)
            counter += 1
        if reference_set is True:
            self.log.info(f'Found {counter} FlyBase {self.datatype} main entities in reference chado instance.')
        else:
            self.log.info(f'Found {counter} FlyBase {self.datatype} main entities in chado.')

    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_cassette_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_entity_props_to_notes('cassette_prop_to_note_mapping')
        # self.map_xrefs()
        self.map_secondary_ids('secondary_identifiers')
        self.map_cassette_associations()

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_cassette_basic(self):
        """Map basic FlyBase transgenic cassette data to the Alliance LinkML object."""
        self.log.info('Map basic cassette info to Alliance object.')
        for cass in self.fb_data_entities.values():
            agr_cass = self.agr_export_type()
            agr_cass.obsolete = cass.chado_obj.is_obsolete
            agr_cass.primary_external_id = f'FB:{cass.uniquename}'
            cass.linkmldto = agr_cass

    def get_datatype_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        # self.get_entity_xrefs(session)
        # Text form From FTA-134 next 3 lines:-
        # feature_relationship, type 'has_reg_region', subject = FBal, object in (FBto, FBsf, FBgn) # from GA30e field
        # feature_relationship, type 'tagged_with', subject = FBal, object in (FBto, FBsf) # from GA30a field
        # feature_relationship, type 'carries_tool' subject = FBal, object in (FBto, FBsf) # from GA30b field

        self.get_entity_relationships(session, 'subject')
        # , rel_type)='has_reg_region',
        #                             entity_type='engineered_region', entity_regex=self.regex['tool'])

    # Elaborate on query_chado_and_export() for the CassetteHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the CassetteHandler."""
        super().query_chado_and_export(session)
        self.generate_export_dict(self.cassette_associations, 'cassette_association_ingest_set')

        return

    def map_cassette_associations(self):
        """Map transgenic cassette associations to Alliance object."""
        self.log.info('Map cassette associations to Alliance object.')
        OBJECT = 1
        SUBJECT = 0
        counter = 0

        cassette_cassette_counter = {}
        for cassette_cassette_key in self.cassette_cassette_rels.keys():
            if self.testing:
                self.log.debug(f'Mapping {cassette_cassette_key} to Alliance object. {self.cassette_cassette_rels[cassette_cassette_key]}')
            try:
                cassette_cassette_counter[cassette_cassette_key[OBJECT]] += 1
            except KeyError:
                cassette_cassette_counter[cassette_cassette_key[OBJECT]] = 1

        for cassette_cassette_key, cassette_cassette_rels in self.cassette_cassette_rels.items():
            object_feature_id = cassette_cassette_key[OBJECT]
            f_object = self.fb_data_entities[object_feature_id]
            object_curie = f'FB:{f_object.uniquename}'
            subject = self.feature_lookup[cassette_cassette_key[SUBJECT]]
            subject_curie = f'FB:{subject["uniquename"]}'
            first_feat_rel = cassette_cassette_rels[0]
            all_pub_ids = []
            for cassette_cassette_rel in cassette_cassette_rels:
                all_pub_ids.extend(cassette_cassette_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            # NOTE: pub 383755 | FlyBase Experimental Tool information Is the only one used
            # for cassettes. But not in lookup pub curies!
            # So pub_curies will be empty.
            pub_curies = self.lookup_pub_curies(all_pub_ids)

            # Adjust allele-gene relation_type as needed.

            rel_type_name = cassette_cassette_rels[0].chado_obj.type.name
            rel_dto = agr_datatypes.CassetteAssociationDTO(
                subject_curie, object_curie,
                pub_curies, False, rel_type_name)
            if f_object.is_obsolete is True or subject['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            first_feat_rel.linkmldto = rel_dto
            self.cassette_associations.append(first_feat_rel)
            counter += 1
        self.log.info(f'Generated {counter} cassette-cassette unique associations.')
        return

    def synthesize_cassette_associations(self):
        """Get cassette relationships."""
        self.log.info('Synthesize cassette.')
        sub_cassette_counter = 0
        obj_cassette_counter = 0
        for cassette in self.fb_data_entities.values():
            relevant_cassette_rels = cassette.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['has_reg_region', 'tagged_with', 'carries_tool'])
            if relevant_cassette_rels:
                sub_cassette_counter += 1
            for cassette_rel in relevant_cassette_rels:
                try:
                    cassette_cassette_key = (cassette_rel.chado_obj.object_id, cassette_rel.chado_obj.subject_id)
                except AttributeError:
                    self.log.error(f"problem {cassette} {cassette_rel}")
                    raise
                try:
                    self.cassette_cassette_rels[cassette_cassette_key].append(cassette_rel)
                except KeyError:
                    self.cassette_cassette_rels[cassette_cassette_key] = [cassette_rel]
                    obj_cassette_counter += 1
        self.log.info(f'Found {obj_cassette_counter} cassettes for {sub_cassette_counter} cassettes.')
        return

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_cassette_associations()
