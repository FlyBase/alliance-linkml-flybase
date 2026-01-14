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

    # NOTE: We have general alleles in here too so we can check we only get the casssettes here.
    #       Also Cassettes are in the Allele code to check we only get these once.
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
        'FBal0322755': 'Mcm3[+tBa]',                 # cassette main type
        'FBal0322754': 'flfl[DeltaRanBD.UAS.Venus]',
        'FBal0296109': 'sSemp1[R41G.UAS]',
        'FBal0193766': 'Gr63a[UAS.cJa]',
        'FBal0239883': 'sd[RNAi.N.UAS]',
        'FBal0000531': 'Amy-p[IX]',
        'FBal0028742': 'Act88F[E334K]',
        'FBal0212171': r'Avic\GFP[UAS.FRT1]',       # in vitro only
        'FBal0290956': 'Csas[21]',                  # curator error: in vitro only in fb_2025_05, fixed for fb_2026_01
                                                    # (now is classical allele reported_as_itself)
        'FBal0392043': r'Avic\GFP[EYFP.3xP3.cUa]',  # in vitro only
        'FBal0028610': 'w[+mC]',                    # Has a secondary identifier to test
        'FBal0045138': 'Sry-delta[SDL1.lacZ]',                 # linked to an FBsf so a str association.
        'FBal0193109': r'Avic\GFP[EGFP.rho.PE.Tag:NLS(tra)]',  # linked to an FBsf so a str association.
        'FBal0250846': r'Scer\GAL4[GMR24E03]',                 # linked to an FBsf so a str association.
        'FBal0041313': r'Ecol\lacZ[eve.1.55] ',                # linked to an FBsf so a str association.
        'dlg1[DeltaSH3.UAS.Tag:FLAG]': 'FBal0083005',  # two refs for same tagged_with (FBrf0099758, FBrf0130114),
                                                       # only one for has_reg_region (tool) (FBrf0099758)
        'cic[Tag:HA]': 'FBal0137284',  # single tagged_with (FBrf0144844, FBrf0180201), single has_reg_region (gene) (FBrf0144844, FBrf0180201)
        'PGRP-LE[UAS.Tag:FLAG]': 'FBal0144698',   # single tagged_with (FBrf0152317, FBrf0212747),
                                                  # two has_reg_region (tool) (UAS = FBrf0212747, UASt = FBrf0152317)
        'wg[PE4.UAS.cCa.Tag:HA]': 'FBal0151333',  # single also_carries (FBrf0173223), single tagged_with (FBrf0167661, FBrf0173223),
                                                  # two has_reg_region (tool) (UAS = FBrf0173223, UASt = FBrf0167661)
    }

    cassette_prop_to_note_mapping = {
        'aminoacid_rep': ('comment', 'note_dtos'),
        'molecular_info': ('comment', 'note_dtos'),
        'nucleotide_sub': ('comment', 'note_dtos'),
        # At the moment, just for code development. (line below)
        # 'internal_notes': ('internal_note', 'note_dtos'),
    }
    cassette_associations = []  # Should delete this one later
    # cassette_component_free_text_associations = []
    cassette_tool_associations = []
    cassette_genomic_entity_associations = []
    cassette_cassette_rels = {}

    def get_general_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_organism_lookup(session)
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
                self.log.debug(f"Cassette entities: {pkey_id}: {result}")
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

    def map_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for a FlyBase entity."""
        self.log.info('Map secondary IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.linkmldto is None:
                continue
            secondary_id_dtos = []
            for secondary_id in fb_data_entity.alt_fb_ids:
                secondary_id_dtos.append(secondary_id)
            sec_id_list = getattr(fb_data_entity.linkmldto, slot_name)
            sec_id_list.extend(secondary_id_dtos)
        return

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
        self.get_entity_relationships(session, 'subject')

    def cassette_dto_type(self, feature):
        """Derive association type from the feature."""
        assoc_type = 'component_free_text'
        # logic to decide which type of Alliance DTO object to use
        if feature["uniquename"].startswith('FBto'):  # cassette component is a TransgenicTool (FBid is a FBto):
            assoc_type = 'tool_association'

        elif feature["uniquename"].startswith('FBsf'):  # cassette component is a seqfeat (FBid is a FBsf):
            assoc_type = 'component_free_text'  # for now, will change to a CassetteGenomicEntityAssociationDTO
            # once we start to submit FBsf features, so keep this loop in place for then even though at the
            # moment its not actually changing the type !

        elif feature["uniquename"].startswith('FBgn'):  # cassette component is a gene (FBid is a FBgn):
            assoc_type = 'genomic_entity_association'
            # if the gene is from another Alliance MOD species AND
            # there is an unambigous mapping of the FBgn to a single Alliance MOD curie:
            #    unless the gene is typically_used_as_tool:
            #    type = 'genomic_entity_association'
            #    # replace the original obj_curie with the single Alliance MOD curie to use
            #      when the genomic_entity_association is made later
            #    obj_curie = the single Alliance MOD curie
            # else:
            #    if the gene is being exported to the Alliance:
            #        unless the gene is being submitted as 'internal':
            #            type = 'genomic_entity_association'
        return assoc_type

    # Elaborate on query_chado_and_export() for the CassetteHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the CassetteHandler."""
        super().query_chado_and_export(session)
        # self.generate_export_dict(self.cassette_component_free_text_associations,
        #                           'cassette_str_association_ingest_set')
        self.generate_export_dict(self.cassette_genomic_entity_associations,
                                  'cassette_genomic_entity_association_ingest_set')
        self.generate_export_dict(self.cassette_tool_associations,
                                  'cassette_transgenic_tool_association_ingest_set')

        return

    def map_cassette_associations(self):
        """Map transgenic cassette associations to Alliance object."""
        self.log.info('Map cassette associations to Alliance object.')
        OBJECT = 1
        SUBJECT = 0
        counter = 0

        map_relationship = {'has_reg_region': 'is_regulated_by',
                            'tagged_with': 'tagged_with',
                            'carries_tool': 'contains'}
        cassette_cassette_counter = {}
        for cassette_cassette_key in self.cassette_cassette_rels.keys():
            if self.testing:
                self.log.debug(f'Mapping {cassette_cassette_key} to Alliance object. {self.cassette_cassette_rels[cassette_cassette_key]}')
            try:
                cassette_cassette_counter[cassette_cassette_key[OBJECT]] += 1
            except KeyError:
                cassette_cassette_counter[cassette_cassette_key[OBJECT]] = 1

        bad_relationship_count = {}
        for cassette_cassette_key, cassette_cassette_rels in self.cassette_cassette_rels.items():
            object_feature_id = cassette_cassette_key[OBJECT]
            f_object = self.fb_data_entities[object_feature_id]
            component_curie = f'FB:{f_object.uniquename}'
            subject = self.feature_lookup[cassette_cassette_key[SUBJECT]]
            cassette_curie = f'FB:{subject["uniquename"]}'

            first_feat_rel = cassette_cassette_rels[0]
            all_pub_ids = []
            for cassette_cassette_rel in cassette_cassette_rels:
                all_pub_ids.extend(cassette_cassette_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            pub_curies = self.lookup_pub_curies(all_pub_ids)

            # Adjust allele-gene relation_type as needed.

            rel_type_name = cassette_cassette_rels[0].chado_obj.type.name
            if rel_type_name in map_relationship:
                rel_type_name = map_relationship[rel_type_name]
            else:
                if rel_type_name not in bad_relationship_count:
                    bad_relationship_count[rel_type_name] = 0
                bad_relationship_count[rel_type_name] += 1
                continue
            assoc_type = self.cassette_dto_type(subject)
            if assoc_type == 'component_free_text':
                # CassetteComponentSlotAnnotationDTO
                if self.testing:
                    print(f"map_cassette_associations: comp:{component_curie} cass:{cassette_curie}")
                # feature = self.feature_lookup[object_feature_id]
                symbol = subject['symbol']
                organism_id = subject['organism_id']
                # pubs = self.lookup_pub_curies(pub_ids)
                taxon_text = self.organism_lookup[organism_id]['full_species_name']
                taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
                rel_dto = agr_datatypes.CassetteComponentSlotAnnotationDTO(
                    rel_type_name, symbol, taxon_curie,
                    taxon_text, pub_curies).dict_export()
                # first_feat_rel.linkmldto = rel_dto
                f_object.linkmldto.cassette_component_dtos.append(rel_dto)
            elif assoc_type == 'tool_association':
                # CassetteTransgenicToolAssociationDTO
                rel_dto = agr_datatypes.CassetteTransgenicToolAssociationDTO(
                    cassette_curie, component_curie,
                    pub_curies, False, rel_type_name)
                first_feat_rel.linkmldto = rel_dto
                self.cassette_tool_associations.append(first_feat_rel)
            elif assoc_type == 'genomic_entity_association':
                # CassetteGenomicEntityAssociationDTO
                rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                    cassette_curie, component_curie,
                    pub_curies, False, rel_type_name)
                first_feat_rel.linkmldto = rel_dto
                self.cassette_genomic_entity_associations.append(first_feat_rel)
            if self.testing:
                self.log.debug(f"{cassette_curie} {component_curie} assoc type is {assoc_type}")
            if f_object.is_obsolete is True or subject['is_obsolete'] is True:
                self.log.error(f"{cassette_curie} {component_curie} should never be obsolete??")
            counter += 1
        for key in bad_relationship_count:
            self.log.error(f'Bad relationship count for {key}: {bad_relationship_count[key]}')
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
        self.log.info(f'Found {obj_cassette_counter} components for {sub_cassette_counter} cassettes.')
        return

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_cassette_associations()
