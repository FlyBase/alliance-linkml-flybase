"""Module:: cassette_handler.

Synopsis:
    A data handler that exports FlyBase data for cassette to Alliance
    TransgenicTool LinkML objects.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

from logging import Logger
import copy
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
        'FBal0198528': 'Pi4KIIalpha[GD9857]',   # Transegnic, superceded by four FBti (superseded_by_transgnc_insertions).
        'FBal0322755': 'Mcm3[+tBa]',                   # cassette main type
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
        'FBal0083005': 'dlg1[DeltaSH3.UAS.Tag:FLAG]',  # two refs for same tagged_with (FBrf0099758, FBrf0130114),
                                                       # only one for has_reg_region (tool) (FBrf0099758)
        'FBal0137284': 'cic[Tag:HA]',  # single tagged_with (FBrf0144844, FBrf0180201), single has_reg_region (gene) (FBrf0144844, FBrf0180201)
        'FBal0144698': 'PGRP-LE[UAS.Tag:FLAG]',   # single tagged_with (FBrf0152317, FBrf0212747),
                                                  # two has_reg_region (tool) (UAS = FBrf0212747, UASt = FBrf0152317)
        'FBal0151333': 'wg[PE4.UAS.cCa.Tag:HA]',  # single also_carries (FBrf0173223), single tagged_with (FBrf0167661, FBrf0173223),
                                                  # two has_reg_region (tool) (UAS = FBrf0173223, UASt = FBrf0167661)
        'FBal0137561': r'Crei\I-CreI[hs.PR]',  #
        'FBal0404843': r'Hsap\CGA[UAS.cLa]',  #
        'FBal0401141': r'Zzzz\VHH[deGradFP.UAS]',  #
        'FBal0051685': r'csw[CS.hs.2sev]',  # Has 2 props
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
        self.build_cvterm_lookup(session)
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
            # sanity check of making sure ALL test data is in the entries.
            unique_names = {}
            # generate dict of uniquename in fb_data_entities
            for entity_id in self.fb_data_entities:
                unique_names[self.fb_data_entities[entity_id].uniquename] = entity_id
            # Are all test examples in the fb_data_entities
            for name, _ in self.test_set.items():
                if name not in unique_names:
                    self.log.error(f"Missing {name} in fb_data_entities. Must be a construct, hopefully.")
            # Are all fb_data_entities in the test set
            for entity_id, entity in self.fb_data_entities.items():
                if entity.uniquename not in self.test_set:
                    self.log.error(f"Missing {entity.uniquename} in test set, could have extras?")

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
        self.get_entity_cvterms(session)
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
            # moment it's not actually changing the type !

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

    def get_comp_type_curies(self, fb_data_entity):
        """Get component_type_curies."""
        component_type_curies = []
        data_key = 'transgenic_product_class'
        if data_key in fb_data_entity.prop_data.keys():
            for bob in fb_data_entity.prop_data[data_key]:
                component_type_curies.append(f"{data_key} {bob['name']}: {bob['type']}:{bob['accession']}")
        return component_type_curies

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

    def lookup_expresses_pub_curies(self, entity):
        """Lookup pub curies for those that uses expresses."""
    #  A. 'expresses' association where is NO encodes_tool feature_relationship
    #    (i.e. the bit that starts "if entity.uniquename not in encoded.keys()"
    #     use this logic in this order.
    #
    #
    # 1. if there are component_type_curies for the cassette FBal,
    #    use the pubs associated with the component_type_curies
    #
    #  1a. in this loop, do an additional check
    #
    #     - if the number of component_type_curies > 1 AND
    #       the total number of pubs associated with all of the component_type_curies > 1
    #       leave the pub_curies blank and write a log message (**see note below)
    #
    #     - otherwise use all the pubs associated with the component_type_curies
    #
    #       **note: There are a handful cases (~5) that have multiple component_type_curies with multiple refs
    #               (when you add up all the refs associated with each component_type_curie).
    #
    #     An example is FBal0051685, csw[CS.hs.2sev]
    #
    #     This has 2 component_type_curies, each with a different reference
    #     (each reference only provides the info to add one of the component_type_curies,
    #      so I can't fix the curation to make the reference the same for both !):
    #
    #     dominant_negative_variant [FBrf0190302]
    #     missense_variant [FBrf0088202]
    #
    #     Since we're combining both component_type_curies into a single 'expresses' association in Alliance,
    #     it would be incorrect to say that both references give evidence for the single 'FBal0051685 (csw[CS.hs.2sev])
    #     expresses FBgn0000382 (csw), dominant_negative_variant, missense_variant' association.
    #
    #     I think its best we just leave the pub_curies blank for these in the linkml code and write a log message
    #     to say which FBal cassette has this issue, rather than trying to do some complex algorithm for this tiny
    #     number of cases (which would involve making two separate 'expresses' annotations
    #      - one for each component_type_curie - to be able to partition the references correctly).
    #        When we switch to using the Alliance as source of truth, curators can use the log message to know which
    #        cassettes have this problem and can fix in the Alliance by making two expresses associations.
    #
    # 2. otherwise, lookup whether the cassette allele has 'molecular_info' featureprop(s) - if it does,
    #    use all the pubs associated with all the 'molecular_info' featureprops to fill in the pub_curies
    #
    #
    # 3. otherwise, if there is a *single reference* for the cassette allele-to-parent gene relationship,
    #    use that as the pub_curies
    #
    # **Note: Its important to do this one as a last resort as its less accurate than using either 1. or 2. above,
    #         because when an allele is merged we often use an FB analysis ref instead of the ref that originally
    #         described the cassette allele, and this FB analysis ref ends up being the single reference in the
    #         cassette allele-to-parent gene relationship.
    #
    # 4. otherwise, leave pub_curies empty (there are currently only ~30 cases that end up here,
    #    because they have nothing in the 'molecular_info' free text note,
    #    so I'll try and fix that by adding the relevant info into chado).
        pubs = set()
        data_key = 'transgenic_product_class'
        if data_key in entity.prop_data.keys():
            for bob in entity.prop_data[data_key]:
                print(f"BOB:{entity.uniquename} {bob}")
                pubs.add(bob['pub'])
        print(f"BOBBY LEPC {entity.uniquename} pubs:{pubs} len:{len(pubs)}")
        pub_curies = list(pubs)
        if len(pub_curies) == 1:  # 1
            print(f"BOBBY: {entity.uniquename} 1 ref")
        elif len(pub_curies) > 1:  # 1 a
            print(f"BOBBY: {entity.uniquename} Multiple comp curie with dif refs {pub_curies}")
            self.log.warning(f"{entity.uniquename} has multiple comp curie with diff refs {pub_curies}")
            pub_curies = []
        elif 'molecular_info' in entity.props_by_type.keys():
            for bob in entity.props_by_type['molecular_info']:
                print(f"BOBBY1: {entity.uniquename} {type(bob)} {dir(bob)}")
            for prop_type, prop_list in entity.props_by_type.items():
                print(f"BOBBY2:  {entity.uniquename} {prop_type}: {len(prop_list)} props")
                for prop in prop_list:
                    print(f"BOBBY3   {entity.uniquename}  - {prop.chado_obj.value[:50] if prop.chado_obj.value else 'None'}...")
        # elif entity.has_molecular_info():  # 2
        #     pass
        # else:  # 3
        #    print(f"BOBBY: {entity.uniquename} 0 refs")
        print(f"BOBBY: Returning {entity.uniquename} {pub_curies}")
        return pub_curies

    def map_cassette_associations(self):
        """Map transgenic cassette-component associations to Alliance object."""
        self.log.info('Map cassette-component associations to Alliance object.')
        CASSETTE = 0
        COMPONENT = 1
        counter = 0

        map_relationship = {'has_reg_region': 'is_regulated_by',
                            'tagged_with': 'tagged_with',
                            'carries_tool': 'contains',
                            'encodes_tool': 'expresses'}
        # cassette_cassette_counter = {}
        for cassette_cassette_key in self.cassette_cassette_rels.keys():
            if self.testing:
                self.log.debug(f'Mapping {cassette_cassette_key} to Alliance object. {self.cassette_cassette_rels[cassette_cassette_key]}')

        bad_relationship_count = {}
        encoded = {}  # dict to store if cassette assoc mapped by encodes_tool
        # go through cassettes and make the cassette-component associations.
        for cassette_cassette_key, cassette_cassette_rels in self.cassette_cassette_rels.items():
            cassette_feature_id = cassette_cassette_key[CASSETTE]
            cassette = self.fb_data_entities[cassette_feature_id]
            cassette_curie = f'FB:{cassette.uniquename}'
            component = self.feature_lookup[cassette_cassette_key[COMPONENT]]
            component_curie = f'FB:{component["uniquename"]}'
            assoc_type = self.cassette_dto_type(component)
            first_feat_rel = cassette_cassette_rels[0]
            all_pub_ids = []
            for cassette_cassette_rel in cassette_cassette_rels:
                all_pub_ids.extend(cassette_cassette_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            pub_curies = self.lookup_pub_curies(all_pub_ids)

            # Adjust cassette-component relation_type as needed.
            rel_type_name = cassette_cassette_rels[0].chado_obj.type.name
            if rel_type_name in map_relationship:
                rel_type_name = map_relationship[rel_type_name]
            else:
                self.log.error(f"Unknown relationship type {rel_type_name}")
                continue
            component_type_curies = []
            if rel_type_name == 'expresses':
                encoded[cassette.uniquename] = 1
                if self.testing:
                    # Cvtermprop type (name) keyed lists of entity_cvterm_ids.
                    for bob in cassette.prop_data.keys():
                        self.log.debug(f"BOBBY: prop_data {cassette.uniquename} {bob} {component_type_curies}")
                component_type_curies = self.get_comp_type_curies(cassette)
            if self.testing:
                self.log.debug(f"BOBBY: comp cur {component_type_curies}")
                self.log.debug(f"\tBOBBY: assoc type->{assoc_type} cass name -> {cassette.uniquename} ctc -> {component_type_curies}")
            if rel_type_name not in bad_relationship_count:
                bad_relationship_count[rel_type_name] = 0
            bad_relationship_count[rel_type_name] += 1
            if assoc_type == 'component_free_text':
                # CassetteComponentSlotAnnotationDTO
                if self.testing:
                    mess = "map_cassette_associations: ComponentSlotAnnotation cass:"
                    mess += f"{cassette_curie} comp:{component_curie} '{rel_type_name}'"
                    self.log.debug(mess)
                symbol = component['symbol']
                organism_id = component['organism_id']
                taxon_text = self.organism_lookup[organism_id]['full_species_name']
                taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
                rel_dto = agr_datatypes.CassetteComponentSlotAnnotationDTO(
                    rel_type_name, symbol, taxon_curie,
                    taxon_text, pub_curies).dict_export()
                cassette.linkmldto.cassette_component_dtos.append(rel_dto)
            elif assoc_type == 'tool_association':
                # CassetteTransgenicToolAssociationDTO
                if self.testing:
                    mess = "map_cassette_associations: TransgenicToolAssociation cass:"
                    mess += f"{cassette_curie} comp:{component_curie} '{rel_type_name}'"
                    self.log.debug(mess)
                rel_dto = agr_datatypes.CassetteTransgenicToolAssociationDTO(
                    cassette_curie, component_curie,
                    pub_curies, False, rel_type_name)
                first_feat_rel.linkmldto = rel_dto
                self.cassette_tool_associations.append(first_feat_rel)
            elif assoc_type == 'genomic_entity_association':
                # CassetteGenomicEntityAssociationDTO
                if self.testing:
                    mess = "map_cassette_associations: GenomicEntityAssociation cass:"
                    mess += f"{cassette_curie} comp:{component_curie} '{rel_type_name}'"
                    mess += f"{'|'.join(component_type_curies)}"
                    self.log.debug(mess)
                rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                    cassette_curie, component_curie,
                    pub_curies, False, rel_type_name,
                    component_type_curies)
                first_feat_rel.linkmldto = rel_dto
                self.cassette_genomic_entity_associations.append(first_feat_rel)
            if cassette.is_obsolete is True or component['is_obsolete'] is True:
                self.log.error(f"{cassette_curie} {component_curie} should never be obsolete??")
            counter += 1
        for key, count in bad_relationship_count.items():
            self.log.error(f'Bad relationship count for {key}: {count}')
        self.log.info(f'Generated {counter} cassette-component unique associations.')

        for entity in self.fb_data_entities.values():
            rels = entity.recall_relationships(
                self.log,
                entity_role='subject',  # 'subject' or 'object'
                rel_types='alleleof',  # str or list of relationship type names
                rel_entity_types='gene'  # (features only) filter by related entity type
            )
            for rel in rels:
                if self.testing:  # just to see more examples look at method always.
                    pub_curies = self.lookup_expresses_pub_curies(entity)
                    print(f"BOBBY TEST: NAME: {entity.uniquename} pub:{pub_curies}")
                if entity.uniquename not in encoded.keys():
                    component_type_curies = self.get_comp_type_curies(entity)
                    pub_curies = self.lookup_expresses_pub_curies(entity)
                    print(f"BOBBY: NAME: {entity.uniquename} pub:{pub_curies}")
                    if self.testing:
                        self.log.debug(f"{entity.uniquename} has parent {rel.chado_obj.object.uniquename}")
                    gene = self.feature_lookup[rel.chado_obj.object.feature_id]
                    assoc_type = self.cassette_dto_type(gene)
                    self.log.debug(f"{entity.uniquename} {gene} has {assoc_type} association")
                    # Always a gene currently BUT might in future have
                    # subset of foreign genes so check now anyway
                    if assoc_type == 'component_free_text':
                        mess = f"cassette {entity.uniquename} has parent {gene.uniquename} "
                        mess += f"BUT assoc_type is {assoc_type} So problem"
                        self.log.error(mess)
                    elif assoc_type == 'genomic_entity_association':
                        # CassetteGenomicEntityAssociationDTO
                        if self.testing:
                            mess = f"map_cassette_associations: GenomicEntityAssociation rel:{rel} cass:"
                            mess += (f"{entity.uniquename} comp:{rel.chado_obj.object.uniquename}"
                                     f" 'expresses' {component_type_curies} ")
                            self.log.debug(mess)
                        rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                            f"FB:{entity.uniquename}",
                            f"FB:{rel.chado_obj.object.uniquename}",
                            pub_curies, False, 'expresses',
                            component_type_curies)  # NEED to add pub_curies still
                        rel.linkmldto = rel_dto
                        self.cassette_genomic_entity_associations.append(rel)
                save_target = False
                for trans in entity.prop_data['transgenic_product_class']:
                    if trans['name'] in ('RNAi_reagent', 'sgRNA', 'antisense'):
                        save_target = True
                if save_target:
                    # Because the relationship is used for both expresses and targets
                    # we want to copy that and not overwrite it.
                    new_rel = copy.copy(rel)  # Create independent copy
                    # CassetteGenomicEntityAssociationDTO
                    if self.testing:
                        mess = "map_cassette_associations: GenomicEntityAssociation "
                        mess += f"rel:{new_rel} cass:{entity.uniquename} comp:{new_rel.chado_obj.object.uniquename} 'targets'"
                        self.log.debug(mess)
                    rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                        f"FB:{entity.uniquename}",
                        f"FB:{new_rel.chado_obj.object.uniquename}",
                        ["NEEDED"], False, 'targets')  # NEED to add pub_curies still
                    new_rel.linkmldto = rel_dto
                    self.cassette_genomic_entity_associations.append(new_rel)
        return

    def synthesize_cassette_associations(self):
        """Get cassette-to-component relationships."""
        self.log.info('Synthesize cassette-to-component relationships.')
        cassette_counter = 0
        component_counter = 0
        for cassette in self.fb_data_entities.values():
            relevant_cassette_rels = cassette.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['has_reg_region', 'tagged_with', 'carries_tool', 'encodes_tool'])
            if relevant_cassette_rels:
                cassette_counter += 1
            # put the data into cassette_cassette_key with the cassette (FBal) first and the component second
            # cassette_rel.chado_obj.subject_id is the feature_of the cassette (chado subject_id)
            # cassette_rel.chado_obj.object_id is the feature_id of the component (chado object_id)
            for cassette_rel in relevant_cassette_rels:
                try:
                    cassette_cassette_key = (cassette_rel.chado_obj.subject_id, cassette_rel.chado_obj.object_id)
                except AttributeError:
                    self.log.error(f"problem {cassette} {cassette_rel}")
                    raise
                try:
                    self.cassette_cassette_rels[cassette_cassette_key].append(cassette_rel)
                except KeyError:
                    self.cassette_cassette_rels[cassette_cassette_key] = [cassette_rel]
                    component_counter += 1
        self.log.info(f'Found {component_counter} components for {cassette_counter} cassettes.')

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_cassette_associations()
