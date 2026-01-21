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
from sqlalchemy.orm import aliased
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler
from harvdev_utils.reporting import (
    Cvterm, Dbxref, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureCvtermPub, FeatureRelationship, FeatureRelationshipPub
)


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
        'FBal0083005': 'dlg1[DeltaSH3.UAS.Tag:FLAG]',  # two refs for same tagged_with (FBrf0099758, FBrf0130114),
                                                       # only one for has_reg_region (tool) (FBrf0099758)
        'FBal0137284': 'cic[Tag:HA]',  # single tagged_with (FBrf0144844, FBrf0180201), single has_reg_region (gene) (FBrf0144844, FBrf0180201)
        'FBal0144698': 'PGRP-LE[UAS.Tag:FLAG]',   # single tagged_with (FBrf0152317, FBrf0212747),
                                                  # two has_reg_region (tool) (UAS = FBrf0212747, UASt = FBrf0152317)
        'FBal0151333': 'wg[PE4.UAS.cCa.Tag:HA]',  # single also_carries (FBrf0173223), single tagged_with (FBrf0167661, FBrf0173223),
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
        self.build_allele_gene_lookup(session)
        self.build_cvterm_lookup(session)

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
        # FTA-137: Get encodes_tool and transgenic_product_class data for cassettes.
        self.get_cassette_encodes_tool_rels(session)
        self.get_cassette_transgenic_product_class(session)
        self.get_cassette_parent_genes(session)

    # FTA-137: Methods for retrieving encodes_tool and transgenic_product_class data.
    def get_cassette_encodes_tool_rels(self, session):
        """Get encodes_tool relationships for cassettes (cassette is the allele subject)."""
        self.log.info('Get encodes_tool relationships for cassettes.')
        component = aliased(Feature, name='component')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['allele']),
            component.is_obsolete.is_(False),
            component.uniquename.op('~')(self.regex['fb_uniquename']),
            Cvterm.name == 'encodes_tool',
        )
        results = session.query(FeatureRelationship).\
            select_from(Feature).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(component, (component.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        # Build pub lookup for relationship results.
        pub_results = session.query(FeatureRelationshipPub).\
            select_from(Feature).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == Feature.feature_id)).\
            join(FeatureRelationshipPub, (FeatureRelationshipPub.feature_relationship_id == FeatureRelationship.feature_relationship_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(component, (component.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        rel_pub_dict = {}
        for pub_result in pub_results:
            try:
                rel_pub_dict[pub_result.feature_relationship_id].append(pub_result.pub_id)
            except KeyError:
                rel_pub_dict[pub_result.feature_relationship_id] = [pub_result.pub_id]
        # Create cassette feature_id-keyed dict of encodes_tool FBRelationship objects.
        cassette_tool_dict = {}
        counter = 0
        for result in results:
            fb_rel = fb_datatypes.FBRelationship(result, 'feature_relationship')
            if result.feature_relationship_id in rel_pub_dict.keys():
                fb_rel.pubs = rel_pub_dict[result.feature_relationship_id]
            try:
                cassette_tool_dict[result.subject_id].append(fb_rel)
                counter += 1
            except KeyError:
                cassette_tool_dict[result.subject_id] = [fb_rel]
                counter += 1
        self.log.info(f'Found {counter} cassette-to-component "encodes_tool" relationships.')
        # Assign to cassette entities.
        cassette_counter = 0
        for cassette in self.fb_data_entities.values():
            cassette_id = cassette.chado_obj.feature_id
            if cassette_id in cassette_tool_dict.keys():
                cassette.encodes_tool_rels.extend(cassette_tool_dict[cassette_id])
                cassette_counter += 1
        self.log.info(f'Assigned encodes_tool relationships to {cassette_counter} cassettes.')
        return

    def get_cassette_transgenic_product_class(self, session):
        """Get transgenic_product_class (GA35) SO terms and their pubs for cassettes."""
        self.log.info('Get transgenic_product_class (GA35) SO terms for cassettes.')
        cvterm = aliased(Cvterm, name='cvterm')
        qualifier = aliased(Cvterm, name='qualifier')
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            cvterm.is_obsolete == 0,
            qualifier.name == 'transgenic_product_class',
        )
        # Main query for GA35 terms.
        results = session.query(Feature, FeatureCvterm, cvterm, Dbxref).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(cvterm, (cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(Dbxref, (Dbxref.dbxref_id == cvterm.dbxref_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(qualifier, (qualifier.cvterm_id == FeatureCvtermprop.type_id)).\
            filter(*filters).\
            distinct()
        # Query for pubs associated with each FeatureCvterm.
        pub_results = session.query(FeatureCvtermPub).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(cvterm, (cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(qualifier, (qualifier.cvterm_id == FeatureCvtermprop.type_id)).\
            join(FeatureCvtermPub, (FeatureCvtermPub.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            filter(*filters).\
            distinct()
        # Build feature_cvterm_id -> pub_ids lookup.
        fc_pub_dict = {}
        for pub_result in pub_results:
            try:
                fc_pub_dict[pub_result.feature_cvterm_id].append(pub_result.pub_id)
            except KeyError:
                fc_pub_dict[pub_result.feature_cvterm_id] = [pub_result.pub_id]
        # Build cassette feature_id -> {term_name: [pub_ids]} and {term_name: curie} lookups.
        cassette_ga35_dict = {}      # feature_id -> {term_name: [pub_ids]}
        cassette_ga35_curie_dict = {}  # feature_id -> {term_name: curie}
        counter = 0
        for result in results:
            feature_id = result.Feature.feature_id
            term_name = result.cvterm.name
            fc_id = result.FeatureCvterm.feature_cvterm_id
            # Build SO curie from dbxref accession (e.g., "SO:0000001").
            so_curie = f'SO:{result.Dbxref.accession}'
            pub_ids = fc_pub_dict.get(fc_id, [])
            if feature_id not in cassette_ga35_dict:
                cassette_ga35_dict[feature_id] = {}
                cassette_ga35_curie_dict[feature_id] = {}
            if term_name not in cassette_ga35_dict[feature_id]:
                cassette_ga35_dict[feature_id][term_name] = pub_ids
                cassette_ga35_curie_dict[feature_id][term_name] = so_curie
                counter += 1
            else:
                # Extend pub list if term already exists (shouldn't happen normally).
                cassette_ga35_dict[feature_id][term_name].extend(pub_ids)
        self.log.info(f'Found {counter} transgenic_product_class (GA35) term annotations.')
        # Assign to cassette entities.
        cassette_counter = 0
        for cassette in self.fb_data_entities.values():
            cassette_id = cassette.chado_obj.feature_id
            if cassette_id in cassette_ga35_dict.keys():
                cassette.transgenic_product_classes = cassette_ga35_dict[cassette_id]
                cassette.transgenic_product_class_curies = cassette_ga35_curie_dict[cassette_id]
                cassette_counter += 1
        self.log.info(f'Assigned transgenic_product_class (GA35) data to {cassette_counter} cassettes.')
        return

    def get_cassette_parent_genes(self, session):
        """Get parent gene for each cassette (via alleleof relationship)."""
        self.log.info('Get parent gene for each cassette.')
        cassette_counter = 0
        for cassette in self.fb_data_entities.values():
            cassette_id = cassette.chado_obj.feature_id
            if cassette_id in self.allele_gene_lookup.keys():
                cassette.parent_gene_id = self.allele_gene_lookup[cassette_id]
                cassette_counter += 1
        self.log.info(f'Assigned parent gene to {cassette_counter} cassettes.')
        return

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
        """Map transgenic cassette-component associations to Alliance object."""
        self.log.info('Map cassette-component associations to Alliance object.')
        CASSETTE = 0
        COMPONENT = 1
        counter = 0

        map_relationship = {'has_reg_region': 'is_regulated_by',
                            'tagged_with': 'tagged_with',
                            'carries_tool': 'contains'}
        # cassette_cassette_counter = {}
        for cassette_cassette_key in self.cassette_cassette_rels.keys():
            if self.testing:
                self.log.debug(f'Mapping {cassette_cassette_key} to Alliance object. {self.cassette_cassette_rels[cassette_cassette_key]}')
            # don't think we need to worry about the count of cassettes to components
            # try:
            #    cassette_cassette_counter[cassette_cassette_key[CASSETTE]] += 1
            # except KeyError:
            #    cassette_cassette_counter[cassette_cassette_key[CASSETTE]] = 1

        bad_relationship_count = {}
        # go through cassettes and make the cassette-component associations.
        for cassette_cassette_key, cassette_cassette_rels in self.cassette_cassette_rels.items():
            cassette_feature_id = cassette_cassette_key[CASSETTE]
            cassette = self.fb_data_entities[cassette_feature_id]
            cassette_curie = f'FB:{cassette.uniquename}'
            component = self.feature_lookup[cassette_cassette_key[COMPONENT]]
            component_curie = f'FB:{component["uniquename"]}'

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
                if rel_type_name not in bad_relationship_count:
                    bad_relationship_count[rel_type_name] = 0
                bad_relationship_count[rel_type_name] += 1
                continue
            assoc_type = self.cassette_dto_type(component)
            if assoc_type == 'component_free_text':
                # CassetteComponentSlotAnnotationDTO
                if self.testing:
                    print(f"map_cassette_associations: cass:{cassette_curie} comp:{component_curie}")
                symbol = component['symbol']
                organism_id = component['organism_id']
                # pubs = self.lookup_pub_curies(pub_ids)
                taxon_text = self.organism_lookup[organism_id]['full_species_name']
                taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
                rel_dto = agr_datatypes.CassetteComponentSlotAnnotationDTO(
                    rel_type_name, symbol, taxon_curie,
                    taxon_text, pub_curies).dict_export()
                # first_feat_rel.linkmldto = rel_dto
                cassette.linkmldto.cassette_component_dtos.append(rel_dto)
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
            if cassette.is_obsolete is True or component['is_obsolete'] is True:
                self.log.error(f"{cassette_curie} {component_curie} should never be obsolete??")
            counter += 1
        for key in bad_relationship_count:
            self.log.error(f'Bad relationship count for {key}: {bad_relationship_count[key]}')
        self.log.info(f'Generated {counter} cassette-component unique associations.')
        return

    def synthesize_cassette_associations(self):
        """Get cassette-to-component relationships."""
        self.log.info('Synthesize cassette-to-component relationships.')
        cassette_counter = 0
        component_counter = 0
        for cassette in self.fb_data_entities.values():
            relevant_cassette_rels = cassette.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['has_reg_region', 'tagged_with', 'carries_tool'])
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
        return

    # FTA-137: Synthesis method for encodes_tool/transgenic_product_class (expresses/targets) associations.
    def synthesize_cassette_encodes_targets(self):
        """Synthesize expresses/targets associations for cassettes based on encodes_tool and GA35 data.

        Logic:
        1. Process encodes_tool relationships -> expresses associations to f_r object
        2. Process GA35 terms -> expresses/targets associations to parent gene

        For GA35 terms:
        - 'RNAi_reagent', 'sgRNA', 'antisense' -> targets relationship
        - Other terms -> expresses relationship (but only if no encodes_tool data for that term)

        Reference comparison is done per-GA35 term to determine if GA35 info can be
        folded into the encodes_tool association or needs a separate parent gene association.
        """
        self.log.info('Synthesize cassette encodes/targets associations (FTA-137).')
        TARGETING_CLASSES = {'RNAi_reagent', 'sgRNA', 'antisense'}

        expresses_counter = 0
        targets_counter = 0
        folded_counter = 0
        unfoldable_counter = 0
        parent_gene_expresses_counter = 0
        parent_gene_targets_counter = 0

        for cassette in self.fb_data_entities.values():
            ga35_classes = set(cassette.transgenic_product_classes.keys())
            has_encodes_tool = bool(cassette.encodes_tool_rels)
            has_ga35 = bool(ga35_classes)

            # === SECTION 1: Process encodes_tool relationships ===
            for encodes_rel in cassette.encodes_tool_rels:
                component_id = encodes_rel.chado_obj.object_id
                rel_pubs = set(encodes_rel.pubs)

                if not has_ga35:
                    # No GA35 -> simple expresses, empty component_type_curies
                    self._create_expresses_association(
                        cassette, component_id, list(rel_pubs),
                        component_type_curies=[])
                    expresses_counter += 1

                elif ga35_classes.intersection(TARGETING_CLASSES):
                    # GA35 contains targeting terms -> expresses with empty curies
                    # (GA35 will be used for targets in section 2)
                    self._create_expresses_association(
                        cassette, component_id, list(rel_pubs),
                        component_type_curies=[])
                    expresses_counter += 1

                else:
                    # GA35 is "something else" -> compare references PER TERM
                    foldable_so_curies = []
                    unfoldable_terms = []

                    for ga35_term, ga35_term_pubs in cassette.transgenic_product_classes.items():
                        if set(ga35_term_pubs) == rel_pubs:
                            # This term's refs match -> can fold into component_type_curies
                            so_curie = cassette.transgenic_product_class_curies.get(ga35_term)
                            if so_curie:
                                foldable_so_curies.append(so_curie)
                                folded_counter += 1
                        else:
                            # This term's refs don't match -> needs separate association
                            unfoldable_terms.append((ga35_term, ga35_term_pubs))
                            unfoldable_counter += 1

                    # Create expresses association to encodes_tool object with foldable curies
                    self._create_expresses_association(
                        cassette, component_id, list(rel_pubs),
                        component_type_curies=foldable_so_curies)
                    expresses_counter += 1

                    # Track unfoldable terms for Section 2
                    cassette.unfoldable_ga35_terms = unfoldable_terms
                    if unfoldable_terms:
                        self.log.info(
                            f"Cassette {cassette.uniquename}: {len(unfoldable_terms)} GA35 terms "
                            f"have mismatched refs - creating separate parent gene associations")

            # === SECTION 2: Process GA35 for parent gene associations ===
            if has_ga35 and cassette.parent_gene_id:
                gene_id = cassette.parent_gene_id

                if ga35_classes.intersection(TARGETING_CLASSES):
                    # Targeting terms -> always make targets association to parent gene
                    # Get the subset of GA35 terms that are targeting terms
                    targeting_terms = ga35_classes.intersection(TARGETING_CLASSES)
                    targeting_so_curies = []
                    targeting_pubs = set()
                    for term in targeting_terms:
                        so_curie = cassette.transgenic_product_class_curies.get(term)
                        if so_curie:
                            targeting_so_curies.append(so_curie)
                        targeting_pubs.update(cassette.transgenic_product_classes.get(term, []))

                    self._create_targets_association(
                        cassette, gene_id, list(targeting_pubs),
                        component_type_curies=targeting_so_curies)
                    targets_counter += 1
                    parent_gene_targets_counter += 1

                elif not has_encodes_tool:
                    # Non-targeting GA35, no encodes_tool -> expresses to parent gene
                    all_so_curies = list(cassette.transgenic_product_class_curies.values())
                    all_ga35_pubs = set()
                    for pubs in cassette.transgenic_product_classes.values():
                        all_ga35_pubs.update(pubs)

                    self._create_expresses_association(
                        cassette, gene_id, list(all_ga35_pubs),
                        component_type_curies=all_so_curies)
                    expresses_counter += 1
                    parent_gene_expresses_counter += 1

                else:
                    # Has encodes_tool - check for unfoldable terms that need separate associations
                    if hasattr(cassette, 'unfoldable_ga35_terms') and cassette.unfoldable_ga35_terms:
                        for ga35_term, ga35_term_pubs in cassette.unfoldable_ga35_terms:
                            so_curie = cassette.transgenic_product_class_curies.get(ga35_term)
                            so_curies = [so_curie] if so_curie else []
                            self._create_expresses_association(
                                cassette, gene_id, list(ga35_term_pubs),
                                component_type_curies=so_curies)
                            expresses_counter += 1
                            parent_gene_expresses_counter += 1
                            self.log.debug(
                                f"Cassette {cassette.uniquename}: created separate parent gene "
                                f"expresses association for GA35 term '{ga35_term}' with mismatched refs")
                    else:
                        # All GA35 terms were folded into encodes_tool association
                        self.log.debug(
                            f"Cassette {cassette.uniquename}: all GA35 terms folded into encodes_tool association")

        self.log.info(f'Created {expresses_counter} expresses associations.')
        self.log.info(f'Created {targets_counter} targets associations.')
        self.log.info(f'Folded {folded_counter} GA35 terms into encodes_tool associations.')
        self.log.info(f'Found {unfoldable_counter} GA35 terms with mismatched refs (separate associations needed).')
        self.log.info(f'Created {parent_gene_expresses_counter} parent gene expresses associations.')
        self.log.info(f'Created {parent_gene_targets_counter} parent gene targets associations.')
        return

    # FTA-137: Helper methods for creating associations.
    def _create_expresses_association(self, cassette, component_id, pub_ids, component_type_curies):
        """Create an expresses association for the cassette.

        Args:
            cassette: The FBCassette object.
            component_id: The feature_id of the component (encodes_tool object or parent gene).
            pub_ids: List of pub_ids for evidence.
            component_type_curies: List of SO term curies for component_type_curies slot.

        """
        self._create_association(cassette, component_id, pub_ids, 'expresses', component_type_curies)

    def _create_targets_association(self, cassette, component_id, pub_ids, component_type_curies):
        """Create a targets association for the cassette.

        Args:
            cassette: The FBCassette object.
            component_id: The feature_id of the component (parent gene).
            pub_ids: List of pub_ids for evidence.
            component_type_curies: List of SO term curies for component_type_curies slot.

        """
        self._create_association(cassette, component_id, pub_ids, 'targets', component_type_curies)

    def _create_association(self, cassette, component_id, pub_ids, relation_name, component_type_curies):
        """Create an association DTO and add it to the appropriate list.

        Args:
            cassette: The FBCassette object.
            component_id: The feature_id of the component.
            pub_ids: List of pub_ids for evidence.
            relation_name: Either 'expresses' or 'targets'.
            component_type_curies: List of SO term curies for component_type_curies slot.

        """
        # Check if component exists in feature_lookup
        if component_id not in self.feature_lookup:
            self.log.warning(f"Component {component_id} not found in feature_lookup for cassette {cassette.uniquename}")
            return

        component = self.feature_lookup[component_id]
        cassette_curie = f'FB:{cassette.uniquename}'
        component_curie = f'FB:{component["uniquename"]}'
        pub_curies = self.lookup_pub_curies(pub_ids)

        # Determine association type based on component type
        assoc_type = self.cassette_dto_type(component)

        if assoc_type == 'genomic_entity_association':
            # CassetteGenomicEntityAssociationDTO - can use component_type_curies
            rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                cassette_curie, component_curie,
                pub_curies, False, relation_name)
            if component_type_curies:
                rel_dto.component_type_curies = component_type_curies
            fb_rel = fb_datatypes.FBExportEntity()
            fb_rel.linkmldto = rel_dto
            self.cassette_genomic_entity_associations.append(fb_rel)

        elif assoc_type == 'tool_association':
            # CassetteTransgenicToolAssociationDTO - cannot use component_type_curies
            rel_dto = agr_datatypes.CassetteTransgenicToolAssociationDTO(
                cassette_curie, component_curie,
                pub_curies, False, relation_name)
            if component_type_curies:
                # Add as notes since component_type_curies not allowed on tool associations
                note_text = f"component_type: {', '.join(component_type_curies)}"
                note_dto = agr_datatypes.NoteDTO('comment', note_text, pub_curies)
                rel_dto.note_dtos = [note_dto.dict_export()]
                self.log.warning(
                    f"Cassette {cassette.uniquename}: component_type_curies added as note for "
                    f"tool association to {component_curie} (not allowed on CassetteTransgenicToolAssociationDTO)")
            fb_rel = fb_datatypes.FBExportEntity()
            fb_rel.linkmldto = rel_dto
            self.cassette_tool_associations.append(fb_rel)

        else:  # component_free_text
            # CassetteComponentSlotAnnotationDTO - inline in cassette
            symbol = component['symbol']
            organism_id = component['organism_id']
            taxon_text = self.organism_lookup[organism_id]['full_species_name']
            taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
            rel_dto = agr_datatypes.CassetteComponentSlotAnnotationDTO(
                relation_name, symbol, taxon_curie,
                taxon_text, pub_curies)
            if component_type_curies:
                # Add as notes since component_type_curies not allowed on free text annotations
                note_text = f"component_type: {', '.join(component_type_curies)}"
                note_dto = agr_datatypes.NoteDTO('comment', note_text, pub_curies)
                rel_dto.note_dtos = [note_dto.dict_export()]
                self.log.warning(
                    f"Cassette {cassette.uniquename}: component_type_curies added as note for "
                    f"free text association to {component_curie} (not allowed on CassetteComponentSlotAnnotationDTO)")
            cassette.linkmldto.cassette_component_dtos.append(rel_dto.dict_export())

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the CassetteHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_cassette_associations()
        # FTA-137: Synthesize encodes_tool/transgenic_product_class associations.
        self.synthesize_cassette_encodes_targets()
