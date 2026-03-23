"""Module:: construct_handler.

Synopsis:
    A data handler that exports FlyBase data for constructs, including their
    associations to other features, to Alliance Construct LinkML objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler
from harvdev_utils.reporting import (
    Cvterm, Feature, FeatureCvterm, FeatureCvtermprop,
    FeatureRelationship, FeatureRelationshipPub
)


class ConstructHandler(FeatureHandler):
    """This object gets, synthesizes and filters construct data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the ConstructHandler object."""
        super().__init__(log, testing)
        self.datatype = 'construct'
        self.fb_export_type = fb_datatypes.FBConstruct
        self.agr_export_type = agr_datatypes.ConstructDTO
        self.primary_export_set = 'construct_ingest_set'

    test_set = {
        'FBtp0008631': 'P{UAS-wg.H.T:HA1}',                       # Expresses FBgn wg, regulated by FBto UASt.
        'FBtp0010648': 'P{wg.FRT.B}',                             # Expresses FBgn wg, regulated by FBgn sev, has FBto FRT.
        'FBtp0145675': 'PBac{UAS-hHTT.ex1.Q97.S13D.mCherry}',     # Expresses FBgn Hsap\HTT, regulated by FBto UAS, tagged with FBto mCherry.
        'FBtp0000074': 'P{ftzG}',                                 # Expresses FBgn ftz, regulated by FBgn ftz.
        'FBtp0000326': 'P{SEV5}',                                 # Expresses FBgn sev, tagged with FBto MYC.
        'FBtp0161516': 'P{lush-GAL4.3}',                          # Expresses FBto GAL4, regulated by FBgn lush.
        'FBtp0057873': 'P{GMR16C10-GAL4}',                        # Expresses FBto GAL4, regulated by GMR16C10 (related to two genes, Brf and lute).
        'FBtp0032215': 'P{GD5007}',                               # Targets FBgn wg, regulated by FBto UASt.
        'FBtp0031452': 'P{GD4157}',                               # Targets FBgn lbe, regulated by FBto UASt.
        'FBtp0145396': 'P{TOE.GS00055}',                          # Targets FBgn wg, regulated by FBto UASt.
        'FBtp0145394': 'P{TKO.GS00469}',                          # Targets FBgn Alp9, Alp10, regulated by FBto UASt.
        'FBtp0000352': 'P{GawB}',                                 # Expresses FBto GAL4, FBgn Scer\GAL4. Report both?
        'FBtp0161256': 'PBac{UAS-G-CEPIA1::TM-2A-TagRFP::TM}',    # 2 FBal; expresses FBto G-CEPIA1, RFP; expresses FBgn Equa\eqFP578, GFP; regulated by UAS.
        'FBtp0051705': 'M{MtnBcDNA-MtnDcDNA.EGFP}',               # has_reg_region MtnB.
        'FBtp0080088': 'P{UAS-Brainbow}',                         # Expresses EBFP2, EGFP, mKO2, has_reg_region UAS; tagged_with HA, MYC, V5; carries lox.
        'FBtp0083738': 'P{GR}',                                   # Is regulated_by FBgn Act5C.
        'FBtp0017594': 'P{UAS(-FRT)ptc.Deltaloop2}',              # Obsolete, has only a non-current symbol synonym - for testing feature lookup.
        'FBtp0001701': 'P{hs-yCDC42.V12}',                        # Expresses Scer CDC42 (SGD:S000004219).
        'FBtp0131348': 'P{UAS-Scer_RCR1.MYC}',                    # Expresses Scer RCR1 (SGD:S000000209).
        'FBtp0001650': 'P{UAS-Cele_ced-3.S}',                     # Expresses Cele ced-3 (WB:WBGene00000417).
        'FBtp0010091': 'P{hs-Drer_nkx2.7.P}',                     # Expresses Drer nkx2.7 (ZFIN:ZDB-GENE-990415-179).
        'FBtp0007421': 'P{hb-Xlh1}',                              # Expresses Xlae h (no XB ID in FB).
        'FBtp0002652': 'P{UAS-mCD8::GFP.L}',                      # Expresses Mmus Cd8a (MGI:88346).
        'FBtp0001429': 'P{UAS-MAP2.A}',                           # Expresses Rnor Map2 (RGD:3044).
        'FBtp0000463': 'P{UAS-MAPT.A}',                           # Expresses Human MAPT (HGNC:6893).
        'FBtp0150381': 'PBac{UAS-SARS-CoV-2-nsp13.B}',            # Expresses SARS-CoV-2 nsp13 (REFSEQ:YP_009725308).
        'FBtp0132292': 'P{U6:2-scw.flySAM2.0}',                   # Exception with `has_transcriptional_unit` as maps to 'FBal0345196'
    }

    # Additional set for export added to the handler.
    construct_associations = []    # Will be a list of FBExportEntity objects (relationships), map to ConstructGenomicEntityAssociationDTO.
    construct_cassette_associations = []    # Will be a list of FBExportEntity objects, map to ConstructCassetteAssociationDTO.
    # Anonymous cassette data.
    anon_cassettes = []
    anon_cassette_tool_associations = []
    anon_cassette_genomic_entity_associations = []

    # Elaborate on get_general_data() for the ConstructHandler.
    def get_general_data(self, session):
        """Extend the method for the ConstructHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session)
        self.build_allele_gene_lookup(session)
        self.build_allele_class_lookup(session)
        self.build_seqfeat_gene_lookup(session)
        self.build_gene_tool_lookup(session)
        return

    # Elaborate on get_datatype_data() for the ConstructHandler.
    def get_allele_encoded_tools(self, session):
        """Get encoded FBto/FBsf objects for the constructs via alleles."""
        self.log.info('Get encoded FBto/FBsf objects for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        component = aliased(Feature, name='component')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            component.is_obsolete.is_(False),
            component.uniquename.op('~')(self.regex['fb_uniquename']),
            Cvterm.name == 'encodes_tool',
        )
        results = session.query(FeatureRelationship).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(component, (component.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        # Build pub lookup for relationship results.
        pub_results = session.query(FeatureRelationshipPub).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
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
        # Create allele feature_id-keyed lists of allele-component "encodes_tool" FBRelationship objects.
        al_tool_dict = {}
        counter = 0
        for result in results:
            fb_rel = fb_datatypes.FBRelationship(result, 'feature_relationship')
            if result.feature_relationship_id in rel_pub_dict.keys():
                fb_rel.pubs = rel_pub_dict[result.feature_relationship_id]
            try:
                al_tool_dict[result.subject_id].append(fb_rel)
                counter += 1
            except KeyError:
                al_tool_dict[result.subject_id] = [fb_rel]
                counter += 1
        self.log.info(f'Found {counter} allele-to-component "encodes_tool" relationships.')
        self.log.info('Now propagate these "encodes_tool" relationships to constructs.')
        cons_counter = 0
        for construct in self.fb_data_entities.values():
            al_cons_rels = construct.recall_relationships(self.log, entity_role='object', rel_types='associated_with', rel_entity_types='allele')
            for al_cons_rel in al_cons_rels:
                allele_id = al_cons_rel.chado_obj.subject_id
                try:
                    construct.al_encodes_tool_rels.extend(al_tool_dict[allele_id])
                    cons_counter += len(al_tool_dict[allele_id])
                except KeyError:
                    pass
        self.log.info(f'Propagated {cons_counter} allele-to-component "encodes_tool" relationships to related constructs.')
        return

    def get_allele_reg_regions(self, session):
        """Get regulatory FBto/FBsf/FBgn objects for the constructs via alleles."""
        self.log.info('Get regulatory FBto/FBsf/FBgn objects for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        component = aliased(Feature, name='component')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            component.is_obsolete.is_(False),
            component.uniquename.op('~')(self.regex['fb_uniquename']),
            Cvterm.name == 'has_reg_region',
        )
        results = session.query(FeatureRelationship).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(component, (component.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        # Build pub lookup for relationship results.
        pub_results = session.query(FeatureRelationshipPub).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
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
        # Create allele feature_id-keyed lists of allele-component "has_reg_region" FeatureRelationship objects.
        al_reg_region_dict = {}
        counter = 0
        for result in results:
            fb_rel = fb_datatypes.FBRelationship(result, 'feature_relationship')
            if result.feature_relationship_id in rel_pub_dict.keys():
                fb_rel.pubs = rel_pub_dict[result.feature_relationship_id]
            try:
                al_reg_region_dict[result.subject_id].append(fb_rel)
                counter += 1
            except KeyError:
                al_reg_region_dict[result.subject_id] = [fb_rel]
                counter += 1
        self.log.info(f'Found {counter} allele-to-component "has_reg_region" relationships.')
        self.log.info('Now propagate "has_reg_region" relationships to constructs.')
        cons_counter = 0
        for construct in self.fb_data_entities.values():
            al_cons_rels = construct.recall_relationships(self.log, entity_role='object', rel_types='associated_with', rel_entity_types='allele')
            for al_cons_rel in al_cons_rels:
                allele_id = al_cons_rel.chado_obj.subject_id
                try:
                    construct.al_reg_region_rels.extend(al_reg_region_dict[allele_id])
                    counter += len(al_reg_region_dict[allele_id])
                except KeyError:
                    pass
        self.log.info(f'Propagated {cons_counter} allele-to-component "has_reg_region" relationships to related constructs.')
        return

    def get_construct_tool_uses(self, session):
        """Get tool_uses FeatureCvterm data for FBtp constructs."""
        self.log.info('Get tool_uses FeatureCvterm data for FBtp constructs.')
        prop_type = aliased(Cvterm, name='prop_type')
        anno_cvterm = aliased(Cvterm, name='anno_cvterm')
        filters = (
            Feature.uniquename.op('~')(self.regex['construct']),
            prop_type.name == 'tool_uses',
        )
        if self.testing:
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(FeatureCvterm, FeatureCvtermprop).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(FeatureCvtermprop, (FeatureCvtermprop.feature_cvterm_id == FeatureCvterm.feature_cvterm_id)).\
            join(prop_type, (prop_type.cvterm_id == FeatureCvtermprop.type_id)).\
            join(anno_cvterm, (anno_cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for fc, fcp in results:
            feature_id = fc.feature_id
            if feature_id not in self.fb_data_entities:
                continue
            construct = self.fb_data_entities[feature_id]
            tool_use_entry = {
                'cvterm_name': fc.cvterm.name,
                'accession': fc.cvterm.dbxref.accession,
                'pub_id': fc.pub_id,
            }
            construct.tool_uses_data.append(tool_use_entry)
            counter += 1
        self.log.info(f'Found {counter} tool_uses annotations for constructs.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the ConstructHandler."""
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
        self.get_allele_encoded_tools(session)
        # marked_with rels already captured by get_entity_relationships(session, 'subject') above.
        self.get_allele_reg_regions(session)
        self.get_construct_tool_uses(session)
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_encoded_tools(self):
        """Synthesize encoded components."""
        self.log.info('Synthesize encoded components.')
        counter = 0
        for construct in self.fb_data_entities.values():
            self.log.debug(f'Assess encoded tools for {construct}.')
            # Reference of related alleles.
            cons_al_rels = construct.recall_relationships(self.log, entity_role='object', rel_types='associated_with', rel_entity_types='allele')
            # self.log.debug(f'{construct} has {len(cons_al_rels)} direct allele relationships.')
            # Direct encodes_tool relationships.
            cons_tool_rels = construct.recall_relationships(self.log, entity_role='subject', rel_types='encodes_tool')
            # self.log.debug(f'{construct} has {len(cons_tool_rels)} direct tool relationships.')
            # construct.expressed_features = {}
            for cons_tool_rel in cons_tool_rels:
                component_id = cons_tool_rel.chado_obj.object_id
                try:
                    construct.expressed_features[component_id].extend(cons_tool_rel.pubs)
                except KeyError:
                    construct.expressed_features[component_id] = cons_tool_rel.pubs
            self.log.debug(f'For {construct}, found {len(construct.expressed_features.keys())} encoded tools via direct relationships.')
            # Indirect encodes_tool relationships.
            # self.log.debug(f'{construct} has {len(construct.al_encodes_tool_rels)} indirect tool relationships via alleles.')
            for al_tool_rel in construct.al_encodes_tool_rels:
                allele_id = al_tool_rel.chado_obj.subject_id
                component_id = al_tool_rel.chado_obj.object_id
                try:
                    construct.expressed_features[component_id].extend(al_tool_rel.pubs)
                except KeyError:
                    construct.expressed_features[component_id] = al_tool_rel.pubs
                # Fold in pubs supporting the construct-allele relationship.
                for cons_al_rel in cons_al_rels:
                    if cons_al_rel.chado_obj.subject_id == allele_id:
                        construct.expressed_features[component_id].extend(cons_al_rel.pubs)
                        # self.log.debug(f'{construct} has these pubs via allele-tool: {cons_al_rel.pubs}')
            counter += len(construct.expressed_features.keys())
        self.log.info(f'Found {counter} encoded tools for constructs via direct and indirect (via allele) relationships.')
        return

    def synthesize_component_genes(self):
        """Synthesize component genes."""
        self.log.info('Synthesize component genes.')
        all_expressed_gene_counter = 0
        all_targeted_gene_counter = 0
        for construct in self.fb_data_entities.values():
            this_expressed_gene_counter = 0
            this_targeted_gene_counter = 0
            # Reference of related alleles.
            cons_al_rels = construct.recall_relationships(self.log, entity_role='object', rel_types='associated_with', rel_entity_types='allele')
            # self.log.debug(f'{construct} has {len(cons_al_rels)} direct allele relationships.')
            for cons_al_rel in cons_al_rels:
                allele_id = cons_al_rel.chado_obj.subject_id
                # Skip obsolete alleles.
                if allele_id not in self.allele_gene_lookup.keys():
                    continue
                gene_id = self.allele_gene_lookup[allele_id]
                # Slot for gene_id depends on the allele class.
                if allele_id in self.transgenic_allele_class_lookup.keys():
                    if set(self.transgenic_allele_class_lookup[allele_id]).intersection({'RNAi_reagent', 'sgRNA', 'antisense'}):
                        gene_slot = getattr(construct, 'targeted_features')
                        this_targeted_gene_counter += 1
                    else:
                        gene_slot = getattr(construct, 'expressed_features')
                        this_expressed_gene_counter += 1
                else:
                    gene_slot = getattr(construct, 'expressed_features')
                    this_expressed_gene_counter += 1
                try:
                    gene_slot[gene_id].extend(cons_al_rel.pubs)
                except KeyError:
                    gene_slot[gene_id] = cons_al_rel.pubs
            all_expressed_gene_counter += this_expressed_gene_counter
            all_targeted_gene_counter += this_targeted_gene_counter
        self.log.info(f'Found {all_expressed_gene_counter} expressed genes and {all_targeted_gene_counter} targeted genes for constructs.')
        return

    def synthesize_reg_regions(self):
        """Synthesize construct reg_region components."""
        self.log.info('Synthesize construct reg_region components.')
        counter = 0
        for construct in self.fb_data_entities.values():
            # Reference of related alleles.
            cons_al_rels = construct.recall_relationships(self.log, entity_role='object', rel_types='associated_with', rel_entity_types='allele')
            # self.log.debug(f'{construct} has {len(cons_al_rels)} direct allele relationships.')
            # Direct has_reg_region relationships (new implementation).
            cons_reg_region_rels = construct.recall_relationships(self.log, entity_role='subject', rel_types='has_reg_region')
            # self.log.debug(f'{construct} has {len(cons_reg_region_rels)} direct reg_region relationships.')
            for cons_reg_region_rel in cons_reg_region_rels:
                reg_region_id = cons_reg_region_rel.chado_obj.object_id
                try:
                    construct.regulating_features[reg_region_id].extend(cons_reg_region_rel.pubs)
                except KeyError:
                    construct.regulating_features[reg_region_id] = cons_reg_region_rel.pubs
            self.log.debug(f'For {construct}, found {len(construct.regulating_features.keys())} encoded reg_regions via direct relationships.')
            # Direct relationships to regulatory_regions (old implementation).
            old_cons_reg_region_rels = construct.recall_relationships(self.log, entity_role='object', rel_entity_types=['region', 'regulatory_region'])
            # self.log.debug(f'{construct} has {len(old_cons_reg_region_rels)} old style direct regulatory_region relationships.')
            for old_cons_reg_region_rel in old_cons_reg_region_rels:
                reg_region_id = old_cons_reg_region_rel.chado_obj.subject_id
                try:
                    construct.regulating_features[reg_region_id].extend(old_cons_reg_region_rel.pubs)
                except KeyError:
                    construct.regulating_features[reg_region_id] = old_cons_reg_region_rel.pubs
            # Indirect encodes_tool relationships.
            for al_reg_region_rel in construct.al_reg_region_rels:
                allele_id = al_reg_region_rel.chado_obj.subject_id
                reg_region_id = al_reg_region_rel.chado_obj.object_id
                try:
                    construct.regulating_features[reg_region_id].extend(al_reg_region_rel.pubs)
                except KeyError:
                    construct.regulating_features[reg_region_id] = al_reg_region_rel.pubs
                # Fold in pubs supporting the construct-allele relationship.
                for cons_al_rel in cons_al_rels:
                    if cons_al_rel.chado_obj.subject_id == allele_id:
                        construct.regulating_features[reg_region_id].extend(cons_al_rel.pubs)
                        # self.log.debug(f'{construct} has these pubs via allele-reg_region: {cons_al_rel.pubs}')
            # Add relationships to genes via seqfeat regulatory_regions.
            for reg_region_id in list(construct.regulating_features.keys()):
                pub_ids = construct.regulating_features[reg_region_id]
                uniquename = self.feature_lookup[reg_region_id]['uniquename']
                feat_type = self.feature_lookup[reg_region_id]['type']
                if uniquename.startswith('FBsf') and feat_type in ['region', 'regulatory_region']:
                    try:
                        gene_ids = self.seqfeat_gene_lookup[reg_region_id]
                    except KeyError:
                        continue
                    for gene_id in gene_ids:
                        try:
                            construct.regulating_features[gene_id].extend(pub_ids)
                        except KeyError:
                            construct.regulating_features[gene_id] = pub_ids
            counter += len(construct.regulating_features.keys())
        self.log.info(f'Found {counter} reg_regions for constructs via direct and indirect (via allele) relationships.')
        return

    def synthesize_redundant_tool_genes(self):
        """For constructs in which a gene and related tool are related, suppress redundant genes."""
        self.log.info('For constructs in which a gene and related tool are related, suppress redundant genes.')
        slot_names = {
            'expressed_features': 'expressed_tool_genes',
            'regulating_features': 'regulating_tool_genes',
        }
        for slot_name, tool_gene_slot_name in slot_names.items():
            self.log.info(f'Prune {slot_name} for constructs.')
            counter = 0
            for construct in self.fb_data_entities.values():
                pruning_list = []
                slot_bin = getattr(construct, slot_name)
                for feature_id in slot_bin.keys():
                    try:
                        all_gene_related_tool_ids = set(self.gene_tool_lookup[feature_id])
                        tool_overlap = all_gene_related_tool_ids.intersection(set(slot_bin.keys()))
                        if tool_overlap:
                            pruning_list.append(feature_id)
                            pruned_gene = f'{self.feature_lookup[feature_id]["name"]} ({self.feature_lookup[feature_id]["uniquename"]})'
                            tool_overlap_str = '|'.join([f'{self.feature_lookup[i]["name"]} ({self.feature_lookup[i]["uniquename"]})' for i in tool_overlap])
                            self.log.debug(f'For {construct}, suppress {pruned_gene} since related tools are more informative: {tool_overlap_str}')
                    except KeyError:
                        pass
                for gene_id in pruning_list:
                    tool_gene_bin = getattr(construct, tool_gene_slot_name)
                    tool_gene_bin.append(gene_id)
                    counter += 1
            self.log.info(f'Will suppress {counter} genes from construct {slot_name} that are better represented as tools.')
        return

    def identify_constructs_needing_anon_cassettes(self):
        """Identify constructs that need anonymous cassettes for direct tool data."""
        self.log.info('Identify constructs needing anonymous cassettes.')
        counter = 0
        excluded_counter = 0
        for construct in self.fb_data_entities.values():
            # Exclude constructs with 'FTA: generic TI construct' internal_notes.
            internal_notes = construct.props_by_type.get('internal_notes', [])
            is_generic_ti = False
            for note in internal_notes:
                if note.chado_obj.value == 'FTA: generic TI construct':
                    is_generic_ti = True
                    break
            if is_generic_ti:
                excluded_counter += 1
                continue
            # Check for direct tool relationships.
            direct_tool_rels = construct.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['encodes_tool', 'has_reg_region', 'tagged_with', 'carries_tool'])
            has_direct_rels = len(direct_tool_rels) > 0
            has_tool_uses = len(construct.tool_uses_data) > 0
            if has_direct_rels or has_tool_uses:
                construct.needs_anon_cassette = True
                counter += 1
        self.log.info(f'Identified {counter} constructs needing anonymous cassettes '
                      f'({excluded_counter} excluded as generic TI constructs).')
        return

    # Elaborate on synthesize_info() for the ConstructHandler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.flag_new_additions_and_obsoletes()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_encoded_tools()
        self.synthesize_component_genes()
        self.synthesize_reg_regions()
        self.synthesize_redundant_tool_genes()
        self.identify_constructs_needing_anon_cassettes()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_construct_basic(self):
        """Map basic FlyBase construct data to the Alliance LinkML object."""
        self.log.info('Map basic construct info to Alliance object.')
        for construct in self.fb_data_entities.values():
            agr_construct = self.agr_export_type()
            agr_construct.obsolete = construct.chado_obj.is_obsolete
            agr_construct.primary_external_id = f'FB:{construct.uniquename}'
            # agr_construct.mod_internal_id = f'FB.feature_id={construct.chado_obj.feature_id}'
            construct.linkmldto = agr_construct
        return

    # Note: There are two distinct mapping types: ConstructComponentSlotAnnotationDTO and ConstructGenomicEntityAssociationDTO.
    #       ConstructComponentSlotAnnotationDTO can list anything (whether or not the component is known to the Alliance).
    #       ConstructGenomicEntityAssociationDTO is restricted to Alliance-submitted objects (currently genes only).
    def map_construct_components(self):
        """Map current construct components to the Alliance LinkML object."""
        self.log.info('Map current construct components to the Alliance LinkML object.')
        component_slots = {
            'expressed_features': 'expresses',
            'targeted_features': 'targets',
            'regulating_features': 'is_regulated_by',
        }
        counter = 0
        for construct in self.fb_data_entities.values():
            for slot_name, rel_type in component_slots.items():
                slot_bin = getattr(construct, slot_name)
                for feature_id, pub_ids in slot_bin.items():
                    # Do not report obsolete components.
                    if self.feature_lookup[feature_id]['is_obsolete'] is True:
                        continue
                    elif slot_name == 'expressed_features' and feature_id in construct.expressed_tool_genes:
                        continue
                    elif slot_name == 'regulating_features' and feature_id in construct.regulating_tool_genes:
                        continue
                    # These are dumped out in the associations
                    if self.feature_lookup[feature_id]['curie'].startswith('FB:FBgn') & \
                            self.organism_lookup[self.feature_lookup[feature_id]['organism_id']]['is_drosophilid'] is True:
                        continue
                    feature = self.feature_lookup[feature_id]
                    symbol = feature['symbol']
                    organism_id = feature['organism_id']
                    pubs = self.lookup_pub_curies(pub_ids)
                    taxon_text = self.organism_lookup[organism_id]['full_species_name']
                    taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
                    component_dto = agr_datatypes.ConstructComponentSlotAnnotationDTO(rel_type, symbol, taxon_curie, taxon_text, pubs).dict_export()
                    construct.linkmldto.construct_component_dtos.append(component_dto)
                    counter += 1
        self.log.info(f'Mapped construct components to {counter} ConstructcomponentDTOs.')
        return

    def map_construct_genomic_associations(self):
        """Map current construct relations to the Alliance LinkML object."""
        self.log.info('Map current construct relations to the Alliance LinkML object.')
        slot_rel_types = {
            'expressed_features': 'expresses',
            'targeted_features': 'targets',
            'regulating_features': 'is_regulated_by',
        }
        counter = 0
        for feature_slot_name, rel_type in slot_rel_types.items():
            self.log.info(f'Sort out Alliance genomic entities from "{feature_slot_name}" to "{rel_type}" associations.')
            for construct in self.fb_data_entities.values():
                component_slot = getattr(construct, feature_slot_name)
                for feature_id, pub_ids in component_slot.items():
                    # Limit reported associations to genes: expand to tools (FBto) and seq features (FBsf) eventually.
                    if self.feature_lookup[feature_id]['type'] != 'gene':
                        continue
                    # Filter out genes that are better represented as tools for this construct.
                    if feature_id in construct.expressed_tool_genes or feature_id in construct.regulating_tool_genes:
                        continue
                    # Limit reported associations to genes lacking a MOD curie (i.e., only FBgn IDs).
                    if not self.feature_lookup[feature_id]['curie'].startswith('FB:FBgn'):
                        continue
                    # Limit reported associations to Drosophilid genes.
                    if self.organism_lookup[self.feature_lookup[feature_id]['organism_id']]['is_drosophilid'] is False:
                        continue
                    cons_curie = f'FB:{construct.uniquename}'
                    obj_curie = self.feature_lookup[feature_id]['curie']
                    pub_curies = self.lookup_pub_curies(pub_ids)
                    fb_rel = fb_datatypes.FBExportEntity()
                    rel_dto = agr_datatypes.ConstructGenomicEntityAssociationDTO(cons_curie, rel_type, obj_curie, pub_curies)
                    if construct.is_obsolete is True or self.feature_lookup[feature_id]['is_obsolete'] is True:
                        rel_dto.obsolete = True
                        rel_dto.internal = True
                    fb_rel.linkmldto = rel_dto
                    self.construct_associations.append(fb_rel)
                    counter += 1
        self.log.info(f'Mapped construct_relationships to {counter} ConstructGenomicEntityAssociationDTOs.')
        return

    def map_construct_cassette_associations(self):
        """Map construct marked_with cassette relations to ConstructCassetteAssociationDTO."""
        self.log.info('Map construct marked_with cassette associations.')
        counter = 0
        has_transcriptional_unit_list = ['FBal0345196', 'FBal0345198', 'FBal0407180',
                                         'FBal0407181', 'FBal0407182', 'FBal0368073']
        for construct in self.fb_data_entities.values():
            marked_with_rels = construct.recall_relationships(
                self.log, entity_role='subject', rel_types='marked_with', rel_entity_types='allele')
            for rel in marked_with_rels:
                cassette_feature_id = rel.chado_obj.object_id
                cassette_uniquename = self.feature_lookup[cassette_feature_id]['uniquename']
                cons_curie = f'FB:{construct.uniquename}'
                cassette_curie = f'FB:{cassette_uniquename}'
                pub_curies = self.lookup_pub_curies(rel.pubs)
                # has_transcriptional_unit_list are special cases per FTA-139.
                if cassette_uniquename in has_transcriptional_unit_list:
                    relation_name = 'has_transcriptional_unit'
                else:
                    relation_name = 'has_selectable_marker'
                fb_rel = fb_datatypes.FBExportEntity()
                rel_dto = agr_datatypes.ConstructCassetteAssociationDTO(
                    cons_curie, relation_name, cassette_curie, pub_curies)
                if construct.is_obsolete is True or self.feature_lookup[cassette_feature_id]['is_obsolete'] is True:
                    rel_dto.obsolete = True
                    rel_dto.internal = True
                fb_rel.linkmldto = rel_dto
                self.construct_cassette_associations.append(fb_rel)
                counter += 1
        self.log.info(f'Mapped {counter} ConstructCassetteAssociationDTOs.')
        return

    def map_anon_cassette_basic(self):
        """Create basic anonymous CassetteDTOs for constructs with direct tool data."""
        self.log.info('Map anonymous cassette basic info.')
        counter = 0
        for construct in self.fb_data_entities.values():
            if not construct.needs_anon_cassette:
                continue
            cassette_id = f'{construct.uniquename}_cas'
            agr_cassette = agr_datatypes.CassetteDTO()
            agr_cassette.placeholder = True
            agr_cassette.primary_external_id = cassette_id
            agr_cassette.obsolete = False
            agr_cassette.cassette_symbol_dto = agr_datatypes.NameSlotAnnotationDTO(
                'nomenclature_symbol', cassette_id, cassette_id, []).dict_export()
            # Build data_provider_dto for the cassette.
            dp_xref = agr_datatypes.CrossReferenceDTO(
                'FB', f'FB:{construct.uniquename}', 'construct', cassette_id).dict_export()
            agr_cassette.data_provider_dto = agr_datatypes.DataProviderDTO(dp_xref).dict_export()
            fb_entity = fb_datatypes.FBExportEntity()
            fb_entity.linkmldto = agr_cassette
            self.anon_cassettes.append(fb_entity)
            construct.anon_cassette_dto = agr_cassette
            counter += 1
        self.log.info(f'Created {counter} anonymous CassetteDTOs.')
        return

    def map_anon_cassette_simple_components(self):
        """Map simple component relationships to anonymous cassette associations."""
        self.log.info('Map anonymous cassette simple components (has_reg_region, tagged_with, carries_tool).')
        rel_type_mapping = {
            'has_reg_region': 'is_regulated_by',
            'tagged_with': 'tagged_with',
            'carries_tool': 'contains',
        }
        counter = 0
        for construct in self.fb_data_entities.values():
            if not construct.needs_anon_cassette:
                continue
            cassette_id = f'{construct.uniquename}_cas'
            direct_rels = construct.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['has_reg_region', 'tagged_with', 'carries_tool'])
            for rel in direct_rels:
                chado_rel_type = rel.chado_obj.type.name
                alliance_rel = rel_type_mapping[chado_rel_type]
                component_id = rel.chado_obj.object_id
                component = self.feature_lookup[component_id]
                component_curie = component['curie']
                pub_curies = self.lookup_pub_curies(rel.pubs)
                if component['uniquename'].startswith('FBto'):
                    fb_rel = fb_datatypes.FBExportEntity()
                    rel_dto = agr_datatypes.CassetteTransgenicToolAssociationDTO(
                        cassette_id, component_curie, pub_curies, False, alliance_rel)
                    fb_rel.linkmldto = rel_dto
                    self.anon_cassette_tool_associations.append(fb_rel)
                elif component['uniquename'].startswith('FBgn'):
                    fb_rel = fb_datatypes.FBExportEntity()
                    rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                        cassette_id, component_curie, pub_curies, False, alliance_rel)
                    fb_rel.linkmldto = rel_dto
                    self.anon_cassette_genomic_entity_associations.append(fb_rel)
                elif component['uniquename'].startswith('FBsf'):
                    symbol = component['symbol']
                    organism_id = component['organism_id']
                    taxon_text = self.organism_lookup[organism_id]['full_species_name']
                    taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
                    comp_dto = agr_datatypes.CassetteComponentSlotAnnotationDTO(
                        alliance_rel, symbol, taxon_curie, taxon_text, pub_curies).dict_export()
                    construct.anon_cassette_dto.cassette_component_dtos.append(comp_dto)
                counter += 1
        self.log.info(f'Mapped {counter} simple component associations for anonymous cassettes.')
        return

    def map_anon_cassette_encodes_tool(self):
        """Map encodes_tool relationships to anonymous cassette associations."""
        self.log.info('Map anonymous cassette encodes_tool relationships.')
        counter = 0
        for construct in self.fb_data_entities.values():
            if not construct.needs_anon_cassette:
                continue
            cassette_id = f'{construct.uniquename}_cas'
            direct_rels = construct.recall_relationships(
                self.log, entity_role='subject', rel_types='encodes_tool')
            for rel in direct_rels:
                component_id = rel.chado_obj.object_id
                component = self.feature_lookup[component_id]
                component_curie = component['curie']
                pub_curies = self.lookup_pub_curies(rel.pubs)
                if component['uniquename'].startswith('FBto'):
                    fb_rel = fb_datatypes.FBExportEntity()
                    rel_dto = agr_datatypes.CassetteTransgenicToolAssociationDTO(
                        cassette_id, component_curie, pub_curies, False, 'expresses')
                    fb_rel.linkmldto = rel_dto
                    self.anon_cassette_tool_associations.append(fb_rel)
                elif component['uniquename'].startswith('FBgn'):
                    fb_rel = fb_datatypes.FBExportEntity()
                    rel_dto = agr_datatypes.CassetteGenomicEntityAssociationDTO(
                        cassette_id, component_curie, pub_curies, False, 'expresses')
                    fb_rel.linkmldto = rel_dto
                    self.anon_cassette_genomic_entity_associations.append(fb_rel)
                elif component['uniquename'].startswith('FBsf'):
                    symbol = component['symbol']
                    organism_id = component['organism_id']
                    taxon_text = self.organism_lookup[organism_id]['full_species_name']
                    taxon_curie = self.organism_lookup[organism_id]['taxon_curie']
                    comp_dto = agr_datatypes.CassetteComponentSlotAnnotationDTO(
                        'expresses', symbol, taxon_curie, taxon_text, pub_curies).dict_export()
                    construct.anon_cassette_dto.cassette_component_dtos.append(comp_dto)
                counter += 1
        self.log.info(f'Mapped {counter} encodes_tool associations for anonymous cassettes.')
        return

    def map_anon_cassette_tool_uses(self):
        """Map tool_uses data to anonymous cassette use DTOs."""
        self.log.info('Map anonymous cassette tool_uses.')
        counter = 0
        for construct in self.fb_data_entities.values():
            if not construct.needs_anon_cassette:
                continue
            if not construct.tool_uses_data:
                continue
            pub_ids = [entry['pub_id'] for entry in construct.tool_uses_data]
            pub_curies = self.lookup_pub_curies(pub_ids)
            use_curies = list(set(
                f'FBcv:{entry["accession"]}' for entry in construct.tool_uses_data))
            slot_dto = agr_datatypes.CassetteUseSlotAnnotationDTO(
                pub_curies, use_curies).dict_export()
            construct.anon_cassette_dto.cassette_use_dtos.append(slot_dto)
            counter += 1
        self.log.info(f'Mapped tool_uses for {counter} anonymous cassettes.')
        return

    def map_anon_cassette_to_construct_association(self):
        """Create ConstructCassetteAssociationDTOs linking anonymous cassettes to constructs."""
        self.log.info('Map anonymous cassette to construct associations.')
        counter = 0
        for construct in self.fb_data_entities.values():
            if not construct.needs_anon_cassette:
                continue
            cassette_id = f'{construct.uniquename}_cas'
            cons_curie = f'FB:{construct.uniquename}'
            # Determine relation_name.
            if construct.tool_uses_data:
                relation_name = 'has_functional_unit'
            else:
                relation_name = 'has_component'
            # Collect ALL pub_ids from all direct MS14 data sources.
            all_pub_ids = []
            direct_rels = construct.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['encodes_tool', 'has_reg_region', 'tagged_with', 'carries_tool'])
            for rel in direct_rels:
                all_pub_ids.extend(rel.pubs)
            for entry in construct.tool_uses_data:
                all_pub_ids.append(entry['pub_id'])
            pub_curies = self.lookup_pub_curies(all_pub_ids)
            unique_pub_curies = list(set(pub_curies))
            # Build the association DTO.
            note_dtos = []
            if len(unique_pub_curies) == 1:
                evidence_curies = unique_pub_curies
            else:
                evidence_curies = []
                if len(unique_pub_curies) > 1:
                    note_dto = agr_datatypes.NoteDTO(
                        'internal_note',
                        'FTA: unable to automatically determine reference for cassette to construct association',
                        []).dict_export()
                    note_dtos.append(note_dto)
            fb_rel = fb_datatypes.FBExportEntity()
            rel_dto = agr_datatypes.ConstructCassetteAssociationDTO(
                cons_curie, relation_name, cassette_id, evidence_curies)
            rel_dto.note_dtos = note_dtos
            fb_rel.linkmldto = rel_dto
            self.construct_cassette_associations.append(fb_rel)
            counter += 1
        self.log.info(f'Created {counter} ConstructCassetteAssociationDTOs for anonymous cassettes.')
        return

    # Elaborate on map_fb_data_to_alliance() for the ConstructHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the ConstructHandler."""
        super().map_fb_data_to_alliance()
        self.map_construct_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        # self.map_xrefs(datatype)    # ConstructDTO lacks cross_reference_dtos attribute.
        self.map_pubs()
        self.map_timestamps()
        # Note - We do not use self.map_secondary_ids('construct_secondary_id_dtos') here.
        #        This is because for reagents, we report only strings, not SecondaryIdSlotAnnotationDTOs.
        for construct in self.fb_data_entities.values():
            construct.linkmldto.secondary_identifiers = construct.alt_fb_ids
        self.map_construct_components()
        self.flag_internal_fb_entities('fb_data_entities')
        self.map_construct_genomic_associations()
        self.flag_internal_fb_entities('construct_associations')
        self.map_construct_cassette_associations()
        self.flag_internal_fb_entities('construct_cassette_associations')
        # Anonymous cassette mapping.
        self.map_anon_cassette_basic()
        self.map_anon_cassette_simple_components()
        self.map_anon_cassette_encodes_tool()
        self.map_anon_cassette_tool_uses()
        self.map_anon_cassette_to_construct_association()
        self.flag_internal_fb_entities('anon_cassettes')
        self.flag_internal_fb_entities('anon_cassette_tool_associations')
        self.flag_internal_fb_entities('anon_cassette_genomic_entity_associations')
        return

    # Elaborate on query_chado_and_export() for the ConstructHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the ConstructHandler."""
        super().query_chado_and_export(session)
        self.flag_unexportable_entities(self.construct_associations, 'construct_genomic_entity_association_ingest_set')
        self.generate_export_dict(self.construct_associations, 'construct_genomic_entity_association_ingest_set')
        self.flag_unexportable_entities(self.construct_cassette_associations, 'construct_cassette_association_ingest_set')
        self.generate_export_dict(self.construct_cassette_associations, 'construct_cassette_association_ingest_set')
        # Anonymous cassette ingest sets.
        self.flag_unexportable_entities(self.anon_cassettes, 'cassette_ingest_set')
        self.generate_export_dict(self.anon_cassettes, 'cassette_ingest_set')
        self.flag_unexportable_entities(
            self.anon_cassette_tool_associations, 'cassette_transgenic_tool_association_ingest_set')
        self.generate_export_dict(
            self.anon_cassette_tool_associations, 'cassette_transgenic_tool_association_ingest_set')
        self.flag_unexportable_entities(
            self.anon_cassette_genomic_entity_associations, 'cassette_genomic_entity_association_ingest_set')
        self.generate_export_dict(
            self.anon_cassette_genomic_entity_associations, 'cassette_genomic_entity_association_ingest_set')
        return
