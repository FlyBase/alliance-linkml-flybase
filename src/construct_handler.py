"""Module:: construct_handler.

Synopsis:
    A data handler that exports FlyBase data for constructs, including their
    associations to other features, to Alliance Construct LinkML objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.production import (
    Cvterm, Feature, FeatureRelationship
)
import fb_datatypes
import agr_datatypes
from feature_handler import FeatureHandler


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
        'FBtp0017594': 'P{UAS(-FRT)ptc.Deltaloop2}'               # Obsolete, has only a non-current symbol synonym - for testing feature lookup.
    }

    # Additional set for export added to the handler.
    construct_associations = []            # Will be a list of FBExportEntity objects (relationships), map to ConstructGenomicEntityAssociationDTO.

    # Lookups needed.
    allele_gene_lookup = {}                # Will be allele feature_id-keyed of a single gene feature_id per allele.
    seqfeat_gene_lookup = {}               # Will be seqfeat feature_id-keyed of a lists of gene feature_ids.
    transgenic_allele_class_lookup = {}    # Will be an allele feature_id-keyed list of "transgenic product class" CV terms.
    gene_tool_lookup = {}                  # Will be gene feature_id-keyed lists of related FBto tools.

    # Elaborate on get_general_data() for the ConstructHandler.
    def get_general_data(self, session):
        """Extend the method for the ConstructHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        self.build_feature_lookup(session)
        self.build_feature_relationship_evidence_lookup(session)
        self.build_allele_class_lookup(session)
        self.build_seqfeat_gene_lookup(session)
        self.build_gene_tool_lookup(session)
        return

    # Elaborate on get_datatype_data() for the ConstructHandler.
    def get_construct_alleles(self, session):
        """Get allele(s) to which constructs belong."""
        self.log.info('Get allele(s) to which constructs belong.')
        self.get_entity_obj_feat_rel_by_type(session, 'parent_allele_rels', rel_type='associated_with', sbj_type='allele', sbj_regex=self.regex['allele'])
        return

    def get_construct_encoded_tools(self, session):
        """Get directly related encoded FBto/FBsf objects for the construct."""
        self.log.info('Get directly related encoded FBto/FBsf objects for the construct.')
        self.get_entity_sbj_feat_rel_by_type(session, 'encodes_tool_rels', rel_type='encodes_tool', obj_regex=self.regex['fb_uniquename'])
        return

    def get_construct_reg_regions(self, session):
        """Get directly related regulatory FBgn/FBto/FBsf objects for the construct."""
        self.log.info('Get directly related regulatory FBgn/FBto/FBsf objects for the construct.')
        self.get_entity_sbj_feat_rel_by_type(session, 'reg_region_rels', rel_type='has_reg_region', obj_regex=self.regex['fb_uniquename'])
        return

    def get_construct_reg_regions_old(self, session):
        """Get directly related regulatory_region FBsf objects, old type of association."""
        self.log.info('Get directly related regulatory_region FBsf objects, old type of association.')
        self.get_entity_obj_feat_rel_by_type(session, 'seqfeat_rels', rel_type='associated_with', sbj_type='regulatory_region', sbj_regex=self.regex['seqfeat'])
        return

    def get_allele_encoded_tools(self, session):
        """Get encoded FBto/FBsf objects for the constructs via alleles."""
        self.log.info('Get encoded FBto/FBsf objects for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        component = aliased(Feature, name='component')
        filters = (
            # allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            # component.is_obsolete.is_(False),
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
        # Create allele feature_id-keyed lists of allele-component "encodes_tool" FeatureRelationship objects.
        al_comp_dict = {}
        counter = 0
        for result in results:
            try:
                al_comp_dict[result.subject_id].append(result)
                counter += 1
            except KeyError:
                al_comp_dict[result.subject_id] = [result]
                counter += 1
        self.log.info(f'Found {counter} allele-to-component "encodes_tool" relationships.')
        self.log.info('Now propagate these "encodes_tool" relationships to constructs.')
        counter = 0
        for construct in self.fb_data_entities.values():
            for allele_rel in construct.parent_allele_rels:
                allele_id = allele_rel.subject_id
                try:
                    construct.al_encodes_tool_rels.extend(al_comp_dict[allele_id])
                    counter += len(construct.al_encodes_tool_rels)
                except KeyError:
                    pass
        self.log.info(f'Propagated {counter} allele-to-component "encodes_tool" relationships to related constructs.')
        return

    def get_allele_reg_regions(self, session):
        """Get regulatory FBto/FBsf/FBgn objects for the constructs via alleles."""
        self.log.info('Get regulatory FBto/FBsf/FBgn objects for the constructs via alleles.')
        allele = aliased(Feature, name='allele')
        component = aliased(Feature, name='component')
        filters = (
            # allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            # component.is_obsolete.is_(False),
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
        # Create allele feature_id-keyed lists of allele-component "has_reg_region" FeatureRelationship objects.
        al_comp_dict = {}
        counter = 0
        for result in results:
            try:
                al_comp_dict[result.subject_id].append(result)
                counter += 1
            except KeyError:
                al_comp_dict[result.subject_id] = [result]
                counter += 1
        self.log.info(f'Found {counter} allele-to-component "has_reg_region" relationships.')
        self.log.info('Now propagate "has_reg_region" relationships to constructs.')
        counter = 0
        for construct in self.fb_data_entities.values():
            for allele_rel in construct.parent_allele_rels:
                allele_id = allele_rel.subject_id
                try:
                    construct.al_reg_region_rels.extend(al_comp_dict[allele_id])
                    counter += len(construct.al_reg_region_rels)
                except KeyError:
                    pass
        self.log.info(f'Propagated {counter} allele-to-component "has_reg_region" relationships to related constructs.')
        return

    def get_cons_genes_via_alleles(self, session):
        """Get genes for the constructs via alleles."""
        self.log.info('Get genes for the constructs via alleles.')
        # First, create allele feature_id-keyed lists of allele-gene "alleleof" FeatureRelationship objects.
        allele = aliased(Feature, name='allele')
        gene = aliased(Feature, name='gene')
        filters = (
            # allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            # gene.is_obsolete.is_(False),
            gene.uniquename.op('~')(self.regex['gene']),
            Cvterm.name == 'alleleof',
        )
        results = session.query(FeatureRelationship).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureRelationship.type_id)).\
            join(gene, (gene.feature_id == FeatureRelationship.object_id)).\
            filter(*filters).\
            distinct()
        al_comp_dict = {}
        counter = 0
        for result in results:
            try:
                al_comp_dict[result.subject_id].append(result)
                counter += 1
            except KeyError:
                al_comp_dict[result.subject_id] = [result]
                counter += 1
        self.log.info(f'Found {counter} allele-to-gene "alleleof" relationships.')
        self.log.info('Now propagate "alleleof" relationships to constructs.')
        counter = 0
        for construct in self.fb_data_entities.values():
            for allele_rel in construct.parent_allele_rels:
                allele_id = allele_rel.subject_id
                try:
                    construct.al_genes.extend(al_comp_dict[allele_id])
                except KeyError:
                    pass
            counter += len(construct.al_genes)
        self.log.info(f'Propagated {counter} allele-to-gene "alleleof" relationships to related constructs.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the ConstructHandler."""
        super().get_datatype_data(session, datatype, fb_export_type, agr_export_type)
        self.get_entities(session, self.datatype, self.fb_export_type)
        self.get_entityprops(session, self.datatype)
        self.get_entity_pubs(session, self.datatype)
        self.get_entity_synonyms(session, self.datatype)
        self.get_entity_fb_xrefs(session, self.datatype)
        self.get_entity_xrefs(session, self.datatype)
        self.get_entity_timestamps(session, self.datatype)
        self.get_construct_alleles(session)
        self.get_construct_encoded_tools(session)
        self.get_construct_reg_regions(session)
        self.get_allele_encoded_tools(session)
        self.get_allele_reg_regions(session)
        self.get_construct_reg_regions_old(session)
        self.get_cons_genes_via_alleles(session)
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_encoded_tools(self):
        """Synthesize encoded components."""
        self.log.info('Synthesize encoded components.')
        counter = 0
        for construct in self.fb_data_entities.values():
            # self.log.debug(f'Assess encoded tools for {construct}.')
            # Direct encodes_tool relationships.
            for rel in construct.encodes_tool_rels:
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pub_ids(rel.feature_relationship_id)
                try:
                    construct.expressed_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.expressed_features[component_id] = pub_ids
            # direct_count = len(construct.expressed_features.keys())
            # self.log.debug(f'For {construct}, found {direct_count} encoded tools via direct relationships.')
            # Indirect encodes_tool relationships.
            for rel in construct.al_encodes_tool_rels:
                allele_id = rel.subject_id
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pub_ids(rel.feature_relationship_id)
                try:
                    construct.expressed_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.expressed_features[component_id] = pub_ids
                # Fold in pubs supporting the construct-allele relationship.
                for al_con_rel in construct.parent_allele_rels:
                    if al_con_rel.subject_id == allele_id:
                        al_con_pub_ids = self.lookup_feat_rel_pub_ids(al_con_rel.feature_relationship_id)
                        construct.expressed_features[component_id].extend(al_con_pub_ids)
            # indirect_count = len(construct.expressed_features.keys()) - direct_count
            # self.log.debug(f'For {construct}, found {indirect_count} encoded tools via indirect allele relationships.')
            counter += len(construct.expressed_features.keys())
        self.log.info(f'Found {counter} encoded tools for constructs via direct and indirect allele relationships.')
        return

    def synthesize_component_genes(self):
        """Synthesize component genes."""
        self.log.info('Synthesize component genes.')
        all_expressed_gene_counter = 0
        all_targeted_gene_counter = 0
        for construct in self.fb_data_entities.values():
            this_expressed_gene_counter = 0
            this_targeted_gene_counter = 0
            # self.log.debug(f'Assess component genes for {construct}.')
            for rel in construct.al_genes:
                allele_id = rel.subject_id
                gene_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pub_ids(rel.feature_relationship_id)
                # Slot for gene_id depends on the allele class.
                try:
                    if set(self.transgenic_allele_class_lookup[allele_id]).intersection({'RNAi_reagent', 'sgRNA', 'antisense'}):
                        gene_slot = getattr(construct, 'targeted_features')
                        this_targeted_gene_counter += 1
                    else:
                        gene_slot = getattr(construct, 'expressed_features')
                        this_expressed_gene_counter += 1
                except KeyError:
                    gene_slot = getattr(construct, 'expressed_features')
                    this_expressed_gene_counter += 1
                try:
                    gene_slot[gene_id].extend(pub_ids)
                except KeyError:
                    gene_slot[gene_id] = pub_ids
                # Fold in pubs supporting the construct-allele relationship.
                for al_con_rel in construct.parent_allele_rels:
                    if al_con_rel.subject_id == allele_id:
                        al_con_pub_ids = self.lookup_feat_rel_pub_ids(al_con_rel.feature_relationship_id)
                        gene_slot[gene_id].extend(al_con_pub_ids)
            # self.log.debug(f'For {construct}, found {this_expressed_gene_counter} expressed genes and {this_targeted_gene_counter} targeted genes.')
            all_expressed_gene_counter += this_expressed_gene_counter
            all_targeted_gene_counter += this_targeted_gene_counter
        self.log.info(f'Found {all_expressed_gene_counter} expressed genes and {all_targeted_gene_counter} targeted genes for constructs.')
        return

    def synthesize_reg_regions(self):
        """Synthesize construct reg_region components."""
        self.log.info('Synthesize construct reg_region components.')
        counter = 0
        for construct in self.fb_data_entities.values():
            # self.log.debug(f'Assess reg_regions for {construct}.')
            # Direct has_reg_region relationships.
            for rel in construct.reg_region_rels:
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pub_ids(rel.feature_relationship_id)
                try:
                    construct.regulating_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.regulating_features[component_id] = pub_ids
            # direct_count = len(construct.regulating_features.keys())
            # self.log.debug(f'For {construct}, found {direct_count} reg_regions via direct relationships.')
            # Direct seqfeat relationships, old_style.
            for rel in construct.seqfeat_rels:
                component_id = rel.subject_id
                pub_ids = self.lookup_feat_rel_pub_ids(rel.feature_relationship_id)
                try:
                    construct.regulating_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.regulating_features[component_id] = pub_ids
            # direct_count_old = len(construct.regulating_features.keys()) - direct_count
            # self.log.debug(f'For {construct}, found {direct_count_old} reg_regions via direct relationships, old style.')
            # Indirect has_reg_region relationships.
            for rel in construct.al_reg_region_rels:
                allele_id = rel.subject_id
                component_id = rel.object_id
                pub_ids = self.lookup_feat_rel_pub_ids(rel.feature_relationship_id)
                try:
                    construct.regulating_features[component_id].extend(pub_ids)
                except KeyError:
                    construct.regulating_features[component_id] = pub_ids
                # Fold in pubs supporting the construct-allele relationship.
                for al_con_rel in construct.parent_allele_rels:
                    if al_con_rel.subject_id == allele_id:
                        al_con_pub_ids = self.lookup_feat_rel_pub_ids(al_con_rel.feature_relationship_id)
                        construct.regulating_features[component_id].extend(al_con_pub_ids)
            # indirect_count = len(construct.regulating_features.keys()) - direct_count - direct_count_old
            # self.log.debug(f'For {construct}, found {indirect_count} reg_regions tools via indirect allele relationships.')
            # Indirect relationships to genes via seqfeats.
            for component_id in list(construct.regulating_features.keys()):
                pub_ids = construct.regulating_features[component_id]
                uniquename = self.feature_lookup[component_id]['uniquename']
                feat_type = self.feature_lookup[component_id]['type']
                if uniquename.startswith('FBsf') and feat_type in ['region', 'regulatory_region']:
                    try:
                        gene_ids = self.seqfeat_gene_lookup[component_id]
                    except KeyError:
                        continue
                    for gene_id in gene_ids:
                        try:
                            construct.regulating_features[gene_id].extend(pub_ids)
                        except KeyError:
                            construct.regulating_features[gene_id] = pub_ids
            # genes_via_seqfeat_count = len(construct.regulating_features.keys()) - direct_count - direct_count_old - indirect_count
            # self.log.debug(f'For {construct}, found an additional {genes_via_seqfeat_count} genes related to seqfeat reg_regions.')
            counter += len(construct.regulating_features.keys())
        self.log.info(f'Found {counter} reg_regions for constructs via direct and indirect allele relationships.')
        return

    def synthesize_redundant_tool_genes(self):
        """For constructs in which a gene and related tool are related, sort out redundant genes."""
        self.log.info('For constructs in which a gene and related tool are related, sort out redundant genes.')
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
                        all_related_tool_ids = set(self.gene_tool_lookup[feature_id])
                        tool_overlap = all_related_tool_ids.intersection(set(slot_bin.keys()))
                        # self.log.debug(f'For {construct}, {self.feature_lookup[feature_id]["name"]} has {len(tool_overlap)} redundantly associated tools.')
                        if tool_overlap:
                            pruning_list.append(feature_id)
                            pruned_gene = f'{self.feature_lookup[feature_id]["name"]} ({self.feature_lookup[feature_id]["uniquename"]})'
                            tool_overlap_str = '|'.join([f'{self.feature_lookup[i]["name"]} ({self.feature_lookup[i]["uniquename"]})' for i in tool_overlap])
                            self.log.debug(f'For {construct}, sort out {pruned_gene} since related tools are more informative: {tool_overlap_str}')
                    except KeyError:
                        pass
                for gene_id in pruning_list:
                    tool_gene_bin = getattr(construct, tool_gene_slot_name)
                    tool_gene_bin.append(gene_id)
                    counter += 1
            self.log.info(f'Pruned {counter} genes from construct {slot_name} that are better represented as tools.')
        return

    def synthesize_construct_genomic_entity_associations(self):
        """Synthesize construct-genomic entity associations."""
        self.log.info('Synthesize construct-genomic entity associations.')
        slot_rel_types = {
            'expressed_features': 'expresses',
            'targeted_features': 'targets',
            'regulating_features': 'is_regulated_by',
        }
        for feature_slot_name, rel_type in slot_rel_types.items():
            self.log.info(f'Sort out Alliance genomic entities from "{feature_slot_name}" to "{rel_type}" associations.')
            counter = 0
            for construct in self.fb_data_entities.values():
                component_slot = getattr(construct, feature_slot_name)
                for feature_id, pub_ids in component_slot.items():
                    if self.feature_lookup[feature_id]['type'] != 'gene' or not self.feature_lookup[feature_id]['uniquename'].startswith('FBgn'):
                        continue
                    feat_rel = fb_datatypes.FBRelationship('feature_relationship', construct.db_primary_id, feature_id, rel_type)
                    feat_rel.pub_ids = pub_ids
                    feat_rel.entity_desc = f'{construct.uniquename} {rel_type} {self.feature_lookup[feature_id]["uniquename"]}'
                    self.construct_associations.append(feat_rel)
                    counter += 1
            self.log.info(f'Synthesized {counter} construct-gene associations.')
        return

    # Elaborate on synthesize_info() for the ConstructHandler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_encoded_tools()
        self.synthesize_component_genes()
        self.synthesize_reg_regions()
        self.synthesize_redundant_tool_genes()
        self.synthesize_construct_genomic_entity_associations()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_construct_basic(self, agr_export_type):
        """Map basic FlyBase construct data to the Alliance LinkML object."""
        self.log.info('Map basic construct info to Alliance object.')
        for construct in self.fb_data_entities.values():
            agr_construct = agr_export_type()
            agr_construct.obsolete = construct.chado_obj.is_obsolete
            agr_construct.mod_entity_id = f'FB:{construct.uniquename}'
            agr_construct.mod_internal_id = str(construct.chado_obj.feature_id)
            construct.linkmldto = agr_construct
        return

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
                    # Do not report genes that are better reported as tools.
                    elif slot_name == 'expressed_features' and feature_id in construct.expressed_tool_genes:
                        continue
                    elif slot_name == 'regulating_features' and feature_id in construct.regulating_tool_genes:
                        continue
                    symbol = self.feature_lookup[feature_id]['symbol']
                    pubs = self.lookup_pub_curies(pub_ids)
                    taxon_text = self.feature_lookup[feature_id]['species']
                    taxon_curie = self.feature_lookup[feature_id]['taxon_id']
                    component_dto = agr_datatypes.ConstructComponentSlotAnnotationDTO(rel_type, symbol, taxon_curie, taxon_text, pubs).dict_export()
                    construct.linkmldto.construct_component_dtos.append(component_dto)
                    counter += 1
        self.log.info(f'Mapped construct components to {counter} ConstructcomponentDTOs.')
        return

    def map_construct_genomic_associations(self):
        """Map current construct relations to the Alliance LinkML object."""
        self.log.info('Map current construct relations to the Alliance LinkML object.')
        counter = 0
        for cons_asso in self.construct_associations:
            cons_curie = f'FB:{self.feature_lookup[cons_asso.subject_id]["uniquename"]}'
            obj_curie = f'FB:{self.feature_lookup[cons_asso.object_id]["uniquename"]}'
            pub_curies = self.lookup_pub_curies(cons_asso.pub_ids)
            rel_dto = agr_datatypes.ConstructGenomicEntityAssociationDTO(cons_curie, cons_asso.rel_type, obj_curie, pub_curies)
            if self.feature_lookup[cons_asso.subject_id]['is_obsolete'] is True or self.feature_lookup[cons_asso.object_id]['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            cons_asso.linkmldto = rel_dto
            counter += 1
        self.log.info(f'Mapped construct_relationships to {counter} ConstructGenomicEntityAssociationDTOs.')
        return

    # Elaborate on map_fb_data_to_alliance() for the ConstructHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the ConstructHandler."""
        super().map_fb_data_to_alliance()
        self.map_construct_basic(agr_export_type)
        self.map_synonyms(datatype, agr_export_type)
        self.map_data_provider_dto(datatype)
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
        return

    # Elaborate on query_chado_and_export() for the ConstructHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the ConstructHandler."""
        super().query_chado_and_export(session, datatype, fb_export_type, agr_export_type)
        self.flag_unexportable_entities(self.construct_associations, 'construct_genomic_entity_association_ingest_set')
        self.generate_export_dict(self.construct_associations, 'construct_genomic_entity_association_ingest_set')
        return
