"""Module:: allele_handlers.

Synopsis:
    A data handler that exports FlyBase data for alleles, insertions and
    aberrations to Alliance Allele LinkML objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
import agr_datatypes
from fb_datatypes import (
    FBAberration, FBAllele, FBBalancer
)
from feature_handler import FeatureHandler
from harvdev_utils.reporting import (
    Cvterm, Feature, FeatureCvterm, FeatureGenotype, FeatureRelationship,
    Genotype, Phenotype, PhenotypeCvterm, Phenstatement, Pub
)
from utils import export_chado_data


class MetaAlleleHandler(FeatureHandler):
    """This objects gets, synthesizes and filters data for various FlyBase features exported as alleles."""
    def __init__(self, log: Logger, testing: bool):
        """Create the generic MetaAlleleHandler."""
        super().__init__(log, testing)
        self.agr_export_type = agr_datatypes.AlleleDTO
        self.primary_export_set = 'allele_ingest_set'

    # Additional sub-methods for map_fb_data_to_alliance().
    def map_metaallele_basic(self):
        """Map basic FlyBase metaallele data to the Alliance LinkML Allele object."""
        self.log.info('Map basic FlyBase metaallele data to the Alliance LinkML Allele object.')
        for metaallele in self.fb_data_entities.values():
            agr_allele = self.agr_export_type()
            agr_allele.obsolete = metaallele.chado_obj.is_obsolete
            agr_allele.primary_external_id = f'FB:{metaallele.uniquename}'
            # agr_allele.mod_internal_id = f'FB.feature_id={metaallele.chado_obj.feature_id}'
            agr_allele.taxon_curie = metaallele.ncbi_taxon_id
            metaallele.linkmldto = agr_allele
        return

    def map_metaallele_database_status(self):
        """Map metaallele database status."""
        self.log.info('Map metaallele database status.')
        counter = 0
        evidence_curies = []
        for metaallele in self.fb_data_entities.values():
            if metaallele.is_obsolete is False:
                db_status = 'approved'
            else:
                db_status = 'deleted'
                counter += 1
            db_status_annotation = agr_datatypes.AlleleDatabaseStatusSlotAnnotationDTO(db_status, evidence_curies)
            metaallele.linkmldto.allele_database_status_dto = db_status_annotation.dict_export()
        self.log.info(f'Marked {counter} {self.datatype}s as having database status "deleted".')
        return

    def map_extinction_info(self):
        """Map extinction info."""
        self.log.info('Map extinction info.')
        counter = 0
        for metaallele in self.fb_data_entities.values():
            if 'availability' in metaallele.props_by_type.keys():
                for prop in metaallele.props_by_type['availability']:
                    if prop.chado_obj.value == 'Stated to be lost.':
                        metaallele.linkmldto.is_extinct = True
            for prop_type in metaallele.props_by_type.keys():
                if prop_type.startswith('derived_stock'):
                    # Stock availability trumps curated extinction comment.
                    if metaallele.linkmldto.is_extinct is True:
                        self.log.warning(f'{metaallele} is stated to be lost but has stocks.')
                    metaallele.linkmldto.is_extinct = False
            if metaallele.linkmldto.is_extinct is True:
                counter += 1
        self.log.info(f'Flagged {counter} {self.datatype}s as extinct.')
        return

    def map_collections(self):
        """Map metaallele collections."""
        self.log.info('Map metaallele collections.')
        for metaallele in self.fb_data_entities.values():
            collections = []
            if metaallele.reagent_colls:
                collections.extend(metaallele.reagent_colls)
            elif metaallele.al_reagent_colls:
                collections.extend(metaallele.al_reagent_colls)
            elif metaallele.ti_reagent_colls:
                collections.extend(metaallele.ti_reagent_colls)
            elif metaallele.tp_reagent_colls:
                collections.extend(metaallele.tp_reagent_colls)
            elif metaallele.sf_reagent_colls:
                collections.extend(metaallele.sf_reagent_colls)
            if collections:
                collections = list(set(collections))
                metaallele.linkmldto.in_collection_name = collections[0].name
                if len(collections) > 1:
                    self.log.warning(f'{metaallele} has many relevant collections: {[i.name for i in collections]}')
        return

    def map_internal_metaallele_status(self):
        """Flag internal metaalleles using metaallele-specific criteria."""
        self.log.info('Flag internal metaalleles using metaallele-specific criteria.')
        internal_gene_counter = 0
        non_dmel_drosophilid_counter = 0
        for metaallele in self.fb_data_entities.values():
            final_org_abbr = self.organism_lookup[metaallele.organism_id]['abbreviation']
            if final_org_abbr != 'Dmel':
                metaallele.linkmldto.internal = True
                metaallele.internal_reasons.append('Non-Dmel')
                non_dmel_drosophilid_counter += 1
            if metaallele.uniquename.startswith('FBal') and metaallele.allele_of_internal_gene is True:
                metaallele.linkmldto.internal = True
                metaallele.internal_reasons.append('Allele of internal type FB gene')
                internal_gene_counter += 1
        if self.datatype == 'allele':
            self.log.info(f'Flagged {internal_gene_counter} alleles of internal-type genes as internal.')
        self.log.info(f'Flagged {non_dmel_drosophilid_counter} non-Dmel Drosophilid alleles as internal.')
        return


class AlleleHandler(MetaAlleleHandler):
    """This object gets, synthesizes and filters allele data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AlleleHandler object."""
        super().__init__(log, testing)
        self.datatype = 'allele'
        self.fb_export_type = FBAllele
        self.at_locus_fbal_fbti_dict = {}      # Will be FBal-feature_id-keyed dict of FBti feature_id lists ("is_represented_at_Alliance_as").
        self.transgenic_fbal_fbti_dict = {}    # Will be FBal-feature_id-keyed dict of FBti feature_id lists (via FBtp to unspecified FBti).
        self.fbti_entities = {}                # Will be feature_id-keyed FBAllele objects generated from FBti insertions.

    test_set = {
        'FBal0137236': 'gukh[142]',             # Insertion allele, part of TI_set_P{hsneo}.BDGP collection.
        'FBal0018482': 'wg[1]',                 # X-ray mutation.
        'FBal0015148': 'Sb[Spi]',               # point mutation.
        'FBal0043981': 'Ecol_lacZ[en-14]',      # Has an allele full name. Relationship to ARG has no pub support.
        'FBal0279489': 'Scer_GAL4[how-24B]',    # Has a 2o ID.
        'FBal0000010': 'alphaTub67C[3]',        # Has semidominant annotation from single locus homozygous genotype.
        'FBal0403680': 'Atg8a[3-31]',           # Has recessive annotation from single locus homozygous genotype.
        'FBal0011189': 'sti[1]',                # Has recessive annotation from single locus hemizygous genotype over a deficiency.
        'FBal0007942': '18w[00053]',            # Has recessive annotation from single locus unspecified zygosity genotype.
        'FBal0015410': 'sei[2]',                # Has codominant annotation from single locus unspecified zygosity genotype.
        'FBal0198096': 'tal[1]',                # Allele of internal type gene tal (gene_with_polycistronic_transcript).
        'FBal0055793': 'wg[UAS.Tag:HA]',        # Allele is directly related to a construct.
        'FBal0048226': 'Dmau_w[a23]',           # Non-Dmel allele related to non-Dmel insertion.
        'FBal0011649': 'Dsim_Lhr[1]',           # Non-Dmel classical allele.
        'FBal0043132': 'Hsap_MAPT[UAS.cAa]',    # Dmel transgenic allele expressing human gene.
        'FBal0062057': 'Scer_CDC42[V12.hs]',    # Dmel transgenic allele expressing yeast gene.
        'FBal0198528': 'CG33269[HMJ22303]',
    }

    # Additional export sets.
    allele_gene_rels = {}            # Will be (allele feature_id, gene feature_id) tuples keying lists of FBRelationships.
    allele_gene_associations = []    # Will be the final list of gene-allele FBRelationships to export (AlleleGeneAssociationDTO under linkmldto attr).

    # Additional reference info.
    allele_class_terms = []          # A list of cvterm_ids for child terms of "allele_class" (FBcv:0000286).
    allele_mutant_type_terms = []    # A list of cvterm_ids for child terms of chromosome_structure_variation or sequence_alteration.
    inheritance_mode_terms = {
        'recessive': 'recessive',
        'dominant': 'dominant',
        'semidominant': 'semi-dominant',
        'codominant': 'codominant'
    }

    # Additional sub-methods for get_general_data().
    def get_key_cvterm_sets(self, session):
        """Get key CV term sets from chado."""
        self.log.info('Get key CV term sets from chado.')
        self.allele_class_terms.extend(self.get_child_cvterms(session, 'allele class', 'FlyBase miscellaneous CV'))
        self.allele_mutant_type_terms.extend(self.get_child_cvterms(session, 'chromosome_structure_variation', 'SO'))
        self.allele_mutant_type_terms.extend(self.get_child_cvterms(session, 'sequence_alteration', 'SO'))
        return

    # Elaborate on get_general_data() for the AlleleHandler.
    def get_general_data(self, session):
        """Extend the method for the AlleleHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.get_key_cvterm_sets(session)
        self.build_feature_lookup(session, feature_types=['construct', 'gene', 'insertion', 'variation', 'transposon'])
        self.get_internal_genes(session)
        return

    # Additional sub-methods for get_datatype_data().
    def find_in_vitro_alleles(self, session):
        """Find in vitro alleles."""
        self.log.info('Find in vitro alleles.')
        cvterm_name_regex = '^in vitro construct'
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            Cvterm.name.op('~')(cvterm_name_regex)
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        ivt_alleles = session.query(Feature).\
            select_from(Feature).\
            join(FeatureCvterm, (FeatureCvterm.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == FeatureCvterm.cvterm_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in ivt_alleles:
            self.fb_data_entities[result.feature_id].in_vitro = True
            counter += 1
        self.log.info(f'Flagged {counter} alleles as "in vitro"')
        return

    def get_phenotypes(self, session):
        """Get phenotypes from single locus genotypes related to alleles."""
        self.log.info('Get phenotypes from single locus genotypes related to alleles.')
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
            Genotype.is_obsolete.is_(False),
            Genotype.description.op('!~')('_'),
            Cvterm.name.in_((self.inheritance_mode_terms.keys())),
            Pub.is_obsolete.is_(False),
            Pub.uniquename.op('~')(self.regex['pub']),
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Feature, Genotype, Phenotype, Cvterm, Pub).\
            join(FeatureGenotype, (FeatureGenotype.feature_id == Feature.feature_id)).\
            join(Genotype, (Genotype.genotype_id == FeatureGenotype.genotype_id)).\
            join(Phenstatement, (Phenstatement.genotype_id == Genotype.genotype_id)).\
            join(Phenotype, (Phenotype.phenotype_id == Phenstatement.phenotype_id)).\
            join(PhenotypeCvterm, (PhenotypeCvterm.phenotype_id == Phenotype.phenotype_id)).\
            join(Cvterm, (Cvterm.cvterm_id == PhenotypeCvterm.cvterm_id)).\
            join(Pub, (Pub.pub_id == Phenstatement.pub_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            self.fb_data_entities[result.Feature.feature_id].phenstatements.append(result)
            counter += 1
        self.log.info(f'Found {counter} allele phenotypes from single locus genotypes.')
        return

    def get_fbal_fbti_replacements(self, session):
        """Find FBal-FBti replacements to make."""
        self.log.info('Build allele-insertion replacement lookup.')
        direct_fbal_fbti_counter = 0
        indirect_fbal_fbti_counter = 0
        # First, get relationships between FBal alleles and at-locus FBti insertions (FBti is inside FBal).
        for allele in self.fb_data_entities.values():
            fbal_fbti_alliance_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='is_represented_at_alliance_as',
                                                                  rel_entity_types=self.feature_subtypes['insertion'])
            for feat_rel in fbal_fbti_alliance_rels:
                insertion = self.feature_lookup[feat_rel.chado_obj.object_id]
                self.log.debug(f'Report {allele} under {insertion["name"]} ({insertion["uniquename"]}).')
                try:
                    self.at_locus_fbal_fbti_dict[allele.db_primary_id].append(insertion["feature_id"])
                    self.log.warning(f'Found another FBti for {allele}, but expected a one-to-one relationship.')
                except KeyError:
                    self.at_locus_fbal_fbti_dict[allele.db_primary_id] = [(insertion["feature_id"])]
                    direct_fbal_fbti_counter += 1
        self.log.info(f'Found {direct_fbal_fbti_counter} FBal alleles to be replaced by at-locus FBti insertions in export file.')
        # Second, get relationships between FBal alleles and transgenic insertions (FBal is inside FBti).
        allele = aliased(Feature, name='allele')
        construct = aliased(Feature, name='construct')
        insertion = aliased(Feature, name='insertion')
        ac_rel = aliased(FeatureRelationship, name='ac_rel')
        ic_rel = aliased(FeatureRelationship, name='ic_rel')
        ac_rel_type = aliased(Cvterm, name='ac_rel_type')
        ic_rel_type = aliased(Cvterm, name='ic_rel_type')
        filters = (
            allele.is_obsolete.is_(False),
            allele.uniquename.op('~')(self.regex['allele']),
            construct.is_obsolete.is_(False),
            construct.uniquename.op('~')(self.regex['construct']),
            insertion.is_obsolete.is_(False),
            insertion.uniquename.op('~')(self.regex['insertion']),
            insertion.name.op('~')('unspecified$'),
            ac_rel_type.name == 'associated_with',
            ic_rel_type.name == 'producedby',
        )
        results = session.query(allele, insertion).\
            select_from(allele).\
            join(ac_rel, (ac_rel.subject_id == allele.feature_id)).\
            join(ac_rel_type, (ac_rel_type.cvterm_id == ac_rel.type_id)).\
            join(construct, (construct.feature_id == ac_rel.object_id)).\
            join(ic_rel, (ic_rel.object_id == construct.feature_id)).\
            join(ic_rel_type, (ic_rel_type.cvterm_id == ic_rel.type_id)).\
            join(insertion, (insertion.feature_id == ic_rel.subject_id)).\
            filter(*filters).\
            distinct()
        for result in results:
            allele_str = f'{result.allele.name} ({result.allele.uniquename})'
            insertion_str = f'{result.insertion.name} ({result.insertion.uniquename})'
            self.log.debug(f'The transgenic allele {allele_str} is related to the generic insertion {insertion_str}.')
            try:
                self.transgenic_fbal_fbti_dict[result.allele.feature_id].append(result.insertion.feature_id)
            except KeyError:
                self.transgenic_fbal_fbti_dict[result.allele.feature_id] = [result.insertion.feature_id]
                indirect_fbal_fbti_counter += 1
        self.log.info(f'Found {indirect_fbal_fbti_counter} FBal alleles to be replaced by transgenic FBti insertions in export file.')
        return

    def get_insertion_entities(self, session):
        """Have the AlleleHandler run the InsertionHandler."""
        separator = 80 * '#'
        self.log.info(f'Have the AlleleHandler run the InsertionHandler.\n{separator}')
        insertion_handler = InsertionHandler(self.log, self.testing)
        export_chado_data(session, self.log, insertion_handler)
        self.fbti_entities = insertion_handler.fb_data_entities
        self.log.info(f'The AlleleHandler obtained {len(self.fbti_entities)} FBti insertions from chado.\n{separator}')
        return

    # Elaborate on get_datatype_data() for the AlleleHandler.
    def get_datatype_data(self, session):
        """Extend the method for the AlleleHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject', rel_type='alleleof', entity_type='gene', entity_regex=self.regex['gene'])
        al_cons_fr_types = ['derived_tp_assoc_alleles', 'associated_with', 'gets_expression_data_from']
        self.get_entity_relationships(session, 'subject', rel_type=al_cons_fr_types, entity_type='construct', entity_regex=self.regex['construct'])
        al_ins_fr_types = ['is_represented_at_alliance_as', 'associated_with']
        self.get_entity_relationships(session, 'subject', rel_type=al_ins_fr_types, entity_type='insertion', entity_regex=self.regex['insertion'])
        self.get_entity_relationships(session, 'object', rel_type='partof', entity_type='variation')
        self.get_entity_cvterms(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_direct_reagent_collections(session)
        self.get_indirect_reagent_collections(session, 'subject', 'associated_with', 'insertion')
        self.get_indirect_reagent_collections(session, 'subject', 'associated_with', 'construct')
        self.get_indirect_reagent_collections(session, 'subject', ['has_reg_region', 'encodes_tool'], 'seqfeat')
        self.get_very_indirect_reagent_collections(session)
        self.get_phenotypes(session)
        self.find_in_vitro_alleles(session)
        self.get_fbal_fbti_replacements(session)
        self.get_insertion_entities(session)
        return

    # Additional sub-methods to be run by synthesize_info() below.
    def add_fbal_to_fbti(self, allele):
        """Add FBal data to a related FBti."""
        fbti_feature_ids = []
        fbti_feature_ids.append(allele.superceded_by_at_locus_insertion)
        fbti_feature_ids.extend(allele.superceded_by_transgnc_insertions)
        lists_to_extend = [
            'dbxrefs',
            'export_warnings',
            'fb_sec_dbxrefs',
            'internal_reasons',
            'new_timestamps',
            'phenstatements',
            'pub_associations',
            'synonyms',
            'timestamps',
        ]
        dicts_of_elements_to_add = [
            'cvt_annos_by_id',
            'rels_by_id',
        ]
        dicts_of_lists_to_add = [
            'props_by_type',
            'cvt_anno_ids_by_cv',
            'cvt_anno_ids_by_prop',
            'cvt_anno_ids_by_term',
            'obj_rel_ids_by_type',
            'sbj_rel_ids_by_type',
        ]

        for fbti_feature_id in fbti_feature_ids:
            if fbti_feature_id not in self.fb_data_entities:
                self.fb_data_entities[fbti_feature_id] = self.fbti_entities[fbti_feature_id]
            insertion = self.fb_data_entities[fbti_feature_id]
            self.log.debug(f'Merge {allele} data into {insertion} data.')
            insertion.alt_fb_ids.append(f'FB:{allele.uniquename}')
            for attr_name in lists_to_extend:
                allele_list = getattr(allele, attr_name)
                insertion_list = getattr(insertion, attr_name)
                self.log.debug(f'For {attr_name}, add {len(allele_list)} allele elements to {len(insertion_list)} insertion elements.')
                insertion_list.extend(allele_list)
            # Add to ID-keyed dict of single chado annotations (key is unique for FBCVtermAnnotation or FBRelationship).
            for attr_name in dicts_of_elements_to_add:
                allele_dict = getattr(allele, attr_name)
                insertion_dict = getattr(insertion, attr_name)
                self.log.debug(f'For {attr_name}, add {len(allele_dict)} allele elements to {len(insertion_dict)} insertion elements.')
                for k, v in allele_dict.items():
                    if k not in insertion_dict.keys():
                        insertion_dict[k] = v
            # Combine lists of annotations.
            for attr_name in dicts_of_lists_to_add:
                allele_dict = getattr(allele, attr_name)
                insertion_dict = getattr(insertion, attr_name)
                self.log.debug(f'For {attr_name}, add {len(allele_dict)} allele lists to {len(insertion_dict)} insertion lists.')
                for k, v in allele_dict.items():
                    try:
                        insertion_dict[k].extend(v)
                    except KeyError:
                        insertion_dict[k] = v
        allele.is_obsolete = False
        allele.for_export = False
        allele.export_warnings.append('Superceded by FBti insertion')
        return

    def merge_fbti_fbal(self):
        """Merge FBal allele info into FBti insertion entities as appropriate."""
        self.log.info('Merge FBal allele info into FBti insertion entities as appropriate.')
        prob_counter = 0
        at_locus_counter = 0
        transgenic_counter = 0
        classical_counter = 0
        for allele in self.fb_data_entities.values():
            if allele.db_primary_id in self.at_locus_fbal_fbti_dict.keys() and allele.db_primary_id in self.transgenic_fbal_fbti_dict.keys():
                self.log.error(f'Allele {allele} unexpectedly has both at-locus and transgenic unspecified FBti insertions.')
                prob_counter += 1
            elif allele.db_primary_id in self.at_locus_fbal_fbti_dict.keys():
                allele.superceded_by_at_locus_insertion = self.at_locus_fbal_fbti_dict[allele.db_primary_id][0]
                self.add_fbal_to_fbti(allele)
                at_locus_counter += 1
            elif allele.db_primary_id in self.transgenic_fbal_fbti_dict.keys():
                allele.superceded_by_transgnc_insertions = self.transgenic_fbal_fbti_dict[allele.db_primary_id]
                self.add_fbal_to_fbti(allele)
                transgenic_counter += 1
            else:
                classical_counter += 1
        self.log.info(f'Found {prob_counter} FBal alleles unexpectedly related to both at-locus and unspecified transgenic FBti insertions.')
        self.log.info(f'Found {at_locus_counter} FBal alleles related to at-locus FBti insertions.')
        self.log.info(f'Found {transgenic_counter} FBal alleles related to unspecified transgenic FBti insertions.')
        self.log.info(f'Found {classical_counter} FBal classical and/or complex alleles to be reported as they are (no FBti replacement).')
        return

    def synthesize_related_features(self):
        """Synthesize allele attributes based on related features."""
        self.log.info('Synthesize allele attributes based on related features.')
        has_construct_counter = 0
        has_dmel_insertion_counter = 0
        has_non_dmel_insertion_counter = 0
        has_args_counter = 0
        for allele in self.fb_data_entities.values():
            if allele.uniquename.startswith('FBti'):
                continue
            # Assess relationships to current constructs.
            relevant_cons_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='derived_tp_assoc_alleles',
                                                             rel_entity_types=self.feature_subtypes['construct'])
            self.log.debug(f'For {allele}, found {len(relevant_cons_rels)} cons rels to review.')
            for cons_rel in relevant_cons_rels:
                construct = self.feature_lookup[cons_rel.chado_obj.object_id]
                if construct['is_obsolete'] is False and construct['uniquename'].startswith('FBtp'):
                    allele.cons_rels.append(cons_rel)
                    has_construct_counter += 1
            # Assess relationships to current insertions (will be used for mutation type assessment of complex alleles).
            relevant_ins_rels = allele.recall_relationships(self.log, entity_role='subject', rel_entity_types=self.feature_subtypes['insertion'])
            # self.log.debug(f'For {allele}, found {len(relevant_ins_rels)} ins rels to review.')
            for ins_rel in relevant_ins_rels:
                insertion = self.feature_lookup[ins_rel.chado_obj.object_id]
                if insertion['is_obsolete'] is False and insertion['uniquename'].startswith('FBti'):
                    if self.organism_lookup[insertion['organism_id']]['abbreviation'] == 'Dmel':
                        # self.log.debug(f'{allele} has Dmel insertion.')
                        allele.dmel_ins_rels.append(ins_rel)
                        has_dmel_insertion_counter += 1
                    else:
                        # self.log.debug(f'{allele} has non-Dmel insertion.')
                        allele.non_dmel_ins_rels.append(ins_rel)
                        has_non_dmel_insertion_counter += 1
            # Assess relationships to ARGs.
            relevant_rels = allele.recall_relationships(self.log, entity_role='object', rel_types='partof', rel_entity_types=self.feature_subtypes['variation'])
            # self.log.debug(f'For {allele}, found {len(relevant_rels)} partof relationships to ARGs.')
            for arg_rel in relevant_rels:
                arg = self.feature_lookup[arg_rel.chado_obj.subject_id]
                if arg['is_obsolete'] is False:
                    allele.arg_rels.append(arg_rel)
                    has_args_counter += 1
        self.log.info(f'Found {has_construct_counter} alleles related to a construct.')
        self.log.info(f'Found {has_dmel_insertion_counter} alleles related to a Dmel insertion.')
        self.log.info(f'Found {has_non_dmel_insertion_counter} alleles related to a non-Dmel insertion.')
        self.log.info(f'Found {has_args_counter} alleles related to a mapped variation feature (ARGs).')
        return

    def synthesize_parent_genes(self):
        """Get allele parent gene IDs."""
        self.log.info('Get allele parent gene IDs.')
        allele_counter = 0
        for allele in self.fb_data_entities.values():
            # Skip transgenic alleles since gene relationships will be reported via Construct/Cassette submissions.
            if allele.cons_rels:
                continue
            parent_gene_ids = []
            relevant_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='alleleof', rel_entity_types='gene')
            self.log.debug(f'For {allele}, found {len(relevant_rels)} alleleof relationships to genes.')
            for allele_gene_rel in relevant_rels:
                parent_gene = self.feature_lookup[allele_gene_rel.chado_obj.object_id]
                if parent_gene['is_obsolete'] is False:
                    parent_gene_ids.append(parent_gene['uniquename'])
            if len(parent_gene_ids) == 1:
                allele.parent_gene_id = parent_gene_ids[0]
                allele_counter += 1
            elif len(parent_gene_ids) == 0 and allele.is_obsolete is False and allele.uniquename.startswith('FBal'):
                self.log.warning(f'Current allele {allele} has no parent gene!')
            elif len(parent_gene_ids) > 1 and allele.is_obsolete is False:
                self.log.warning(f'{allele} has many parent genes!')
        self.log.info(f'Found parental gene for {allele_counter} alleles.')
        return

    def flag_alleles_of_internal_genes(self):
        """Flag alleles of internal genes."""
        self.log.info('Flag alleles of internal genes.')
        counter = 0
        for allele in self.fb_data_entities.values():
            if allele.parent_gene_id in self.internal_gene_ids:
                allele.allele_of_internal_gene = True
                counter += 1
        self.log.info(f'Flagged {counter} alleles as related to internal genes.')
        return

    def adjust_allele_organism(self):
        """Adjust organism for classical non-Dmel alleles."""
        self.log.info('Adjust organism for classical non-Dmel alleles.')
        counter = 0
        for allele in self.fb_data_entities.values():
            # Ignore for insertions.
            if allele.uniquename.startswith('FBti'):
                continue
            # SET ALL ALLELES TO HAVE DEFAULT "Dmel" ORGANISM.
            allele.organism_id = 1
            # Skip alleles that exist in Dmel (certainly or most likely).
            if allele.org_abbr == 'Dmel':
                continue
            elif allele.dmel_ins_rels:
                continue
            # Mapped mutations (ARGs) are exclusively Dmel genome mapped features for Dmel alleles.
            elif allele.arg_rels is True:
                continue
            # The assumption here is that transgenic alleles for non-Dmel species are inserted into the Dmel genome.
            # This may not always be true, but there way to ascertain the host species in the database, so we assume this is so.
            elif allele.cons_rels is True:
                continue
            # Find clear evidence that allele is non-Dmel classical allele.
            is_non_dmel_classical = False
            # A non-Dmel insertion is good evidence that the allele exists in a non-Dmel species.
            if allele.non_dmel_ins_rels is True:
                is_non_dmel_classical = True
            # If there is no indication that the allele is transgenic, we assume that it exists in a non-Dmel species.
            elif self.organism_lookup[allele.organism_id]['is_drosophilid'] is True and allele.in_vitro is False:
                is_non_dmel_classical = True
            # For truly non-Dmel alleles, revert organism_id to that of the related allele chado object.
            if is_non_dmel_classical is True:
                allele.organism_id = allele.chado_obj.organism_id
                adj_org_abbr = self.organism_lookup[allele.organism_id]['abbreviation']
                self.log.debug(f'Non-Dmel allele: id={allele.uniquename}, name={allele.name}, org_abbr={adj_org_abbr}')
                counter += 1
        self.log.info(f'Adjusted organism to be "non-Dmel" for {counter} alleles.')
        return

    def synthesize_allele_gene_associations(self):
        """Synthesize allele-to-gene associations."""
        self.log.info('Synthesize allele-to-gene associations.')
        gene_counter = 0
        allele_counter = 0
        # Need to code for the rare possibility that gene-allele is represented by many feature_relationships.
        for allele in self.fb_data_entities.values():
            if allele.superceded_by_at_locus_insertion or allele.superceded_by_transgnc_insertion:
                continue
            # Skip transgenic alleles.
            elif allele.uniquename.startswith('FBti') and allele.name.endswith('unspecified'):
                continue
            relevant_gene_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='alleleof', rel_entity_types='gene')
            if relevant_gene_rels:
                allele_counter += 1
            # self.log.debug(f'For {gene}, found {len(relevant_allele_rels)} allele rels to review.')
            for gene_rel in relevant_gene_rels:
                gene_feature_id = gene_rel.chado_obj.object_id
                gene = self.feature_lookup[gene_feature_id]
                # Suppress allele-gene associations involving non-Drosophilid genes (which are not exported to the Alliance).
                if self.organism_lookup[gene['organism_id']]['is_drosophilid'] is False:
                    continue
                allele_gene_key = (allele.db_primary_id, gene_feature_id)
                try:
                    self.allele_gene_rels[allele_gene_key].append(gene_rel)
                except KeyError:
                    self.allele_gene_rels[allele_gene_key] = [gene_rel]
                    gene_counter += 1
        self.log.info(f'Found {gene_counter} genes for {allele_counter} alleles.')
        return

    # Elaborate on synthesize_info() for the AlleleHandler.
    def synthesize_info(self):
        """Extend the method for the AlleleHandler."""
        super().synthesize_info()
        self.merge_fbti_fbal()
        self.flag_new_additions_and_obsoletes()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_related_features()
        self.synthesize_parent_genes()
        self.flag_alleles_of_internal_genes()
        self.adjust_allele_organism()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_allele_gene_associations()
        return

    # Additional methods to be run by map_fb_data_to_alliance() below.
    def map_insertion_mutation_types(self):
        """Map insertion mutation types."""
        self.log.info('Map insertion mutation types.')
        # Collect evidence for both types of possible values: mobile_element_insertion (SO:0001837) or transgenic_insertion (SO:0001218).
        # Then, pick the best option and report relevant pubs.
        counter = 0
        te_insertion_subtypes = [
            'natTE_isolate',
            'natTE_isolate_named',
            'natTE_partial_named',
            'natTE_sequenced_strain_1',
        ]
        for insertion in self.fb_data_entities.values():
            if insertion.uniquename.startswith('FBal'):
                continue
            mutation_type_curie = None
            te_pub_ids = []
            tp_pub_ids = []
            # First look at producedby relationships. Every FBti should have a single current FBtp or FBte associated in this way.
            inserted_element_rels = insertion.recall_relationships(self.log, entity_role='subject', rel_types='producedby')
            for inserted_element_rel in inserted_element_rels:
                inserted_element_id = inserted_element_rel.chado_obj.object_id
                # Code for data quirks.
                if inserted_element_id not in self.feature_lookup.keys():
                    continue
                inserted_element = self.feature_lookup[inserted_element_id]
                if inserted_element['is_obsolete'] is True:
                    continue
                if inserted_element['uniquename'].startswith('FBte'):
                    mutation_type_curie = 'SO:0001837'    # mobile_element_insertion
                    te_pub_ids.extend(inserted_element_rel.pubs)
                elif inserted_element['uniquename'].startswith('FBtp'):
                    mutation_type_curie = 'SO:0001218'    # transgenic_insertion
                    tp_pub_ids.extend(inserted_element_rel.pubs)
            # Then look at TI_subtype annotations. Opportunity to change transgenic_insertion to mobile_element_insertion.
            mutation_type_annotations = insertion.recall_cvterm_annotations(self.log, cv_names='TI_subtype')
            for mutation_type_annotation in mutation_type_annotations:
                insertion_subtype_term = self.cvterm_lookup[mutation_type_annotation.chado_obj.cvterm_id]['name']
                if insertion_subtype_term in te_insertion_subtypes:
                    mutation_type_curie = 'SO:0001837'    # mobile_element_insertion
                    te_pub_ids.append(mutation_type_annotation.chado_obj.pub_id)
                else:
                    tp_pub_ids.append(mutation_type_annotation.chado_obj.pub_id)
                    if mutation_type_curie is None:
                        mutation_type_curie = 'SO:0001218'    # transgenic_insertion
            # Pick the mutation type and relevant pubs.
            if mutation_type_curie is None:
                if insertion.is_obsolete is False:
                    self.log.error(f'Could not determine mutation_type for {insertion}')
                continue
            elif mutation_type_curie == 'SO:0001837':
                pub_curies = self.lookup_pub_curies(te_pub_ids)
            else:
                pub_curies = self.lookup_pub_curies(tp_pub_ids)
            mutant_type_annotation = agr_datatypes.AlleleMutationTypeSlotAnnotationDTO(mutation_type_curie, pub_curies)
            insertion.linkmldto.allele_mutation_type_dtos.append(mutant_type_annotation.dict_export())
            counter += 1
        self.log.info(f'Mapped {counter} mutation type annotations for insertions.')
        return

    def map_allele_mutation_types(self):
        """Map allele mutation types."""
        self.log.info('Map allele mutation types.')
        mutation_type_conversion = {
            'transposable_element_insertion_site': 'SO:0001218',    # transgenic_insertion (from feature type)
            'insertion_site': 'SO:0001218',                         # transgenic_insertion (from feature type)
            'transposable_element': 'SO:0001837',                   # mobile_element_insertion (from feature type)
        }
        counter = 0
        for allele in self.fb_data_entities.values():
            if allele.uniquename.startswith('FBti'):
                continue
            mutation_types = {}    # Will be a dict of mutation type curies and supporting pub ids.
            relevant_ins_rels = []
            relevant_ins_rels.extend(allele.dmel_ins_rels)
            relevant_ins_rels.extend(allele.non_dmel_ins_rels)
            for arg_rel in allele.arg_rels:
                mutation_type_curie = None
                fb_feat_type_id = arg_rel.chado_obj.subject.type_id
                fb_feat_type_name = self.cvterm_lookup[fb_feat_type_id]['name']
                if fb_feat_type_id in self.allele_mutant_type_terms:
                    mutation_type_curie = self.cvterm_lookup[fb_feat_type_id]['curie']
                if mutation_type_curie is None:
                    continue
                if mutation_type_curie in mutation_types.keys():
                    mutation_types[mutation_type_curie].extend(arg_rel.pubs)
                else:
                    mutation_types[mutation_type_curie] = arg_rel.pubs
            for ins_rel in relevant_ins_rels:
                mutation_type_curie = None
                fb_feat_type_id = ins_rel.chado_obj.object.type_id
                fb_feat_type_name = self.cvterm_lookup[fb_feat_type_id]['name']
                if fb_feat_type_name in mutation_type_conversion.keys():
                    mutation_type_curie = mutation_type_conversion[fb_feat_type_name]
                if mutation_type_curie is None:
                    continue
                if mutation_type_curie in mutation_types.keys():
                    mutation_types[mutation_type_curie].extend(ins_rel.pubs)
                else:
                    mutation_types[mutation_type_curie] = ins_rel.pubs
            for mutation_type_curie, pub_ids in mutation_types.items():
                pub_curies = self.lookup_pub_curies(pub_ids)
                mutant_type_annotation = agr_datatypes.AlleleMutationTypeSlotAnnotationDTO(mutation_type_curie, pub_curies)
                allele.linkmldto.allele_mutation_type_dtos.append(mutant_type_annotation.dict_export())
                counter += 1
        self.log.info(f'Mapped {counter} mutation type annotations.')
        return

    def map_inheritance_modes(self):
        """Map allele inheritance modes."""
        self.log.info('Map allele inheritance modes.')
        INHERITANCE_MODE_NAME = 0
        PHENOTYPE_TERM_CURIE = 1
        PHENOTYPE_STATEMENT = 2
        for allele in self.fb_data_entities.values():
            # For each allele, gather phenotype data relevant to inheritance mode as a dict.
            # Keys will be tuple of (inheritance_mode, phenotype_term_curie, phenotype_statement)
            # Value for each key will be a list of pub_ids in support of the annotation.
            inheritance_data = {}
            for phenstmt in allele.phenstatements:
                genotype_uniquename = phenstmt.Genotype.uniquename
                genotype_description = phenstmt.Genotype.description
                fb_inheritance_mode_name = phenstmt.Cvterm.name
                inheritance_mode_name = self.inheritance_mode_terms[fb_inheritance_mode_name]
                phenotype_term_curie = self.cvterm_lookup[phenstmt.Phenotype.cvalue_id]['curie']
                phenotype_statement = phenstmt.Phenotype.uniquename
                pheno_key = (inheritance_mode_name, phenotype_term_curie, phenotype_statement)
                anno_desc = f'genotype={genotype_uniquename} ({genotype_description}); '
                anno_desc += f'inheritance={inheritance_mode_name}; '
                anno_desc += f'pheno_curie={phenotype_term_curie}; '
                anno_desc += f'pheno_stmt={phenotype_statement}.'
                # self.log.debug(f'Assess this phenotype annotation: {anno_desc}')
                # Filter out single locus genotypes having many distinct alleles.
                single_allele_genotype = False
                if '|' not in genotype_description:
                    # self.log.debug('GENOTYPE HAS ONLY ONE ALLELE LISTED, UNSPECIFIED ZYGOSITY.')
                    single_allele_genotype = True
                elif genotype_description.startswith('FBab'):
                    # self.log.debug('GENOTYPE IS HEMIZYGOUS: ALLELE/DEFICIENCY.')
                    single_allele_genotype = True
                else:
                    features = genotype_description.split('|')
                    if features[0] == features[1]:
                        # self.log.debug('GENOTYPE IS HOMOZYGOUS.')
                        single_allele_genotype = True
                    else:
                        for feature in features:
                            if feature == '+':
                                single_allele_genotype = True
                            elif feature.endswith('[+]'):
                                single_allele_genotype = True
                        # if single_allele_genotype is True:
                            # self.log.debug('GENOTYPE IS HOMOZYGOUS.')
                if single_allele_genotype is False:
                    # self.log.debug('GENOTYPE IS HETEROZYGOUS FOR TWO NON-WT ALLELES - SKIP IT.')
                    continue
                try:
                    inheritance_data[pheno_key].append(phenstmt.Pub.pub_id)
                except KeyError:
                    inheritance_data[pheno_key] = [phenstmt.Pub.pub_id]
            # Convert data into Alliance slot annotations.
            for pheno_key, pub_ids in inheritance_data.items():
                pub_curies = self.lookup_pub_curies(pub_ids)
                agr_anno_dto = agr_datatypes.AlleleInheritanceModeSlotAnnotationDTO(pheno_key[INHERITANCE_MODE_NAME], pheno_key[PHENOTYPE_TERM_CURIE],
                                                                                    pheno_key[PHENOTYPE_STATEMENT], pub_curies)
                allele.linkmldto.allele_inheritance_mode_dtos.append(agr_anno_dto.dict_export())
        return

    def map_allele_gene_associations(self):
        """Map gene-allele associations to Alliance object."""
        self.log.info('Map gene-allele associations to Alliance object.')
        ALLELE = 0
        GENE = 1
        counter = 0
        for allele_gene_key, allele_gene_rels in self.allele_gene_rels.items():
            allele = self.fb_data_entities[allele_gene_key[ALLELE]]
            allele_curie = f'FB:{allele.uniquename}'
            gene = self.feature_lookup[allele_gene_key[GENE]]
            gene_curie = f'FB:{gene["uniquename"]}'
            first_feat_rel = allele_gene_rels[0]
            all_pub_ids = []
            for allele_gene_rel in allele_gene_rels:
                all_pub_ids.extend(allele_gene_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            pub_curies = self.lookup_pub_curies(all_pub_ids)
            rel_dto = agr_datatypes.AlleleGeneAssociationDTO(allele_curie, 'is_allele_of', gene_curie, pub_curies)
            if allele.is_obsolete is True or gene['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            first_feat_rel.linkmldto = rel_dto
            self.allele_gene_associations.append(first_feat_rel)
            counter += 1
        self.log.info(f'Generated {counter} allele-gene unique associations.')
        return

    # Elaborate on map_fb_data_to_alliance() for the AlleleHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AlleleHandler."""
        super().map_fb_data_to_alliance()
        self.map_metaallele_basic()
        self.map_metaallele_database_status()
        self.map_internal_metaallele_status()
        self.map_insertion_mutation_types()
        self.map_allele_mutation_types()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_extinction_info()
        self.map_collections()
        self.map_inheritance_modes()
        self.map_pubs()    # Suppress if load times are slow.
        self.map_timestamps()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        self.map_allele_gene_associations()
        self.flag_internal_fb_entities('allele_gene_associations')
        return

    # Elaborate on query_chado_and_export() for the AlleleHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the AlleleHandler."""
        super().query_chado_and_export(session)
        self.flag_unexportable_entities(self.allele_gene_associations, 'allele_gene_association_ingest_set')
        self.generate_export_dict(self.allele_gene_associations, 'allele_gene_association_ingest_set')
        return


class InsertionHandler(MetaAlleleHandler):
    """This object gets, synthesizes and filters insertion data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the InsertionHandler object."""
        super().__init__(log, testing)
        self.datatype = 'insertion'
        self.fb_export_type = FBAllele

    # Types: 228747 transposable_element_insertion_site; 7726 insertion_site; 5753 transposable element; 3573 match (internal).
    # Relationships: 234754 FBti(producedby)FBtp; 64920 FBal(associated_with)FBti.
    test_set = {
        'FBti0000040': 'P{hsneo}Xrp1[142]',         # type=transposable_element_insertion_site. Location trap. FBal-(associated_with)->FBti-(producedby)->FBtp
        'FBti0151770': 'P{UAS-stnB.M}vl',           # type=transposable_element_insertion_site. FBti-(producedby)->FBtp<-(associated_with)-FBal.
        'FBti0167947': 'TI{TI}wg[GFP]',             # type=insertion_site. FBti-(producedby)->FBtp<-(associated_with)-FBal.
        'FBti0018862': '17.6{}804',                 # type=17.6{}804; this insertion shares its uniquename with two internal "match" features.
        'FBti0016979': 'P{PZ}Vha44[06072b]',        # type=transposable_element_insertion_site. Direct "associated_with" association to a gene.
        'FBti0186374': 'P{TOE.GS00088}attP40',      # type=transposable_element_insertion_site. Related to TRiP-OE-VPR collection via FBtp0116301.
        'FBti0178263': 'TI{TI}Rab1[EYFP]',          # type=insertion_site. Related to YRab collection via FBal0314192.
        'FBti0164639': 'P{TRiP.HMJ22303}attP40',    # type=transposable_element_insertion_site. Related to TRiP-3 collection via FBtp0097015-FBsf0000443916.
        'FBti0009227': 'P{PZ}Ubx[Ubx-Plac61]',      # type=transposable_element_insertion_site. Lacks TI_subtype, producedby FBtp0000210 (P{PZ}, no pubs).
        'FBti0168248': 'TI{TI}sano[KO1]',           # type=insertion_site. Lacks TI_subtype, producedby FBtp0099201 (TI{TI}).
        'FBti0074148': 'blood{}cl[1]',              # type=transposable_element_insertion_site. Lacks TI_subtype, producedby FBte0000279 (blood).
        'FBti0018906': 'Doc2{}650',                 # type=transposable_element. Has synTE_insertion TI_subtype, producedby FBte0000103 (Doc2-element).
        'FBti0186554': 'gypsy{5}y-TDmh1',           # type=transposable_element_insertion_site, natTE_partial_named TI_subtype, producedby FBtp0012975 gypsy{5'}
        'FBti0000005': 'P{hsneo}102',               # type=transposable_element_insertion_site. Has synTE_insertion subtype, producedby FBtp0000078 (P{hsneo}).
        'FBti0145331': 'P{PTGAL}26',                # type=transposable_element_insertion_site. Has synTE_insertion TI_subtype, no current producedby.
        'FBti0128384': 'FBti0128384',               # Had null mutation_type_curies bug.
    }

    # Elaborate on get_general_data() for the InsertionHandler.
    def get_general_data(self, session):
        """Suppress the method for the InsertionHandler."""
        self.log.info('DO NOT GET FLYBASE INSERTION DATA FROM CHADO via InsertionHandler; use AlleleHandler.')
        return

    # Elaborate on get_datatype_data() for the InsertionHandler.
    def get_datatype_data(self, session):
        """Extend the method for the InsertionHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject')
        self.get_entity_relationships(session, 'object')
        self.get_entity_cvterms(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_direct_reagent_collections(session)
        self.get_indirect_reagent_collections(session, 'object', 'associated_with', 'allele')
        # self.get_indirect_reagent_collections(session, 'subject', 'producedby', 'construct')
        # self.get_very_indirect_reagent_collections(session)    # Suppressed because it's slow and perhaps too indirect.
        return

    # Elaborate on synthesize_info() for the InsertionHandler.
    def synthesize_info(self):
        """Suppress the method for the InsertionHandler."""
        self.log.info('DO NOT SYNTHESIZE FLYBASE ALLELE DATA FROM CHADO via the InsertionHandler; use AlleleHandler.')
        return

    # Elaborate on map_fb_data_to_alliance() for the InsertionHandler.
    def map_fb_data_to_alliance(self):
        """Suppress the method for the InsertionHandler."""
        self.log.info('DO NOT MAP FLYBASE DATA TO ALLIANCE DATA via the InsertionHandler; use AlleleHandler.')
        return

    # Elaborate on query_chado_and_export() for the InsertionHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the InsertionHandler."""
        super().query_chado_and_export(session)
        return


class AberrationHandler(MetaAlleleHandler):
    """This object gets, synthesizes and filters aberration data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AberrationHandler object."""
        super().__init__(log, testing)
        self.datatype = 'aberration'
        self.fb_export_type = FBAberration

    test_set = {
        'FBab0000001': 'Df(2R)03072',           # Random selection.
        'FBab0000006': 'Df(3L)ZN47',            # Has many genes associated in many ways.
        'FBab0024587': 'Dp(1;f)8D',             # Unusual feature type: "free duplication".
        'FBab0005448': 'In(3LR)P88',            # Many distinct "wt_aberr" type CV term annotations.
        'FBab0038557': 'Dmau_Int(3)46.22',      # Unusual annotation: introgressed_chromosome_region (SO:0000664). Assign 'NCBITaxon:32644' (unidentified).
        'FBab0038658': 'Dsim_Int(2L)S',         # Unusual annotation: introgressed_chromosome_region (SO:0000664). Assign 'NCBITaxon:32644' (unidentified).
        'FBab0047489': 'Dp(3;3)NA18',           # Unusual annotation: direct_tandem_duplication (SO:1000039).
        'FBab0010504': 'T(2;3)G16DTE35B-3P',    # Unusual annotation: assortment_derived_deficiency_plus_duplication (SO:0000801).
    }

    # Additional export sets.
    aberration_gene_rels = {}            # Will be (aberration feature_id, gene feature_id, fr type_id) tuples keying lists of FBRelationships.
    aberration_gene_associations = []    # Will be the final list of gene-aberration FBRelationships to export (AlleleGeneAssociationDTO under linkmldto attr).

    # Additional reference info.
    chr_str_var_terms = []    # A list of cvterm_ids for child terms of "chromosome_structure_variation" (SO:0000240).
    seq_alt_terms = []        # A list of cvterm_ids for child terms of "sequence_alteration" (SO:0001059).
    str_var_terms = []        # A list of cvterm_ids for child terms of "structural_variant" (SO:0001537).

    # Additional sub-methods for get_general_data().
    def get_key_cvterm_sets_for_aberrations(self, session):
        """Get key CV term sets for aberrations from chado."""
        self.log.info('Get key CV term sets for aberrations from chado.')
        self.chr_str_var_terms.extend(self.get_child_cvterms(session, 'chromosome_structure_variation', 'SO'))
        self.seq_alt_terms.extend(self.get_child_cvterms(session, 'sequence_alteration', 'SO'))
        self.str_var_terms.extend(self.get_child_cvterms(session, 'structural_variant', 'SO'))
        return

    # Elaborate on get_general_data() for the AberrationHandler.
    def get_general_data(self, session):
        """Extend the method for the AberrationHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['gene'])
        self.build_allele_gene_lookup(session)
        self.get_key_cvterm_sets_for_aberrations(session)
        return

    # Additional sub-methods for get_datatype_data().
    # Placeholder.

    # Elaborate on get_datatype_data() for the AberrationHandler.
    def get_datatype_data(self, session):
        """Extend the method for the AberrationHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject')
        self.get_entity_relationships(session, 'object')
        self.get_entity_cvterms(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_direct_reagent_collections(session)
        return

    # Additional sub-methods to be run by synthesize_info() below.
    def adjust_aberration_organism(self):
        """Adjust organism for aberrations, if needed."""
        self.log.info('Adjust organism for aberrations, if needed.')
        counter = 0
        for aberration in self.fb_data_entities.values():
            if 'introgressed_chromosome_region' in aberration.cvt_anno_ids_by_term.keys():
                aberration.organism_id = 509   # Set organism_id corresponding to "Unknown".
                counter += 1
        self.log.info(f'Adjusted organism to be "unspecified" for {counter} aberrations.')
        return

    def qc_aberration_mutation_types(self):
        """Run QC checks on aberration mutation type annotations."""
        self.log.info('Run QC checks on aberration mutation type annotations.')
        prop_types = ['wt_class', 'aberr_class']
        for aberration in self.fb_data_entities.values():
            # self.log.debug(f'Assess {aberration} mutation types.')
            for prop_type_name in prop_types:
                annotations = aberration.recall_cvterm_annotations(self.log, prop_type_names=prop_type_name)
                # self.log.debug(f'Retrieved {len(annotations)} annotations with prop_type_name={prop_type_name}')
                for annotation in annotations:
                    cvterm_id = annotation.chado_obj.cvterm_id
                    if cvterm_id not in self.chr_str_var_terms:
                        cvterm_name = self.cvterm_lookup[cvterm_id]['name']
                        cvterm_curie = self.cvterm_lookup[cvterm_id]['curie']
                        pub_curie = annotation.chado_obj.pub.uniquename
                        in_seq_alt = False
                        if cvterm_id in self.seq_alt_terms:
                            in_seq_alt = True
                        in_str_var = False
                        if cvterm_id in self.str_var_terms:
                            in_str_var = True
                        report_str = f'{aberration.uniquename}\t{aberration.name}\t{prop_type_name}'
                        report_str += f'\tin_seq_alt={in_seq_alt}\tin_str_var={in_str_var}'
                        report_str += f'\t{cvterm_name}\t{cvterm_curie}\t{pub_curie}'
                        self.log.debug(f'QCCHK:{report_str}')
        return

    def synthesize_aberration_gene_associations(self):
        """Synthesize aberration-to-gene associations."""
        self.log.info('Synthesize aberration-to-gene associations.')
        aberration_counter = 0
        gene_rel_counter = 0
        for aberration in self.fb_data_entities.values():
            relevant_sbj_gene_rels = aberration.recall_relationships(self.log, entity_role='subject', rel_entity_types='gene')
            relevant_obj_gene_rels = aberration.recall_relationships(self.log, entity_role='object', rel_entity_types='gene')
            if relevant_sbj_gene_rels or relevant_obj_gene_rels:
                aberration_counter += 1
            for sbj_gene_rel in relevant_sbj_gene_rels:
                rel_key = (aberration.db_primary_id, sbj_gene_rel.chado_obj.object_id, sbj_gene_rel.chado_obj.type_id)
                try:
                    self.aberration_gene_rels[rel_key].append(sbj_gene_rel)
                except KeyError:
                    self.aberration_gene_rels[rel_key] = [sbj_gene_rel]
                    gene_rel_counter += 1
            for obj_gene_rel in relevant_obj_gene_rels:
                rel_key = (aberration.db_primary_id, obj_gene_rel.chado_obj.subject_id, obj_gene_rel.chado_obj.type_id)
                try:
                    self.aberration_gene_rels[rel_key].append(obj_gene_rel)
                except KeyError:
                    self.aberration_gene_rels[rel_key] = [obj_gene_rel]
                    gene_rel_counter += 1
        self.log.info(f'Found {gene_rel_counter} aberration-gene relationships for {aberration_counter} aberrations.')
        return

    # Elaborate on synthesize_info() for the AberrationHandler.
    def synthesize_info(self):
        """Extend the method for the AberrationHandler."""
        super().synthesize_info()
        self.flag_new_additions_and_obsoletes()
        self.adjust_aberration_organism()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_aberration_gene_associations()
        self.synthesize_ncbi_taxon_id()
        self.qc_aberration_mutation_types()
        return

    # Additional methods to be run by map_fb_data_to_alliance() below.
    def map_aberration_mutation_types(self):
        """Map aberration mutation types."""
        self.log.info('Map aberration mutation types.')
        counter = 0
        for aberration in self.fb_data_entities.values():
            mutation_types = {}    # curie-keyed dict of pub_ids.
            mutation_type_annotations = aberration.recall_cvterm_annotations(self.log, cv_names='SO', prop_type_names='wt_class')
            for mutation_type_annotation in mutation_type_annotations:
                mutation_type_curie = self.cvterm_lookup[mutation_type_annotation.chado_obj.cvterm_id]['curie']
                if mutation_type_curie in mutation_types.keys():
                    mutation_types[mutation_type_curie].append(mutation_type_annotation.chado_obj.pub_id)
                else:
                    mutation_types[mutation_type_curie] = [mutation_type_annotation.chado_obj.pub_id]
            for mutation_type_curie, pub_ids in mutation_types.items():
                pub_curies = self.lookup_pub_curies(pub_ids)
                mutant_type_annotation = agr_datatypes.AlleleMutationTypeSlotAnnotationDTO(mutation_type_curie, pub_curies)
                aberration.linkmldto.allele_mutation_type_dtos.append(mutant_type_annotation.dict_export())
                counter += 1
        self.log.info(f'Mapped {counter} mutation type annotations for aberrations.')
        return

    def map_aberration_gene_associations(self):
        """Map aberration-gene associations to Alliance object."""
        self.log.info('Map aberration-gene associations to Alliance object.')
        ABERRATION = 0
        GENE = 1
        REL_TYPE = 2
        fb_agr_aberr_rel_mapping = {
            'deletes': 'full_deletion',
            'molec_deletes': 'full_deletion',
            'part_deletes': 'partial_deletion',
            'molec_partdeletes': 'partial_deletion',
            'nondeletes': 'mutation_does_not_delete',
            'molec_nondeletes': 'mutation_does_not_delete',
            'duplicates': 'full_duplication',
            'molec_dups': 'full_duplication',
            'part_duplicates': 'partial_duplication',
            'molec_partdups': 'partial_duplication',
            'nonduplicates': 'mutation_does_not_duplicate',
            'molec_nondups': 'mutation_does_not_duplicate',
        }
        counter = 0
        for rel_key, aberration_gene_rels in self.aberration_gene_rels.items():
            aberration = self.fb_data_entities[rel_key[ABERRATION]]
            aberration_curie = f'FB:{aberration.uniquename}'
            gene = self.feature_lookup[rel_key[GENE]]
            gene_curie = f'FB:{gene["uniquename"]}'
            fb_rel_type_name = self.cvterm_lookup[rel_key[REL_TYPE]]['name']
            if fb_rel_type_name not in fb_agr_aberr_rel_mapping.keys():
                continue
            agr_rel_type_name = fb_agr_aberr_rel_mapping[fb_rel_type_name]
            first_feat_rel = aberration_gene_rels[0]
            all_pub_ids = []
            for aberration_gene_rel in aberration_gene_rels:
                all_pub_ids.extend(aberration_gene_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            pub_curies = self.lookup_pub_curies(all_pub_ids)
            rel_dto = agr_datatypes.AlleleGeneAssociationDTO(aberration_curie, agr_rel_type_name, gene_curie, pub_curies)
            if aberration.is_obsolete is True or gene['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            first_feat_rel.linkmldto = rel_dto
            self.aberration_gene_associations.append(first_feat_rel)
            counter += 1
        self.log.info(f'Generated {counter} aberration-gene unique associations.')
        return

    # Elaborate on map_fb_data_to_alliance() for the AberrationHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AberrationHandler."""
        super().map_fb_data_to_alliance()
        self.map_metaallele_basic()
        self.map_metaallele_database_status()
        self.map_internal_metaallele_status()
        self.map_aberration_mutation_types()
        self.map_aberration_gene_associations()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_extinction_info()
        self.map_collections()
        self.map_pubs()    # Suppress if load times are slow.
        self.map_timestamps()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        # self.flag_internal_fb_entities('aberration_gene_associations')
        return

    # Elaborate on query_chado_and_export() for the AberrationHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the AberrationHandler."""
        super().query_chado_and_export(session)
        # self.flag_unexportable_entities(self.aberration_gene_associations, 'allele_gene_association_ingest_set')
        # self.generate_export_dict(self.aberration_gene_associations, 'allele_gene_association_ingest_set')
        return


class BalancerHandler(MetaAlleleHandler):
    """This object gets, synthesizes and filters balancer data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the BalancerHandler object."""
        super().__init__(log, testing)
        self.datatype = 'balancer'
        self.fb_export_type = FBBalancer

    test_set = {
        'FBba0000001': 'CyO-Df(2R)B80',    # A random selection.
    }

    # Additional sub-methods for get_general_data().
    # Placeholder.

    # Elaborate on get_general_data() for the BalancerHandler.
    def get_general_data(self, session):
        """Extend the method for the BalancerHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        return

    # Additional sub-methods for get_datatype_data().
    # Placeholder.

    # Elaborate on get_datatype_data() for the BalancerHandler.
    def get_datatype_data(self, session):
        """Extend the method for the BalancerHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject')
        self.get_entity_relationships(session, 'object')
        self.get_entity_cvterms(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_direct_reagent_collections(session)
        return

    # Additional sub-methods to be run by synthesize_info() below.
    # Placeholder.

    # Elaborate on synthesize_info() for the BalancerHandler.
    def synthesize_info(self):
        """Extend the method for the BalancerHandler."""
        super().synthesize_info()
        self.flag_new_additions_and_obsoletes()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_ncbi_taxon_id()
        return

    # Additional methods to be run by map_fb_data_to_alliance() below.
    def map_balancer_mutation_types(self):
        """Map balancer mutation types."""
        self.log.info('Map balancer mutation types.')
        counter = 0
        for balancer in self.fb_data_entities.values():
            mutation_type_curie = 'SO:1000183'    # chromosome_structure_variation
            pub_curies = []
            mutant_type_annotation = agr_datatypes.AlleleMutationTypeSlotAnnotationDTO(mutation_type_curie, pub_curies)
            balancer.linkmldto.allele_mutation_type_dtos.append(mutant_type_annotation.dict_export())
            counter += 1
        self.log.info(f'Mapped {counter} mutation type annotations for balancers.')
        return

    # Elaborate on map_fb_data_to_alliance() for the BalancerHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the BalancerHandler."""
        super().map_fb_data_to_alliance()
        self.map_metaallele_basic()
        self.map_metaallele_database_status()
        self.map_internal_metaallele_status()
        self.map_balancer_mutation_types()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_extinction_info()
        self.map_collections()
        self.map_pubs()    # Suppress if load times are slow.
        self.map_timestamps()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Elaborate on query_chado_and_export() for the BalancerHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the BalancerHandler."""
        super().query_chado_and_export(session)
        return
