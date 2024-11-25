"""Module:: allele_handler.

Synopsis:
    A data handler that exports FlyBase data for alleles to Alliance Gene LinkML
    objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
import agr_datatypes
from fb_datatypes import FBAllele
from feature_handler import FeatureHandler
from harvdev_utils.production import (
    Cvterm, Feature, FeatureCvterm, FeatureGenotype, FeatureRelationship,
    Genotype, Library, LibraryFeature, LibraryFeatureprop, Phenotype,
    PhenotypeCvterm, Phenstatement, Pub
)


class AlleleHandler(FeatureHandler):
    """This object gets, synthesizes and filters allele data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AlleleHandler object."""
        super().__init__(log, testing)
        self.datatype = 'allele'
        self.fb_export_type = FBAllele
        self.agr_export_type = agr_datatypes.AlleleDTO
        self.primary_export_set = 'allele_ingest_set'

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
    }

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
        self.get_key_cvterm_sets(session)
        self.build_ncbi_taxon_lookup(session)
        self.get_drosophilid_organisms(session)
        self.build_feature_lookup(session, feature_types=['construct', 'gene', 'insertion', 'variation'])
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

    def get_direct_collections(self, session):
        """Find collections directly related to alleles."""
        self.log.info('Find collections directly related to alleles.')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        filters = (
            Feature.uniquename.op('~')(self.regex['allele']),
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
        self.log.info(f'Found {counter} direct allele-collection associations.')
        return

    def get_indirect_collections(self, session):
        """Find collections indirectly related to alleles via insertions or constructs."""
        self.log.info('Find collections indirectly related to alleles via insertions or constructs.')
        allele = aliased(Feature, name='allele')
        feature = aliased(Feature, name='feature')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        featreltype = aliased(Cvterm, name='featreltype')
        filters = (
            allele.uniquename.op('~')(self.regex['allele']),
            feature.uniquename.op('~')(self.regex['consins']),
            feature.is_obsolete.is_(False),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.regex['library']),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection',
            featreltype.name == 'associated_with'
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (allele.uniquename.in_((self.test_set.keys())), )
        indirect_collections = session.query(allele, feature, Library).\
            select_from(allele).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == allele.feature_id)).\
            join(featreltype, (featreltype.cvterm_id == FeatureRelationship.type_id)).\
            join(feature, (feature.feature_id == FeatureRelationship.object_id)).\
            join(LibraryFeature, (LibraryFeature.feature_id == feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        fbti_counter = 0
        fbtp_counter = 0
        for result in indirect_collections:
            if result.feature.uniquename.startswith('FBti'):
                self.fb_data_entities[result.allele.feature_id].ins_colls.append(result.Library)
                fbti_counter += 1
            elif result.feature.uniquename.startswith('FBtp'):
                self.fb_data_entities[result.allele.feature_id].cons_colls.append(result.Library)
                fbtp_counter += 1
        self.log.info(f'Found {fbti_counter} insertion-mediated allele-collection associations.')
        self.log.info(f'Found {fbtp_counter} construct-mediated allele-collection associations.')
        return

    def get_sf_collections(self, session):
        """Find collections indirectly related to alleles via sequence features."""
        self.log.info('Find collections indirectly related to alleles via sequence features.')
        allele = aliased(Feature, name='allele')
        construct = aliased(Feature, name='construct')
        seqfeat = aliased(Feature, name='seqfeat')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        featreltype = aliased(Cvterm, name='featreltype')
        allele_construct = aliased(FeatureRelationship, name='allele_construct')
        seqfeat_construct = aliased(FeatureRelationship, name='seqfeat_construct')
        filters = (
            allele.uniquename.op('~')(self.regex['allele']),
            construct.uniquename.op('~')(self.regex['construct']),
            seqfeat.uniquename.op('~')(self.regex['seqfeat']),
            construct.is_obsolete.is_(False),
            seqfeat.is_obsolete.is_(False),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.regex['library']),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection',
            featreltype.name == 'associated_with'
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (allele.uniquename.in_((self.test_set.keys())), )
        sf_collections = session.query(allele, Library).\
            select_from(allele).\
            join(allele_construct, (allele_construct.subject_id == allele.feature_id)).\
            join(construct, (construct.feature_id == allele_construct.object_id)).\
            join(featreltype, (featreltype.cvterm_id == allele_construct.type_id)).\
            join(seqfeat_construct, (seqfeat_construct.object_id == construct.feature_id)).\
            join(seqfeat, (seqfeat.feature_id == seqfeat_construct.subject_id)).\
            join(LibraryFeature, (LibraryFeature.feature_id == seqfeat.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in sf_collections:
            self.fb_data_entities[result.allele.feature_id].sf_colls.append(result.Library)
            counter += 1
        self.log.info(f'Found {counter} sequence feature-mediated allele-library associations.')
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

    # Elaborate on get_datatype_data() for the AlleleHandler.
    def get_datatype_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entity_relationships(session, 'subject', rel_type='alleleof', entity_type='gene', entity_regex=self.regex['gene'])
        self.get_entity_relationships(session, 'subject', rel_type='derived_tp_assoc_alleles', entity_type='construct', entity_regex=self.regex['construct'])
        self.get_entity_relationships(session, 'subject', rel_type='associated_with', entity_type='insertion', entity_regex=self.regex['insertion'])
        self.get_entity_relationships(session, 'object', rel_type='partof', entity_type='variation')
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_phenotypes(session)
        self.find_in_vitro_alleles(session)
        self.get_direct_collections(session)
        self.get_indirect_collections(session)
        self.get_sf_collections(session)
        return

    # Add sub-methods to be run by synthesize_info() below.
    def synthesize_parent_genes(self):
        """Get allele parent gene IDs."""
        self.log.info('Get allele parent gene IDs.')
        allele_counter = 0
        for allele in self.fb_data_entities.values():
            parent_gene_ids = []
            relevant_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='alleleof', rel_entity_types='gene')
            # self.log.debug(f'For {allele}, found {len(relevant_rels)} alleleof relationships to genes.')
            for allele_gene_rel in relevant_rels:
                parent_gene = self.feature_lookup[allele_gene_rel.chado_obj.object_id]
                if parent_gene['is_obsolete'] is False:
                    parent_gene_ids.append(parent_gene['uniquename'])
            if len(parent_gene_ids) == 1:
                allele.parent_gene_id = parent_gene_ids[0]
                allele_counter += 1
            elif len(parent_gene_ids) == 0 and allele.is_obsolete is False:
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

    def synthesize_related_features(self):
        """Synthesize allele attributes based on related features."""
        self.log.info('Synthesize allele attributes based on related features.')
        has_construct_counter = 0
        has_dmel_insertion_counter = 0
        has_non_dmel_insertion_counter = 0
        has_args_counter = 0
        for allele in self.fb_data_entities.values():
            # Assess relationships to current constructs.
            relevant_cons_rels = allele.recall_relationships(self.log, entity_role='subject', rel_types='derived_tp_assoc_alleles',
                                                             rel_entity_types=self.feature_subtypes['construct'])
            # self.log.debug(f'For {allele}, found {len(relevant_cons_rels)} cons rels to review.')
            for cons_rel in relevant_cons_rels:
                construct = self.feature_lookup[cons_rel.chado_obj.object_id]
                if construct['is_obsolete'] is False and construct['uniquename'].startswith('FBtp'):
                    allele.cons_rels.append(cons_rel)
                    has_construct_counter += 1
            # Assess relationships to current insertions.
            relevant_ins_rels = allele.recall_relationships(self.log, entity_role='subject', rel_entity_types=self.feature_subtypes['insertion'])
            # self.log.debug(f'For {allele}, found {len(relevant_ins_rels)} ins rels to review.')
            for ins_rel in relevant_ins_rels:
                insertion = self.feature_lookup[ins_rel.chado_obj.object_id]
                if insertion['is_obsolete'] is False and insertion['uniquename'].startswith('FBti'):
                    if insertion['species'] == 'Drosophila melanogaster':
                        allele.dmel_ins_rels.append(ins_rel)
                        has_dmel_insertion_counter += 1
                    else:
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

    def adjust_allele_organism(self):
        """Adjust organism for classical non-Dmel alleles."""
        self.log.info('Adjust organism for classical non-Dmel alleles.')
        # The default allele.adj_org_abbr is Dmel. Change this if needed.
        counter = 0
        for allele in self.fb_data_entities.values():
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
            elif allele.organism_id in self.drosophilid_list and allele.in_vitro is False:
                is_non_dmel_classical = True
            # Make the adjustment.
            if is_non_dmel_classical is True:
                allele.adj_org_abbr = allele.org_abbr
                self.log.debug(f'Non-Dmel allele: id={allele.uniquename}, name={allele.name}, org_abbr={allele.org_abbr}')
                counter += 1
        self.log.info(f'Adjusted organism to be "non-Dmel" for {counter} alleles.')
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
        self.flag_alleles_of_internal_genes()
        self.synthesize_related_features()
        self.adjust_allele_organism()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_allele_basic(self):
        """Map basic FlyBase allele data to the Alliance LinkML object."""
        self.log.info('Map basic allele info to Alliance object.')
        for allele in self.fb_data_entities.values():
            agr_allele = self.agr_export_type()
            agr_allele.obsolete = allele.chado_obj.is_obsolete
            agr_allele.mod_entity_id = f'FB:{allele.uniquename}'
            agr_allele.mod_internal_id = str(allele.chado_obj.feature_id)
            agr_allele.taxon_curie = allele.ncbi_taxon_id
            allele.linkmldto = agr_allele
        return

    def map_allele_database_status(self):
        """Map allele database status."""
        self.log.info('Map allele database status.')
        evidence_curies = []
        for allele in self.fb_data_entities.values():
            if allele.is_obsolete is False:
                db_status = 'approved'
            else:
                db_status = 'deleted'
            db_status_annotation = agr_datatypes.AlleleDatabaseStatusSlotAnnotationDTO(db_status, evidence_curies)
            allele.linkmldto.allele_database_status_dto = db_status_annotation.dict_export()
        return

    def map_collections(self):
        """Map allele collections."""
        self.log.info('Map allele collections.')
        for allele in self.fb_data_entities.values():
            collections = []
            if allele.direct_colls:
                collections.extend(allele.direct_colls)
            elif allele.ins_colls:
                collections.extend(allele.ins_colls)
            elif allele.cons_colls:
                collections.extend(allele.cons_colls)
            elif allele.sf_colls:
                collections.extend(allele.sf_colls)
            if collections:
                collections = list(set(collections))
                allele.linkmldto.in_collection_name = collections[0].name
                if len(collections) > 1:
                    self.log.warning(f'{allele} has many relevant collections: {[i.name for i in collections]}')
        return

    def map_extinction_info(self):
        """Map extinction info."""
        self.log.info('Map extinction info.')
        counter = 0
        for allele in self.fb_data_entities.values():
            if 'availability' in allele.props_by_type.keys():
                for prop in allele.props_by_type['availability']:
                    if prop.chado_obj.value == 'Stated to be lost.':
                        allele.linkmldto.is_extinct = True
            for prop_type in allele.props_by_type.keys():
                if prop_type.startswith('derived_stock'):
                    # Stock availability trumps curated extinction comment.
                    if allele.linkmldto.is_extinct is True:
                        self.log.warning(f'{allele} is stated to be lost but has stocks.')
                    allele.linkmldto.is_extinct = False
            if allele.linkmldto.is_extinct is True:
                counter += 1
        self.log.info(f'Flagged {counter} alleles as extinct.')
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

    def map_internal_allele_status(self):
        """Flag internal alleles using allele-specific criteria."""
        self.log.info('Flag internal alleles using allele-specific criteria.')
        internal_gene_counter = 0
        non_dmel_drosophilid_counter = 0
        for allele in self.fb_data_entities.values():
            if allele.allele_of_internal_gene is True:
                allele.linkmldto.internal = True
                allele.internal_reasons.append('Allele of internal type FB gene.')
                internal_gene_counter += 1
            if allele.adj_org_abbr != 'Dmel':
                allele.linkmldto.internal = True
                allele.internal_reasons.append('Allele of non-Dmel Drosophilid.')
                non_dmel_drosophilid_counter += 1
        self.log.info(f'Flagged {internal_gene_counter} alleles of internal-type genes as internal.')
        self.log.info(f'Flagged {non_dmel_drosophilid_counter} non-Dmel Drosophilid alleles as internal.')
        return

    def map_mutation_types(self):
        """Map mutation types."""
        self.log.info('Map mutation types.')
        insertion_conversion = {
            'transposable_element_insertion_site': 'SO:0001218',    # transgenic_insertion
            'transposable_element': 'SO:0001837',                   # mobile_element_insertion
            'insertion': 'SO:0000667'                               # insertion
        }
        counter = 0
        for allele in self.fb_data_entities.values():
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
                if fb_feat_type_name in insertion_conversion.keys():
                    mutation_type_curie = insertion_conversion[fb_feat_type_name]
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

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_allele_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_internal_allele_status()
        self.map_extinction_info()
        self.map_inheritance_modes()
        self.map_collections()
        self.map_allele_database_status()
        self.map_mutation_types()
        # self.map_pubs()    # TEMPORARILY SUPPRESS UNTIL LOAD SPEED IMPROVES
        self.map_timestamps()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Elaborate on query_chado_and_export() for the GeneHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the GeneHandler."""
        super().query_chado_and_export(session)
        return
