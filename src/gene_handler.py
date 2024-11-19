"""Module:: gene_handler.

Synopsis:
    A data handler that exports FlyBase data for genes to Alliance Gene LinkML
    objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import csv
import re
from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class GeneHandler(FeatureHandler):
    """This object gets, synthesizes and filters gene data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the GeneHandler object."""
        super().__init__(log, testing)
        self.datatype = 'gene'
        self.fb_export_type = fb_datatypes.FBGene
        self.agr_export_type = agr_datatypes.GeneDTO
        self.primary_export_set = 'gene_ingest_set'

    test_set = {
        'FBgn0284084': 'wg',                # Current annotated nuclear protein_coding gene.
        'FBgn0004009': 'wg',                # Obsolete annotated nuclear protein_coding gene.
        'FBgn0013687': 'mt:ori',            # Current localized but unannotated mitochondrial gene.
        'FBgn0013678': 'mt:Cyt-b',          # Current annotated mitochondrial protein_coding gene.
        'FBgn0019661': 'lncRNA:roX1',       # Current annotated nuclear ncRNA gene.
        'FBgn0262451': 'mir-ban',           # Current annotated nuclear miRNA gene.
        'FBgn0034365': 'CG5335',            # Current annotated gene with CG symbol.
        'FBgn0003884': 'alphaTub84B',       # Current annotated gene with non-ASCII char in symbol.
        'FBgn0263477': 'scaRNA:PsiU1-6',    # Current annotated gene needs systematic synonym dto generated from annotation ID.
        'FBgn0030179': 'CG12094',           # Obsolete unannotated gene, should not get systematic name but needs symbol.
        'FBgn0108495': 'Dere\\GG16260',     # Current unannotated non-Dmel with systematic name.
        'FBgn0031087': 'CG12656',           # Current withdrawn gene.
        'FBgn0000154': 'Bar',               # Current unannotated gene.
        'FBgn0001200': 'His4',              # Current unannotated gene family.
        'FBgn0087003': 'tal',               # Current unannotated oddball.
        'FBgn0015267': 'Mmus\\Abl1',        # Current mouse gene with MGI xref.
    }

    # Additional reference info.
    gene_allele_rels = {}            # Will be (allele feature_id, gene feature_id) tuples keying lists of FBRelationships.
    gene_allele_associations = []    # Will be the final list of gene-allele FBRelationships to export.
    pthr_dict = {}                   # Will be an 1:1 FBgn_ID-PTHR xref dict.

    # Elaborate on get_general_data() for the GeneHandler.
    def get_general_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        self.build_feature_lookup(session, feature_types=['gene', 'allele'])
        self.get_internal_genes(session)
        self.get_chr_info(session)
        return

    # Elaborate on get_datatype_data() for the GeneHandler.
    def get_panther_info(self):
        """Extract panther information from file."""
        self.log.info('Extract panther information from file.')
        filepath = '/src/input/PTHR18.0_fruit_fly'
        tsv_file = open(filepath, "r")
        tsvin = csv.reader(tsv_file, delimiter='\t')
        FB = 0
        PTHR = 3
        counter = 0
        gene_regex = r'FBgn[0-9]{7}'
        for row in tsvin:
            fields = len(row)
            if fields:  # Ignore blank lines
                if re.search(gene_regex, row[FB]) and re.search(self.regex['panther'], row[PTHR]):
                    self.pthr_dict[re.search(gene_regex, row[FB]).group(0)] = re.search(self.regex['panther'], row[PTHR]).group(0)
                    counter += 1
        self.log.info(f'Processed {counter} lines from the panther orthology file.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        # self.get_entity_relationships(session, 'subject')
        self.get_entity_relationships(session, 'object', rel_type='alleleof', entity_type='allele', entity_regex=self.regex['allele'])
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        self.get_panther_info()
        self.get_annotation_ids(session)
        self.get_chr_featurelocs(session)
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_gene_type(self):
        """Synthesize gene type."""
        self.log.info('Synthesize gene type.')
        for gene in self.fb_data_entities.values():
            if 'promoted_gene_type' not in gene.props_by_type.keys():
                continue
            gene_types = gene.props_by_type['promoted_gene_type']
            if len(gene_types) == 1:
                prop_value = gene_types[0].chado_obj.value
                gene.gene_type_name = prop_value[11:-1]
                gene.gene_type_id = prop_value[1:10].replace('SO', 'SO:')
            elif len(gene_types) > 1:
                self.log.warning(f'{gene} has many promoted gene types.')
        return

    def synthesize_gene_alleles(self):
        """Synthesize gene allele relationships."""
        self.log.info('Synthesize gene allele relationships.')
        gene_counter = 0
        allele_counter = 0
        # Need to code for the rare possibility that gene-allele is represented by many feature_relationships.
        for gene in self.fb_data_entities.values():
            relevant_allele_rels = gene.recall_relationships(self.log, entity_role='object', rel_types='alleleof', rel_entity_types='allele')
            if relevant_allele_rels:
                gene_counter += 1
            # self.log.debug(f'For {gene}, found {len(relevant_allele_rels)} allele rels to review.')
            for allele_rel in relevant_allele_rels:
                allele_feature_id = allele_rel.chado_obj.subject_id
                allele_gene_key = (allele_feature_id, gene.db_primary_id)
                try:
                    self.gene_allele_rels[allele_gene_key].extend(allele_rel)
                except KeyError:
                    self.gene_allele_rels[allele_gene_key] = [allele_rel]
                    allele_counter += 1
        self.log.info(f'Found {allele_counter} alleles for {gene_counter} genes.')
        return

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_info(self):
        """Extend the method for the GeneHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_gene_type()
        self.synthesize_anno_ids()
        self.synthesize_gene_alleles()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_gene_basic(self):
        """Map basic FlyBase gene data to the Alliance LinkML object."""
        self.log.info('Map basic gene info to Alliance object.')
        for gene in self.fb_data_entities.values():
            agr_gene = self.agr_export_type()
            agr_gene.obsolete = gene.chado_obj.is_obsolete
            agr_gene.mod_entity_id = f'FB:{gene.uniquename}'
            agr_gene.mod_internal_id = str(gene.chado_obj.feature_id)
            agr_gene.taxon_curie = gene.ncbi_taxon_id
            gene.linkmldto = agr_gene
        return

    def map_gene_type(self):
        """Map gene type."""
        self.log.info('Map gene type to Alliance object.')
        for gene in self.fb_data_entities.values():
            gene.linkmldto.gene_type_curie = gene.gene_type_id
        return

    def map_gene_snapshot(self):
        """Map gene snapshot."""
        self.log.info('Map gene snapshot to Alliance object.')
        for gene in self.fb_data_entities.values():
            if 'gene_summary_text' not in gene.props_by_type.keys():
                continue
            gene_snapshots = gene.props_by_type['gene_summary_text']
            if len(gene_snapshots) == 1:
                note_type_name = 'MOD_provided_gene_description'
                free_text = gene_snapshots[0].chado_obj.value.replace('@', '')
                pub_curies = ['FB:FBrf0232436']
                snapshot_note_dto = agr_datatypes.NoteDTO(note_type_name, free_text, pub_curies).dict_export()
                gene.linkmldto.related_notes.append(snapshot_note_dto)
            elif len(gene_snapshots) > 1:
                self.log.warning(f'{gene} has many gene snapshots.')
        return

    def map_gene_panther_xrefs(self):
        """Add panther xrefs."""
        self.log.info('Map panther xrefs to Alliance object.')
        for gene in self.fb_data_entities.values():
            if gene.uniquename not in self.pthr_dict.keys():
                return
            # Build Alliance xref DTO
            prefix = 'PANTHER'
            page_area = 'FB'
            curie = f'{prefix}:{self.pthr_dict[gene.uniquename]}'
            display_name = curie
            xref_dto = agr_datatypes.CrossReferenceDTO(prefix, curie, page_area, display_name).dict_export()
            gene.linkmldto.cross_reference_dtos.append(xref_dto)
        return

    def map_gene_allele_associations(self):
        """Map gene-allele associations to Alliance object."""
        self.log.info('Map gene-allele associations to Alliance object.')
        ALLELE = 0
        GENE = 1
        counter = 0
        for allele_gene_key, allele_gene_rels in self.gene_allele_rels.items():
            allele_curie = f'FB:{self.feature_lookup[allele_gene_key[ALLELE]]["uniquename"]}'
            gene_curie = f'FB:{self.feature_lookup[allele_gene_key[GENE]]["uniquename"]}'
            first_feat_rel = allele_gene_rels[0]
            if not allele_curie.startswith('FB:FBal'):
                first_feat_rel.for_export = False
                first_feat_rel.export_warnings.append('"alleleof" relationship includes non-allele')
            all_pub_ids = []
            for allele_gene_rel in allele_gene_rels:
                all_pub_ids.extend(allele_gene_rel.pubs)
            first_feat_rel.pubs = all_pub_ids
            pub_curies = self.lookup_pub_curies(all_pub_ids)
            rel_dto = agr_datatypes.AlleleGeneAssociationDTO(allele_curie, 'is_allele_of', gene_curie, pub_curies)
            if self.feature_lookup[allele_gene_key[ALLELE]]['is_obsolete'] is True or self.feature_lookup[allele_gene_key[GENE]]['is_obsolete'] is True:
                rel_dto.obsolete = True
                rel_dto.internal = True
            first_feat_rel.linkmldto = rel_dto
            self.gene_allele_associations.append(first_feat_rel)
            counter += 1
        self.log.info(f'Generated {counter} allele-gene unique associations.')
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_gene_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        # self.map_pubs()    # Suppress until LinkML Gene gets reference_curies slot.
        self.map_timestamps()
        self.map_secondary_ids('gene_secondary_id_dtos')
        # self.map_gene_snapshot()
        self.map_gene_type()
        self.map_gene_panther_xrefs()
        self.map_anno_ids_to_secondary_ids('gene_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        self.map_gene_allele_associations()
        self.flag_internal_fb_entities('gene_allele_associations')
        return

    # Elaborate on query_chado_and_export() for the GeneHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the GeneHandler."""
        super().query_chado_and_export(session)
        self.flag_unexportable_entities(self.gene_allele_associations, 'allele_gene_association_ingest_set')
        self.generate_export_dict(self.gene_allele_associations, 'allele_gene_association_ingest_set')
        return
