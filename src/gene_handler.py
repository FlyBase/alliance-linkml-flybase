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
        'FBgn0284084': 'wg',                  # Current annotated nuclear protein_coding gene.
        'FBgn0004009': 'wg',                  # Obsolete annotated nuclear protein_coding gene.
        'FBgn0044027': 'Ori66Dbeta',          # Current unannotated gene representing origin_of_replication.
        'FBgn0013687': 'mt:ori',              # Current localized but unannotated mitochondrial gene.
        'FBgn0013678': 'mt:Cyt-b',            # Current annotated mitochondrial protein_coding gene.
        'FBgn0019661': 'lncRNA:roX1',         # Current annotated nuclear ncRNA gene.
        'FBgn0262451': 'mir-ban',             # Current annotated nuclear miRNA gene.
        'FBgn0034365': 'CG5335',              # Current annotated gene with CG symbol.
        'FBgn0003884': 'alphaTub84B',         # Current annotated gene with non-ASCII char in symbol.
        'FBgn0263477': 'scaRNA:PsiU1-6',      # Current annotated gene needs systematic synonym dto generated from annotation ID.
        'FBgn0030179': 'CG12094',             # Obsolete unannotated gene, should not get systematic name but needs symbol.
        'FBgn0066164': 'Dsim_Hmr',            # Current Dsim gene.
        'FBgn0016335': 'Dsim_HeT-A_gag',      # Current Dsim gene for TE CDS.
        'FBgn0108495': 'Dere_GG16260',        # Obsolete unannotated non-Dmel obsolete gene with systematic name.
        'FBgn0031087': 'CG12656',             # Current withdrawn gene.
        'FBgn0000154': 'Bar',                 # Current unannotated gene.
        'FBgn0001200': 'His4',                # Current unannotated gene family.
        'FBgn0087003': 'tal',                 # Current unannotated oddball.
        'FBgn0015267': 'Mmus_Abl1',           # Current mouse gene, MGI:87859.
        'FBgn0026367': 'Scer_GAL80',          # Current yeast gene, SGD:S000004515.
        'FBgn0287889': 'SARS-CoV-2_ORF3a',    # Current SARS-CoV2 gene, REFSEQ:YP_009724391.
    }

    # Additional reference info.
    pthr_dict = {}                   # Will be an 1:1 FBgn_ID-PTHR xref dict.

    # Elaborate on get_general_data() for the GeneHandler.
    def get_general_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
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

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_info(self):
        """Extend the method for the GeneHandler."""
        super().synthesize_info()
        self.flag_new_additions_and_obsoletes()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.synthesize_gene_type()
        self.synthesize_anno_ids()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_gene_basic(self):
        """Map basic FlyBase gene data to the Alliance LinkML object."""
        self.log.info('Map basic gene info to Alliance object.')
        for gene in self.fb_data_entities.values():
            agr_gene = self.agr_export_type()
            agr_gene.obsolete = gene.chado_obj.is_obsolete
            agr_gene.primary_external_id = f'FB:{gene.uniquename}'
            # agr_gene.mod_internal_id = f'FB.feature_id={gene.db_primary_id}'
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

    def flag_unexportable_genes(self):
        """Flag unexportable genes."""
        self.log.info('Flag unexportable genes.')
        counter = 0
        for gene in self.fb_data_entities.values():
            if self.organism_lookup[gene.organism_id]['is_drosophilid'] is False:
                gene.for_export = False
                gene.export_warnings.append('Non-Drosophilid gene')
                counter += 1
        self.log.info(f'Flagged {counter} genes as unexportable.')
        return

    def flag_internal_genes(self):
        """Flag internal genes."""
        self.log.info('Flag internal genes.')
        counter = 0
        for gene in self.fb_data_entities.values():
            if self.organism_lookup[gene.organism_id]['abbreviation'] != 'Dmel':
                gene.linkmldto.internal = True
                gene.internal_reasons.append('Non-Dmel')
            if gene.uniquename in self.internal_gene_ids:
                gene.linkmldto.internal = True
                gene.internal_reasons.append('Internal gene type.')
            if gene.linkmldto.internal is True:
                counter += 1
        self.log.info(f'Flagged {counter} genes as internal for gene-specific reasons')
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
        self.flag_internal_genes()
        self.flag_internal_fb_entities('fb_data_entities')
        self.flag_unexportable_genes()
        return

    # Elaborate on query_chado_and_export() for the GeneHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the GeneHandler."""
        super().query_chado_and_export(session)
        return
