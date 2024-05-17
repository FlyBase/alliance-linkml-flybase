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
from feature_handler import FeatureHandler


class GeneHandler(FeatureHandler):
    """This object gets, synthesizes and filters gene data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the GeneHandler object."""
        super().__init__(log, fb_data_type, testing)
        self.pthr_dict = {}    # Will be an 1:1 FBgn_ID-PTHR xref dict.

    test_set = {
        'FBgn0284084': 'wg',                # Current annotated nuclear protein_coding gene.
        'FBgn0004009': 'wg',                # Obsolete annotated nuclear protein_coding gene.
        'FBgn0013687': 'mt:ori',            # Current localized but unannotated mitochondrial gene.
        'FBgn0013678': 'mt:Cyt-b',          # Current annotated mitochondrial protein_coding gene.
        'FBgn0019661': 'lncRNA:roX1',       # Current annotated nuclear ncRNA gene.
        'FBgn0262451': 'mir-ban',           # Current annotated nuclear miRNA gene.
        'FBgn0034365': 'CG5335',            # Current annotated gene with CG symbol.
        'FBgn0003884': 'alphaTub84B',       # Current annotated gene with non-ASCII char in symbol.
        'FBgn0263477': 'scaRNA:PsiU1-6',    # Current annotated gene needs systematic synonym dto.
        'FBgn0030179': 'CG12094',           # Obsolete unannotated gene, should not get systematic name but needs symbol.
        'FBgn0108495': 'Dere\\GG16260',     # Current unannotated non-Dmel with systematic name.
        'FBgn0031087': 'CG12656',           # Current withdrawn gene.
        'FBgn0000154': 'Bar',               # Current unannotated gene.
        'FBgn0001200': 'His4',              # Current unannotated gene family.
        'FBgn0087003': 'tal',               # Current unannotated oddball.
        'FBgn0015267': 'Mmus\\Abl1',        # Current mouse gene with MGI xref.
    }
    # Elaborate on export filters for GeneHandler.
    required_fields = {
        'gene_ingest_set': [
            'curie',
            'data_provider_dto',
            'gene_symbol_dto',
            'internal',
            'taxon_curie',
        ],
    }
    output_fields = {
        'gene_ingest_set': [
            'created_by_curie',
            'cross_reference_dtos',
            'curie',
            'data_provider_dto',
            'date_created',
            'date_updated',
            'gene_full_name_dto',
            'gene_symbol_dto',
            'gene_synonym_dtos',
            'gene_systematic_name_dto',
            'gene_secondary_id_dtos',
            'gene_type_curie',
            'internal',
            'obsolete',
            # 'related_notes',    # Not present in GeneDTO.
            'taxon_curie',
            'updated_by_curie',
        ],
    }
    fb_agr_db_dict = {
        'EntrezGene': 'NCBI_Gene',
        'RNAcentral': 'RNAcentral',
        # 'UniProt/GCRP': 'UniProt/GCRP',
        'UniProt/Swiss-Prot': 'UniProtKB',
        'UniProt/TrEMBL': 'UniProtKB',
        'SGD': 'SGD',
        'WormBase': 'WB',
        'ZFIN': 'ZFIN',
        'Xenbase': 'Xenbase',
        'RGD': 'RGD',
        'MGD': 'MGI',
        'MGI': 'MGI'
    }
    # Reference dicts.
    internal_gene_types = [
        'engineered_fusion_gene',
        'engineered_region',
        'gene_group',
        'gene_with_polycistronic_transcript',
        'insulator',
        'mitochondrial_sequence',
        'origin_of_replication',
        'region',
        'regulatory_region',
        'repeat_region',
        'satellite_DNA',
        'transposable_element_gene'
    ]

    # Elaborate on get_general_data() for the GeneHandler.
    def get_general_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_general_data(session)
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

    def get_gene_snapshots(self, session):
        """Get human-written gene summaries."""
        self.log.info('Get human-written gene summaries.')
        results = self.get_featureprops_by_type(session, 'gene_summary_text')
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].gene_snapshots.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} gene snapshots for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} gene snapshots for {self.fb_data_type} entities.')
        return

    def get_gene_types(self, session):
        """Get promoted_gene_type for genes."""
        self.log.info('Get promoted_gene_type for genes.')
        results = self.get_featureprops_by_type(session, 'promoted_gene_type')
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].gene_type_names.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} promoted_gene_types for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} promoted_gene_types for {self.fb_data_type} entities.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session)
        self.get_panther_info()
        self.get_annotation_ids(session)
        self.get_chr_featurelocs(session)
        self.get_gene_snapshots(session)
        self.get_gene_types(session)
        return

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_gene_type(self):
        """Synthesize gene type."""
        self.log.info('Synthesize gene type.')
        for gene in self.fb_data_entities.values():
            if len(gene.gene_type_names) == 1:
                prop_value = gene.gene_type_names[0].value
                gene.gene_type_name = prop_value[11:-1]
                gene.gene_type_id = prop_value[1:10].replace('SO', 'SO:')
            elif len(gene.gene_type_names) > 1:
                self.log.warning(f'{gene} has many promoted gene types.')
        return

    def synthesize_info(self):
        """Extend the method for the GeneHandler."""
        super().synthesize_info()
        self.synthesize_gene_type()
        self.synthesize_anno_ids()
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_gene_basic(self):
        """Map basic FlyBase gene data to the Alliance LinkML object."""
        self.log.info('Map basic gene info to Alliance object.')
        for gene in self.fb_data_entities.values():
            agr_gene = agr_datatypes.GeneDTO()
            agr_gene.obsolete = gene.chado_obj.is_obsolete
            agr_gene.curie = f'FB:{gene.uniquename}'
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
            if len(gene.gene_snapshots) == 1:
                note_type_name = 'MOD_provided_gene_description'
                free_text = gene.gene_snapshots[0].value.replace('@', '')
                pub_curies = ['FB:FBrf0232436']
                snapshot_note_dto = agr_datatypes.NoteDTO(note_type_name, free_text, pub_curies).dict_export()
                gene.linkmldto.related_notes.append(snapshot_note_dto)
            elif len(gene.gene_snapshots) > 1:
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

    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_gene_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_pubs()
        self.map_xrefs()
        self.map_timestamps()
        self.map_secondary_ids('gene_secondary_id_dtos')
        self.map_gene_snapshot()
        self.map_gene_type()
        self.map_gene_panther_xrefs()
        self.map_anno_ids_to_secondary_ids('gene_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return
