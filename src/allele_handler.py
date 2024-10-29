"""Module:: allele_handler.

Synopsis:
    A data handler that exports FlyBase data for alleles to Alliance Gene LinkML
    objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from agr_datatypes import AlleleDTO
from fb_datatypes import FBAllele
from feature_handler import FeatureHandler


class AlleleHandler(FeatureHandler):
    """This object gets, synthesizes and filters allele data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AlleleHandler object."""
        super().__init__(log, testing)

    datatype = 'allele'
    fb_export_type = FBAllele
    agr_export_type = AlleleDTO
    primary_agr_ingest_type = 'allele_ingest_set'

    test_set = {
        'FBal0137236': 'gukh[142]',    # P{hsneo}Xrp1142 insertion allele.
    }

    # Elaborate on get_general_data() for the AlleleHandler.
    def get_general_data(self, session):
        """Extend the method for the AlleleHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        self.build_feature_lookup(session)
        self.build_feature_relationship_evidence_lookup(session)
        return

    def get_parent_genes(self, session):
        """Get parent genes for alleles."""
        self.log.info('Get parent genes for alleles.')
        self.get_entity_sbj_feat_rel_by_type(session, 'allele_rels', rel_type='alleleof', obj_type='gene', obj_regex=self.regex['gene'])
        return

    def get_datatype_data(self, session, datatype, fb_export_type, agr_export_type):
        """Extend the method for the GeneHandler."""
        super().get_datatype_data(session, datatype, fb_export_type, agr_export_type)
        self.get_entities(session, self.datatype, self.fb_export_type)
        self.get_entity_pubs(session, self.datatype)
        self.get_entity_synonyms(session, self.datatype)
        self.get_entity_fb_xrefs(session, self.datatype)
        self.get_entity_xrefs(session, self.datatype)
        self.get_entity_timestamps(session, self.datatype)
        self.get_parent_genes(session)
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_parent_genes(self):
        """Get allele parent gene IDs."""
        self.log.info('Get allele parent gene IDs.')
        allele_counter = 0
        for allele in self.fb_data_entities.values():
            parent_gene_ids = []
            for gene_feature_id in allele.parent_gene_rels:
                parent_gene = self.feature_lookup[gene_feature_id]
                if parent_gene['is_obsolete'] is False:
                    parent_gene_ids.append(parent_gene['uniquename'])
            if len(parent_gene_ids) == 1:
                allele.parent_gene_id = parent_gene_ids[0]
                allele_counter += 1
            elif len(parent_gene_ids) == 0:
                self.log.warning(f'{allele} has no parent gene!')
            else:
                self.log.warning(f'{allele} has many parent genes!')
        self.log.info(f'Found parental gene for {allele_counter} alleles.')
        return

    # Elaborate on synthesize_info() for the GeneHandler.
    def synthesize_info(self, datatype, fb_export_type, agr_export_type):
        """Extend the method for the GeneHandler."""
        super().synthesize_info(datatype, fb_export_type, agr_export_type)
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_props()
        self.synthesize_pubs()
        self.synthesize_parent_genes()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_allele_basic(self, agr_export_type):
        """Map basic FlyBase allele data to the Alliance LinkML object."""
        self.log.info('Map basic allele info to Alliance object.')
        for allele in self.fb_data_entities.values():
            agr_allele = agr_export_type()
            agr_allele.obsolete = allele.chado_obj.is_obsolete
            agr_allele.mod_entity_id = f'FB:{allele.uniquename}'
            agr_allele.mod_internal_id = str(allele.chado_obj.feature_id)
            agr_allele.taxon_curie = allele.ncbi_taxon_id
            allele.linkmldto = agr_allele
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneHandler.
    def map_fb_data_to_alliance(self, datatype, fb_export_type, agr_export_type):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance(datatype, fb_export_type, agr_export_type)
        self.map_allele_basic(agr_export_type)
        self.map_synonyms()
        self.map_data_provider_dto()
        # self.map_pubs()    # TEMPORARILY SUPPRESS UNTIL LOAD SPEED IMPROVES
        self.map_xrefs()
        self.map_timestamps()
        self.map_secondary_ids('allele_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Elaborate on query_chado_and_export() for the GeneHandler.
    def query_chado_and_export(self, session, datatype, fb_export_type, agr_export_type):
        """Elaborate on query_chado_and_export method for the GeneHandler."""
        super().query_chado_and_export(session, datatype, fb_export_type, agr_export_type)
        return
