"""Module:: gene_seqfeat_handler.

Synopsis:
    A data handler that exports FlyBase sequence features (FBsf) of selected types
    to Alliance Gene LinkML objects.

    The exported FBsf types are declared in DataHandler.feature_subtypes['gene_seqfeat'],
    so the generic PrimaryEntityHandler driver query picks them up without any bespoke SQL
    here. The feature.type SO term is exported as the gene_type_curie. No types are enabled
    yet (FTA-284); types such as TSS may be added later.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class GeneSeqFeatHandler(FeatureHandler):
    """This object gets, synthesizes and filters FBsf data for export as Alliance genes."""
    def __init__(self, log: Logger, testing: bool):
        """Create the GeneSeqFeatHandler object."""
        super().__init__(log, testing)
        self.datatype = 'gene_seqfeat'
        self.fb_export_type = fb_datatypes.FBGeneSeqFeat
        self.agr_export_type = agr_datatypes.GeneDTO
        self.primary_export_set = 'gene_ingest_set'

    # Add examples here once a type is enabled in feature_subtypes['gene_seqfeat'].
    test_set = {}

    # FlyBase props exported as Alliance notes.
    gene_seqfeat_prop_to_note_mapping = {
        'comment': ('comment', 'note_dtos'),    # SF6
    }

    # Add methods to be run by get_general_data() below.
    def get_general_data(self, session):
        """Extend the method for the GeneSeqFeatHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_organism_lookup(session)
        self.build_cvterm_lookup(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_datatype_data(self, session):
        """Extend the method for the GeneSeqFeatHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_info(self):
        """Extend the method for the GeneSeqFeatHandler."""
        super().synthesize_info()
        self.flag_new_additions_and_obsoletes()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_gene_seqfeat_basic(self):
        """Map basic FlyBase FBsf data to the Alliance LinkML object."""
        self.log.info('Map basic FBsf info to Alliance gene object.')
        for seqfeat in self.fb_data_entities.values():
            agr_gene = self.agr_export_type()
            agr_gene.obsolete = seqfeat.chado_obj.is_obsolete
            agr_gene.primary_external_id = f'FB:{seqfeat.uniquename}'
            agr_gene.taxon_curie = seqfeat.ncbi_taxon_id
            seqfeat.linkmldto = agr_gene
        return

    def map_gene_seqfeat_type(self):
        """Map the FBsf feature.type SO term to the gene_type_curie."""
        self.log.info('Map FBsf feature type to Alliance gene type.')
        counter = 0
        for seqfeat in self.fb_data_entities.values():
            if seqfeat.linkmldto is None:
                continue
            seqfeat.linkmldto.gene_type_curie = self.cvterm_lookup[seqfeat.chado_obj.type_id]['curie']
            counter += 1
        self.log.info(f'Mapped {counter} FBsf feature types to Alliance gene types.')
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneSeqFeatHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneSeqFeatHandler."""
        super().map_fb_data_to_alliance()
        self.map_gene_seqfeat_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_timestamps()
        self.map_secondary_ids('gene_secondary_id_dtos')
        self.map_gene_seqfeat_type()
        self.map_entity_props_to_notes('gene_seqfeat_prop_to_note_mapping')
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Elaborate on query_chado_and_export() for the GeneSeqFeatHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the GeneSeqFeatHandler."""
        super().query_chado_and_export(session)
        return
