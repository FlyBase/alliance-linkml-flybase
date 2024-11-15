"""Module:: agm_handlers.

Synopsis:
    Data handlers that export FlyBase data for strains and genotypes to
    Alliance AffectedGenomicModel (AGM) LinkML objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
import agr_datatypes
import fb_datatypes
from entity_handler import PrimaryEntityHandler


class StrainHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters strain data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the StrainHandler object."""
        super().__init__(log, testing)
        self.datatype = 'strain'
        self.fb_export_type = fb_datatypes.FBStrain
        self.agr_export_type = agr_datatypes.AffectedGenomicModelDTO
        self.primary_export_set = 'agm_ingest_set'

    test_set = {
        'FBsn0000001': 'Oregon-R-modENCODE',
        'FBsn0000091': 'DGRP-373',
        'FBsn0000272': 'iso-1',
        'FBsn0001072': 'DSPR-B1-019',
        'FBsn0000283': 'MV2-25',
        'FBsn0000284': 'DGRP_Flyland',
    }

    # Elaborate on get_general_data() for the StrainHandler.
    def get_general_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        return

    # Elaborate on get_datatype_data() for the StrainHandler.
    def get_datatype_data(self, session):
        """Extend the method for the StrainHandler."""
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
        return

    # Elaborate on synthesize_info() for the StrainHandler.
    def synthesize_info(self):
        """Extend the method for the StrainHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        return

    # Elaborate on map_fb_data_to_alliance() for the StrainHandler.
    def map_strain_basic(self):
        """Map basic FlyBase strain data to the Alliance object."""
        self.log.info('Map basic strain info.')
        bob_counter = 0
        billy_counter = 0
        bad_counter = 0
        for strain in self.fb_data_entities.values():
            # BILLYBOB: Test recall_relationships() for non-feature entities.
            bob_rels = strain.recall_relationships(self.log, entity_role='subject', rel_types='derived_from')
            self.log.info(f'For {strain}, found {len(bob_rels)} relationships where strain is the subject and the rel_types="derived_from".')
            bob_counter += len(bob_rels)
            billy_rels = strain.recall_relationships(self.log, entity_role='object', rel_types='derived_from')
            self.log.info(f'For {strain}, found {len(billy_rels)} relationships where strain is the object and the rel_types="derived_from".')
            billy_counter += len(billy_rels)
            bad_rels = strain.recall_relationships(self.log, entity_role='object', rel_types='billybob')
            self.log.info(f'For {strain}, found {len(bad_rels)} relationships where strain is the object and the rel_types="billybob".')
            bad_counter += len(bad_rels)
            agr_strain = self.agr_export_type()
            agr_strain.obsolete = strain.chado_obj.is_obsolete
            agr_strain.mod_entity_id = f'FB:{strain.uniquename}'
            agr_strain.mod_internal_id = str(strain.chado_obj.strain_id)
            agr_strain.taxon_curie = strain.ncbi_taxon_id
            agr_strain.name = strain.name
            agr_strain.subtype_name = 'strain'
            strain.linkmldto = agr_strain
        self.log.info(f'BILLYBOB: bob_counter={bob_counter}')
        self.log.info(f'BILLYBOB: billy_counter={billy_counter}')
        self.log.info(f'BILLYBOB: bad_counter={bad_counter}')
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the StrainHandler."""
        super().map_fb_data_to_alliance()
        self.map_strain_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_pubs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return
