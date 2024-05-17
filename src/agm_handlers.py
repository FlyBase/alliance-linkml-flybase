"""Module:: agm_handlers.

Synopsis:
    Data handlers that export FlyBase data for strains and genotypes to
    Alliance AffectedGenomicModel (AGM) LinkML objects.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import csv
import datetime
import json
import re
import strict_rfc3339
from logging import Logger
from sqlalchemy.orm import aliased, Session
# from sqlalchemy.inspection import inspect
from harvdev_utils.char_conversions import sub_sup_sgml_to_html
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Organism, OrganismDbxref, Pub, PubDbxref, Featureloc,
    Strain, StrainPub, StrainSynonym, StrainDbxref, Strainprop, StrainpropPub, StrainCvterm, StrainCvtermprop,
    Feature, FeaturePub, FeatureSynonym, FeatureDbxref, Featureprop, FeaturepropPub, FeatureCvterm, FeatureCvtermprop,
    FeatureRelationship, FeatureRelationshipPub, Synonym
)
import fb_datatypes
import agr_datatypes
from handler import PrimaryEntityHandler


class StrainHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters strain data for export."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the StrainHandler object."""
        super().__init__(log, fb_data_type, testing)

    test_set = {
        'FBsn0000001': 'Oregon-R-modENCODE',
        'FBsn0000091': 'DGRP-373',
        'FBsn0000272': 'iso-1',
        'FBsn0001072': 'DSPR-B1-019',
        'FBsn0000283': 'MV2-25',
        'FBsn0000284': 'DGRP_Flyland',
    }
    # Elaborate on export filters for StrainHandler.
    required_fields = {
        'agm_ingest_set': [
            'curie',
            'data_provider_dto',
            'internal',
            'subtype_name',
            'taxon_curie',
        ],
    }
    output_fields = {
        'agm_ingest_set': [
            'agm_secondary_id_dtos',
            'created_by_curie',
            'cross_reference_dtos',
            'curie',
            'data_provider_dto',
            'date_created',
            'date_updated',
            'internal',
            'updated_by_curie',
            'name',
            'obsolete',
            'reference_curies',
            'subtype_name',
            'taxon_curie'
        ],
    }

    # Elaborate on get_general_data() for the StrainHandler.
    def get_general_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() for the StrainHandler.
    def get_datatype_data(self, session):
        """Extend the method for the StrainHandler."""
        super().get_datatype_data(session)
        return

    # Elaborate on synthesize_info() for the StrainHandler.
    def synthesize_info(self):
        """Extend the method for the StrainHandler."""
        super().synthesize_info()
        return

    # Elaborate on map_fb_data_to_alliance() for the StrainHandler.
    def map_strain_basic(self):
        """Map basic FlyBase strain data to the Alliance object."""
        self.log.info('Map basic strain info.')
        for strain in self.fb_data_entities.values():
            agr_strain = agr_datatypes.AffectedGenomicModelDTO()
            agr_strain.obsolete = strain.chado_obj.is_obsolete
            agr_strain.curie = f'FB:{strain.uniquename}'
            agr_strain.taxon_curie = strain.ncbi_taxon_id
            agr_strain.name = strain.name
            agr_strain.subtype_name = 'strain'
            strain.linkmldto = agr_strain
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the StrainHandler."""
        super().map_fb_data_to_alliance()
        self.map_strain_basic()
        self.map_data_provider_dto()
        self.map_pubs()
        self.map_xrefs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return
