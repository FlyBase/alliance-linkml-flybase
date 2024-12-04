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
from harvdev_utils.production import (
    Db, Dbxref, GenotypeDbxref
)


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
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'balancer', 'insertion', 'gene'])
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
        self.get_entity_cvterms(session)
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
        for strain in self.fb_data_entities.values():
            agr_strain = self.agr_export_type()
            agr_strain.obsolete = strain.chado_obj.is_obsolete
            agr_strain.mod_entity_id = f'FB:{strain.uniquename}'
            agr_strain.mod_internal_id = str(strain.db_primary_id)
            agr_strain.taxon_curie = strain.ncbi_taxon_id
            agr_strain.name = strain.name
            agr_strain.subtype_name = 'strain'
            strain.linkmldto = agr_strain
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the StrainHandler."""
        super().map_fb_data_to_alliance()
        self.map_strain_basic()
        # self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_pubs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return


class GenotypeHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters genotype data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the GenotypeHandler object."""
        super().__init__(log, testing)
        self.datatype = 'genotype'
        self.fb_export_type = fb_datatypes.FBGenotype
        self.agr_export_type = agr_datatypes.AffectedGenomicModelDTO
        self.primary_export_set = 'agm_ingest_set'

    test_set = {
        2: 'Ab(1)ZWD16 | FBab0027942 | FBgo0000002',    # The first genotype in the table.
    }

    # Elaborate on get_general_data() for the GenotypeHandler.
    def get_general_data(self, session):
        """Extend the method for the GenotypeHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_ncbi_taxon_lookup(session)
        return

    # Additional sub-methods for get_datatype_data().
    # Note that for genotypes, the "uniquename" is not the FBgo accession, so need to get these another way.
    # Also note that the FBgo IDs are currently internal.
    def get_genotype_fb_curies(self, session):
        """Get FlyBase curies for genotypes."""
        self.log.info('Get FlyBase curies for genotypes.')
        filters = (
            GenotypeDbxref.is_current.is_(True),
            Dbxref.accession.op('~')(self.regex['genotype']),
            Db.name == 'FlyBase',
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (GenotypeDbxref.genotype_id.in_((self.test_set.keys())), )
        results = session.query(GenotypeDbxref, Dbxref).\
            select_from(GenotypeDbxref).\
            join(Dbxref, (Dbxref.dbxref_id == GenotypeDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            try:
                self.fb_data_entities[result.GenotypeDbxref.genotype_id].fb_curie = result.Dbxref.accession
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} FBgo IDs for {self.datatype}s.')
        return

    # Elaborate on get_datatype_data() for the GenotypeHandler.
    def get_datatype_data(self, session):
        """Extend the method for the GenotypeHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_genotype_fb_curies(session)
        self.get_entity_cvterms(session)
        self.get_entityprops(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        return

    # Elaborate on synthesize_info() for the GenotypeHandler.
    def synthesize_info(self):
        """Extend the method for the GenotypeHandler."""
        super().synthesize_info()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        return

    # Elaborate on map_fb_data_to_alliance() for the GenotypeHandler.
    def map_genotype_basic(self):
        """Map basic FlyBase genotype data to the Alliance object."""
        self.log.info('Map basic FlyBase genotype data to the Alliance object.')
        counter = 0
        for genotype in self.fb_data_entities.values():
            # Skip genotypes associated with stock import for now.
            if genotype.fb_curie is None:
                continue
            agr_genotype = self.agr_export_type()
            agr_genotype.obsolete = genotype.chado_obj.is_obsolete
            agr_genotype.mod_entity_id = f'FB:{genotype.fb_curie}'
            agr_genotype.mod_internal_id = str(genotype.db_primary_id)
            agr_genotype.taxon_curie = genotype.ncbi_taxon_id
            agr_genotype.name = genotype.uniquename
            agr_genotype.subtype_name = 'genotype'
            genotype.linkmldto = agr_genotype
            counter += 1
        self.log.info(f'{counter} genotypes could be mapped to the Alliance LinkML model.')
        return

    def map_fb_data_to_alliance(self):
        """Extend the method for the GenotypeHandler."""
        super().map_fb_data_to_alliance()
        self.map_genotype_basic()
        # self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_pubs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return
