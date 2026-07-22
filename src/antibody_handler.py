"""Module:: antibody_handler.

Synopsis:
    A data handler that exports FlyBase gene-level antibody data to Alliance
    Antibody LinkML objects. FlyBase has no first class antibody entities; each
    exported antibody is synthesized from data hung off a current gene (FBgn):
      1. "Lab-generated" antibodies, from a gene's "reported_antibod_gen"
         Featureprop (clonality in the prop value) and its FeaturepropPub pub(s).
      2. "Commercial" antibodies, from a gene FeatureDbxref to the "DSHB" or
         "Cell Signaling Technology" db, with clonality in a "linkout" Dbxrefprop
         value, and the antibody ID in the Dbxref accession.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy import func
from sqlalchemy.orm import aliased
from harvdev_utils.reporting import (
    Cvterm, Db, Dbxref, Dbxrefprop, Feature, FeatureDbxref,
    Featureprop, FeaturepropPub, Pub
)
import agr_datatypes
import fb_datatypes
from handler import DataHandler


class AntibodyHandler(DataHandler):
    """This object gets, synthesizes and filters gene-level antibody data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the AntibodyHandler object."""
        super().__init__(log, testing)
        self.datatype = 'antibody'
        self.fb_export_type = fb_datatypes.FBAntibody
        self.agr_export_type = agr_datatypes.AntibodyDTO
        self.primary_export_set = 'antibody_ingest_set'

    # Map of commercial antibody chado db.name to the short "source" label used in the antibody name.
    commercial_antibody_dbs = {
        'DSHB': 'DSHB',
        'Cell Signaling Technology': 'CST',
    }
    # Antibody clonality values that FlyBase reports.
    valid_clonalities = ['monoclonal', 'polyclonal']

    def normalize_clonality(self, raw_value):
        """Return "monoclonal" or "polyclonal" from a raw prop value, or None if neither is present."""
        if raw_value is None:
            return None
        value = raw_value.strip().lower()
        for clonality in self.valid_clonalities:
            if clonality in value:
                return clonality
        return None

    # Elaborate on get_general_data() for the AntibodyHandler.
    def get_general_data(self, session):
        """Extend the method for the AntibodyHandler."""
        super().get_general_data(session)
        self.build_organism_lookup(session)
        self.build_bibliography(session)
        return

    # Elaborate on get_datatype_data() for the AntibodyHandler.
    def get_datatype_data(self, session):
        """Extend the method for the AntibodyHandler."""
        super().get_datatype_data(session)
        self.get_lab_antibodies(session)
        self.get_commercial_antibodies(session)
        return

    # Add methods to be run by get_datatype_data() below.
    def get_lab_antibodies(self, session):
        """Get lab-generated antibodies from gene "reported_antibod_gen" featureprops."""
        self.log.info('Get lab-generated antibodies from gene "reported_antibod_gen" featureprops.')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['gene']),
            Cvterm.name == 'reported_antibod_gen',
            Pub.is_obsolete.is_(False),
        )
        results = session.query(Feature, Featureprop, Pub).\
            select_from(Feature).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(Cvterm, (Cvterm.cvterm_id == Featureprop.type_id)).\
            join(FeaturepropPub, (FeaturepropPub.featureprop_id == Featureprop.featureprop_id)).\
            join(Pub, (Pub.pub_id == FeaturepropPub.pub_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        skipped_counter = 0
        for result in results:
            clonality = self.normalize_clonality(result.Featureprop.value)
            if clonality is None:
                skipped_counter += 1
                continue
            gene = result.Feature
            pub = result.Pub
            uniq_key = f'{gene.uniquename}_{clonality}_{pub.uniquename}'
            if uniq_key in self.fb_data_entities.keys():
                continue
            antibody = self.fb_export_type(uniq_key, gene.uniquename, gene.name, gene.organism_id, clonality)
            antibody.pub_id = pub.pub_id
            self.fb_data_entities[uniq_key] = antibody
            counter += 1
        self.log.info(f'Found {counter} lab-generated antibodies ({skipped_counter} skipped for unrecognized clonality).')
        return

    def get_commercial_antibodies(self, session):
        """Get commercial antibodies from gene DSHB/CST feature_dbxrefs."""
        self.log.info('Get commercial antibodies from gene DSHB/CST feature_dbxrefs.')
        prop_type = aliased(Cvterm, name='prop_type')
        filters = (
            Feature.is_obsolete.is_(False),
            Feature.uniquename.op('~')(self.regex['gene']),
            FeatureDbxref.is_current.is_(True),
            Db.name.in_((list(self.commercial_antibody_dbs.keys()))),
            prop_type.name == 'linkout',
            func.lower(Dbxrefprop.value).in_((self.valid_clonalities)),
        )
        results = session.query(Feature, Db, Dbxref, Dbxrefprop).\
            select_from(Feature).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            join(Dbxrefprop, (Dbxrefprop.dbxref_id == Dbxref.dbxref_id)).\
            join(prop_type, (prop_type.cvterm_id == Dbxrefprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            clonality = result.Dbxrefprop.value.strip().lower()
            gene = result.Feature
            source = self.commercial_antibody_dbs[result.Db.name]
            accession = result.Dbxref.accession
            uniq_key = f'{gene.uniquename}_{clonality}_{source}_{accession}'
            if uniq_key in self.fb_data_entities.keys():
                continue
            antibody = self.fb_export_type(uniq_key, gene.uniquename, gene.name, gene.organism_id, clonality)
            antibody.source = source
            antibody.accession = accession
            self.fb_data_entities[uniq_key] = antibody
            counter += 1
        self.log.info(f'Found {counter} commercial antibodies.')
        return

    # Elaborate on synthesize_info() for the AntibodyHandler.
    def synthesize_info(self):
        """Extend the method for the AntibodyHandler."""
        super().synthesize_info()
        self.synthesize_antibody_names()
        self.synthesize_antibody_references()
        self.synthesize_antigen_taxon()
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_antibody_names(self):
        """Build the made-up name for each antibody, joining components by underscores."""
        self.log.info('Synthesize antibody names.')
        for antibody in self.fb_data_entities.values():
            if antibody.source is not None:
                # Commercial: gene name, clonality, source (DSHB/CST), antibody ID.
                components = [antibody.gene_name, antibody.clonality, antibody.source, antibody.accession]
            else:
                # Lab-generated: gene name, clonality, pub curie (PMID if available, else FB:FBrf).
                pub_curie = self.lookup_single_pub_curie(antibody.pub_id)
                components = [antibody.gene_name, antibody.clonality, pub_curie]
            antibody.antibody_name = '_'.join([str(component) for component in components])
        return

    def synthesize_antibody_references(self):
        """Look up the supporting pub curie for each lab-generated antibody."""
        self.log.info('Look up supporting pub curies for lab-generated antibodies.')
        counter = 0
        for antibody in self.fb_data_entities.values():
            if antibody.pub_id is None:
                continue
            pub_curie = self.lookup_single_pub_curie(antibody.pub_id)
            if pub_curie is not None and pub_curie != 'FB:unattributed':
                antibody.reference_curie = pub_curie
                counter += 1
        self.log.info(f'Found supporting pub curies for {counter} lab-generated antibodies.')
        return

    def synthesize_antigen_taxon(self):
        """Determine the antigen taxon for each antibody from its target gene organism."""
        self.log.info('Determine the antigen taxon for each antibody from the target gene organism.')
        for antibody in self.fb_data_entities.values():
            try:
                antibody.antigen_taxon_curie = self.organism_lookup[antibody.gene_organism_id]['taxon_curie']
            except KeyError:
                self.log.warning(f'No taxon found for organism_id={antibody.gene_organism_id} ({antibody}).')
        return

    # Elaborate on map_fb_data_to_alliance() for the AntibodyHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the AntibodyHandler."""
        super().map_fb_data_to_alliance()
        self.map_antibody_basic()
        self.map_data_provider_dto()
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_antibody_basic(self):
        """Map basic FlyBase antibody data to the Alliance LinkML object."""
        self.log.info('Map basic antibody info to the Alliance object.')
        for antibody in self.fb_data_entities.values():
            agr_antibody = self.agr_export_type()
            agr_antibody.name = antibody.antibody_name
            agr_antibody.clonality_name = antibody.clonality
            agr_antibody.antigen_taxon_curie = antibody.antigen_taxon_curie
            agr_antibody.antibody_target_gene_identifiers = [f'FB:{antibody.gene_uniquename}']
            if antibody.reference_curie is not None:
                agr_antibody.reference_curies = [antibody.reference_curie]
                agr_antibody.original_reference_curie = antibody.reference_curie
            antibody.linkmldto = agr_antibody
        return

    def map_data_provider_dto(self):
        """Map the data provider to the Alliance object.

        Antibodies have no dedicated FB report page (and no FB curie), so the
        DataProviderDTO carries only the "FB" source organization, with no
        cross-reference.
        """
        self.log.info('Map data provider to the Alliance object.')
        for antibody in self.fb_data_entities.values():
            if antibody.linkmldto is None:
                continue
            antibody.linkmldto.data_provider_dto = agr_datatypes.DataProviderDTO(None).dict_export()
        return
