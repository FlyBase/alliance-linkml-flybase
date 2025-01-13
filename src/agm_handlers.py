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
from harvdev_utils.reporting import (
    Db, Dbxref, Feature, FeatureGenotype, Genotype, GenotypeDbxref
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
        self.build_organism_lookup(session)
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
            # agr_strain.mod_internal_id = f'FB.strain_id={strain.db_primary_id}'
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
        2: 'Ab(1)ZWD16 | FBab0027942',                                           # The first genotype in the table.
        294012: 'P{CH1226-43A10} lz<up>L</up>',                                  # Has FBal and FBtp associated.
        223641: 'Dp1<up>EP2422</up> P{hsp26-pt-T}39C-12',                        # Has FBal and FBti associated.
        452205: 'wg<up>1</up>/wg<up>GBM</up>',                                   # Transheterozygous wg[1]/wg[GBM].
        450391: 'wg<up>l-8</up>/wg<up>l-8</up>',                                 # Homozygous wg[l-8] allele.
        367896: 'Tak1<up>2</up>/Tak1<up>+</up>',                                 # Heterozygous Tak1[2] over wt allele.
        515567: 'Sdc<up>12</up>/Sdc<up>unspecified</up>',                        # Sdc[12]/Sdc[unspecified].
        166899: 'Df(2R)173/PCNA<up>D-292</up>',                                  # PCNA[D-292]/Df(2R)173.
        168332: 'Df(3L)Ez7 hay<up>nc2.tMa</up>',                                 # Df(3L)Ez7 + hay[nc2.tMa] (diff complementation groups).
        166704: 'shi<up>EM33</up> shi<up>t15</up>',                              # shi[EM33] + shi[t15] (allele + rescue construct).
        219912: 'dpp<up>s4</up> wg<up>l-17</up>',                                # dpp[s4], wg[l-17].
        199449: 'Hsap_MAPT<up>UAS.cAa</up> Scer_GAL4<up>smid-C161</up>',         # GAL4 drives UAS human construct.
        32369: 'Dmau_w<up>a23</up>',                                             # Single non-Dmel classical allele.
        466842: 'Dsim_mir-983<up>KO</up>',                                       # Single non-Dmel classical allele.
        340500: 'Dsim_Lhr<up>1</up> Dsim_Lhr<up>2+16aa.Tag:HA</up>',             # non-Dmel classical allele + transgene.
        167097: 'Df(3R)CA3 Scer_GAL4<up>GMR.PF</up> pim<up>UAS.Tag:MYC</up>',    # Df, GAL4 and UAS.
        515875: 'Dp(1;Y)y<up>+</up>/tyn<up>1</up>',                              # tyn[1] / Dp(1:Y).
        202134: 'Tp(3;1)P115/pad<up>1</up>',                                     # pad[1] / Tp(3;1)P115.
        182943: 'Df(3L)emc/+',                                                   # Heterozygous deficiency.
        334079: 'Dsim_Int(2L)D/+ Dsim_Int(2L)S/+ Nup160<up>EP372</up>',          # Dsim x Dmel hybrid (?).
        1105: 'Df(2R)Dark2',                                                     # One of only four production genotypes associated with a stock: FBst0007156.
        171479: 'Df(1)52 P{w<up>+</up>4&Dgr;4.3} lncRNA:roX1<up>ex6</up> lncRNA:roX2<up>Hsp83.PH</up> | FBab0029971_FBal0099841_FBal0127187_FBtp0016778',
        525357: 'w[*]; betaTub60D[2] Kr[If-1]|CyO',                              # Genotype from stock; genotype_id here is for FB2024_06 only.
    }

    # Elaborate on get_general_data() for the GenotypeHandler.
    def get_general_data(self, session):
        """Extend the method for the GenotypeHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'balancer', 'construct', 'insertion'])
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

    def get_feature_genotypes(self, session):
        """Get genotype components."""
        self.log.info('Get genotype components.')
        # Limit to only current genotypes and their current features.
        filters = (
            Genotype.is_obsolete.is_(False),
            Feature.is_obsolete.is_(False),
        )
        if self.testing:
            filters += (Genotype.genotype_id.in_(self.test_set.keys()), )
        results = session.query(FeatureGenotype).\
            select_from(Genotype).\
            join(FeatureGenotype, (FeatureGenotype.genotype_id == Genotype.genotype_id)).\
            join(Feature, (Feature.feature_id == FeatureGenotype.feature_id)).\
            filter(*filters).\
            distinct()
        fg_counter = 0
        for result in results:
            if result.cgroup in self.fb_data_entities[result.genotype_id].feature_genotypes.keys():
                self.fb_data_entities[result.genotype_id].feature_genotypes[result.cgroup].append(result)
                fg_counter += 1
            else:
                self.fb_data_entities[result.genotype_id].feature_genotypes[result.cgroup] = [result]
                fg_counter += 1
        self.log.info(f'Found {fg_counter} feature_genotype entries for genotypes.')
        return

    # Elaborate on get_datatype_data() for the GenotypeHandler.
    def get_datatype_data(self, session):
        """Extend the method for the GenotypeHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_genotype_fb_curies(session)
        self.get_feature_genotypes(session)
        self.get_entity_cvterms(session)
        self.get_entityprops(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        return

    # Additional sub-methods for synthesize_info().
    def prune_stock_only_genotypes(self):
        """Prune genotypes related only to fly stocks."""
        self.log.info('Prune genotypes related only to fly stocks.')
        # Skip genotypes associated only with fly stocks (these genotypes have no FBgo, no description, but have name).
        # Fly stocks and related genotype info will be exported by some other means.
        stock_only_genotype_ids = []
        for genotype in self.fb_data_entities.values():
            if genotype.fb_curie is None:
                stock_only_genotype_ids.append(genotype.db_primary_id)
        for genotype_id in stock_only_genotype_ids:
            del(self.fb_data_entities[genotype_id])
        self.log.info(f'Pruned {len(stock_only_genotype_ids)} stock-only genotypes from the dataset.')
        return

    def synthesize_genotype_components(self):
        """Determine genotype feature components to report and their zygosity."""
        self.log.info('Determine genotype feature components to report and their zygosity.')
        component_counter = 0
        genotype_counter = 0
        for genotype in self.fb_data_entities.values():
            genotype.component_features = {
                'unspecified zygosity': [],
                'homozygous': [],
                'heterozygous': [],
                'hemizygous': [],
            }
            # Handle genotypes from production, which are properly structured.
            if genotype.name is None:
                for cgroup, fg_list in genotype.feature_genotypes.items():
                    cgroup_feature_list = [i.feature_id for i in fg_list]
                    if len(cgroup_feature_list) > 2:
                        self.log.warning(f'{genotype}, cgroup={cgroup} has too many components!')
                    cgroup_feature_set = set(cgroup_feature_list)
                    if len(cgroup_feature_list) == 1:
                        zygosity = 'unspecified zygosity'
                    elif len(cgroup_feature_list) == 2 and len(cgroup_feature_set) == 1:
                        zygosity = 'homozygous'
                    elif len(cgroup_feature_list) == 2 and len(cgroup_feature_set) == 2:
                        zygosity = 'heterozygous'
                    for feature_id in cgroup_feature_set:
                        # Filter out reporting of "bogus symbols" (internal feature placeholders for wild-type alleles).
                        if feature_id in self.feature_lookup.keys():
                            genotype.component_features[zygosity].append(feature_id)
            # Handle stock-only genotypes, which are not properly structured (all feature_genotype associations have cgroup=0 and rank=0).
            else:
                feature_list = [i.feature_id for i in genotype.feature_genotypes[0]]
                feature_count = {}
                for feature_id in feature_list:
                    if feature_id in feature_count.keys():
                        feature_count[feature_id] += 1
                    else:
                        feature_count[feature_id] = 1
                for feature_id, feature_count in feature_count.items():
                    if feature_count == 1:
                        zygosity = 'unspecified zygosity'
                    else:
                        zygosity = 'homozygous'
                    # Filter out reporting of "bogus symbols" (internal feature placeholders for wild-type alleles).
                    if feature_id in self.feature_lookup.keys():
                        genotype.component_features[zygosity].append(feature_id)
                        component_counter += 1
            if genotype.component_features:
                genotype_counter += 1
        self.log.info(f'Found {component_counter} components for {genotype_counter} genotypes for export.')
        return

    def synthesize_genotype_ncbi_taxon_id(self):
        """Determine the NCBITaxon ID for FB genotypes."""
        self.log.info('Determine the NCBITaxon ID for FB genotypes.')
        # Need to handle many different cases.
        # 1. Genotypes with features - need to find Drosophilids (or hybrids of them).
        # 2. Genotypes for non-Dros stocks.
        # for genotype in self.fb_data_entities.values():
        #     feature_id_list = []
        #     for cgroup_feature_genotype_list in genotype.feature_genotypes.values():
        #         cgroup_feature_id_list = [i.feature_id for i in cgroup_feature_genotype_list]
        #         feature_id_list.extend(cgroup_feature_id_list)
        #     for feature_id in feature_id_list:
        #         if feature_id in self.feature_lookup.keys():
        #             feature = self.feature_lookup[feature_id]
        #             organism = self.organism_lookup[feature['organism_id']]
        #             if organism['abbreviation'] == 'Dmel':
        #                 continue
        #             if feature['uniquename'].startswith('FBti')
        #     # Catch cases where the FB data entity has no organism_id: e.g., genotype.
        #     # These datatypes will have special handling in the datatype-specific handlers.
        #     try:
        #         organism_id = fb_data_entity.chado_obj.organism_id
        #     except AttributeError:
        #         self.log.warning(f'No organism_id for {fb_data_entity}.')
        #         return
        #     # Catch cases where the FB data entity has no corresponding NCBITaxon ID.
        #     fb_data_entity.ncbi_taxon_id = self.organism_lookup[organism_id]['taxon_curie']
        #     if fb_data_entity.ncbi_taxon_id == 'NCBITaxon:32644':
        #         self.log.warning(f'{fb_data_entity} has "unidentified" NCBITaxon ID.')
        return

    # Elaborate on synthesize_info() for the GenotypeHandler.
    def synthesize_info(self):
        """Extend the method for the GenotypeHandler."""
        super().synthesize_info()
        self.prune_stock_only_genotypes()
        self.synthesize_genotype_components()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        return

    # Additional sub-methods for map_fb_data_to_alliance().
    def map_genotype_basic(self):
        """Map basic FlyBase genotype data to the Alliance object."""
        self.log.info('Map basic FlyBase genotype data to the Alliance object.')
        counter = 0
        for genotype in self.fb_data_entities.values():
            agr_genotype = self.agr_export_type()
            agr_genotype.obsolete = genotype.chado_obj.is_obsolete
            agr_genotype.mod_entity_id = f'FB:{genotype.fb_curie}'
            agr_genotype.taxon_curie = genotype.ncbi_taxon_id
            agr_genotype.name = genotype.uniquename
            agr_genotype.subtype_name = 'genotype'
            genotype.linkmldto = agr_genotype
            counter += 1
        self.log.info(f'{counter} genotypes could be mapped to the Alliance LinkML model.')
        return

    def map_genotype_components(self):
        """Map genotype components."""
        self.log.info('Map genotype components.')
        for genotype in self.fb_data_entities.values():
            for zygosity, component_feature_id_list in genotype.component_features.items():
                for feature_id in component_feature_id_list:
                    component_curie = self.feature_lookup[feature_id]['curie']
                    agm_component = agr_datatypes.AffectedGenomicModelComponentDTO(component_curie, zygosity)
                    genotype.linkmldto.component_dtos.append(agm_component.dict_export())
        return

    # Elaborate on map_fb_data_to_alliance() for the GenotypeHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GenotypeHandler."""
        super().map_fb_data_to_alliance()
        self.map_genotype_basic()
        self.map_genotype_components()
        # self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_pubs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        return
