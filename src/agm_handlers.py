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
from harvdev_utils.char_conversions import (
    sgml_to_plain_text, sub_sup_sgml_to_plain_text
)
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
        self.flag_new_additions_and_obsoletes()
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
            agr_strain.primary_external_id = f'FB:{strain.uniquename}'
            agr_strain.taxon_curie = strain.ncbi_taxon_id
            agr_strain.name = strain.name
            agr_strain.subtype_name = 'strain'
            if strain.ncbi_taxon_id != 'NCBITaxon:7227':
                agr_strain.internal = True
                strain.internal_reasons.append('Non-Dmel')
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
        2: 'Ab(1)ZWD16 | FBab0027942',                             # The first genotype in the table.
        452205: 'wg[1]/wg[GBM]',                                   # Compound heterozygous wg[1]/wg[GBM].
        450391: 'wg[l-8]/wg[l-8]',                                 # Homozygous wg[l-8] allele.
        367896: 'Tak1[2]/Tak1[+]',                                 # Simple heterozygous Tak1[2] over wt allele.
        515567: 'Sdc[12]/Sdc[unspecified]',                        # Sdc[12]/Sdc[unspecified].
        167742: 'Psn[-]',                                          # A single unknown mutant allele represented by a bogus symbol; has pheno_comp data.
        416752: 'crb[Y10A]/crb[-]',                                # crb allele with unknown mutant crb (bogus symbol); has pheno_comp data.
        169086: 'dl<up>-</up> dl<up>S70A</up>',                    # A transgenic dl (FBal0095462) with unknown mutant dl (bogus symbol); has pheno_comp data.
        166899: 'Df(2R)173/PCNA[D-292]',                           # PCNA[D-292]/Df(2R)173.
        168332: 'Df(3L)Ez7 hay[nc2.tMa]',                          # Df(3L)Ez7 + hay[nc2.tMa] (diff complementation groups).
        166704: 'shi[EM33] shi[t15]',                              # shi[EM33] + shi[t15] (allele + rescue construct).
        219912: 'dpp[s4] wg[l-17]',                                # dpp[s4], wg[l-17].
        199449: 'Hsap_MAPT[UAS.cAa] Scer_GAL4[smid-C161]',         # GAL4 drives UAS human construct.
        32369: 'Dmau_w[a23]',                                      # Single non-Dmel classical allele -> non-Dmel genotype.
        466842: 'Dsim_mir-983[KO]',                                # Single non-Dmel classical allele -> non-Dmel genotype.
        340500: 'Dsim_Lhr[1] Dsim_Lhr[2+16aa.Tag:HA]',             # non-Dmel classical allele + transgene -> non-Dmel genotype.
        525097: 'Dvir_tra[tra.WT]',                                # Genotype where Dvir tra replaces Dmel tra in Dmel; Dvir allele related to Dmel FBti.
        167097: 'Df(3R)CA3 Scer_GAL4[GMR.PF] pim[UAS.Tag:MYC]',    # Df, GAL4 and UAS.
        515875: 'Dp(1;Y)y[+]/tyn[1]',                              # tyn[1] / Dp(1:Y).
        202134: 'Tp(3;1)P115/pad[1]',                              # pad[1] / Tp(3;1)P115.
        182943: 'Df(3L)emc/+',                                     # Heterozygous deficiency.
        334079: 'Dsim_Int(2L)D/+ Dsim_Int(2L)S/+ Nup160[EP372]',   # Dsim FBabs that are actually Dsim chr parts introgressed into Dmel.
        1105: 'Df(2R)Dark2',                                       # One of only four production genotypes associated with a stock: FBst0007156.
        365272: 'Dsim_Cyp6g1[UAS.cHa]',                            # Genotype carrying one allele of UAS Dsim gene (in Dmel).
        365273: 'Dsim_Cyp6g1[UAS.cHa] Scer_GAL4[Cyp6g1.HR]',       # Genotype carrying one allele of UAS Dsim gene plus GAL4 (in Dmel).
        371290: 'Hsap_MAPT[UAS.cAa] Scer_GAL4[GMR.PU]',            # Genotype of GAL4-driven Hsap construct (in Dmel).
        171479: 'Df(1)52 P{w[+]4&Dgr;4.3} lncRNA:roX1[ex6] lncRNA:roX2[Hsp83.PH] | FBab0029971_FBal0099841_FBal0127187_FBtp0016778',    # Has FBtp - no export.
        294012: 'P{CH1226-43A10} lz[L]',                           # Has FBal and FBtp associated - no export.
        223641: 'Dp1[EP2422] P{hsp26-pt-T}39C-12',                 # Has FBal and FBti associated - no export.
        169272: 'P{wA}4-4 brm[2]',                                 # Has FBti directly related - no export.
        526093: 'daw[d05680]',                                     # A genotype new to FB2025_01 (not in FB2024_06 reference) - test incremental export.
        510779: 'OBSOLETE510779:w[<up>]Ecol_lexA.110]',            # A genotype newly obsoleted in FB2025_01 - test incremental export.

        # 525357: 'w[*]; betaTub60D[2] Kr[If-1]|CyO',                              # Genotype from stock; genotype_id here is for FB2024_06 only.
    }

    # Additional export sets.
    agm_component_associations = []    # A list of AgmAlleleAssociationDTOs.

    # Additional reference info.
    dmel_insertion_allele_ids = []    # feature_ids for alleles related to FBti insertions (associated_with/progenitor).
    transgenic_allele_ids = []        # feature_ids for alleles related to FBtp constructs.
    in_vitro_allele_ids = []          # feature_ids for alleles having "in vitro construct" CV term annotations.
    introgressed_aberr_ids = []       # feature_ids for aberrations of type "introgressed_chromosome".

    # Elaborate on get_general_data() for the GenotypeHandler.
    def get_general_data(self, session):
        """Extend the method for the GenotypeHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_cvterm_lookup(session)
        self.build_organism_lookup(session)
        self.build_feature_lookup(session, feature_types=['aberration', 'allele', 'balancer', 'construct', 'insertion', 'bogus symbol'])
        self.get_dmel_insertion_allele_ids(session)
        self.get_transgenic_allele_ids(session)
        self.get_in_vitro_allele_ids(session)
        self.get_introgressed_aberration_ids(session)
        return

    # Additional sub-methods for get_datatype_data().
    # Note that for genotypes, the "uniquename" is not the FBgo accession, so need to get these from genotype_dbxref table.
    # Also note that the FBgo IDs are currently internal (not exposed publicly).
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
            del self.fb_data_entities[genotype_id]
        self.log.info(f'Pruned {len(stock_only_genotype_ids)} stock-only genotypes from the dataset.')
        return

    def clean_up_genotype_name(self):
        """Clean up genotype name."""
        self.log.info('Clean up genotype name.')
        for genotype in self.fb_data_entities.values():
            genotype.name = sub_sup_sgml_to_plain_text(genotype.name)
            genotype.name = sgml_to_plain_text(genotype.name)
            genotype.entity_desc = f'{genotype.name} (genotype_id={genotype.db_primary_id}, {genotype.fb_curie})'
        return

    def synthesize_genotype_ncbi_taxon_id(self):
        """Determine the NCBITaxon ID for FB genotypes."""
        self.log.info('Determine the NCBITaxon ID for FB genotypes.')
        non_dmel_genotype_counter = 0
        unspecified_species_genotype_counter = 0
        for genotype in self.fb_data_entities.values():
            feature_id_list = []
            organism_id_list = []
            for cgroup_feature_genotype_list in genotype.feature_genotypes.values():
                cgroup_feature_id_list = [i.feature_id for i in cgroup_feature_genotype_list]
                feature_id_list.extend(cgroup_feature_id_list)
            feature_id_set = set(feature_id_list)
            # For each feature, determine the relevant organism.
            for feature_id in feature_id_set:
                feature = self.feature_lookup[feature_id]
                org_id = self.feature_lookup[feature_id]['organism_id']
                is_of_drosophilid = self.organism_lookup[org_id]['is_drosophilid']
                # Skip bogus symbols (relevant features should be in the feature_lookup dict).
                if feature_id not in self.feature_lookup.keys() or self.feature_lookup[feature_id]['type'] == 'bogus symbol':
                    continue
                # Skip transgenic constructs.
                elif feature['uniquename'].startswith('FBtp'):
                    continue
                # Record FBti organism.
                elif feature['uniquename'].startswith('FBti'):
                    organism_id_list.append(org_id)
                # Record FBab organism (unless introgressed).
                elif feature['uniquename'].startswith('FBab'):
                    if feature_id in self.introgressed_aberr_ids:
                        continue
                    else:
                        organism_id_list.append(org_id)
                # Record FBba organism.
                elif feature['uniquename'].startswith('FBba'):
                    organism_id_list.append(org_id)
                # Process alleles.
                elif feature['uniquename'].startswith('FBal'):
                    # Skip transgenic alleles: related to FBtp constructs.
                    if feature_id in self.transgenic_allele_ids:
                        continue
                    # Skip transgenic alleles: annotated as "in vitro construct".
                    elif feature_id in self.in_vitro_allele_ids:
                        continue
                    # Skip non-Drosophilid alleles.
                    elif is_of_drosophilid is False:
                        continue
                    # Record Dmel for Dmel insertions.
                    elif feature_id in self.dmel_insertion_allele_ids:
                        organism_id_list.append(1)    # i.e., append organism_id=1 (represents Dmel).
                    # Record allele organism: by process of elimination, allele must be non-Dmel Drosophilid classical/insertional.
                    else:
                        organism_id_list.append(org_id)
            # Analyze the set of organism_ids left.
            organism_id_list = list(set(organism_id_list))
            # If there are no components that provide definitive organism info, the default Dmel persists.
            if len(organism_id_list) == 0:
                pass    # The default.
            # If there is non-ambiguous definitive organism info, use that for the genotype taxon.
            elif len(organism_id_list) == 1:
                org_id = organism_id_list[0]
                genotype.ncbi_taxon_id = self.organism_lookup[org_id]['taxon_curie']
                if org_id != 1:
                    non_dmel_genotype_counter += 1
            # If organism info is ambiguous, the taxon is given as "unidentified".
            else:
                genotype.ncbi_taxon_id = 'NCBITaxon:32644'
                unspecified_species_genotype_counter += 1
        self.log.info(f'Found {non_dmel_genotype_counter} non-Dmel Drosophilid genotypes.')
        self.log.info(f'For {unspecified_species_genotype_counter} genotypes, taxon is unidentified.')
        return

    def synthesize_genotype_components(self):
        """Determine genotype feature components to report and their zygosity."""
        self.log.info('Determine genotype feature components to report and their zygosity.')
        component_counter = 0
        genotype_counter = 0
        for genotype in self.fb_data_entities.values():
            # Fail-safe filter against processing stock-only genotypes.
            if genotype.fb_curie is None:
                continue
            genotype_counter += 1
            genotype.component_features = {
                'unspecified zygosity': [],
                'homozygous': [],
                'simple heterozygous': [],
                'compound heterozygous': [],
                # 'hemizygous': [],
            }
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
                    zygosity = 'compound heterozygous'
                    for feature_id in cgroup_feature_set:
                        if self.feature_lookup[feature_id]['type'] == 'bogus symbol' and '+' in self.feature_lookup[feature_id]['name']:
                            zygosity = 'simple heterozygous'
                for feature_id in cgroup_feature_set:
                    # Filter out reporting of "bogus symbols" (internal feature placeholders for wild-type alleles).
                    if feature_id in self.feature_lookup.keys() and self.feature_lookup[feature_id]['type'] != 'bogus symbol':
                        genotype.component_features[zygosity].append(feature_id)
                        component_counter += 1
        self.log.info(f'Found {component_counter} components for {genotype_counter} genotypes for export.')
        return

    def flag_non_compliant_genotypes(self):
        """Flag non-Alliance-compliant genotypes."""
        self.log.info('Flag non-Alliance-compliant genotypes.')
        compliant_counter = 0
        non_compliant_counter = 0
        for genotype in self.fb_data_entities.values():
            if genotype.is_obsolete is True:
                genotype.export_warnings.append('Suppress obsolete genotypes from Alliance export')
                genotype.for_export = False
            if 'alliance_compliant' in genotype.cvt_anno_ids_by_term.keys():
                genotype.is_alliance_compliant = True
                compliant_counter += 1
            else:
                genotype.is_alliance_compliant = False
                non_compliant_counter += 1
        self.log.info(f'Found {compliant_counter} Alliance-compliant genotypes.')
        self.log.info(f'Found {non_compliant_counter} Alliance-compliant genotypes.')
        return

    # Elaborate on synthesize_info() for the GenotypeHandler.
    def synthesize_info(self):
        """Extend the method for the GenotypeHandler."""
        super().synthesize_info()
        self.prune_stock_only_genotypes()
        self.flag_new_additions_and_obsoletes()
        self.clean_up_genotype_name()
        self.synthesize_genotype_ncbi_taxon_id()
        self.synthesize_genotype_components()
        self.synthesize_secondary_ids()
        self.synthesize_synonyms()
        self.synthesize_pubs()
        self.flag_non_compliant_genotypes()
        return

    # Additional sub-methods for map_fb_data_to_alliance().
    def map_genotype_basic(self):
        """Map basic FlyBase genotype data to the Alliance object."""
        self.log.info('Map basic FlyBase genotype data to the Alliance object.')
        counter = 0
        for genotype in self.fb_data_entities.values():
            agr_genotype = self.agr_export_type()
            agr_genotype.obsolete = genotype.chado_obj.is_obsolete
            agr_genotype.primary_external_id = f'FB:{genotype.fb_curie}'
            agr_genotype.taxon_curie = genotype.ncbi_taxon_id
            agr_genotype.name = genotype.name
            agr_genotype.subtype_name = 'genotype'
            if genotype.ncbi_taxon_id != 'NCBITaxon:7227':
                agr_genotype.internal = True
                genotype.internal_reasons.append('Non-Dmel')
            if genotype.is_alliance_compliant is False:
                agr_genotype.internal = True
                genotype.internal_reasons.append('Non-Alliance-compliant')
            genotype.linkmldto = agr_genotype
            counter += 1
        self.log.info(f'{counter} genotypes could be mapped to the Alliance LinkML model.')
        return

    def map_genotype_components(self):
        """Map genotype components."""
        self.log.info('Map genotype components.')
        for genotype in self.fb_data_entities.values():
            if genotype.for_export is False:
                continue
            for zygosity, component_feature_id_list in genotype.component_features.items():
                for feature_id in component_feature_id_list:
                    geno_allele_rel = fb_datatypes.FBExportEntity()
                    component_curie = self.feature_lookup[feature_id]['curie']
                    # Do not report FBtp construct or FBba balancer components for genotypes.
                    # FBal alleles, FBti insertions, and FBab aberrations are ok.
                    if component_curie.startswith('FB:FBtp') or component_curie.startswith('FB:FBba'):
                        continue
                    geno_allele_rel.linkmldto = agr_datatypes.AgmAlleleAssociationDTO(genotype.linkmldto.primary_external_id,
                                                                                      component_curie, zygosity)
                    self.agm_component_associations.append(geno_allele_rel)
        return

    # Elaborate on map_fb_data_to_alliance() for the GenotypeHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GenotypeHandler."""
        super().map_fb_data_to_alliance()
        self.map_genotype_basic()
        self.map_genotype_components()
        # self.map_synonyms()    # Suppressed until AGM has proper support for synonyms.
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_pubs()
        self.map_timestamps()
        self.map_secondary_ids('agm_secondary_id_dtos')
        self.flag_internal_fb_entities('fb_data_entities')
        self.flag_internal_fb_entities('agm_component_associations')
        return

    # Elaborate on query_chado_and_export() for the GenotypeHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the GenotypeHandler."""
        super().query_chado_and_export(session)
        self.flag_unexportable_entities(self.agm_component_associations, 'agm_allele_association_ingest_set')
        self.generate_export_dict(self.agm_component_associations, 'agm_allele_association_ingest_set')
        return
