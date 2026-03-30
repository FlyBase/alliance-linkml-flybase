"""Module:: gene_group_handler.

Synopsis:
    A data handler that exports FlyBase gene group (FBgg) data to Alliance
    FunctionalGeneSet LinkML objects.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

from logging import Logger
from os import environ
from sqlalchemy.orm import aliased
from harvdev_utils.reporting import (
    Cvterm, Feature, Grp, Grpmember, GrpmemberPub, FeatureGrpmember
)
import agr_datatypes
import fb_datatypes
from entity_handler import PrimaryEntityHandler


class GeneGroupHandler(PrimaryEntityHandler):
    """This object gets, synthesizes and filters gene group data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the GeneGroupHandler object."""
        super().__init__(log, testing)
        self.datatype = 'grp'
        self.fb_export_type = fb_datatypes.FBGeneGroup
        self.agr_export_type = agr_datatypes.FunctionalGeneSetDTO
        self.primary_export_set = 'functional_gene_set_ingest_set'
        self.gene_member_associations = []
        self.gene_member_data = {}
        self.export_data_for_tsv = []

    test_set = {
        'FBgg0000113': 'ARP',
        'FBgg0000960': 'ALP',
        'FBgg0000162': 'WNT',
    }

    gene_group_prop_to_note_mapping = {
        'gg_description': ('summary', 'note_dtos'),
    }

    # Relationship type mappings: grp_relationship type → output field.
    parent_rel_types = ['parent_grp']
    related_rel_types = ['undefined_grp']

    # Gene-group-specific xref db mappings (FTA-171).
    gene_group_db_dict = {
        'HGNC-GG1': 'HGNC_Group',
        'WB-GG': 'WBGeneClass',
        'ComplexPortal': 'ComplexPortal',
    }
    gene_group_page_area = {
        'HGNC_Group': 'gene_group',
        'WBGeneClass': 'gene_class',
        'ComplexPortal': 'default',
    }

    # Elaborate on get_general_data() for the GeneGroupHandler.
    def get_general_data(self, session):
        """Extend the method for the GeneGroupHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_organism_lookup(session)
        return

    # Elaborate on get_datatype_data() for the GeneGroupHandler.
    def get_datatype_data(self, session):
        """Extend the method for the GeneGroupHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_cvterms(session)
        self.get_entity_relationships(session, 'subject')
        self.get_entity_relationships(session, 'object')
        self.get_entity_timestamps(session)
        if environ.get('EXPORT_GG_XREFS') == 'YES':
            self.get_entity_xrefs(session)
        else:
            self.log.info('Skipping gene group xrefs (EXPORT_GG_XREFS not set to "YES").')
        self.get_gene_members(session)
        return

    def get_gene_members(self, session):
        """Get gene members of gene groups via FeatureGrpmember + Grpmember."""
        self.log.info('Get gene members for gene groups.')
        grp_ids = list(self.fb_data_entities.keys())
        if not grp_ids:
            self.log.info('No gene groups loaded; skipping gene member query.')
            return

        gene_type = aliased(Cvterm, name='gene_type')
        results = session.query(
            Feature.uniquename,
            Grp.uniquename,
            Grpmember.grpmember_id).\
            select_from(FeatureGrpmember).\
            join(Grpmember, Grpmember.grpmember_id == FeatureGrpmember.grpmember_id).\
            join(Grp, Grp.grp_id == Grpmember.grp_id).\
            join(Feature, Feature.feature_id == FeatureGrpmember.feature_id).\
            join(gene_type, gene_type.cvterm_id == Feature.type_id).\
            filter(
                Grpmember.grp_id.in_(grp_ids),
                Feature.uniquename.op('~')(self.regex['gene']),
                gene_type.name == 'gene',
                Feature.is_obsolete.is_(False),
                Grp.is_obsolete.is_(False)).\
            distinct()

        # Build grpmember_id → (gene_uniquename, grp_uniquename) mapping.
        grpmember_ids = []
        # Key: (gene_uniquename, grp_uniquename) → set of pub_ids.
        self.gene_member_data = {}
        grpmember_to_key = {}
        for gene_uname, grp_uname, grpmember_id in results:
            key = (gene_uname, grp_uname)
            if key not in self.gene_member_data:
                self.gene_member_data[key] = set()
            grpmember_ids.append(grpmember_id)
            grpmember_to_key[grpmember_id] = key

        # Get pubs for each membership.
        if grpmember_ids:
            pub_results = session.query(
                GrpmemberPub.grpmember_id,
                GrpmemberPub.pub_id
            ).filter(
                GrpmemberPub.grpmember_id.in_(grpmember_ids)
            ).distinct()
            for grpmember_id, pub_id in pub_results:
                if grpmember_id in grpmember_to_key:
                    self.gene_member_data[grpmember_to_key[grpmember_id]].add(pub_id)

        self.log.info(f'Found {len(self.gene_member_data)} gene-to-gene-group memberships.')
        return

    # Elaborate on synthesize_info() for the GeneGroupHandler.
    def synthesize_info(self):
        """Extend the method for the GeneGroupHandler."""
        super().synthesize_info()
        # Note: synthesize_ncbi_taxon_id() is not called here because the grp
        # table has no organism_id. Taxon is hardcoded in map_gene_group_basic().
        self.synthesize_synonyms()
        self.synthesize_pubs()
        return

    def process_for_tsv_export(self):
        """Process gene group data for export to TSV."""
        self.log.info('Process gene groups for TSV export.')
        counter = 0
        for gene_group in self.fb_data_entities.values():
            if gene_group.linkmldto is None:
                continue
            dto = gene_group.linkmldto
            # Collect synonyms.
            synonyms = []
            for syno in dto.set_synonym_dtos:
                synonyms.append(syno.get('format_text', ''))
            # Collect GO terms.
            go_mf = ' | '.join(dto.set_go_mf_term_curies) if dto.set_go_mf_term_curies else ''
            go_bp = ' | '.join(dto.set_go_bp_term_curies) if dto.set_go_bp_term_curies else ''
            go_cc = ' | '.join(dto.set_go_cc_term_curies) if dto.set_go_cc_term_curies else ''
            # Collect parent/related groups.
            parents = ' | '.join(dto.parent_set_identifiers) if dto.parent_set_identifiers else ''
            related = ' | '.join(dto.related_set_identifiers) if dto.related_set_identifiers else ''
            # Collect gene members.
            gene_members = []
            for (gene_uname, grp_uname) in self.gene_member_data.keys():
                if grp_uname == gene_group.uniquename:
                    gene_members.append(f'FB:{gene_uname}')
            # Collect notes.
            notes = []
            for note in dto.note_dtos:
                notes.append(note.get('free_text', ''))
            tsv_row = {
                'gene_group_id': f'FB:{gene_group.uniquename}',
                'symbol': dto.symbol or '',
                'full_name': dto.full_name or '',
                'synonyms': ' | '.join(synonyms),
                'description': ' | '.join(notes),
                'go_molecular_function': go_mf,
                'go_biological_process': go_bp,
                'go_cellular_component': go_cc,
                'parent_groups': parents,
                'related_groups': related,
                'gene_members': ' | '.join(gene_members),
            }
            self.export_data_for_tsv.append(tsv_row)
            counter += 1
        self.log.info(f'Generated {counter} gene group TSV rows.')
        return

    # Mapping methods.
    def map_gene_group_basic(self):
        """Map basic FlyBase gene group data to the Alliance FunctionalGeneSetDTO."""
        self.log.info('Map basic gene group info.')
        for gene_group in self.fb_data_entities.values():
            agr_gg = self.agr_export_type()
            agr_gg.obsolete = gene_group.chado_obj.is_obsolete
            agr_gg.primary_external_id = f'FB:{gene_group.uniquename}'
            # Gene groups are all Drosophila melanogaster.
            agr_gg.taxon_curie = 'NCBITaxon:7227'
            # Set symbol and full_name as plain strings.
            if gene_group.curr_fb_symbol:
                agr_gg.symbol = gene_group.curr_fb_symbol
            else:
                agr_gg.symbol = gene_group.name
            if gene_group.curr_fb_fullname:
                agr_gg.full_name = gene_group.curr_fb_fullname
            gene_group.linkmldto = agr_gg
        return

    def map_gene_group_synonyms(self):
        """Map gene group synonyms to set_synonym_dtos."""
        self.log.info('Map gene group synonyms.')
        for gene_group in self.fb_data_entities.values():
            if gene_group.linkmldto is None:
                continue
            synonym_dtos = []
            for syno_dict in gene_group.synonym_dict.values():
                # Skip current symbol and current fullname (they go into
                # symbol/full_name string fields, not set_synonym_dtos).
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] == 'nomenclature_symbol':
                    if syno_dict['format_text'] == gene_group.name:
                        continue
                if syno_dict['is_current'] is True and syno_dict['name_type_name'] == 'full_name':
                    continue
                name_dto = agr_datatypes.NameSlotAnnotationDTO(
                    syno_dict['name_type_name'],
                    syno_dict['format_text'],
                    syno_dict['display_text'],
                    syno_dict['pub_curies']
                ).dict_export()
                name_dto['internal'] = syno_dict['is_internal']
                synonym_dtos.append(name_dto)
            gene_group.linkmldto.set_synonym_dtos = synonym_dtos
        return

    def map_gene_group_notes(self):
        """Map gene group props to Alliance notes."""
        self.log.info('Map gene group notes.')
        self.map_entity_props_to_notes('gene_group_prop_to_note_mapping')
        return

    def map_gene_group_go_terms(self):
        """Map gene group GO term annotations to set_go_*_term_curies."""
        self.log.info('Map gene group GO terms.')
        go_cv_to_slot = {
            'molecular_function': 'set_go_mf_term_curies',
            'biological_process': 'set_go_bp_term_curies',
            'cellular_component': 'set_go_cc_term_curies',
        }
        for gene_group in self.fb_data_entities.values():
            if gene_group.linkmldto is None:
                continue
            for cv_name, slot_name in go_cv_to_slot.items():
                go_curies = []
                cvt_anno_ids = gene_group.cvt_anno_ids_by_cv.get(cv_name, [])
                for cvt_anno_id in cvt_anno_ids:
                    cvt_anno = gene_group.cvt_annos_by_id[cvt_anno_id]
                    db_name = cvt_anno.chado_obj.cvterm.dbxref.db.name
                    accession = cvt_anno.chado_obj.cvterm.dbxref.accession
                    go_curie = f'{db_name}:{accession}'
                    if go_curie not in go_curies:
                        go_curies.append(go_curie)
                setattr(gene_group.linkmldto, slot_name, go_curies)
        return

    def map_gene_group_relationships(self):
        """Map gene group relationships to parent_set_identifiers and related_set_identifiers."""
        self.log.info('Map gene group relationships.')
        # Build grp_id → uniquename lookup from loaded entities.
        grp_id_lookup = {}
        for grp_id, entity in self.fb_data_entities.items():
            grp_id_lookup[grp_id] = entity.uniquename
        for gene_group in self.fb_data_entities.values():
            if gene_group.linkmldto is None:
                continue
            parent_ids = []
            related_ids = []
            # Check subject relationships (this group is the subject).
            for rel_type, rel_id_list in gene_group.sbj_rel_ids_by_type.items():
                for rel_id in rel_id_list:
                    rel = gene_group.rels_by_id[rel_id]
                    object_grp_id = rel.chado_obj.object_id
                    if object_grp_id not in grp_id_lookup:
                        self.log.warning(
                            f'Related grp_id {object_grp_id} for {gene_group} '
                            f'(rel_type={rel_type}) not in loaded gene groups.')
                        continue
                    related_curie = f'FB:{grp_id_lookup[object_grp_id]}'
                    if rel_type in self.parent_rel_types:
                        if related_curie not in parent_ids:
                            parent_ids.append(related_curie)
                    elif rel_type in self.related_rel_types:
                        if related_curie not in related_ids:
                            related_ids.append(related_curie)
                    else:
                        self.log.debug(
                            f'Unhandled grp_relationship type "{rel_type}" '
                            f'for {gene_group}.')
            gene_group.linkmldto.parent_set_identifiers = parent_ids
            gene_group.linkmldto.related_set_identifiers = related_ids
        return

    def map_gene_group_xrefs(self):
        """Map gene group cross-references to cross_reference_dtos (FTA-171)."""
        self.log.info('Map gene group cross-references.')
        for gene_group in self.fb_data_entities.values():
            if gene_group.linkmldto is None:
                continue
            cross_reference_dtos = []
            # Add the FB self-xref.
            curie = f'FB:{gene_group.uniquename}'
            page_area = self.page_area_conversion.get(self.datatype, self.datatype)
            fb_xref = agr_datatypes.CrossReferenceDTO(
                'FB', curie, page_area, curie
            ).dict_export()
            cross_reference_dtos.append(fb_xref)
            # Add external xrefs from gene-group-specific databases.
            for xref in gene_group.dbxrefs:
                db_name = xref.dbxref.db.name
                if db_name not in self.gene_group_db_dict:
                    self.log.debug(
                        f'Skipping unrecognized xref db "{db_name}" '
                        f'for {gene_group}.')
                    continue
                prefix = self.gene_group_db_dict[db_name]
                xref_page_area = self.gene_group_page_area.get(prefix, page_area)
                accession = xref.dbxref.accession
                # Strip version/suffix after '#' (e.g., WB-GG "arp#01--10" → "arp").
                if '#' in accession:
                    accession = accession.split('#')[0]
                xref_curie = f'{prefix}:{accession}'
                display_name = xref.dbxref.description if xref.dbxref.description else xref_curie
                xref_dto = agr_datatypes.CrossReferenceDTO(
                    prefix, xref_curie, xref_page_area, display_name
                ).dict_export()
                cross_reference_dtos.append(xref_dto)
            gene_group.linkmldto.cross_reference_dtos = cross_reference_dtos
        return

    def map_gene_member_associations(self):
        """Map gene-to-gene-group memberships to GeneFunctionalGeneSetAssociationDTOs."""
        self.log.info('Map gene member associations.')
        self.gene_member_associations = []
        if not hasattr(self, 'gene_member_data'):
            return
        for (gene_uname, grp_uname), pub_ids in self.gene_member_data.items():
            dto = agr_datatypes.GeneFunctionalGeneSetAssociationDTO()
            dto.gene_identifier = f'FB:{gene_uname}'
            dto.functional_gene_set_identifier = f'FB:{grp_uname}'
            dto.relation_name = 'is_member_of'
            dto.evidence_curies = self.lookup_pub_curies(list(pub_ids))
            # Create a wrapper FBExportEntity for the export pipeline.
            wrapper = fb_datatypes.FBExportEntity()
            wrapper.linkmldto = dto
            wrapper.entity_desc = f'{gene_uname} is_member_of {grp_uname}'
            self.gene_member_associations.append(wrapper)
        self.log.info(f'Created {len(self.gene_member_associations)} gene member association DTOs.')
        return

    # Elaborate on map_fb_data_to_alliance() for the GeneGroupHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneGroupHandler."""
        super().map_fb_data_to_alliance()
        self.map_gene_group_basic()
        self.map_gene_group_synonyms()
        self.map_data_provider_dto()
        self.map_timestamps()
        self.map_gene_group_notes()
        self.map_gene_group_go_terms()
        self.map_gene_group_relationships()
        if environ.get('EXPORT_GG_XREFS') == 'YES':
            self.map_gene_group_xrefs()
        else:
            self.log.info('Skipping gene group xref mapping (EXPORT_GG_XREFS not set to "YES").')
        self.map_gene_member_associations()
        self.flag_internal_fb_entities('fb_data_entities')
        return

    # Elaborate on query_chado_and_export() for the GeneGroupHandler.
    def query_chado_and_export(self, session):
        """Extend the method for the GeneGroupHandler."""
        super().query_chado_and_export(session)
        self.flag_unexportable_entities(
            self.gene_member_associations,
            'gene_functional_gene_set_association_ingest_set')
        self.generate_export_dict(
            self.gene_member_associations,
            'gene_functional_gene_set_association_ingest_set')
        return
