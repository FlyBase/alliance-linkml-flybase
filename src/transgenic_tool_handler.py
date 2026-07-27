"""Module:: tool_handler.

Synopsis:
    A data handler that exports FlyBase data for experimental tools to Alliance
    TransgenicTool LinkML objects.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

# import csv
# import re
from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class ExperimentalToolHandler(FeatureHandler):
    """This object gets, synthesizes and filters exp tool data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the ExperimentalToolHandler object."""
        super().__init__(log, testing)
        self.datatype = 'tool'
        self.fb_export_type = fb_datatypes.FBTool
        self.agr_export_type = agr_datatypes.TransgenicToolDTO
        self.primary_export_set = 'transgenic_tool_ingest_set'
        self.tool_associations = []    # FBRelationship carriers, one per exported tool-tool association.
        self.tool_tool_rels = {}       # Lists of FBRelationship objects keyed by (subject_id, object_id).

    test_set = {
        'FBto0000001': 'C-Cerulean',  # First one
        'FBto0000027': 'EGFP',
        'FBto0000417': 'sgGFP',
        'FBto0000921': 'Sapphire',
        'FBto0000606': 'AflIII',    # Has UniProtKB:E3VX96
    }

    transgenic_tool_prop_to_note_mapping = {
        'description': ('summary', 'note_dtos'),
        'misc': ('comment', 'note_dtos'),
        # At the moment, just for code development. (line below)
        # 'internal_notes': ('internal_note', 'note_dtos'),
    }

    def get_general_data(self, session):
        """Extend the method for the AlleleHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_feature_lookup(session, feature_types=['tool'])
        return

    def get_datatype_data(self, session):
        """Extend the method for the ExperimentalToolHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        # Query both directions: a compatible_tool relationship is only ever attached to the tool
        # playing the queried role, so a single-role query leaves the other tool of the pair without
        # the relationship (FTA-132).
        for role in ('subject', 'object'):
            self.get_entity_relationships(session, role, rel_type='compatible_tool',
                                          entity_type='engineered_region', entity_regex=self.regex['tool'])

        # self.build_feature_lookup(session)

    def map_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for a FlyBase entity."""
        self.log.info('Map secondary IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.linkmldto is None:
                continue
            secondary_id_dtos = []
            for secondary_id in fb_data_entity.alt_fb_ids:
                secondary_id_dtos.append(secondary_id)
            sec_id_list = getattr(fb_data_entity.linkmldto, slot_name)
            sec_id_list.extend(secondary_id_dtos)
        return

    # Elaborate on map_fb_data_to_alliance() for the ExpToolHandler.
    def map_fb_data_to_alliance(self):
        """Extend the method for the GeneHandler."""
        super().map_fb_data_to_alliance()
        self.map_tool_basic()
        self.map_synonyms()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_entity_props_to_notes('transgenic_tool_prop_to_note_mapping')
        self.map_secondary_ids('secondary_identifiers')
        self.map_tool_associations()
        # Cascade chado-obsolete -> internal=True (matches every other handler).
        self.flag_internal_fb_entities('fb_data_entities')

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_tool_basic(self):
        """Map basic FlyBase transgenic tool data to the Alliance LinkML object."""
        self.log.info('Map basic transgenic tool info to Alliance object.')
        for tool in self.fb_data_entities.values():
            agr_tool = self.agr_export_type()
            agr_tool.obsolete = tool.chado_obj.is_obsolete
            agr_tool.primary_external_id = f'FB:{tool.uniquename}'
            tool.linkmldto = agr_tool
        return

    def synthesize_tool_associations(self):
        """Collect compatible_tool relationships from both directions and group them by tool pair."""
        self.log.info('Synthesize transgenic tool associations.')
        # Each feature_relationship is fetched once per role, producing two distinct FBRelationship
        # objects for the same chado row, so dedupe on the feature_relationship_id.
        rels_by_id = {}
        rels_per_role = {'subject': 0, 'object': 0}
        tools_with_rels = set()
        for tool in self.fb_data_entities.values():
            for role in rels_per_role.keys():
                role_rels = tool.recall_relationships(self.log, entity_role=role, rel_types='compatible_tool')
                rels_per_role[role] += len(role_rels)
                if role_rels:
                    tools_with_rels.add(tool.db_primary_id)
                for tool_rel in role_rels:
                    rels_by_id[tool_rel.db_primary_id] = tool_rel
        for role, role_count in rels_per_role.items():
            self.log.info(f'Recalled {role_count} compatible_tool relationships where the tool is the {role}.')
        for tool_rel in rels_by_id.values():
            tool_tool_key = (tool_rel.chado_obj.subject_id, tool_rel.chado_obj.object_id)
            try:
                self.tool_tool_rels[tool_tool_key].append(tool_rel)
            except KeyError:
                self.tool_tool_rels[tool_tool_key] = [tool_rel]
        msg = f'Found {len(rels_by_id)} distinct compatible_tool feature_relationships in '
        msg += f'{len(self.tool_tool_rels)} subject-object combinations, for {len(tools_with_rels)} tools.'
        self.log.info(msg)
        return

    # Elaborate on synthesize_info() for the Handler.
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_tool_associations()

    def build_tool_association(self, chado_obj, subject_curie, object_curie, pub_ids, pub_curies, is_obsolete):
        """Build an FBRelationship carrier holding one TransgenicToolAssociationDTO.

        Each exported association needs its own carrier object, because generate_export_dict() reads
        linkmldto (and for_export/internal_reasons/export_warnings) off one object per output row.

        Args:
            chado_obj (FeatureRelationship): The chado relationship the association derives from.
            subject_curie (str): The FB curie of the subject tool.
            object_curie (str): The FB curie of the object tool.
            pub_ids (list): The chado pub_ids supporting the association.
            pub_curies (list): The FB curies of the supporting pubs.
            is_obsolete (bool): True if either tool of the pair is obsolete in chado.

        Returns:
            An FBRelationship object with its linkmldto set.

        """
        carrier = fb_datatypes.FBRelationship(chado_obj, 'feature_relationship')
        carrier.pubs = pub_ids
        rel_dto = agr_datatypes.TransgenicToolAssociationDTO(
            subject_curie, object_curie,
            pub_curies, False, 'compatible_tool')
        if is_obsolete is True:
            rel_dto.obsolete = True
            rel_dto.internal = True
        carrier.linkmldto = rel_dto
        return carrier

    def map_tool_associations(self):
        """Map transgenic tool associations to Alliance objects, in both directions."""
        self.log.info('Map tool associations to Alliance object.')
        # Merge the chado rows into one record per unordered tool pair, so that a pair stored in both
        # directions contributes a single set of pubs rather than two competing ones.
        pair_records = {}
        for tool_tool_key, tool_tool_rels in self.tool_tool_rels.items():
            pair_key = tuple(sorted(tool_tool_key))
            try:
                pair_record = pair_records[pair_key]
            except KeyError:
                pair_record = {'rels': [], 'pub_ids': [], 'chado_directions': set()}
                pair_records[pair_key] = pair_record
            pair_record['rels'].extend(tool_tool_rels)
            pair_record['chado_directions'].add(tool_tool_key)
            for tool_tool_rel in tool_tool_rels:
                pair_record['pub_ids'].extend(tool_tool_rel.pubs)

        exported_pairs = set()    # (subject_curie, object_curie) tuples already exported.
        both_ways_in_chado = 0
        reciprocal_counter = 0
        skipped_counter = 0
        for pair_key, pair_record in pair_records.items():
            if self.testing:
                self.log.debug(f'Mapping {pair_key} to Alliance object. {pair_record["rels"]}')
            try:
                feat_a = self.feature_lookup[pair_key[0]]
                feat_b = self.feature_lookup[pair_key[1]]
            except KeyError:
                self.log.warning(f'Skip compatible_tool pair {pair_key}: a feature is missing from the tool lookup.')
                skipped_counter += 1
                continue
            curie_a = f'FB:{feat_a["uniquename"]}'
            curie_b = f'FB:{feat_b["uniquename"]}'
            is_obsolete = feat_a['is_obsolete'] is True or feat_b['is_obsolete'] is True
            all_pub_ids = pair_record['pub_ids']
            # NOTE: pub 383755 | FlyBase Experimental Tool information Is the only one used
            # for tools. But not in lookup pub curies!
            # So pub_curies will be empty.
            pub_curies = self.lookup_pub_curies(all_pub_ids)
            if len(pair_record['chado_directions']) == 2:
                both_ways_in_chado += 1
            elif curie_a != curie_b:
                reciprocal_counter += 1
            chado_obj = pair_record['rels'][0].chado_obj
            # Export both directions: the Alliance model treats the association as directed, so a
            # tool only shows its compatible partner if it is named as the subject (FTA-132).
            for subject_curie, object_curie in ((curie_a, curie_b), (curie_b, curie_a)):
                if (subject_curie, object_curie) in exported_pairs:
                    continue
                self.tool_associations.append(
                    self.build_tool_association(chado_obj, subject_curie, object_curie,
                                                all_pub_ids, pub_curies, is_obsolete))
                exported_pairs.add((subject_curie, object_curie))
        self.log.info(f'Found {len(pair_records)} unique tool-tool pairs, {both_ways_in_chado} of which chado '
                      'stores in both directions.')
        self.log.info(f'Synthesized {reciprocal_counter} reciprocal associations not present in chado.')
        if skipped_counter:
            self.log.warning(f'Skipped {skipped_counter} tool-tool pairs missing from the tool lookup.')
        self.log.info(f'Generated {len(self.tool_associations)} tool-tool unique associations.')
        return

    # Elaborate on query_chado_and_export() for the TransgenicToolHandler.
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the TransgenicToolHandler."""
        super().query_chado_and_export(session)
        self.generate_export_dict(self.tool_associations, 'tool_association_ingest_set')

        return
