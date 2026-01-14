# Map data form FB to Alliance
    Convert the fb_data_entities got DTO structured data (json)

     Main config is in the handler

    # Cassette_handler.py
    class CassetteHandler(FeatureHandler):
    """This object gets, synthesizes and filters cassette data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the CassetteHandler object."""
        super().__init__(log, testing)
        self.datatype = 'cassette'
        self.fb_export_type = fb_datatypes.FBCassette
        self.agr_export_type = agr_datatypes.CassetteDTO
        self.primary_export_set = 'cassette_ingest_set'
        self.incremental_update = False

#### Cassette_handler.py
        def map_fb_data_to_alliance(self):
            """Extend the method for the GeneHandler."""
            self.map_cassette_basic()
            self.map_synonyms()
            self.map_data_provider_dto()
            self.map_entity_props_to_notes('cassette_prop_to_note_mapping')
            self.map_secondary_ids('secondary_identifiers')
            self.map_cassette_associations()

#### map_cassette_basic
    Note: self.agr_export_type() is CassetteDTO(ReagentDTO)
          So we are creating a DTO object for each fb_data_entitie
          and storing it in linkmldto

        def map_cassette_basic(self):
        """Map basic FlyBase transgenic cassette data to the Alliance LinkML object."""
        self.log.info('Map basic cassette info to Alliance object.')
        for cass in self.fb_data_entities.values():
            agr_cass = self.agr_export_type()
            agr_cass.obsolete = cass.chado_obj.is_obsolete
            agr_cass.primary_external_id = f'FB:{cass.uniquename}'
            cass.linkmldto = agr_cass

#### entity_handler.py (So general for all)
    1) So add synonyms to fb_data_entities in linkmldto.synonym_bin

        def map_synonyms(self):
            ...
            setattr(fb_data_entity.linkmldto, linkml_synonym_slots['synonym_bin'], linkml_synonym_bins['synonym_bin'])

    2) Add data_provider_dto
       
        def map_data_provider_dto(self):
            """Return the DataProviderDTO for the FB data entity."""
            # Note - this method is depends on previous determination of fb_data_entity.curr_fb_symbol by synthesize_synonyms(), if applicable.
            for fb_data_entity in self.fb_data_entities.values():
                referenced_curie = f'FB:{fb_data_entity.uniquename}'
                display_name = fb_data_entity.curr_fb_symbol
                page_area = self.datatype
                dp_xref = agr_datatypes.CrossReferenceDTO('FB', referenced_curie, page_area, display_name).dict_export()
                fb_data_entity.linkmldto.data_provider_dto = agr_datatypes.DataProviderDTO(dp_xref).dict_export()

    3) map_entity_props_to_notes('cassette_prop_to_note_mapping')
        
      mapping_dict_name = 'cassette_prop_to_note_mapping'
      for fb_prop_type, note_specs in mapping_dict.items():
          agr_note_type_name = note_specs[NOTE_TYPE_NAME]
          agr_slot_name = note_specs[NOTE_SLOT_NAME]
          agr_notes = self.convert_prop_to_note(entity, fb_prop_type, agr_note_type_name)
          agr_note_slot = getattr(entity.linkmldto, agr_slot_name)
          agr_note_slot.extend(agr_notes)

#### cassette_handler.py
    Add 1) 'secondary_identifiers' to fb_data_entity.linkmldto
    and 2) associations to :- 
        self.cassette_tool_associations or
        self.cassette_genomic_entity_associations
    depending on the type.

    1) map_secondary_ids(self, slot_name)
        def map_secondary_ids(self, slot_name):
            for fb_data_entity in self.fb_data_entities.values():
                secondary_id_dtos = []
                for secondary_id in fb_data_entity.alt_fb_ids:
                    secondary_id_dtos.append(secondary_id)
                sec_id_list = getattr(fb_data_entity.linkmldto, slot_name)
                sec_id_list.extend(secondary_id_dtos)
        return

    2) map_cassette_associations(self)
        def map_cassette_associations(self):
            for cassette_cassette_key, cassette_cassette_rels in self.cassette_cassette_rels.items():
                ... example of tool_association
                rel_dto = agr_datatypes.CassetteTransgenicToolAssociationDTO(
                    object_curie, subject_curie,
                    pub_curies, False, rel_type_name)
                first_feat_rel.linkmldto = rel_dto
                self.cassette_tool_associations.append(first_feat_rel)
