# Synthesize data

## cassette_handler.py
    def synthesize_info(self):
        """Extend the method for the ConstructHandler."""
        super().synthesize_info()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_cassette_associations()

### synthesize_synonyms
        for fb_data_entity in self.fb_data_entities.values():
            # For each entity, gather synonym_id-keyed dict of synonym info.
            for feat_syno in fb_data_entity.synonyms:
                try:
                    fb_data_entity.synonym_dict[feat_syno.synonym_id]['is_current'].append(feat_syno.is_current)
                    fb_data_entity.synonym_dict[feat_syno.synonym_id]['is_internal'].append(feat_syno.is_internal)
                    fb_data_entity.synonym_dict[feat_syno.synonym_id]['pub_ids'].append(feat_syno.pub_id)

### synthesize secondary  ids
        for fb_data_entity in self.fb_data_entities.values():
            secondary_ids = []
            for xref in fb_data_entity.fb_sec_dbxrefs:
                secondary_ids.append(f'FB:{xref.dbxref.accession}')
            fb_data_entity.alt_fb_ids = list(set(fb_data_entity.alt_fb_ids).union(set(secondary_ids)))
        
### synthesize_cassette_associations

          for cassette in self.fb_data_entities.values():
            relevant_cassette_rels = cassette.recall_relationships(
                self.log, entity_role='subject',
                rel_types=['has_reg_region', 'tagged_with', 'carries_tool'])
            for cassette_rel in relevant_cassette_rels:
                    cassette_cassette_key = (cassette_rel.chado_obj.object_id, cassette_rel.chado_obj.subject_id)
                    self.cassette_cassette_rels[cassette_cassette_key].append(cassette_rel)
