# Get the features data of interest.

### cassette_handler.py
    def get_datatype_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_relationships(session, 'subject')

### get entities:
        joins FeatureRelationship, CVterm, relationFeature
        BUt later we call get_entity_relationships ??
            results = session.query(chado_table).\
            select_from(chado_table).\
            join(Cvterm, (Cvterm.cvterm_id == chado_table.type_id)).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == chado_table.feature_id)).\
            join(feat_object, (feat_object.feature_id == FeatureRelationship.object_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
            for result in results:
               self.fb_data_entities[pkey_id] = self.fb_export_type(result)
               # creates for each entity:-
                 self.has_reg_region = []
                 self.tagged_with = []
                 self.carries_tool = []

### get_entityprops
    self.fb_data_entities[subject_id].props_by_type[prop.chado_obj.type.name] = [prop]

### get_entity_pubs
    self.fb_data_entities[entity_pkey_id].pub_associations.append(result)

### get_entity_synonyms
    self.fb_data_entities[entity_pkey_id].synonyms.append(result)

### get_entity_fb_xrefs
    self.fb_data_entities[entity_pkey_id].synonyms.append(result)

### get_entity_relationships
#### For cassettes 'subject' chosen
    for rel_result in rel_results:
        rel_id = getattr(rel_result, f'{chado_type}_relationship_id')
        rel_dict[rel_id] = fb_datatypes.FBRelationship(rel_result, f'{chado_type}_relationship')
    for rel_id, rel in rel_dict.items():
        self.fb_data_entities[entity_id].rels_by_id[rel_id] = rel
