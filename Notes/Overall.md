## So in this example we are using cassettes.
## This is an overview
# AGR_data_retrieval_curation_casste.py
    def main():

        # Get the data and process it.
        cassette_handler = CassetteHandler(log, testing)
        if reference_session:
            export_chado_data

# utils.py:
    def export_chado_data(session: Session, log: Logger, object_to_execute: DataHandler, **kwargs):
        object_to_execute.query_chado_and_export(session)

# cassette_handler.py
    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the CassetteHandler."""
        super().query_chado_and_export(session)
        # self.generate_export_dict(self.cassette_component_free_text_associations,
        #                           'cassette_str_association_ingest_set')
        self.generate_export_dict(self.cassette_genomic_entity_associations,
                                  'cassette_genomic_entity_association_ingest_set')
        self.generate_export_dict(self.cassette_tool_associations,
                                  'cassette_transgenic_tool_association_ingest_set')
# handler.py (super)
    def query_chado_and_export(self, session):
        """Wrapper that runs all methods within an SQLAlchemy session."""
        self.log.info(f'RUN MAIN {type(self)} QUERY_CHADO_AND_EXPORT() METHOD.')
        self.get_general_data(session)
        self.get_datatype_data(session)
        self.synthesize_info()
        self.map_fb_data_to_alliance()

        self.flag_unexportable_entities(self.fb_data_entities.values(), self.primary_export_set)
        self.generate_export_dict(self.fb_data_entities.values(), self.primary_export_set)
        return