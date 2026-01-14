# Export the json files

#### cassette_handler.py
    Default dump main json for features.
    Then dump the associations which are in there own dicts.

    def query_chado_and_export(self, session):
        """Elaborate on query_chado_and_export method for the CassetteHandler."""
        super().query_chado_and_export(session)
                          'cassette_str_association_ingest_set')
        self.generate_export_dict(self.cassette_genomic_entity_associations,
                                  'cassette_genomic_entity_association_ingest_set')
        self.generate_export_dict(self.cassette_tool_associations,
                                  'cassette_transgenic_tool_association_ingest_set')

        return

#### handler.py (super)
    Important bit is:-
       for attr in i.linkmldto.__dict__.keys()
    So anything tacked on linkmldtp will be outputted


    def generate_export_dict(self, input_list: list, output_set_name: str)
            for attr in i.linkmldto.__dict__.keys():
                # self.log.debug(f'Assess this attr: {attr}')
                if attr in i.linkmldto.internal_fields:
                    # self.log.debug(f'Skip this field: {attr}')
                    continue
                elif attr in i.linkmldto.required_fields:
                    # self.log.debug(f'Export required field: {attr}')
                    export_agr_dict[attr] = getattr(i.linkmldto, attr)
                elif getattr(i.linkmldto, attr) is not None and getattr(i.linkmldto, attr) != []:
                    # self.log.debug(f'Export optional non-empty field: {attr}')
                    export_agr_dict[attr] = getattr(i.linkmldto, attr)
                else:
                    # self.log.debug(f'Empty value for this attr: {attr}')
                    pass
            if self.incremental_update is False:
                self.export_data[output_set_name].append(export_agr_dict)
            elif i.is_new_addition is True or i.is_new_obsolete is True:
                self.export_data[output_set_name].append(export_agr_dict)
                incremental_count += 1
        