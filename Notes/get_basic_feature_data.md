## What data is obtained from chado and how is it stored.
# methods are trimmed to shoe bare minimum that is needed to explain


# cassette_handler.py
    def get_general_data(self, session):
        """Extend the method for the CassetteHandler."""
        super().get_general_data(session)

        self.build_bibliography(session)
            # Gets list of FBrf snf pmcs
            self.fbrf_bibliography[uniquename] = pub
            self.bibliography[pub_id] = f'PMID:{xref.Dbxref.accession}'
            self.bibliography[pub_id] = 'FB:uniquename'

        self.build_organism_lookup(session)
             # Gets organism info
             self.organism_lookup[organism_id]['taxon_curie'] = NCBI id
                                              ['is_drosophilid'] = True/ None
                                              ['official_db'] = value ?

        self.build_feature_lookup(session, feature_types=['cassette', 'construct', 'allele', 'tool', 'gene', 'seqfeat'])

# handler.py
     def build_feature_lookup(self, session, **kwargs):
          # Get feature objects
          ... Get all features with specific regex
          feat_results = session.query(..
                join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
                filter(*feat_filters).\
                distinct()
          for result in feat_results:
                feat_dict = {
                    'feature_id': result[FEATURE_ID],
                    'uniquename': result[UNIQUENAME],
                    'curie': f'FB:{result[UNIQUENAME]}',    # Replaced by MOD curies as applicable down below.
                    'is_obsolete': result[OBSOLETE],
                    'type': result[TYPE],
                    'organism_id': result[ORG_ID],
                    'name': result[NAME],
                    'symbol': result[NAME],
                    'exported': self.feat_type_export[feat_type],
                }
                self.feature_lookup[result[FEATURE_ID]] = feat_dict
                feat_counter += 1
           for result in syno_results:
               self.feature_lookup[result[FEATURE_ID]]['symbol'] = sub_sup_sgml_to_html(result[SYMBOL])
