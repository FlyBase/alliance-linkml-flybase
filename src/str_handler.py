"""Module:: str_handler.

Synopsis:
    A data handler that exports FlyBase sequence targeting reagent (STR) data to
    Alliance SequenceTargetingReagent LinkML objects.

    The STR data class is the subset of FlyBase sequence features (FBsf) whose
    feature.type is "RNAi_reagent" or "sgRNA" (FTA-224). That subset is declared in
    DataHandler.feature_subtypes['str'], so the generic PrimaryEntityHandler driver
    query picks it up without any bespoke SQL here.

Author(s):
    Ian Longden ilongden@morgan.harvard.edu

"""

from logging import Logger
import agr_datatypes
import fb_datatypes
from feature_handler import FeatureHandler


class SequenceTargetingReagentHandler(FeatureHandler):
    """This object gets, synthesizes and filters STR data for export."""
    def __init__(self, log: Logger, testing: bool):
        """Create the SequenceTargetingReagentHandler object."""
        super().__init__(log, testing)
        self.datatype = 'str'
        self.fb_export_type = fb_datatypes.FBStr
        self.agr_export_type = agr_datatypes.SequenceTargetingReagentDTO
        self.primary_export_set = 'str_ingest_set'

    test_set = {
        'FBsf0000443916': 'dsRNA-HMJ22303',     # RNAi_reagent, TRiP-3 collection (see allele_handlers FBti0164639).
        'FBsf0000078690': 'dsRNA-GD9857',       # RNAi_reagent, "encodes_tool" component of cassette FBal0198528.
        'FBsf0000000002': 'DRSC00004',          # RNAi_reagent with a secondary FB ID.
        'FBsf0000003641': 'DRSC03748',          # RNAi_reagent with two non-current symbols (HDC03394, HDC03690) and a 2o FB ID.
        'FBsf0000437305': 'dsRNA-HMS01879',     # The only obsolete STR in fb_2026_02: must export as internal.
        'FBsf0000910691': 'sgRNA-GS00469.1',    # sgRNA, "encodes_tool" component of cassette FBal0365029.
        'FBsf0000910692': 'sgRNA-GS00469.2',    # sgRNA, "encodes_tool" component of cassette FBal0365030.
        'FBsf0000910567': 'sgRNA-GS00055.1',    # sgRNA, one of two components of cassette FBal0365146.
        # FTA-240: note-bearing entities.
        'FBsf0000032837': 'DRSC37925',          # "comment": reagent no longer exists at the DRSC.
        'FBsf0000407135': 'BKN28870',           # "comment" carrying both an SGML entity (&ggr;) and @symbol@ markup.
        'FBsf0000077885': 'dsRNA-GD8893',       # "internalnotes": must export as an internal_note.
    }

    # FTA-240: FlyBase props exported as Alliance notes. Same two chado prop types, and the same
    # Alliance note types, as construct_prop_to_note_mapping - "internalnotes" is in
    # convert_prop_to_note()'s internal_note_types, so those notes are flagged internal for us.
    str_prop_to_note_mapping = {
        'comment': ('comment', 'note_dtos'),               # SF6
        'internalnotes': ('internal_note', 'note_dtos'),   # SF8 (marked internal automatically)
    }

    # Add methods to be run by get_general_data() below.
    def get_general_data(self, session):
        """Extend the method for the SequenceTargetingReagentHandler."""
        super().get_general_data(session)
        self.build_bibliography(session)
        self.build_organism_lookup(session)
        self.build_cvterm_lookup(session)
        # NB - no build_feature_lookup(): the minimal STR submission makes no
        # associations to other features, so nothing needs the lookup.
        return

    # Add methods to be run by get_datatype_data() below.
    def get_datatype_data(self, session):
        """Extend the method for the SequenceTargetingReagentHandler."""
        super().get_datatype_data(session)
        self.get_entities(session)
        self.get_entityprops(session)    # FTA-240: "comment" and "internalnotes" props become notes.
        self.get_entity_pubs(session)
        self.get_entity_synonyms(session)
        self.get_entity_fb_xrefs(session)
        self.get_entity_xrefs(session)
        self.get_entity_timestamps(session)
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_info(self):
        """Extend the method for the SequenceTargetingReagentHandler."""
        super().synthesize_info()
        self.synthesize_ncbi_taxon_id()
        self.synthesize_synonyms()
        self.synthesize_secondary_ids()
        self.synthesize_pubs()
        self.flag_new_additions_and_obsoletes()
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_fb_data_to_alliance(self):
        """Extend the method for the SequenceTargetingReagentHandler."""
        super().map_fb_data_to_alliance()
        self.map_str_basic()
        self.map_str_names()
        self.map_timestamps()
        self.map_pubs()
        self.map_data_provider_dto()
        self.map_xrefs()
        self.map_secondary_ids('secondary_identifiers')
        self.map_entity_props_to_notes('str_prop_to_note_mapping')
        # Cascade chado-obsolete -> internal=True (matches every other handler).
        self.flag_internal_fb_entities('fb_data_entities')
        return

    def map_str_basic(self):
        """Map basic FlyBase STR data to the Alliance LinkML object."""
        self.log.info('Map basic STR info to Alliance object.')
        for str_entity in self.fb_data_entities.values():
            agr_str = self.agr_export_type()
            agr_str.obsolete = str_entity.chado_obj.is_obsolete
            agr_str.primary_external_id = f'FB:{str_entity.uniquename}'
            agr_str.taxon_curie = str_entity.ncbi_taxon_id
            str_entity.linkmldto = agr_str
        return

    def map_str_names(self):
        """Map the STR symbol and synonyms to the Alliance LinkML object (FTA-225, FTA-226).

        SequenceTargetingReagentDTO holds "name" as a plain string and "synonyms" as a plain
        list of strings, so PrimaryEntityHandler.map_synonyms() does not apply: it only fills
        "*_symbol_dto"/"*_synonym_dtos" slots with NameSlotAnnotationDTOs. The consequence is
        that neither pub attribution nor is_current status can be carried for a synonym; that
        loss is inherent to the LinkML slot, not a shortcut taken here.

        The format_text (plain) form of each name is used rather than the display_text
        (SGML-to-HTML) form, since the target slots are bare strings.
        """
        self.log.info('Map STR names and synonyms to Alliance object.')
        name_counter = 0
        fallback_counter = 0
        synonym_counter = 0
        symbol_type_names = ['nomenclature_symbol', 'systematic_name']
        for str_entity in self.fb_data_entities.values():
            if str_entity.linkmldto is None:
                continue
            # 1. The one current symbol becomes "name". Every STR in fb_2026_02 has exactly one,
            #    but prefer the symbol matching the feature's own name if the data gives several.
            current_symbols = [i['format_text'] for i in str_entity.synonym_dict.values()
                               if i['is_current'] is True and i['name_type_name'] in symbol_type_names]
            if not current_symbols:
                self.log.warning(f'No current symbol found for {str_entity}: using feature.name instead.')
                str_entity.linkmldto.name = str_entity.name
                fallback_counter += 1
            else:
                chosen_symbol = current_symbols[0]
                if str_entity.name in current_symbols:
                    chosen_symbol = str_entity.name
                if len(current_symbols) > 1:
                    other_symbols = ', '.join([i for i in current_symbols if i != chosen_symbol])
                    self.log.warning(f'Found many current symbols for {str_entity}: exporting "{chosen_symbol}" '
                                     f'and keeping these as synonyms: {other_symbols}')
                str_entity.linkmldto.name = chosen_symbol
                name_counter += 1
            # 2. Every other name becomes an (unattributed) entry in "synonyms".
            synonyms = {i['format_text'] for i in str_entity.synonym_dict.values()}
            synonyms.discard(str_entity.linkmldto.name)
            str_entity.linkmldto.synonyms = sorted(synonyms)
            synonym_counter += len(synonyms)
        self.log.info(f'Mapped {name_counter} STR names from a current symbol, and {fallback_counter} from feature.name.')
        self.log.info(f'Mapped {synonym_counter} STR synonyms.')
        return

    def map_secondary_ids(self, slot_name):
        """Map secondary FB IDs to the Alliance LinkML object.

        Overrides PrimaryEntityHandler.map_secondary_ids(), which wraps each ID in a
        SecondaryIdSlotAnnotationDTO. The "secondary_identifiers" slot on
        SequenceTargetingReagentDTO is a plain list of curies, so the raw "FB:FBsf..."
        strings are used instead (same approach as CassetteHandler).
        """
        self.log.info('Map secondary IDs to Alliance object.')
        counter = 0
        for fb_data_entity in self.fb_data_entities.values():
            if fb_data_entity.linkmldto is None:
                continue
            sec_id_list = getattr(fb_data_entity.linkmldto, slot_name)
            sec_id_list.extend(fb_data_entity.alt_fb_ids)
            counter += len(fb_data_entity.alt_fb_ids)
        self.log.info(f'Mapped {counter} secondary FB IDs for {self.datatype} entities.')
        return
