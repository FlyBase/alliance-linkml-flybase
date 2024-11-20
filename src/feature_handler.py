"""Module:: feature_handler.

Synopsis:
    A generic data handler that runs basic processes common to the mapping of
    all FlyBase features in chado to Alliance LinkML objects. Built upon these
    core handlers are more specific handlers to specific features (see other
    handler files in this repo).

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Featureloc, Feature, FeatureDbxref, Library, LibraryFeature, LibraryFeatureprop
)
import agr_datatypes
from entity_handler import PrimaryEntityHandler


class FeatureHandler(PrimaryEntityHandler):
    """A generic, abstract handler for that gets data for FlyBase features and maps it to the Alliance LinkML model."""
    def __init__(self, log: Logger, testing: bool):
        """Create the FeatureHandler object."""
        super().__init__(log, testing)

    # Add methods to be run by get_general_data() below.
    # Placeholder.

    # Add methods to be run by get_datatype_data() below.
    def get_annotation_ids(self, session):
        """Get annotation IDs (current and non-current)."""
        self.log.info('Get annotation IDs (current and non-current).')
        filters = (
            Feature.uniquename.op('~')(self.regex[self.datatype]),
            Feature.is_analysis.is_(False),
            Db.name == 'FlyBase Annotation IDs'
        )
        if self.datatype in self.feature_subtypes.keys():
            self.log.info(f'Filter main table by these feature_subtypes: {self.feature_subtypes[self.datatype]}')
            filters += (Cvterm.name.in_((self.feature_subtypes[self.datatype])), )
        if self.testing:
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(FeatureDbxref).\
            select_from(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(FeatureDbxref, (FeatureDbxref.feature_id == Feature.feature_id)).\
            join(Dbxref, (Dbxref.dbxref_id == FeatureDbxref.dbxref_id)).\
            join(Db, (Db.db_id == Dbxref.db_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].fb_anno_dbxrefs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} annotation IDs for {self.datatype} entities.')
        self.log.info(f'Ignored {pass_counter} annotation IDs for {self.datatype} entities.')
        return

    def get_chr_featurelocs(self, session):
        """Get chromosomal featureloc data."""
        self.log.info('Get chromosomal featureloc data.')
        filters = (
            Feature.is_analysis.is_(False),
            Featureloc.srcfeature_id.in_((self.chr_dict.keys()))
        )
        if self.datatype in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.datatype]}')
            filters += (Feature.uniquename.op('~')(self.regex[self.datatype]), )
        if self.datatype in self.feature_subtypes.keys():
            self.log.info(f'Filter main table by these feature_subtypes: {self.feature_subtypes[self.datatype]}')
            filters += (Cvterm.name.in_((self.feature_subtypes[self.datatype])), )
        if self.testing:
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Featureloc).\
            select_from(Feature).\
            join(Cvterm, (Cvterm.cvterm_id == Feature.type_id)).\
            join(Featureloc, (Featureloc.feature_id == Feature.feature_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        pass_counter = 0
        for result in results:
            entity_pkey_id = result.feature_id
            try:
                self.fb_data_entities[entity_pkey_id].chr_flocs.append(result)
                counter += 1
            except KeyError:
                pass_counter += 1
        self.log.info(f'Found {counter} chromosomal featurelocs for {self.datatype} entities.')
        self.log.info(f'Ignored {pass_counter} chromosomal featurelocs for {self.datatype} entities.')
        return

    def get_direct_reagent_collections(self, session):
        """Find reagent collections directly related to primary features."""
        self.log.info('Find reagent collections directly related to primary features.')
        libtype = aliased(Cvterm, name='libtype')
        libfeattype = aliased(Cvterm, name='libfeattype')
        filters = (
            Feature.uniquename.op('~')(self.regex[self.datatype]),
            Library.is_obsolete.is_(False),
            Library.uniquename.op('~')(self.regex['library']),
            libtype.name == 'reagent collection',
            libfeattype.name == 'member_of_reagent_collection'
        )
        if self.testing:
            self.log.info(f'TESTING: limit to these entities: {self.test_set}')
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        collections = session.query(Feature, Library).\
            select_from(Feature).\
            join(LibraryFeature, (LibraryFeature.feature_id == Feature.feature_id)).\
            join(Library, (Library.library_id == LibraryFeature.library_id)).\
            join(LibraryFeatureprop, (LibraryFeatureprop.library_feature_id == LibraryFeature.library_feature_id)).\
            join(libtype, (libtype.cvterm_id == Library.type_id)).\
            join(libfeattype, (libfeattype.cvterm_id == LibraryFeatureprop.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in collections:
            self.fb_data_entities[result.Feature.feature_id].reagent_colls.append(result.Library)
            counter += 1
        self.log.info(f'Found {counter} direct allele-collection associations.')
        return

    # Add methods to be run by synthesize_info() below.
    def synthesize_anno_ids(self):
        """Synthesize annotation IDs."""
        self.log.info('Synthesize annotation IDs.')
        for fb_data_entity in self.fb_data_entities.values():
            current_anno_ids = []
            alt_anno_ids = []
            for xref in fb_data_entity.fb_anno_dbxrefs:
                if xref.is_current is True:
                    current_anno_ids.append(xref.dbxref.accession)
                else:
                    alt_anno_ids.append(xref.dbxref.accession)
            # Get the one current annotation ID.
            if len(current_anno_ids) == 1:
                fb_data_entity.curr_anno_id = current_anno_ids[0]
            elif len(current_anno_ids) > 1:
                self.log.warning(f'{fb_data_entity} has {len(current_anno_ids)} current annotations IDs.')
                fb_data_entity.alt_anno_ids.extend(current_anno_ids)
            # Record old annotation IDs.
            fb_data_entity.alt_anno_ids.extend(alt_anno_ids)
        return

    # Add methods to be run by map_fb_data_to_alliance() below.
    def map_anno_ids_to_secondary_ids(self, slot_name):
        """Return a list of Alliance SecondaryIdSlotAnnotationDTOs for annotation IDs."""
        self.log.info('Map annotation IDs to Alliance object.')
        for fb_data_entity in self.fb_data_entities.values():
            anno_ids = []
            if fb_data_entity.curr_anno_id:
                anno_ids.append(fb_data_entity.curr_anno_id)
            if fb_data_entity.alt_anno_ids:
                anno_ids.extend(fb_data_entity.alt_anno_ids)
            anno_secondary_id_dtos = []
            for anno_id in anno_ids:
                sec_dto = agr_datatypes.SecondaryIdSlotAnnotationDTO(f'FB:{anno_id}', []).dict_export()
                anno_secondary_id_dtos.append(sec_dto)
            curr_sec_id_dtos = getattr(fb_data_entity.linkmldto, slot_name)
            curr_sec_id_dtos.extend(anno_secondary_id_dtos)
        return anno_secondary_id_dtos
