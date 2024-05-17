"""Module:: feature_handler.

Synopsis:
    Core data handler objects that run basic processes common to the mapping of
    all FlyBase features in chado to Alliance LinkML objects. Built upon these
    core handlers are more specific handlers to specific features (see other
    handler files in this repo).

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

from logging import Logger
from sqlalchemy.orm import aliased
from harvdev_utils.production import (
    Cvterm, Db, Dbxref, Featureloc, Feature, FeatureDbxref, Featureprop,
    FeatureRelationship
)
import agr_datatypes
from handler import PrimaryEntityHandler


class FeatureHandler(PrimaryEntityHandler):
    """A generic, abstract handler for that gets data for FlyBase features and maps it to the Alliance LinkML model."""
    def __init__(self, log: Logger, fb_data_type: str, testing: bool):
        """Create the FeatureHandler object."""
        super().__init__(log, fb_data_type, testing)

    def get_general_data(self, session):
        """Extend the method for the FeatureHandler."""
        super().get_general_data(session)
        return

    # Elaborate on get_datatype_data() for the FeatureHandler; sub-methods might only be used in some more specific Datahandlers.
    def get_annotation_ids(self, session):
        """Get annotation IDs (current and non-current)."""
        self.log.info('Get annotation IDs (current and non-current).')
        filters = (
            Feature.uniquename.op('~')(self.regex[self.fb_data_type]),
            Feature.is_analysis.is_(False),
            Db.name == 'FlyBase Annotation IDs'
        )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (Cvterm.name.in_((self.subtypes[self.fb_data_type])), )
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
        self.log.info(f'Found {counter} annotation IDs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} annotation IDs for {self.fb_data_type} entities.')
        return

    def get_chr_featurelocs(self, session):
        """Get chromosomal featureloc data."""
        self.log.info('Get chromosomal featureloc data.')
        filters = (
            Feature.is_analysis.is_(False),
            Featureloc.srcfeature_id.in_((self.chr_dict.keys()))
        )
        if self.fb_data_type in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.fb_data_type]}')
            filters += (Feature.uniquename.op('~')(self.regex[self.fb_data_type]), )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (Cvterm.name.in_((self.subtypes[self.fb_data_type])), )
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
        self.log.info(f'Found {counter} chromosomal featurelocs for {self.fb_data_type} entities.')
        self.log.info(f'Ignored {pass_counter} chromosomal featurelocs for {self.fb_data_type} entities.')
        return

    def get_featureprops_by_type(self, session, fprop_type):
        """Return a list of featureprops of a given type."""
        self.log.info(f'Get featureprops of type {fprop_type} for {self.fb_data_type} entities.')
        feat_type = aliased(Cvterm, name='feat_type')
        prop_type = aliased(Cvterm, name='prop_type')
        filters = (
            Feature.is_analysis.is_(False),
            prop_type.name == fprop_type
        )
        if self.fb_data_type in self.regex.keys():
            self.log.info(f'Use this regex: {self.regex[self.fb_data_type]}')
            filters += (Feature.uniquename.op('~')(self.regex[self.fb_data_type]), )
        if self.fb_data_type in self.subtypes.keys():
            self.log.info(f'Filter main table by these subtypes: {self.subtypes[self.fb_data_type]}')
            filters += (feat_type.name.in_((self.subtypes[self.fb_data_type])), )
        if self.testing:
            filters += (Feature.uniquename.in_((self.test_set.keys())), )
        results = session.query(Featureprop).\
            select_from(Feature).\
            join(feat_type, (feat_type.cvterm_id == Feature.type_id)).\
            join(Featureprop, (Featureprop.feature_id == Feature.feature_id)).\
            join(prop_type, (prop_type.cvterm_id == Featureprop.type_id)).\
            filter(*filters).\
            distinct()
        return results

    def get_entity_sbj_feat_rel_by_type(self, session, slot_name, **kwargs):
        """Get a list of FeatureRelationships for primary feature entities (subject).

        Get FeatureRelationships where the primary feature entity is the subject,
        restricted by relationship type and related object.

        Args:
            session (Session): SQLAlchemy session for db queries.
            slot_name (str): The name of the FB data entity attribute name to which results are appended.

        Keyword Args:
            rel_type (str): The CV term name for the feature_relationship of interest. If none given, any rel_type allowed.
            obj_type (list): A list of CV terms for the object feature types. If none given, any object feature type allowed.
            obj_regex (str): The regex for the object feature uniquename. If none given, any object uniquename allowed.

        """
        self.log.info(f'Add feature_relationships to "{slot_name}" with these criteria: {kwargs}')
        object = aliased(Feature, name='object')
        rel_type = aliased(Cvterm, name='rel_type')
        obj_type = aliased(Cvterm, name='obj_type')
        filters = ()
        try:
            filters += (rel_type.name == kwargs['rel_type'], )
        except KeyError:
            pass
        try:
            filters += (obj_type.name == kwargs['obj_type'], )
        except KeyError:
            pass
        try:
            filters += (object.uniquename.op('~')(kwargs['obj_regex']), )
        except KeyError:
            pass
        results = session.query(FeatureRelationship).\
            select_from(object).\
            join(FeatureRelationship, (FeatureRelationship.object_id == object.feature_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            join(obj_type, (obj_type.cvterm_id == object.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.fb_data_entities[result.subject_id].__dict__[slot_name].append(result)
                counter += 1
            except KeyError:
                pass
        self.log.info(f'Added {counter} feature_relationship results to "{slot_name}" list.')
        return

    def get_entity_obj_feat_rel_by_type(self, session, slot_name, **kwargs):
        """Get a list of FeatureRelationships for primary feature entities (object).

        Get FeatureRelationships where the primary feature entity is the object,
        restricted by relationship type and related subject.

        Args:
            session (Session): SQLAlchemy session for db queries.
            slot_name (str): The name of the FB data entity attribute name to which results are appended.

        Keyword Args:
            rel_type (str): The CV term name for the feature_relationship of interest. If none given, any rel_type allowed.
            sbj_type (str): The CV term for the subject feature types. If none given, any subject feature type allowed.
            sbj_regex (str): The regex for the subject feature uniquename. If none given, any subject uniquename allowed.

        """
        self.log.info(f'Add feature_relationships to "{slot_name}" with these criteria: {kwargs}')
        subject = aliased(Feature, name='subject')
        # object = aliased(Feature, name='object')
        rel_type = aliased(Cvterm, name='rel_type')
        sbj_type = aliased(Cvterm, name='sbj_type')
        # filters = (
        #     object.feature_id.in_(self.fb_data_entities.keys()),
        # )
        filters = ()
        try:
            filters += (rel_type.name == kwargs['rel_type'], )
        except KeyError:
            pass
        try:
            filters += (sbj_type.name == kwargs['sbj_type'], )
        except KeyError:
            pass
        try:
            filters += (subject.uniquename.op('~')(kwargs['sbj_regex']), )
        except KeyError:
            pass
        results = session.query(FeatureRelationship).\
            select_from(subject).\
            join(sbj_type, (sbj_type.cvterm_id == subject.type_id)).\
            join(FeatureRelationship, (FeatureRelationship.subject_id == subject.feature_id)).\
            join(rel_type, (rel_type.cvterm_id == FeatureRelationship.type_id)).\
            filter(*filters).\
            distinct()
        counter = 0
        for result in results:
            try:
                self.fb_data_entities[result.object_id].__dict__[slot_name].append(result)
                counter += 1
            except KeyError:
                pass
        self.log.info(f'Added {counter} feature_relationship results to "{slot_name}" list.')
        return

    def get_datatype_data(self, session):
        """Extend the method for the FeatureHandler."""
        super().get_datatype_data(session)
        return

    # Elaborate on synthesize_info() for the FeatureHandler; sub-methods might only be used in some more specific DataHandlers.
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

    def synthesize_info(self):
        """Extend the method for the FeatureHandler."""
        super().synthesize_info()
        return

    # Elaborate on map_fb_data_to_alliance() for the FeatureHandler; sub-methods might only be used in some more specific DataHandlers.
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

    def map_fb_data_to_alliance(self):
        """Extend the method for the FeatureHandler."""
        super().map_fb_data_to_alliance()
        return
