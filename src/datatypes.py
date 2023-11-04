"""Module:: datatypes.

Synopsis:
    Datatype objects representing FlyBaseData retrieval utils from FlyBase chado for export to the Alliance in
    LinkML-compliant JSON files.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

"""

import json
from logging import Logger
from sqlalchemy.orm import Session
from harvdev_utils.production import (
    Cv, Cvterm, Db, Dbxref, Organism, OrganismDbxref, Pub, PubDbxref
)


# FlyBase Classes
class FlyBaseDataEntity(object):
    """A generic FlyBase data entity with all it related data."""
    def __init__(self, chado_obj):
        """Create the generic FlyBase data entity object from the main db entry/entries."""
        # Associated data.
        self.chado_obj = None     # The primary SQLAlchemy chado object (or tuple of objects).
        self.db_uniq_id = None    # The chado table primary key (or concatenation of primary keys).
        self.fb_uname = None      # The FlyBase uniquename, if applicable.
        self.uniq_key = None      # A string derived from the uniquely defining properties of the entity.
        self.linkmldto = None     # The Alliance LinkML data transfer object model to which data will be mapped.
        # Export status and messages.
        self.for_export = True        # Change to False if object should be excluded from export.
        self.internal_reasons = []    # Reasons for marking an object as internal in the export file.
        self.export_warnings = []     # Reasons for suppressing an object from the export file.


# Alliance Classes for FlyBase data.
class AuditedObjectDTO(object):
    """Base Alliance class."""
    def __init__(self):
        """Create base AuditedObjectDTO for FlyBase objects."""
        self.internal = False
        self.obsolete = False
        self.date_created = None
        self.date_updated = None
        self.created_by_curie = 'FB:FB_curator'
        self.updated_by_curie = 'FB:FB_curator'
        self.test_list = []


class DataProviderDTO(AuditedObjectDTO):
    """DataProvider class."""
    def __init__(self):
        """Create simple DataProviderDTO for FlyBase objects."""
        self.source_organization_abbreviation = 'FB'
        self.cross_reference_dto = {'prefix': 'FB', 'internal': False}

