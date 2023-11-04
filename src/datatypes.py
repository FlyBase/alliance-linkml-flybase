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
        self.chado_obj = None    # The primary SQLAlchemy chado object (or tuple of objects).

    db_uniq_id = None    # The chado table primary key (or concatenation of primary keys).
    fb_uname = None      # The FlyBase uniquename, if applicable.
    uniq_key = None      # A string derived from the uniquely defining properties of the entity.
    linkmldto = None     # The Alliance LinkML data transfer object model to which data will be mapped.


# Alliance Classes for FlyBase data.
class AuditedObjectDTO(object):
    """Base Alliance class."""
    def __init__(self):
        """Create base AuditedObjectDTO for FlyBase objects."""

    internal = False
    obsolete = False
    date_created = None
    date_updated = None
    created_by_curie = 'FB:FB_curator'
    updated_by_curie = 'FB:FB_curator'
    test_list = []


class DataProviderDTO(AuditedObjectDTO):
    """DataProvider class."""
    def __init__(self):
        """Create simple DataProviderDTO for FlyBase objects."""

    source_organization_abbreviation = 'FB'
    cross_reference_dto = {'prefix': 'FB', 'internal': False}
