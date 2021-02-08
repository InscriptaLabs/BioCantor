"""
Enumeration utilities.
"""
from enum import Enum


class HasMemberMixin(Enum):
    """Adds a `has_value()` and `has_name()` convenience method to enumerations."""

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_

    @classmethod
    def has_name(cls, name):
        return name in cls.__members__
