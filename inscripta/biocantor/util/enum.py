"""
Enumeration utilities.
"""
from enum import Enum


class HasMemberMixin(Enum):
    """Adds a `has_value()` convenience method to enumerations."""

    @classmethod
    def has_value(cls, value):
        return value in cls._value2member_map_
