"""
The :class:`Parent` class defines a hierarchical relationship between objects. :class:`Location` objects,
:class:`Sequence` objects, and :class:`Parent` objects can have a :class:`Parent`. A :class:`Parent` object defines its
relationship with its child object through a variety of optional constructor parameters.
"""

from functools import singledispatch

from inscripta.biocantor.parent.parent import Parent, SequenceType  # noqa: F401


@singledispatch
def make_parent(obj) -> Parent:
    raise TypeError("{} not supported".format(type(obj)))


@make_parent.register(str)
def _(obj) -> Parent:
    return Parent(id=obj)


@make_parent.register(Parent)
def _(obj) -> Parent:
    return obj
