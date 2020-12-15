"""
:class:`Location` objects represent features defined with respect to a coordinate system. The :class:`Location` API
includes rich feature arithmetic methods. Additionally, a :class:`Location` object can define its relationship with
a :class:`Parent` object, situating it within a potentially arbitrary hierarchy of coordinate systems; the
:class:`Location` API provides rich coordinate and location conversion methods.
"""

from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import make_parent, Parent


@make_parent.register(Location)
def _(obj) -> Parent:
    return Parent(location=obj)


@make_parent.register(Strand)
def _(obj) -> Parent:
    return Parent(strand=obj)
