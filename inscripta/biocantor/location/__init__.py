from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import make_parent, Parent


@make_parent.register(Location)
def _(obj) -> Parent:
    return Parent(location=obj)


@make_parent.register(Strand)
def _(obj) -> Parent:
    return Parent(strand=obj)
