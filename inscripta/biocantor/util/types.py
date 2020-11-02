from typing import Union

from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.location.location import Location


ParentInputType = Union[Sequence, str, Location, Strand, Parent]
