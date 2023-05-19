from typing import Union

from biocantor.location.strand import Strand
from biocantor.parent import Parent
from biocantor.sequence import Sequence
from biocantor.location.location import Location


ParentInputType = Union[Sequence, str, Location, Strand, Parent]
