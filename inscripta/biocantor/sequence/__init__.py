"""
The :class:`Sequence` class defines a sequence with an :class:`Alphabet`. :class:`Sequence` integrates into the
:class:`Location` and :class:`Parent` paradigm, so that a :class:`Sequence` can be situated with respect to a parent
and can include child features.
"""

from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.parent import make_parent, Parent


@make_parent.register(Sequence)
def _(obj) -> Parent:
    return Parent(sequence=obj)
