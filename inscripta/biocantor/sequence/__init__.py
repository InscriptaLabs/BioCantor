from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.parent import make_parent, Parent


@make_parent.register(Sequence)
def _(obj) -> Parent:
    return Parent(sequence=obj)
