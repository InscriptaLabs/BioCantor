from functools import singledispatch

from inscripta.biocantor.parent.parent import Parent


@singledispatch
def make_parent(obj) -> Parent:
    raise TypeError("{} not supported".format(type(obj)))


@make_parent.register(str)
def _(obj) -> Parent:
    return Parent(id=obj)


@make_parent.register(Parent)
def _(obj) -> Parent:
    return obj
