from enum import Enum
from functools import total_ordering
from inscripta.biocantor.exc import UnsupportedOperationException, InvalidStrandException


@total_ordering
class Strand(Enum):
    PLUS = 1
    MINUS = -1
    UNSTRANDED = 0

    def __str__(self):
        return str(self.to_symbol())

    @staticmethod
    def from_symbol(value: str):
        """Converts string representation of a strand to a Strand"""
        if value == "+":
            return Strand.PLUS
        if value == "-":
            return Strand.MINUS
        raise ValueError("{} is not a valid string representation of a strand".format(value))

    def to_symbol(self) -> str:
        if self == Strand.PLUS:
            return "+"
        if self == Strand.MINUS:
            return "-"
        return "."

    @staticmethod
    def from_int(value: int):
        """Converts integer representation of a strand to a Strand"""
        return Strand(value)  # Raises ValueError for invalid int

    @staticmethod
    def _order():
        return {Strand.PLUS: 1, Strand.MINUS: 2, Strand.UNSTRANDED: 3}

    def __lt__(self, other):
        if not type(other) is Strand:
            raise ValueError("Cannot compare {} to {}".format(type(self).__name__, type(other).__name__))
        order = Strand._order()
        return order[self] < order[other]

    def reverse(self):
        """Returns the opposite of this Strand"""
        if self == Strand.PLUS:
            return Strand.MINUS
        if self == Strand.MINUS:
            return Strand.PLUS
        if self == Strand.UNSTRANDED:
            return Strand.UNSTRANDED
        raise UnsupportedOperationException("Not implemented for {}".format(self))

    def relative_to(self, other: "Strand") -> "Strand":
        """Returns the orientation of this strand relative to the given other strand.
        Note: this operator is commutative.
        """
        if self == Strand.UNSTRANDED or other == Strand.UNSTRANDED:
            return Strand.UNSTRANDED
        if self == other:
            return Strand.PLUS
        return Strand.MINUS

    def assert_directional(self):
        """Raises InvalidStrandException if this Strand does not have a defined direction (plus or minus)"""
        if self not in [Strand.PLUS, Strand.MINUS]:
            raise InvalidStrandException("Strand {} does not have a defined direction".format(self))
