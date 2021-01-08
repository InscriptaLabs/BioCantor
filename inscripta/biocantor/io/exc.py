"""
I/O exceptions.
"""
from inscripta.biocantor.exc import BioCantorException


class BioCantorIOException(BioCantorException):
    pass


class InvalidInputError(BioCantorIOException):
    pass


class StrandViolationWarning(UserWarning):
    pass
