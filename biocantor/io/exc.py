"""
I/O exceptions.
"""
from biocantor.exc import BioCantorException


class BioCantorIOException(BioCantorException):
    pass


class InvalidInputError(BioCantorIOException):
    pass


class StrandViolationWarning(UserWarning):
    pass


class InvalidIntervalWarning(UserWarning):
    pass


class InvalidCDSIntervalWarning(InvalidIntervalWarning):
    pass


class DuplicateFeatureWarning(UserWarning):
    pass


class DuplicateTranscriptWarning(UserWarning):
    pass


class DuplicateSequenceException(BioCantorIOException):
    pass
