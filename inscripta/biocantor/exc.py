class InvalidStrandException(Exception):
    pass


class InvalidPositionException(Exception):
    pass


class UnsupportedOperationException(TypeError):
    pass


class LocationException(ValueError):
    pass


class LocationOverlapException(LocationException):
    pass


class EmptyLocationException(LocationException):
    pass


class AlphabetError(ValueError):
    pass


class NoSuchAncestorException(Exception):
    pass


class ValidationError(Exception):
    pass


class InvalidCDSIntervalError(Exception):
    pass


class EmptySequenceFastaError(Exception):
    pass
