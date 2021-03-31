class BioCantorException(Exception):
    """
    Base exception class for BioCantor.
    """

    pass


class InvalidStrandException(BioCantorException):
    """
    Raised when an operation is performed on an invalid strand -- usually this is when an operation
    is strand-specific but the provided Strand is Unstranded.
    """

    pass


class InvalidPositionException(BioCantorException):
    """
    Raised when a position is outside of a valid range for the operation being performed.
    """

    pass


class UnsupportedOperationException(BioCantorException):
    """
    Raised when an object is being used for a comparison like intersection or union in a way that is unsupported.
    """

    pass


class LocationException(BioCantorException):
    """
    Raised when a Location constructor is given invalid inputs, such as an unequal number of start/end positions.
    """

    pass


class LocationOverlapException(LocationException):
    """
    Raised when a Location operation that requires overlap is being performed on Locations that do not overlap.
    """

    pass


class EmptyLocationException(LocationException):
    """
    Raised when a Location operation that requires a non-empty location is being performed on an EmptyLocation.
    """

    pass


class AlphabetError(BioCantorException):
    """
    Raised when an operation on an Alphabet is unsupported for the provided Alphabet.
    """

    pass


class NoSuchAncestorException(BioCantorException):
    """
    Raised when the provided Location or Parent does not have an Ancestor of the specified type.
    """

    pass


class ValidationException(BioCantorException):
    """
    Raised when object constructors are given invalid inputs that are not LocationExceptions.
    """

    pass


class InvalidCDSIntervalError(ValidationException):
    """
    Raised when CDS-specific constructor items are invalid.
    """

    pass


class EmptySequenceFastaError(BioCantorException):
    """
    Raised when FASTA export is attempted on an empty Sequence object.
    """

    pass


class NoncodingTranscriptError(BioCantorException):
    """
    Raised when operations that require a coding transcript are performed on a non-coding transcript.
    """

    pass


class InvalidAnnotationError(BioCantorException):
    """
    Raised when Collection objects have invalid arguments.
    """

    pass


class InvalidQueryError(BioCantorException):
    """
    Raised when range queries on Collection objects are invalid.
    """

    pass


class ParentException(BioCantorException):
    """
    Generic exception involving Parent objects.
    """

    pass


class MismatchedParentException(ParentException):
    """
    Raised when operations that require comparing Parent objects are performed with incompatible Parent objects.
    """

    pass


class NullParentException(ParentException):
    """
    Raised when operations that require a Parent are performed on a parentless Location.
    """

    pass


class NullSequenceException(ParentException):
    """
    Raised when a Location operation is performed that requires a Sequence and there is no Sequence.
    """

    pass


class MismatchedFrameException(BioCantorException):
    """
    Raised when the frames list of a CDS does not match the location object.
    """

    pass
