from biocantor.io.exc import BioCantorIOException, InvalidInputError


class GenBankParserError(InvalidInputError):
    """
    Raised when there is an error parsing a genbank file.
    """

    pass


class GenBankLocusTagError(GenBankParserError):
    """
    Raised when a locus tag error is encountered in a genbank file.
    """


class GenBankMultipleTranscriptFeatureError(GenBankParserError):
    """
    Raised when more than one transcript-level feature maps to a gene in a genbank file.
    """


class EmptyGenBankError(GenBankParserError):
    """
    Raised when a genbank parsing event returns zero genes and zero features.
    """


class GenBankExportError(BioCantorIOException):
    """
    Raised when an error is encountered exporting a genbank file.
    """

    pass


class GenBankLocationException(GenBankParserError):
    """
    Raised when there is an issue with locations (usually BioPython makes them None)
    """

    pass


class GenBankNullStrandException(GenBankParserError):
    """
    Raised when BioPython constructs a feature with a null Strand.
    """

    pass


class GenBankUnknownFeatureWarning(Warning):
    pass


class GenBankEmptyGeneWarning(Warning):
    pass


class GenBankDuplicateLocusTagWarning(Warning):
    pass


class UnknownGenBankFeatureWarning(Warning):
    pass
