from inscripta.biocantor.io.exc import InvalidInputError


class GenBankValidationError(InvalidInputError):
    pass


class GenBankParserError(InvalidInputError):
    pass


class GenBankExportError(InvalidInputError):
    pass
