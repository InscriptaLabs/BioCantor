from inscripta.biocantor.exc import BioCantorException


class BEDExportException(BioCantorException):
    pass


class BEDMissingSequenceNameError(BEDExportException):
    pass
