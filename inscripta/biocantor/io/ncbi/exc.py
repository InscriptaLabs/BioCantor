from inscripta.biocantor.io.exc import BioCantorIOException


class TblExportException(BioCantorIOException):
    pass


class LocusTagException(TblExportException):
    pass
