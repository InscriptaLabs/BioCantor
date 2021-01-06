"""
Test exporting to NCBI TBL format.

All of these TBL files have been validated to pass the tbl2asn error validator tool. This tool is not easy
to acquire and so is not packaged for these unit tests.
"""
import pytest
from inscripta.biocantor.io.ncbi.tbl_writer import collection_to_tbl, LocusTagException
from inscripta.biocantor.io.gff3.parser import parse_gff3_embedded_fasta
from inscripta.biocantor.io.parser import ParsedAnnotationRecord


@pytest.mark.parametrize("gff3,expected_tbl",
                         [
                             (
                                "Inscripta_BL21_embedded_extra_contig.gff3", "Inscripta_BL21_embedded_extra_contig.tbl"
                             ),
                             (
                                 "insO_frameshift.gff3", "insO_frameshift.tbl"
                             ),
                             (
                                 "PEG10_offset_gff3_fasta.gff3", "PEG10_offset_gff3_fasta.tbl"
                             )
                         ])
def test_tbl_export_from_gff3(test_data_dir, tmp_path, gff3, expected_tbl):
    gff3 = test_data_dir / gff3
    recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3)))
    tmp = tmp_path / "tmp.tbl"
    with open(tmp, "w") as fh:
        collection_to_tbl(recs, fh, locus_tag_prefix="test")
    with open(tmp) as fh1, open(test_data_dir / expected_tbl) as fh2:
        assert fh1.read() == fh2.read()


def test_export_exception(test_data_dir, tmp_path):
    """Files without locus tags on all features raise an exception if no locus tag is provided"""
    gff3 = test_data_dir / "insO_frameshift.gff3"
    recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3)))
    tmp = tmp_path / "tmp.tbl"
    with pytest.raises(LocusTagException):
        with open(tmp, "w") as fh:
            collection_to_tbl(recs, fh)
