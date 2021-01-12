"""
Test writing GFF3 files.
"""
import pytest

from inscripta.biocantor.io.genbank.parser import parse_genbank
from inscripta.biocantor.io.gff3.parser import parse_standard_gff3
from inscripta.biocantor.io.gff3.writer import collection_to_gff3
from inscripta.biocantor.io.parser import ParsedAnnotationRecord


class TestGff3Writer:
    """Test GFF3 writing"""

    @pytest.mark.parametrize(
        "gbk,gff3,add_sequences",
        [("INSC1006_chrI.gbff", "INSC1006_chrI.gff3", True), ("INSC1003.gbk", "INSC1003.gff3", False)],
    )
    def test_genbank_to_gff(self, test_data_dir, tmp_path, gbk, gff3, add_sequences):
        """
        INSC1006_chrI.gff3 and INSC1003.gff3 were created from INSC1006_chrI.gbff and INSC1003.gbk
        respectively, so we can compare to the source file.
        """
        gbk = test_data_dir / gbk
        with open(gbk, "r") as fh:
            parsed = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

        tmp_gff = tmp_path / "tmp.gff"
        with open(tmp_gff, "w") as fh:
            collection_to_gff3(parsed, fh, add_sequences=add_sequences)

        for l1, l2 in zip(open(tmp_gff), open(test_data_dir / gff3)):
            assert l1 == l2

    def test_to_gff(self, test_data_dir, tmp_path):
        """Parse GFF, write to disk, parse, compare"""
        gff = test_data_dir / "INSC1006_chrI.gff3"
        parsed = list(parse_standard_gff3(gff))
        a = [x.to_annotation_collection() for x in parsed]

        tmp_gff = tmp_path / "tmp.gff"
        with open(tmp_gff, "w") as fh:
            collection_to_gff3(a, fh)
        gff2 = list(parse_standard_gff3(tmp_gff))
        a2 = [x.to_annotation_collection() for x in gff2]

        for c1, c2 in zip(a, a2):
            assert c1.to_dict() == c2.to_dict()
