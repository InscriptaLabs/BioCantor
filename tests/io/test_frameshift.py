"""
Prove that we can handle -1 frameshifts properly when modeled in the input data.
"""
import pytest
from inscripta.biocantor.exc import MismatchedFrameException
from inscripta.biocantor.io.gff3.parser import parse_gff3_embedded_fasta
from inscripta.biocantor.io.genbank.parser import parse_genbank
from inscripta.biocantor.io.parser import ParsedAnnotationRecord


class TestParseFrameshifts:
    def test_parse_inso(self, test_data_dir):
        """This proves we handle frame and phase"""
        gbk = test_data_dir / "insO_frameshift.gbk"
        gff3 = test_data_dir / "insO_frameshift.gff3"

        with open(gbk, "r") as fh:
            gbk_rec = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        gff3_rec = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3)))[0]

        expected_protein = (
            "MKKRNFSAEFKRESAQLVVDQKYTVADAAKAMDVGLSTMTRWVKQLRDERQGKTPKASPITPEQIEIRKLRKKLQRIEMENEILKKNRP"
            "EKPDGRRAVLRSQVLELHGISHGSAGARSIATMATRRGYQMGRWLAGRLMKELGLVSCQQPTHRYKRGGHEHVAIPNYLERQFAVTEPNQV"
            "WCGDVTYIWTGKRWAYLAVVLDLFARKPVGWAMSFSPDSRLTMKALEMAWETRGKPVGVMFQSDQGSHYTSRQFRQLLWRYRIRQSMSRR"
            "GNCWDNSPMERFFRSLKNEWVPATGYVSFSDAAHAITDYIVGYYSALRPHEYNGGLPPNESENRYWKNSNAEASFS*"
        )
        assert (
            str(gbk_rec.genes[0].get_primary_protein())
            == str(gff3_rec.genes[0].get_primary_protein())
            == expected_protein
        )

    def test_broken_frameshift(self, test_data_dir):
        """If I merge the transcript, the frames list no longer matches the location and an exception is raised."""
        gbk = test_data_dir / "insO_frameshift.gbk"
        with open(gbk, "r") as fh:
            gbk_rec = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))[0]

        cds = gbk_rec.genes[0].get_primary_transcript().cds
        cds._location = cds._location.merge_overlapping()
        with pytest.raises(MismatchedFrameException):
            _ = cds.translate()

    def test_parse_peg10(self, test_data_dir):
        """PEG10 is a human gene with a -1 frameshift"""
        gff3 = test_data_dir / "PEG10_offset_gff3_fasta.gff3"
        gff3_rec = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_gff3_embedded_fasta(gff3)))[0]
        tx = gff3_rec.genes[0].transcripts[0]
        assert not tx.has_in_frame_stop
