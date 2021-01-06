"""
Test exporting FASTA files from various kinds of collections.
"""
from inscripta.biocantor.io.fasta.fasta import collection_to_fasta
from inscripta.biocantor.io.genbank.parser import parse_genbank
from inscripta.biocantor.io.parser import ParsedAnnotationRecord


def test_collection_to_fasta_from_genbank(test_data_dir, tmp_path):
    """This FASTA export matches exactly because there are no FASTA comments."""
    gbk = test_data_dir / "INSC1006_chrI.gbff"
    with open(gbk, "r") as fh:
        parsed = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

    tmp_fasta = tmp_path / "tmp.fasta"
    with open(tmp_fasta, "w") as fh:
        collection_to_fasta(parsed, fh)

    with open(tmp_fasta, "r") as fh1, open(test_data_dir / "INSC1006_chrI.fa", "r") as fh2:
        assert fh1.read() == fh2.read()


def test_collection_to_fasta_from_genbank_fasta_header(test_data_dir, tmp_path):
    """Inscripta_BL21.fa has FASTA comments, and so the sequence will match but the comments will be lost."""
    gbk = test_data_dir / "Inscripta_BL21.gbk"
    with open(gbk, "r") as fh:
        parsed = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

    tmp_fasta = tmp_path / "tmp.fasta"
    with open(tmp_fasta, "w") as fh:
        collection_to_fasta(parsed, fh)

    with open(tmp_fasta, "r") as fh1, open(test_data_dir / "Inscripta_BL21.fa", "r") as fh2:
        f1 = fh1.readlines()
        f2 = fh2.readlines()
        assert f1[1:] == f2[1:]
        assert f1[0] != f2[0]
        assert f1[0].split()[0] == f2[0].split()[0]


def test_records_to_fasta_from_genbank_fasta_header(test_data_dir, tmp_path):
    """Because we are exporting from a ParsedAnnotationRecord directly the FASTA comments are retained."""
    gbk = test_data_dir / "Inscripta_BL21.gbk"
    with open(gbk, "r") as fh:
        parsed = parse_genbank(fh)
        tmp_fasta = tmp_path / "tmp.fasta"
        with open(tmp_fasta, "w") as ofh:
            for rec in parsed:
                rec.to_fasta(ofh)

    with open(tmp_fasta, "r") as fh1, open(test_data_dir / "Inscripta_BL21.fa", "r") as fh2:
        assert fh1.read() == fh2.read()


def test_records_to_fasta_from_genbank(test_data_dir, tmp_path):
    """INSC1006_chrI.gbff will have FASTA comments when exported by BioPython,
    and so the sequence will match but the comments will be lost."""
    gbk = test_data_dir / "INSC1006_chrI.gbff"
    with open(gbk, "r") as fh:
        parsed = parse_genbank(fh)
        tmp_fasta = tmp_path / "tmp.fasta"
        with open(tmp_fasta, "w") as ofh:
            for rec in parsed:
                rec.to_fasta(ofh)

    with open(tmp_fasta, "r") as fh1, open(test_data_dir / "INSC1006_chrI.fa", "r") as fh2:
        f1 = fh1.readlines()
        f2 = fh2.readlines()
        assert f1[1:] == f2[1:]
        assert f1[0] != f2[0]
        assert f1[0].split()[0] == f2[0].split()[0]
