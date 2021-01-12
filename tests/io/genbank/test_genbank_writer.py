"""
Test GenBank writing.
"""
import pytest
from inscripta.biocantor.io.genbank.parser import parse_genbank
from inscripta.biocantor.io.genbank.writer import collection_to_genbank, GenbankFlavor
from inscripta.biocantor.io.parser import ParsedAnnotationRecord


@pytest.mark.parametrize(
    "gbk,flavor", [("INSC1006_chrI.gbff", GenbankFlavor.EUKARYOTIC), ("INSC1003.gbk", GenbankFlavor.PROKARYOTIC)]
)
def test_sequences(test_data_dir, tmp_path, gbk, flavor):
    """Parse a genbank, write it to disk, then parse it again and compare."""
    gbk = test_data_dir / gbk
    with open(gbk, "r") as fh:
        collections = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

    tmp_gbk = tmp_path / "tmp.gbk"
    with open(tmp_gbk, "w") as fh:
        collection_to_genbank(collections, fh, flavor)

    with open(tmp_gbk, "r") as fh:
        new_collection = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(fh)))

    assert len(collections[0].genes) == len(new_collection[0].genes)

    for gene_a, gene_b in zip(collections[0], new_collection[0]):
        assert gene_a.location == gene_b.location
        tx_a = gene_a.transcripts[0]
        tx_b = gene_b.transcripts[0]
        assert tx_a.location == tx_b.location
        if tx_a.is_coding:
            assert tx_a.cds.location == tx_b.cds.location
            assert tx_a.get_protein_sequence() == tx_b.get_protein_sequence()
        assert tx_a.get_transcript_sequence() == tx_b.get_transcript_sequence()
