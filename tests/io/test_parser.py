import pytest
from biocantor.location.location_impl import SingleInterval
from biocantor.parent import Parent, SequenceType
from biocantor.sequence.alphabet import Alphabet
from biocantor.sequence.sequence import Sequence, Strand
from biocantor.io.parser import seq_chunk_to_parent, seq_to_parent


def test_seq_to_parent():
    seq = "ATGCATGC"
    seq_id = "TestSeq"
    obs = seq_to_parent(seq, seq_id=seq_id)
    assert obs == Parent(
        sequence=Sequence(seq, Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME, id=seq_id),
        location=SingleInterval(0, len(seq), Strand.PLUS),
    )


@pytest.mark.parametrize("strand", [Strand.PLUS, Strand.MINUS])
def test_seq_chunk_to_parent(strand):
    obs = seq_chunk_to_parent("ATGCATGC", "TestSeq", 200, 208, strand=strand)
    assert obs == Parent(
        id="TestSeq:200-208",
        sequence_type=SequenceType.SEQUENCE_CHUNK,
        strand=None,
        location=None,
        sequence=Sequence(
            data="ATGCATGC",
            id="TestSeq:200-208",
            alphabet=Alphabet.NT_EXTENDED_GAPPED,
            type=SequenceType.SEQUENCE_CHUNK,
            parent=Parent(
                id="TestSeq",
                sequence_type=SequenceType.CHROMOSOME,
                strand=strand,
                location=SingleInterval(
                    200,
                    208,
                    strand,
                    parent=Parent(
                        id="TestSeq",
                        sequence_type=SequenceType.CHROMOSOME,
                        strand=strand,
                        location=SingleInterval(200, 208, strand),
                        sequence=None,
                        parent=None,
                    ),
                ),
                sequence=None,
                parent=Parent(
                    id="TestSeq",
                    sequence_type=SequenceType.CHROMOSOME,
                    strand=strand,
                    location=SingleInterval(200, 208, strand, parent=None),
                    sequence=None,
                    parent=None,
                ),
            ),
        ),
    )
    assert obs.has_ancestor_of_type(SequenceType.CHROMOSOME)
    assert obs.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK)
