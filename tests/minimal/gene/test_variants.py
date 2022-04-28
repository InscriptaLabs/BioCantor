import pytest

from inscripta.biocantor.exc import LocationOverlapException
from inscripta.biocantor.gene.variants import VariantInterval, VariantIntervalCollection
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence.sequence import SequenceType, Sequence, Alphabet
from inscripta.biocantor.location import SingleInterval, Strand

snp_1 = VariantInterval(start=1, end=2, sequence="G", variant_type="SNV")
insertion_5 = VariantInterval(start=5, end=6, sequence="GGC", variant_type="insertion")
# this deletion is VCF-style and contains a left-anchoring base
deletion_11_13 = VariantInterval(start=10, end=13, sequence="T", variant_type="deletion")
# these deletions are not VCF-style, and have an empty ALT
deletion_12_13 = VariantInterval(start=12, end=13, sequence="", variant_type="deletion")
deletion_13_15 = VariantInterval(start=13, end=15, sequence="", variant_type="deletion")


seq = Parent(sequence=Sequence("ACTCTCTCTATCTCATCCAC", Alphabet.NT_EXTENDED_GAPPED))
snp_1_seq = VariantInterval(start=1, end=2, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=seq)
insertion_5_seq = VariantInterval(
    start=5, end=6, sequence="GGC", variant_type="insertion", parent_or_seq_chunk_parent=seq
)
# this deletion is VCF-style and contains a left-anchoring base
deletion_11_13_seq = VariantInterval(
    start=10, end=13, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=seq
)
# these deletions are not VCF-style, and have an empty ALT
deletion_12_13_seq = VariantInterval(
    start=12, end=13, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=seq
)
deletion_13_15_seq = VariantInterval(
    start=13, end=15, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=seq
)

chunk_seq = Parent(
    id="chunk_seq_offset_100",
    sequence=Sequence(
        str(seq.sequence),
        Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            location=SingleInterval(
                100,
                100 + len(str(seq.sequence)),
                Strand.PLUS,
                parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
            )
        ),
    ),
)
snp_1_chunk = VariantInterval(
    start=101, end=102, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=chunk_seq
)
insertion_5_chunk = VariantInterval(
    start=105, end=106, sequence="GGC", variant_type="insertion", parent_or_seq_chunk_parent=chunk_seq
)
# this deletion is VCF-style and contains a left-anchoring base
deletion_11_13_chunk = VariantInterval(
    start=110, end=113, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=chunk_seq
)
# these deletions are not VCF-style, and have an empty ALT
deletion_12_13_chunk = VariantInterval(
    start=112, end=113, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=chunk_seq
)
deletion_13_15_chunk = VariantInterval(
    start=113, end=115, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=chunk_seq
)


class TestVariants:
    @pytest.mark.parametrize(
        "variant,expected",
        [
            [snp_1_seq, "AGTCTCTCTATCTCATCCAC"],
            [snp_1_chunk, "AGTCTCTCTATCTCATCCAC"],
            [insertion_5_seq, "ACTCTGGCTCTATCTCATCCAC"],
            [insertion_5_chunk, "ACTCTGGCTCTATCTCATCCAC"],
            [deletion_11_13_seq, "ACTCTCTCTATCATCCAC"],
            [deletion_11_13_chunk, "ACTCTCTCTATCATCCAC"],
            [deletion_13_15_seq, "ACTCTCTCTATCTTCCAC"],
            [deletion_13_15_chunk, "ACTCTCTCTATCTTCCAC"],
        ],
    )
    def test_alternative_sequence(self, variant, expected):
        assert str(variant.alternative_genomic_sequence) == expected


class TestVariantCollections:
    @pytest.mark.parametrize(
        "variant_collection,expected",
        [
            [
                VariantIntervalCollection(
                    [snp_1_seq, insertion_5_seq, deletion_11_13_seq], parent_or_seq_chunk_parent=seq
                ),
                "AGTCTGGCTCTATCATCCAC",
            ],
            [
                VariantIntervalCollection(
                    [snp_1_chunk, insertion_5_chunk, deletion_11_13_chunk], parent_or_seq_chunk_parent=chunk_seq
                ),
                "AGTCTGGCTCTATCATCCAC",
            ],
            [
                VariantIntervalCollection(
                    [snp_1_seq, insertion_5_seq, deletion_11_13_seq, deletion_13_15_seq], parent_or_seq_chunk_parent=seq
                ),
                "AGTCTGGCTCTATTCCAC",
            ],
            [
                VariantIntervalCollection(
                    [snp_1_chunk, insertion_5_chunk, deletion_11_13_chunk, deletion_13_15_chunk],
                    parent_or_seq_chunk_parent=chunk_seq,
                ),
                "AGTCTGGCTCTATTCCAC",
            ],
        ],
    )
    def test_alternative_sequence(self, variant_collection, expected):
        assert str(variant_collection.alternative_genomic_sequence) == expected

    @pytest.mark.parametrize(
        "variants,expected_exception",
        [
            [
                [snp_1, snp_1],
                LocationOverlapException,
            ],
            [
                [deletion_11_13, deletion_12_13],
                LocationOverlapException,
            ],
        ],
    )
    def test_exceptions(self, variants, expected_exception):
        with pytest.raises(expected_exception):
            _ = VariantIntervalCollection(variants)
