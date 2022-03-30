import pytest

from inscripta.biocantor.exc import LocationOverlapException
from inscripta.biocantor.gene.variants import VariantInterval, VariantIntervalCollection
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence, Alphabet

seq = Parent(sequence=Sequence("ACTCTCTCTATCTCATCCAC", Alphabet.NT_EXTENDED_GAPPED))
snp_1 = VariantInterval(start=1, end=2, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=seq)
insertion_5 = VariantInterval(start=5, end=5, sequence="GG", variant_type="insertion", parent_or_seq_chunk_parent=seq)
deletion_10_13 = VariantInterval(
    start=10, end=13, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=seq
)
deletion_12_13 = VariantInterval(start=12, end=13, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=seq)
deletion_13_15 = VariantInterval(start=13, end=15, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=seq)


class TestVariants:
    @pytest.mark.parametrize(
        "variant,expected",
        [
            [snp_1, "AGTCTCTCTATCTCATCCAC"],
            [insertion_5, "ACTCTGGCTCTATCTCATCCAC"],
            [deletion_10_13, "ACTCTCTCTATCATCCAC"],
            [deletion_13_15, "ACTCTCTCTATCTTCCAC"],
        ],
    )
    def test_edited_sequence(self, variant, expected):
        assert str(variant.edited_genomic_sequence) == expected


class TestVariantCollections:
    @pytest.mark.parametrize(
        "variant_collection,expected",
        [
            [
                VariantIntervalCollection([snp_1, insertion_5, deletion_10_13], parent_or_seq_chunk_parent=seq),
                "AGTCTGGCTCTATCATCCAC",
            ],
            [
                VariantIntervalCollection(
                    [snp_1, insertion_5, deletion_10_13, deletion_13_15], parent_or_seq_chunk_parent=seq
                ),
                "AGTCTGGCTCTATTCCAC",
            ],
        ],
    )
    def test_edited_sequence(self, variant_collection, expected):
        assert str(variant_collection.edited_genomic_sequence) == expected

    @pytest.mark.parametrize(
        "variants,expected_exception",
        [
            [
                [snp_1, snp_1],
                LocationOverlapException,
            ],
            [
                [deletion_10_13, deletion_12_13],
                LocationOverlapException,
            ],
        ],
    )
    def test_exceptions(self, variants, expected_exception):
        with pytest.raises(expected_exception):
            _ = VariantIntervalCollection(variants)
