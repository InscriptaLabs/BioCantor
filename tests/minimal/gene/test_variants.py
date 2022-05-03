import pytest

from inscripta.biocantor.exc import LocationOverlapException, NullSequenceException
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

    @pytest.mark.parametrize(
        "variant,expected",
        [
            [
                insertion_5_seq,
                Parent(
                    sequence=Sequence(
                        "ACTCTGGCTCTATCTCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                    ),
                    location=SingleInterval(0, 22, Strand.PLUS),
                ),
            ],
            [
                insertion_5_chunk,
                Parent(
                    id="chunk_seq_offset_100:100-122",
                    sequence=Sequence(
                        "ACTCTGGCTCTATCTCATCCAC",
                        Alphabet.NT_EXTENDED_GAPPED,
                        type=SequenceType.SEQUENCE_CHUNK,
                        id="chunk_seq_offset_100:100-122",
                        parent=Parent(
                            location=SingleInterval(
                                100,
                                122,
                                Strand.PLUS,
                                parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
            ],
            [
                deletion_11_13_seq,
                Parent(
                    sequence=Sequence(
                        "ACTCTCTCTATCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                    ),
                    location=SingleInterval(0, 18, Strand.PLUS),
                ),
            ],
            [
                deletion_11_13_chunk,
                Parent(
                    id="chunk_seq_offset_100:100-118",
                    sequence=Sequence(
                        "ACTCTCTCTATCATCCAC",
                        Alphabet.NT_EXTENDED_GAPPED,
                        type=SequenceType.SEQUENCE_CHUNK,
                        id="chunk_seq_offset_100:100-118",
                        parent=Parent(
                            location=SingleInterval(
                                100,
                                118,
                                Strand.PLUS,
                                parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
            ],
        ],
    )
    def test_parent_with_alternative_sequence(self, variant, expected):
        assert variant.parent_with_alternative_sequence == expected

    @pytest.mark.parametrize(
        "variant,location,expected",
        [
            [
                insertion_5_seq,
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(
                    0,
                    12,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence(
                            "ACTCTGGCTCTATCTCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                        )
                    ),
                ),
            ],
            # lifting 100-110 in chromosome coordinates results in 0-12 in chunk coordinates
            # because the VariantInterval is in chunk coordinates
            [
                insertion_5_chunk,
                SingleInterval(100, 110, Strand.PLUS),
                SingleInterval(
                    0,
                    12,
                    Strand.PLUS,
                    parent=Parent(
                        id="chunk_seq_offset_100:100-122",
                        sequence=Sequence(
                            "ACTCTGGCTCTATCTCATCCAC",
                            Alphabet.NT_EXTENDED_GAPPED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            id="chunk_seq_offset_100:100-122",
                            parent=Parent(
                                location=SingleInterval(
                                    100,
                                    122,
                                    Strand.PLUS,
                                    parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ],
            # lifting 0-10 in chunk coordinates results in 0-12 in chunk coordinates
            [
                insertion_5_chunk,
                SingleInterval(0, 10, Strand.PLUS, parent=chunk_seq),
                SingleInterval(
                    0,
                    12,
                    Strand.PLUS,
                    parent=Parent(
                        id="chunk_seq_offset_100:100-122",
                        sequence=Sequence(
                            "ACTCTGGCTCTATCTCATCCAC",
                            Alphabet.NT_EXTENDED_GAPPED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            id="chunk_seq_offset_100:100-122",
                            parent=Parent(
                                location=SingleInterval(
                                    100,
                                    122,
                                    Strand.PLUS,
                                    parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ],
            [
                deletion_11_13_seq,
                SingleInterval(0, 15, Strand.PLUS),
                SingleInterval(
                    0,
                    13,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence(
                            "ACTCTCTCTATCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                        )
                    ),
                ),
            ],
            [
                deletion_11_13_chunk,
                SingleInterval(100, 115, Strand.PLUS),
                SingleInterval(
                    0,
                    13,
                    Strand.PLUS,
                    parent=Parent(
                        id="chunk_seq_offset_100:100-118",
                        sequence=Sequence(
                            "ACTCTCTCTATCATCCAC",
                            Alphabet.NT_EXTENDED_GAPPED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            id="chunk_seq_offset_100:100-118",
                            parent=Parent(
                                location=SingleInterval(
                                    100,
                                    118,
                                    Strand.PLUS,
                                    parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ],
            [
                deletion_11_13_chunk,
                SingleInterval(0, 15, Strand.PLUS, parent=chunk_seq),
                SingleInterval(
                    0,
                    13,
                    Strand.PLUS,
                    parent=Parent(
                        id="chunk_seq_offset_100:100-118",
                        sequence=Sequence(
                            "ACTCTCTCTATCATCCAC",
                            Alphabet.NT_EXTENDED_GAPPED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            id="chunk_seq_offset_100:100-118",
                            parent=Parent(
                                location=SingleInterval(
                                    100,
                                    118,
                                    Strand.PLUS,
                                    parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ],
            # sequenceless liftover has no Parent
            [deletion_11_13, SingleInterval(0, 15, Strand.PLUS), SingleInterval(0, 13, Strand.PLUS)],
        ],
    )
    def test_lift_over_location(self, variant, location, expected):
        assert variant.lift_over_location(location) == expected

    def test_sequenceless_exceptions(self):
        with pytest.raises(NullSequenceException):
            _ = snp_1.alternative_genomic_sequence
        with pytest.raises(NullSequenceException):
            _ = snp_1.parent_with_alternative_sequence


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
        "variant_collection,expected",
        [
            [
                VariantIntervalCollection(
                    [snp_1_seq, insertion_5_seq, deletion_11_13_seq], parent_or_seq_chunk_parent=seq
                ),
                Parent(
                    sequence=Sequence(
                        "AGTCTGGCTCTATCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                    ),
                    location=SingleInterval(0, 20, Strand.PLUS),
                ),
            ],
            [
                VariantIntervalCollection(
                    [snp_1_chunk, insertion_5_chunk, deletion_11_13_chunk], parent_or_seq_chunk_parent=chunk_seq
                ),
                Parent(
                    id="chunk_seq_offset_100:100-120",
                    sequence=Sequence(
                        "AGTCTGGCTCTATCATCCAC",
                        Alphabet.NT_EXTENDED_GAPPED,
                        type=SequenceType.SEQUENCE_CHUNK,
                        id="chunk_seq_offset_100:100-120",
                        parent=Parent(
                            location=SingleInterval(
                                100,
                                120,
                                Strand.PLUS,
                                parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
            ],
            [
                VariantIntervalCollection(
                    [snp_1_seq, insertion_5_seq, deletion_11_13_seq, deletion_13_15_seq], parent_or_seq_chunk_parent=seq
                ),
                Parent(
                    sequence=Sequence(
                        "AGTCTGGCTCTATTCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                    ),
                    location=SingleInterval(0, 18, Strand.PLUS),
                ),
            ],
            [
                VariantIntervalCollection(
                    [snp_1_chunk, insertion_5_chunk, deletion_11_13_chunk, deletion_13_15_chunk],
                    parent_or_seq_chunk_parent=chunk_seq,
                ),
                Parent(
                    id="chunk_seq_offset_100:100-118",
                    sequence=Sequence(
                        "AGTCTGGCTCTATTCCAC",
                        Alphabet.NT_EXTENDED_GAPPED,
                        type=SequenceType.SEQUENCE_CHUNK,
                        id="chunk_seq_offset_100:100-118",
                        parent=Parent(
                            location=SingleInterval(
                                100,
                                118,
                                Strand.PLUS,
                                parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
            ],
        ],
    )
    def test_parent_with_alternative_sequence(self, variant_collection, expected):
        assert variant_collection.parent_with_alternative_sequence == expected

    @pytest.mark.parametrize(
        "variant_collection,location,expected",
        [
            [
                VariantIntervalCollection(
                    [snp_1_seq, insertion_5_seq, deletion_11_13_seq], parent_or_seq_chunk_parent=seq
                ),
                SingleInterval(0, 16, Strand.PLUS),
                SingleInterval(
                    0,
                    16,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence(
                            "AGTCTGGCTCTATCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                        )
                    ),
                ),
            ],
            [
                VariantIntervalCollection(
                    [snp_1_chunk, insertion_5_chunk, deletion_11_13_chunk], parent_or_seq_chunk_parent=chunk_seq
                ),
                SingleInterval(100, 116, Strand.PLUS),
                SingleInterval(
                    0,
                    16,
                    Strand.PLUS,
                    Parent(
                        id="chunk_seq_offset_100:100-120",
                        sequence=Sequence(
                            "AGTCTGGCTCTATCATCCAC",
                            Alphabet.NT_EXTENDED_GAPPED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            id="chunk_seq_offset_100:100-120",
                            parent=Parent(
                                location=SingleInterval(
                                    100,
                                    120,
                                    Strand.PLUS,
                                    parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ],
            [
                VariantIntervalCollection(
                    [snp_1_seq, insertion_5_seq, deletion_11_13_seq, deletion_13_15_seq], parent_or_seq_chunk_parent=seq
                ),
                SingleInterval(0, 16, Strand.PLUS),
                SingleInterval(
                    0,
                    14,
                    Strand.PLUS,
                    parent=Parent(
                        sequence=Sequence(
                            "AGTCTGGCTCTATTCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                        )
                    ),
                ),
            ],
            [
                VariantIntervalCollection(
                    [snp_1_chunk, insertion_5_chunk, deletion_11_13_chunk, deletion_13_15_chunk],
                    parent_or_seq_chunk_parent=chunk_seq,
                ),
                SingleInterval(100, 116, Strand.PLUS),
                SingleInterval(
                    0,
                    14,
                    Strand.PLUS,
                    Parent(
                        id="chunk_seq_offset_100:100-118",
                        sequence=Sequence(
                            "AGTCTGGCTCTATTCCAC",
                            Alphabet.NT_EXTENDED_GAPPED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            id="chunk_seq_offset_100:100-118",
                            parent=Parent(
                                location=SingleInterval(
                                    100,
                                    118,
                                    Strand.PLUS,
                                    parent=Parent(id="chunk_seq_offset_100", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ],
        ],
    )
    def test_lift_over_location(self, variant_collection, location, expected):
        assert variant_collection.lift_over_location(location) == expected

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
    def test_constructor_exceptions(self, variants, expected_exception):
        with pytest.raises(expected_exception):
            _ = VariantIntervalCollection(variants)

    def test_sequenceless_exceptions(self):
        collection = VariantIntervalCollection([snp_1, deletion_11_13])
        with pytest.raises(NullSequenceException):
            _ = collection.alternative_genomic_sequence
        with pytest.raises(NullSequenceException):
            _ = collection.parent_with_alternative_sequence
