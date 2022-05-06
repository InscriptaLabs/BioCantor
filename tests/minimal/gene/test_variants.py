from uuid import UUID

import pytest

from inscripta.biocantor.exc import LocationOverlapException, NullSequenceException
from inscripta.biocantor.gene import (
    GeneInterval,
    AnnotationCollection,
    TranscriptInterval,
    FeatureInterval,
    FeatureIntervalCollection,
    CDSInterval,
    CDSFrame,
)
from inscripta.biocantor.gene.variants import VariantInterval, VariantIntervalCollection
from inscripta.biocantor.location import SingleInterval, Strand
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence.sequence import SequenceType, Sequence, Alphabet

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

snp1_ins5_del11to13_seq = VariantIntervalCollection(
    [
        VariantInterval(start=1, end=2, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=seq),
        VariantInterval(start=5, end=6, sequence="GGC", variant_type="insertion", parent_or_seq_chunk_parent=seq),
        VariantInterval(start=10, end=13, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=seq),
    ],
    parent_or_seq_chunk_parent=seq,
)
snp1_ins5_del11to13_seq_chunk = VariantIntervalCollection(
    [
        VariantInterval(start=101, end=102, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=chunk_seq),
        VariantInterval(
            start=105, end=106, sequence="GGC", variant_type="insertion", parent_or_seq_chunk_parent=chunk_seq
        ),
        VariantInterval(
            start=110, end=113, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=chunk_seq
        ),
    ],
    parent_or_seq_chunk_parent=chunk_seq,
)

snp1_ins5_del11to15_seq = VariantIntervalCollection(
    [
        VariantInterval(start=1, end=2, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=seq),
        VariantInterval(start=5, end=6, sequence="GGC", variant_type="insertion", parent_or_seq_chunk_parent=seq),
        VariantInterval(start=10, end=13, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=seq),
        VariantInterval(start=13, end=15, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=seq),
    ],
    parent_or_seq_chunk_parent=seq,
)

snp1_ins5_del11to15_seq_chunk = VariantIntervalCollection(
    [
        VariantInterval(start=101, end=102, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=chunk_seq),
        VariantInterval(
            start=105, end=106, sequence="GGC", variant_type="insertion", parent_or_seq_chunk_parent=chunk_seq
        ),
        VariantInterval(
            start=110, end=113, sequence="T", variant_type="deletion", parent_or_seq_chunk_parent=chunk_seq
        ),
        VariantInterval(start=113, end=115, sequence="", variant_type="deletion", parent_or_seq_chunk_parent=chunk_seq),
    ],
    parent_or_seq_chunk_parent=chunk_seq,
)


class TestVariants:
    @pytest.mark.parametrize(
        "variant,expected",
        [
            # SNP
            # ref: ACTCTCTCTATCTCATCCAC
            # alt: AGTCTCTCTATCTCATCCAC
            [snp_1_seq, "AGTCTCTCTATCTCATCCAC"],
            [snp_1_chunk, "AGTCTCTCTATCTCATCCAC"],
            # GG insertion at 5
            # ref: ACTCT  CTCTATCTCATCCAC
            # alt: ACTCTGGCTCTATCTCATCCAC
            [insertion_5_seq, "ACTCTGGCTCTATCTCATCCAC"],
            [insertion_5_chunk, "ACTCTGGCTCTATCTCATCCAC"],
            # deletion from 11-13
            # ref: ACTCTCTCTATCTCATCCAC
            # alt: ACTCTCTCTAT  CATCCAC
            [deletion_11_13_seq, "ACTCTCTCTATCATCCAC"],
            [deletion_11_13_chunk, "ACTCTCTCTATCATCCAC"],
            # deletion at 13-15
            # ref: ACTCTCTCTATCTCATCCAC
            # alt: ACTCTCTCTATCT  TCCAC
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
        snp = VariantInterval(start=1, end=2, sequence="G", variant_type="SNV")
        with pytest.raises(NullSequenceException):
            _ = snp.alternative_genomic_sequence
        with pytest.raises(NullSequenceException):
            _ = snp.parent_with_alternative_sequence


class TestVariantCollections:
    @pytest.mark.parametrize(
        "variant_collection,expected",
        [
            # SNP + insertion + deletion
            #      0 1 2 3 4     5 6 7 8 9 10111213141516171819
            # ref: A C T C T     C T C T A T C T C A T C C A C
            #
            # alt: A G T C T G G C T C T A T     C A T C C A C
            #        ^
            [
                snp1_ins5_del11to13_seq,
                "AGTCTGGCTCTATCATCCAC",
            ],
            [
                snp1_ins5_del11to13_seq_chunk,
                "AGTCTGGCTCTATCATCCAC",
            ],
            # SNP + insertion + deletion + deletion (abutting deletions)
            #      0 1 2 3 4     5 6 7 8 9 10111213141516171819
            # ref: A C T C T     C T C T A T C T C A T C C A C
            #
            # alt: A G T C T G G C T C T A T         T C C A C
            #        ^
            [
                snp1_ins5_del11to15_seq,
                "AGTCTGGCTCTATTCCAC",
            ],
            [
                snp1_ins5_del11to15_seq_chunk,
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
                snp1_ins5_del11to13_seq,
                Parent(
                    sequence=Sequence(
                        "AGTCTGGCTCTATCATCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                    ),
                    location=SingleInterval(0, 20, Strand.PLUS),
                ),
            ],
            [
                snp1_ins5_del11to13_seq_chunk,
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
                snp1_ins5_del11to15_seq,
                Parent(
                    sequence=Sequence(
                        "AGTCTGGCTCTATTCCAC", type=SequenceType.CHROMOSOME, alphabet=Alphabet.NT_EXTENDED_GAPPED
                    ),
                    location=SingleInterval(0, 18, Strand.PLUS),
                ),
            ],
            [
                snp1_ins5_del11to15_seq_chunk,
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
                snp1_ins5_del11to13_seq,
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
                snp1_ins5_del11to13_seq_chunk,
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
                snp1_ins5_del11to15_seq,
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
                snp1_ins5_del11to15_seq_chunk,
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
                [
                    VariantInterval(start=1, end=2, sequence="G", variant_type="SNV"),
                    VariantInterval(start=1, end=2, sequence="G", variant_type="SNV"),
                ],
                LocationOverlapException,
            ],
            [
                [
                    VariantInterval(start=10, end=13, sequence="T", variant_type="deletion"),
                    VariantInterval(start=12, end=13, sequence="", variant_type="deletion"),
                ],
                LocationOverlapException,
            ],
        ],
    )
    def test_constructor_exceptions(self, variants, expected_exception):
        with pytest.raises(expected_exception):
            _ = VariantIntervalCollection(variants)

    def test_sequenceless_exceptions(self):
        collection = VariantIntervalCollection(
            [
                VariantInterval(start=1, end=2, sequence="G", variant_type="SNV"),
                VariantInterval(start=10, end=13, sequence="T", variant_type="deletion"),
            ]
        )
        with pytest.raises(NullSequenceException):
            _ = collection.alternative_genomic_sequence
        with pytest.raises(NullSequenceException):
            _ = collection.parent_with_alternative_sequence


class TestVariantIncorporation:
    """
    Test incorporating variants into cds/transcript/feature/genes/etc
    """

    @pytest.mark.parametrize(
        "variant,cds,expected_sequence",
        [
            # SNP
            # ref: ACTCTCTCTATCTCA
            # alt: AGTCTCTCTATCTCA
            (
                snp_1_seq,
                CDSInterval([0], [15], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=seq),
                "AGTCTCTCTATCTCA",
            ),
            (
                snp_1_chunk,
                CDSInterval([100], [115], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chunk_seq),
                "AGTCTCTCTATCTCA",
            ),
            # insertion
            # the trailing CA in the reference is not in the final sequence because it is no longer
            # a multiple of three
            # genomic (copied from above):
            # GG insertion at 5
            # ref: ACTCT  CTCTATCTCATCCAC
            # alt: ACTCTGGCTCTATCTCATCCAC
            #
            # CDS:
            # ref: ACTCT  CTCTATCTCA
            # alt: ACTCTGGCTCTATCTCA
            #      ^^^VVV^^^VVV^^^
            (
                insertion_5_seq,
                CDSInterval([0], [15], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=seq),
                "ACTCTGGCTCTATCT",
            ),
            (
                insertion_5_chunk,
                CDSInterval([100], [115], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chunk_seq),
                "ACTCTGGCTCTATCT",
            ),
            # deletion
            # the trailing A in the reference is not in the final sequence because it is no longer
            # a multiple of three
            #
            # genomic (copied from above):
            # deletion from 11-13
            # ref: ACTCTCTCTATCTCATCCAC
            # alt: ACTCTCTCTAT  CATCCAC
            #
            # CDS:
            # ref: ACTCTCTCTATCTCA
            # alt: ACTCTCTCTAT  C
            #      ^^^VVV^^^VV  V
            (
                deletion_11_13_seq,
                CDSInterval([0], [15], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=seq),
                "ACTCTCTCTATC",
            ),
            (
                deletion_11_13_chunk,
                CDSInterval([100], [115], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chunk_seq),
                "ACTCTCTCTATC",
            ),
        ],
    )
    def test_cds_variant(self, variant, cds, expected_sequence):
        assert str(cds.incorporate_variants(variant).extract_sequence()) == expected_sequence

    @pytest.mark.parametrize(
        "variant,cds,expected_sequence",
        [
            # SNP + insertion + deletion + deletion (abutting deletions)
            #      0 1 2 3 4     5 6 7 8 9 10111213141516171819
            # ref: A C T C T     C T C T A T C T C A T C C A C
            #
            # alt: A G T C T G G C T C T A T         T C C A C
            #        ^
            # ORF: ^ ^ ^ V V V ^ ^ ^ V V V
            (
                snp1_ins5_del11to15_seq,
                CDSInterval([0], [15], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=seq),
                "AGTCTGGCTCTA",
            ),
            (
                snp1_ins5_del11to15_seq_chunk,
                CDSInterval([100], [115], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chunk_seq),
                "AGTCTGGCTCTA",
            ),
        ],
    )
    def test_cds_variant_collection(self, variant, cds, expected_sequence):
        assert str(cds.incorporate_variants(variant).extract_sequence()) == expected_sequence

    @pytest.mark.parametrize(
        "variant,transcript,expected_sequence",
        [
            (
                snp_1_seq,
                TranscriptInterval([0], [15], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                "AGTCTCTCTATCTCA",
            ),
            (
                snp_1_chunk,
                TranscriptInterval([100], [115], Strand.PLUS, parent_or_seq_chunk_parent=chunk_seq),
                "AGTCTCTCTATCTCA",
            ),
            (
                deletion_11_13_seq,
                TranscriptInterval([0], [15], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                "ACTCTCTCTATCA",
            ),
            (
                deletion_11_13_chunk,
                TranscriptInterval([100], [115], Strand.PLUS, parent_or_seq_chunk_parent=chunk_seq),
                "ACTCTCTCTATCA",
            ),
        ],
    )
    def test_transcript_variant(self, variant, transcript, expected_sequence):
        assert str(transcript.incorporate_variants(variant).get_genomic_sequence()) == expected_sequence

    @pytest.mark.parametrize(
        "variant,transcript,expected_sequence",
        [
            # SNP + insertion + deletion + deletion (abutting deletions)
            #      0 1 2 3 4     5 6 7 8 9 10111213141516171819
            # ref: A C T C T     C T C T A T C T C A T C C A C
            #
            # alt: A G T C T G G C T C T A T
            (
                snp1_ins5_del11to15_seq,
                TranscriptInterval([0], [15], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                "AGTCTGGCTCTAT",
            ),
            (
                snp1_ins5_del11to15_seq_chunk,
                TranscriptInterval([100], [115], Strand.PLUS, parent_or_seq_chunk_parent=chunk_seq),
                "AGTCTGGCTCTAT",
            ),
        ],
    )
    def test_transcript_variant_collection(self, variant, transcript, expected_sequence):
        assert str(transcript.incorporate_variants(variant).get_genomic_sequence()) == expected_sequence

    @pytest.mark.parametrize(
        "variant,gene,expected_sequences",
        [
            # SNP + insertion + deletion + deletion (abutting deletions)
            #      0 1 2 3 4     5 6 7 8 9 10111213141516171819
            # ref: A C T C T     C T C T A T C T C A T C C A C
            # spliced isoform 1 (0-3:+, 5-10:+, 12-18:+)
            #      A C T         C T C T A T     C A T C C
            #
            #      0 1 2 3 4 5 6 7 8 9 101112        1314151617
            # alt: A G T C T G G C T C T A T         T C C A C
            # spliced isoform 1: (0-3:+, 5-11:+, 12-16:+)
            #      A G T     G G C T C     T         T C C
            (
                snp1_ins5_del11to15_seq,
                GeneInterval(
                    [
                        TranscriptInterval([0], [15], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                        TranscriptInterval([0], [15], Strand.MINUS, parent_or_seq_chunk_parent=seq),
                        TranscriptInterval([0, 5, 12], [3, 10, 18], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                    ],
                    parent_or_seq_chunk_parent=seq,
                ),
                ["AGTCTGGCTCTAT", "ATAGAGCCAGACT", "AGTGGCTCTTTCC"],
            )
        ],
    )
    def test_gene_variant_collection(self, variant, gene, expected_sequences):
        new_gene = gene.incorporate_variants(variant)
        for tx, expected_sequence in zip(new_gene.transcripts, expected_sequences):
            assert str(tx.get_spliced_sequence()) == expected_sequence

    @pytest.mark.parametrize(
        "variant,annotation_collection,expected_sequences",
        [
            (
                snp1_ins5_del11to15_seq,
                AnnotationCollection(
                    [
                        FeatureIntervalCollection(
                            [FeatureInterval([0], [15], Strand.PLUS, parent_or_seq_chunk_parent=seq)],
                            parent_or_seq_chunk_parent=seq,
                        )
                    ],
                    [
                        GeneInterval(
                            [
                                TranscriptInterval([0], [15], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                                TranscriptInterval([0], [15], Strand.MINUS, parent_or_seq_chunk_parent=seq),
                                TranscriptInterval(
                                    [0, 5, 12], [3, 10, 18], Strand.PLUS, parent_or_seq_chunk_parent=seq
                                ),
                            ],
                            parent_or_seq_chunk_parent=seq,
                        )
                    ],
                    parent_or_seq_chunk_parent=seq,
                ),
                ["AGTCTGGCTCTAT", "ATAGAGCCAGACT", "AGTGGCTCTTTCC", "AGTCTGGCTCTAT"],
            ),
            (
                snp1_ins5_del11to15_seq_chunk,
                AnnotationCollection(
                    [
                        FeatureIntervalCollection(
                            [FeatureInterval([100], [115], Strand.PLUS, parent_or_seq_chunk_parent=chunk_seq)],
                            parent_or_seq_chunk_parent=chunk_seq,
                        ),
                    ],
                    [
                        GeneInterval(
                            [
                                TranscriptInterval([100], [115], Strand.PLUS, parent_or_seq_chunk_parent=chunk_seq),
                                TranscriptInterval([100], [115], Strand.MINUS, parent_or_seq_chunk_parent=chunk_seq),
                                TranscriptInterval(
                                    [100, 105, 112], [103, 110, 118], Strand.PLUS, parent_or_seq_chunk_parent=chunk_seq
                                ),
                            ],
                            parent_or_seq_chunk_parent=chunk_seq,
                        )
                    ],
                    parent_or_seq_chunk_parent=chunk_seq,
                ),
                ["AGTCTGGCTCTAT", "ATAGAGCCAGACT", "AGTGGCTCTTTCC", "AGTCTGGCTCTAT"],
            ),
        ],
    )
    def test_annotation_collection(self, variant, annotation_collection, expected_sequences):
        new_collection = annotation_collection.incorporate_variants(variant)
        i = 0
        for child in new_collection:
            for leaf in child:
                assert str(leaf.get_spliced_sequence()) == expected_sequences[i]
                i += 1


class TestAnnotationCollection:
    ac1 = AnnotationCollection(
        [FeatureIntervalCollection([FeatureInterval([0], [15], Strand.PLUS)])],
        [
            GeneInterval(
                [
                    # this transcript does not overlap the variant, however ends up in the map
                    # because the gene does
                    TranscriptInterval([17], [20], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                    TranscriptInterval([0, 5, 12], [3, 10, 18], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                ],
                parent_or_seq_chunk_parent=seq,
            ),
            # this gene does not overlap, does not end up in the map
            GeneInterval(
                [
                    # this transcript does not overlap the variant, however ends up in the map
                    # because the gene does
                    TranscriptInterval([17], [20], Strand.PLUS, parent_or_seq_chunk_parent=seq),
                ],
                parent_or_seq_chunk_parent=seq,
            ),
        ],
        [
            VariantIntervalCollection(
                [
                    VariantInterval(start=1, end=2, sequence="G", variant_type="SNV", parent_or_seq_chunk_parent=seq),
                    VariantInterval(
                        start=5,
                        end=6,
                        sequence="GGC",
                        variant_type="insertion",
                        parent_or_seq_chunk_parent=seq,
                    ),
                    VariantInterval(
                        start=10,
                        end=13,
                        sequence="T",
                        variant_type="deletion",
                        parent_or_seq_chunk_parent=seq,
                    ),
                    VariantInterval(
                        start=13,
                        end=15,
                        sequence="",
                        variant_type="deletion",
                        parent_or_seq_chunk_parent=seq,
                    ),
                ],
                guid=UUID("d972d0aa-ca79-9f2e-e985-86c16f23797e"),
                parent_or_seq_chunk_parent=seq,
            )
        ],
        parent_or_seq_chunk_parent=seq,
    )

    @pytest.mark.parametrize(
        "annotation_collection_with_variant,expected_map",
        [
            (
                ac1,
                "{UUID('d972d0aa-ca79-9f2e-e985-86c16f23797e'): "
                "[GeneInterval(identifiers=set(), Intervals:TranscriptInterval((15-18:+), cds=[None], symbol=None),"
                "TranscriptInterval((0-3:+, 5-11:+, 12-16:+), cds=[None], symbol=None)), "
                "FeatureIntervalCollection(identifiers=set(), Intervals:FeatureInterval((0-13:+), name=None))]}",
            )
        ],
    )
    def test_annotation_collection(self, annotation_collection_with_variant, expected_map):
        assert str(annotation_collection_with_variant.alternative_haplotype_mapping) == expected_map
