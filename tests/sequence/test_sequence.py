import pytest
from Bio.Seq import Seq

from inscripta.biocantor.exc import (
    AlphabetError,
    NoSuchAncestorException,
    EmptySequenceFastaError,
    InvalidStrandException,
)
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.location.strand import Strand


class TestAlphabet:
    def test_is_nucleotide(self):
        assert Alphabet.NT_STRICT.is_nucleotide_alphabet()
        assert Alphabet.NT_STRICT_GAPPED.is_nucleotide_alphabet()
        assert Alphabet.NT_EXTENDED.is_nucleotide_alphabet()
        assert Alphabet.NT_EXTENDED_GAPPED.is_nucleotide_alphabet()
        assert not Alphabet.AA.is_nucleotide_alphabet()
        assert not Alphabet.GENERIC.is_nucleotide_alphabet()


class TestSequence:
    def test_init(self):
        sequence = Sequence(
            "ACTG",
            Alphabet.NT_STRICT,
            id="id",
            type="seqtype_1",
            parent=Parent(
                id="parent",
                sequence_type="seqtype_2",
                location=SingleInterval(5, 9, Strand.MINUS, parent="parent"),
            ),
        )
        # Sequence data
        assert sequence.sequence == Seq("ACTG")
        assert str(sequence) == "ACTG"
        # Alphabet
        assert sequence.alphabet == Alphabet.NT_STRICT
        # ID
        assert sequence.id == "id"
        # Sequence type
        assert sequence.sequence_type == "seqtype_1"
        # Parent ID
        assert sequence.parent_id == "parent"
        assert Sequence("A", Alphabet.NT_STRICT, parent="parent").parent_id == "parent"
        assert (
            Sequence(
                "A",
                Alphabet.NT_STRICT,
                parent=Parent(location=SingleInterval(5, 6, Strand.MINUS, parent="parent")),
            ).parent_id
            == "parent"
        )
        assert Sequence("A", Alphabet.NT_STRICT).parent_id is None
        # Parent type
        assert sequence.parent_type == "seqtype_2"
        # Parent strand
        assert sequence.parent_strand == Strand.MINUS
        assert Sequence("A", Alphabet.NT_STRICT, parent=Strand.UNSTRANDED).parent_strand == Strand.UNSTRANDED
        # Location on parent
        assert (
            Sequence(
                "A",
                Alphabet.NT_STRICT,
                parent=SingleInterval(3, 4, Strand.UNSTRANDED),
            ).parent_strand
            == Strand.UNSTRANDED
        )
        assert Sequence("A", Alphabet.NT_STRICT).parent_strand is None
        assert sequence.location_on_parent == SingleInterval(5, 9, Strand.MINUS, parent="parent")
        # No alphabet validation
        Sequence("xxx", Alphabet.NT_STRICT, validate_alphabet=False)

    @pytest.mark.parametrize(
        "data,alphabet,parent_id,parent_type,parent_strand,location_on_parent,expected_exception",
        [
            ("A-C", Alphabet.NT_STRICT, None, None, None, None, AlphabetError),
            (
                "ACG",
                Alphabet.NT_STRICT,
                None,
                None,
                None,
                SingleInterval(0, 4, Strand.PLUS),
                ValueError,
            ),
            (
                "ATT",
                Alphabet.NT_STRICT,
                "parent1",
                None,
                None,
                SingleInterval(0, 3, Strand.PLUS, parent="parent2"),
                ValueError,
            ),
            (
                "GGG",
                Alphabet.NT_STRICT,
                None,
                None,
                Strand.MINUS,
                SingleInterval(0, 3, Strand.PLUS),
                ValueError,
            ),
            (
                "GGG",
                Alphabet.NT_STRICT,
                None,
                "seqtype_2",
                None,
                SingleInterval(
                    0,
                    3,
                    Strand.PLUS,
                    parent=Parent(sequence_type="seqtype_3"),
                ),
                ValueError,
            ),
        ],
    )
    def test_init_invalid_params(
        self,
        data,
        alphabet,
        parent_id,
        parent_type,
        parent_strand,
        location_on_parent,
        expected_exception,
    ):
        with pytest.raises(expected_exception):
            Sequence(
                data,
                alphabet,
                parent=Parent(
                    id=parent_id,
                    sequence_type=parent_type,
                    strand=parent_strand,
                    location=location_on_parent,
                ),
            )

    @pytest.mark.parametrize(
        "sequence,other,expected",
        [
            (
                Sequence("AAAA", Alphabet.NT_STRICT, validate_alphabet=False),
                "AAAA",
                False,
            ),
            (
                Sequence("AAAA", Alphabet.NT_STRICT, validate_alphabet=False),
                Sequence("AAAA", Alphabet.NT_STRICT, validate_alphabet=True),
                True,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                True,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                Sequence(
                    "AAAa",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                False,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                Sequence(
                    "AAAA",
                    Alphabet.NT_EXTENDED,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                False,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq2",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                False,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype_1",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype_2",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                False,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.PLUS, None),
                    validate_alphabet=False,
                ),
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype",
                    parent=SingleInterval(0, 4, Strand.UNSTRANDED, None),
                    validate_alphabet=False,
                ),
                False,
            ),
            (
                Sequence("AAAA", Alphabet.NT_STRICT, parent="parent1"),
                Sequence("AAAA", Alphabet.NT_STRICT, parent="parent2"),
                False,
            ),
            (
                Sequence("AAAA", Alphabet.NT_STRICT, parent="parent"),
                Sequence("AAAA", Alphabet.NT_STRICT),
                False,
            ),
            (
                Sequence("AAAA", Alphabet.NT_STRICT, parent=Strand.UNSTRANDED),
                Sequence("AAAA", Alphabet.NT_STRICT, parent=Strand.PLUS),
                False,
            ),
            (
                Sequence("AAAA", Alphabet.NT_STRICT, parent=Strand.UNSTRANDED),
                Sequence("AAAA", Alphabet.NT_STRICT),
                False,
            ),
            (
                Sequence(
                    "AAAA",
                    Alphabet.NT_STRICT,
                    parent="seqtype",
                ),
                Sequence("AAAA", Alphabet.NT_STRICT),
                False,
            ),
        ],
    )
    def test_equals(self, sequence, other, expected):
        assert (sequence == other) is expected
        assert (other == sequence) is expected

    def test_str(self):
        assert str(Sequence("AAAAt", Alphabet.NT_EXTENDED_GAPPED)) == "AAAAt"

    def test_len(self):
        assert len(Sequence("", Alphabet.NT_EXTENDED_GAPPED)) == 0
        assert len(Sequence("AAAAt", Alphabet.NT_EXTENDED_GAPPED)) == 5

    @pytest.mark.parametrize(
        "seq,key,exp",
        [
            # No parent
            (Sequence("acgtacgt", Alphabet.NT_STRICT), 3, Sequence("t", Alphabet.NT_STRICT)),
            (Sequence("acgtacgt", Alphabet.NT_STRICT), slice(3, 6), Sequence("tac", Alphabet.NT_STRICT)),
            (Sequence("acgtacgt", Alphabet.NT_STRICT), slice(3, 10), Sequence("tacgt", Alphabet.NT_STRICT)),
            # Parent with location; slice
            (
                Sequence("actgactg", Alphabet.NT_STRICT, parent=SingleInterval(0, 8, Strand.PLUS)),
                slice(3, 6),
                Sequence("gac", Alphabet.NT_STRICT, parent=SingleInterval(3, 6, Strand.PLUS)),
            ),
            (
                Sequence("actgactg", Alphabet.NT_STRICT, parent=SingleInterval(0, 8, Strand.MINUS)),
                slice(3, 6),
                Sequence("gac", Alphabet.NT_STRICT, parent=SingleInterval(2, 5, Strand.MINUS)),
            ),
            # Parent with location; single position
            (
                Sequence("actgactg", Alphabet.NT_STRICT, parent=SingleInterval(0, 8, Strand.PLUS)),
                3,
                Sequence("g", Alphabet.NT_STRICT, parent=SingleInterval(3, 4, Strand.PLUS)),
            ),
            (
                Sequence("actgactg", Alphabet.NT_STRICT, parent=SingleInterval(0, 8, Strand.MINUS)),
                3,
                Sequence("g", Alphabet.NT_STRICT, parent=SingleInterval(4, 5, Strand.MINUS)),
            ),
            # Parent without full location
            (
                Sequence("actgactg", Alphabet.NT_STRICT, parent="parent"),
                slice(3, 6),
                Sequence("gac", Alphabet.NT_STRICT, parent="parent"),
            ),
            (
                Sequence("actgactg", Alphabet.NT_STRICT, parent=Strand.UNSTRANDED),
                slice(3, 6),
                Sequence("gac", Alphabet.NT_STRICT, parent=Strand.UNSTRANDED),
            ),
        ],
    )
    def test_getitem(self, seq, key, exp):
        assert seq[key] == exp

    def test_getitem_error(self):
        with pytest.raises(InvalidStrandException):
            Sequence("actgactg", Alphabet.NT_STRICT, parent=SingleInterval(0, 8, Strand.UNSTRANDED))[3:6]

    @pytest.mark.parametrize(
        "sequence,alphabet,validate_alphabet",
        [
            ("", Alphabet.NT_STRICT, True),
            ("acgtACGT", Alphabet.NT_STRICT, True),
            ("N", Alphabet.NT_STRICT, False),
            ("acNNNw", Alphabet.NT_EXTENDED, True),
            ("AN-", Alphabet.NT_EXTENDED, False),
            ("GG--AAA", Alphabet.NT_STRICT_GAPPED, True),
            ("AN-", Alphabet.NT_STRICT_GAPPED, False),
            ("nnAAw-cg", Alphabet.NT_EXTENDED_GAPPED, True),
            ("xxx", Alphabet.NT_EXTENDED_GAPPED, False),
            ("MWT*", Alphabet.AA, True),
            ("T*-", Alphabet.AA, False),
            ("ABCDE-", Alphabet.GENERIC, True),
            ("*", Alphabet.GENERIC, False),
        ],
    )
    def test_validate_alphabet(self, sequence, alphabet, validate_alphabet):
        Sequence(sequence, alphabet, validate_alphabet=validate_alphabet)

    @pytest.mark.parametrize(
        "sequence,alphabet",
        [
            ("N", Alphabet.NT_STRICT),
            ("A-", Alphabet.NT_EXTENDED),
            ("AN-", Alphabet.NT_STRICT_GAPPED),
            ("E", Alphabet.NT_EXTENDED_GAPPED),
            ("R-", Alphabet.AA),
            ("?", Alphabet.GENERIC),
        ],
    )
    def test_validate_alphabet_error(self, sequence, alphabet):
        with pytest.raises(AlphabetError):
            Sequence(sequence, alphabet, validate_alphabet=True)

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            (Sequence("A", Alphabet.NT_STRICT, parent="parent"), "parent"),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(id="parent", location=SingleInterval(10, 11, Strand.UNSTRANDED)),
                ),
                "parent",
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        location=SingleInterval(10, 11, Strand.UNSTRANDED, parent="parent"),
                    ),
                ),
                "parent",
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(location=SingleInterval(10, 11, Strand.UNSTRANDED, parent="parent")),
                ),
                "parent",
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(10, 11, Strand.UNSTRANDED),
                ),
                None,
            ),
            (Sequence("A", Alphabet.NT_STRICT), None),
        ],
    )
    def test_parent_id(self, sequence, expected):
        assert sequence.parent_id == expected

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            (Sequence("A", Alphabet.NT_STRICT), None),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(sequence_type="seqtype"),
                ),
                "seqtype",
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        location=SingleInterval(
                            0,
                            1,
                            Strand.PLUS,
                            parent=Parent(sequence_type="seqtype"),
                        )
                    ),
                ),
                "seqtype",
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        location=SingleInterval(
                            0,
                            1,
                            Strand.PLUS,
                            parent=Sequence("AA", Alphabet.NT_STRICT, type="seqtype"),
                        )
                    ),
                ),
                "seqtype",
            ),
        ],
    )
    def test_parent_type(self, sequence, expected):
        assert sequence.parent_type == expected

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            (Sequence("A", Alphabet.NT_STRICT), None),
            (
                Sequence("A", Alphabet.NT_STRICT, parent=Strand.MINUS),
                Strand.MINUS,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(10, 11, Strand.MINUS),
                ),
                Strand.MINUS,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(10, 11, Strand.MINUS),
                ),
                Strand.MINUS,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        strand=Strand.MINUS,
                        location=SingleInterval(10, 11, Strand.MINUS),
                    ),
                ),
                Strand.MINUS,
            ),
        ],
    )
    def test_parent_strand(self, sequence, expected):
        assert sequence.parent_strand == expected

    @pytest.mark.parametrize(
        "sequence,new_id,new_type,expected",
        [
            (
                Sequence("", Alphabet.NT_STRICT),
                None,
                None,
                Sequence("", Alphabet.NT_STRICT),
            ),
            (
                Sequence("ACGtacgT", Alphabet.NT_STRICT),
                None,
                None,
                Sequence("AcgtaCGT", Alphabet.NT_STRICT),
            ),
            (
                Sequence("ATUGCYRSWKMBdhvnNVHDbmkwsrycguta", Alphabet.NT_EXTENDED),
                None,
                None,
                Sequence("taacgryswmkvHDBNnbdhVKMWSYRGCAAT", Alphabet.NT_EXTENDED),
            ),
            (
                Sequence("--A-CGta", Alphabet.NT_STRICT_GAPPED),
                None,
                None,
                Sequence("taCG-T--", Alphabet.NT_STRICT_GAPPED),
            ),
            (
                Sequence("AtUC-N-", Alphabet.NT_EXTENDED_GAPPED),
                None,
                None,
                Sequence("-N-GAaT", Alphabet.NT_EXTENDED_GAPPED),
            ),
            (
                Sequence("ACGta", Alphabet.NT_STRICT),
                "new_id",
                "seqtype",
                Sequence(
                    "taCGT",
                    Alphabet.NT_STRICT,
                    id="new_id",
                    type="seqtype",
                ),
            ),
            (
                Sequence("ACGta", Alphabet.NT_STRICT, parent=Strand.PLUS),
                None,
                None,
                Sequence("taCGT", Alphabet.NT_STRICT, parent=Strand.MINUS),
            ),
            (
                Sequence(
                    "ACGta",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(5, 10, Strand.PLUS),
                ),
                None,
                None,
                Sequence(
                    "taCGT",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(5, 10, Strand.MINUS),
                ),
            ),
        ],
    )
    def test_reverse_complement(self, sequence, new_id, new_type, expected):
        assert sequence.reverse_complement(new_id=new_id, new_type=new_type) == expected

    @pytest.mark.parametrize(
        "sequence",
        [
            Sequence("AAA", Alphabet.AA),
            Sequence("AAA", Alphabet.GENERIC),
            Sequence("xxx", Alphabet.NT_STRICT, validate_alphabet=False),
        ],
    )
    def test_reverse_complement_error(self, sequence):
        with pytest.raises(AlphabetError):
            sequence.reverse_complement()

    @pytest.mark.parametrize(
        "seq1,seq2,new_id,data_only,expected",
        [
            (
                Sequence("", Alphabet.NT_STRICT, parent="parent1"),
                Sequence("", Alphabet.NT_STRICT, parent="parent2"),
                "new_id",
                True,
                Sequence("", Alphabet.NT_STRICT, id="new_id"),
            ),
            (
                Sequence("AA", Alphabet.NT_STRICT, parent="parent1"),
                Sequence("TT", Alphabet.NT_STRICT, parent="parent2"),
                "new_id",
                True,
                Sequence("AATT", Alphabet.NT_STRICT, id="new_id"),
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    type="seqtype_1",
                    parent=Parent(
                        id="parent1",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 7, Strand.PLUS),
                    ),
                ),
                Sequence(
                    "TT",
                    Alphabet.NT_STRICT,
                    type="seqtype_2",
                    parent=Parent(
                        id="parent1",
                        strand=Strand.MINUS,
                        location=SingleInterval(20, 22, Strand.MINUS),
                    ),
                ),
                None,
                True,
                Sequence("AATT", Alphabet.NT_STRICT),
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(
                        id="parent1",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 7, Strand.PLUS),
                    ),
                ),
                Sequence(
                    "TT",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(id="parent1", strand=Strand.PLUS),
                ),
                None,
                False,
                Sequence(
                    "AATT",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent="parent1",
                ),
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(id="parent1", strand=Strand.PLUS),
                ),
                Sequence(
                    "TT",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(
                        id="parent1",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 7, Strand.PLUS),
                    ),
                ),
                None,
                False,
                Sequence(
                    "AATT",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent="parent1",
                ),
            ),
            (
                Sequence(
                    "CC",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(id="parent", location=SingleInterval(3, 5, Strand.PLUS)),
                ),
                Sequence(
                    "TT",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(id="parent", location=SingleInterval(10, 12, Strand.PLUS)),
                ),
                "new_id",
                True,
                Sequence("CCTT", Alphabet.NT_STRICT, id="new_id"),
            ),
            (
                Sequence(
                    "CC",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(id="parent", location=SingleInterval(3, 5, Strand.PLUS)),
                ),
                Sequence(
                    "TT",
                    Alphabet.NT_STRICT,
                    type="seqtype",
                    parent=Parent(id="parent", location=SingleInterval(0, 2, Strand.PLUS)),
                ),
                "new_id",
                True,
                Sequence("CCTT", Alphabet.NT_STRICT, id="new_id"),
            ),
            (
                Sequence("AA", Alphabet.NT_STRICT, id="seq1", parent="parent"),
                Sequence("", Alphabet.NT_STRICT, id="seq2", parent="parent"),
                None,
                False,
                Sequence("AA", Alphabet.NT_STRICT, parent="parent"),
            ),
            (
                Sequence(
                    "ACT",
                    Alphabet.NT_STRICT,
                    parent=Parent(id="parent", location=SingleInterval(2, 5, Strand.PLUS)),
                ),
                Sequence(
                    "GGA",
                    Alphabet.NT_STRICT,
                    parent=Parent(id="parent", location=SingleInterval(8, 11, Strand.PLUS)),
                ),
                "new_id",
                False,
                Sequence(
                    "ACTGGA",
                    Alphabet.NT_STRICT,
                    id="new_id",
                    parent=Parent(
                        id="parent",
                        location=CompoundInterval.from_single_intervals(
                            [
                                SingleInterval(2, 5, Strand.PLUS),
                                SingleInterval(8, 11, Strand.PLUS),
                            ]
                        ),
                    ),
                ),
            ),
            (
                Sequence(
                    "ACT",
                    Alphabet.NT_STRICT,
                    parent=Parent(id="parent", location=SingleInterval(8, 11, Strand.MINUS)),
                ),
                Sequence(
                    "GGA",
                    Alphabet.NT_STRICT,
                    parent=Parent(id="parent", location=SingleInterval(2, 5, Strand.MINUS)),
                ),
                "new_id",
                False,
                Sequence(
                    "ACTGGA",
                    Alphabet.NT_STRICT,
                    id="new_id",
                    parent=Parent(
                        id="parent",
                        location=CompoundInterval.from_single_intervals(
                            [
                                SingleInterval(2, 5, Strand.MINUS),
                                SingleInterval(8, 11, Strand.MINUS),
                            ]
                        ),
                    ),
                ),
            ),
            (
                Sequence(
                    "ACT",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype_1",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        strand=Strand.PLUS,
                        location=SingleInterval(5, 8, Strand.PLUS),
                    ),
                ),
                Sequence(
                    "GCG",
                    Alphabet.NT_STRICT,
                    id="seq1",
                    type="seqtype_1",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        strand=Strand.PLUS,
                        location=SingleInterval(15, 18, Strand.PLUS),
                    ),
                ),
                None,
                False,
                Sequence(
                    "ACTGCG",
                    Alphabet.NT_STRICT,
                    type="seqtype_1",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        strand=Strand.PLUS,
                        location=CompoundInterval.from_single_intervals(
                            [
                                SingleInterval(5, 8, Strand.PLUS),
                                SingleInterval(15, 18, Strand.PLUS),
                            ]
                        ),
                    ),
                ),
            ),
        ],
    )
    def test_append(self, seq1, seq2, new_id, data_only, expected):
        assert seq1.append(seq2, new_id, data_only) == expected

    @pytest.mark.parametrize(
        "seq1,seq2,new_id,data_only",
        [
            (
                Sequence("AA", Alphabet.NT_STRICT),
                Sequence("TT", Alphabet.NT_EXTENDED),
                None,
                True,
            ),
            (
                Sequence("AA", Alphabet.NT_STRICT, type="seqtype_1"),
                Sequence("AA", Alphabet.NT_STRICT, type="seqtype_2"),
                None,
                False,
            ),
            (
                Sequence("AA", Alphabet.NT_STRICT, parent="parent1"),
                Sequence("AA", Alphabet.NT_STRICT, parent="parent2"),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(10, 12, Strand.PLUS),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(5, 7, Strand.PLUS),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(2, 4, Strand.MINUS),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(5, 7, Strand.MINUS),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(10, 12, Strand.PLUS),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(11, 13, Strand.PLUS),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(10, 12, Strand.UNSTRANDED),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=SingleInterval(15, 17, Strand.UNSTRANDED),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent1",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent2",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_1",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype",
                        strand=Strand.MINUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAAT", Alphabet.NT_STRICT),
                        parent="grandparent",
                    ),
                ),
                None,
                False,
            ),
            (
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent1",
                    ),
                ),
                Sequence(
                    "AA",
                    Alphabet.NT_STRICT,
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype",
                        strand=Strand.PLUS,
                        sequence=Sequence("AAA", Alphabet.NT_STRICT),
                        parent="grandparent2",
                    ),
                ),
                None,
                False,
            ),
        ],
    )
    def test_append_error(self, seq1, seq2, new_id, data_only):
        with pytest.raises(ValueError):
            seq1.append(seq2, new_id, data_only)

    @pytest.mark.parametrize(
        "sequence,sequence_type,include_self,expected",
        [
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                True,
                Parent(
                    sequence=Sequence(
                        "A",
                        Alphabet.NT_STRICT,
                        id="self",
                        type="seqtype",
                        parent=Parent(id="parent", sequence_type="seqtype"),
                    )
                ),
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                False,
                Parent(id="parent", sequence_type="seqtype"),
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype_1",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                    ),
                ),
                "seqtype_2",
                True,
                Parent(
                    id="parent",
                    sequence_type="seqtype_2",
                    parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                ),
            ),
        ],
    )
    def test_first_ancestor_of_type(self, sequence, sequence_type, include_self, expected):
        assert sequence.first_ancestor_of_type(sequence_type, include_self=include_self) == expected

    @pytest.mark.parametrize(
        "sequence,sequence_type,include_self",
        [
            (
                Sequence("A", Alphabet.NT_STRICT, id="self"),
                "seqtype_1",
                True,
            ),
            (
                Sequence("A", Alphabet.NT_STRICT, id="self", parent="parent"),
                "seqtype_1",
                True,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype_2",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_1",
                        parent=Parent(id="grandparent", sequence_type="seqtype_1"),
                    ),
                ),
                "seqtype_3",
                True,
            ),
        ],
    )
    def test_first_ancestor_of_type_error(self, sequence, sequence_type, include_self):
        with pytest.raises(NoSuchAncestorException):
            sequence.first_ancestor_of_type(sequence_type, include_self=include_self)

    @pytest.mark.parametrize(
        "sequence,sequence_type,include_self,expected",
        [
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                True,
                True,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype",
                    parent=Parent(id="parent", sequence_type="seqtype"),
                ),
                "seqtype",
                False,
                True,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype_1",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_2",
                        parent=Parent(id="grandparent", sequence_type="seqtype_2"),
                    ),
                ),
                "seqtype_2",
                True,
                True,
            ),
            (
                Sequence("A", Alphabet.NT_STRICT, id="self"),
                "seqtype_1",
                True,
                False,
            ),
            (
                Sequence("A", Alphabet.NT_STRICT, id="self", parent="parent"),
                "seqtype_1",
                True,
                False,
            ),
            (
                Sequence(
                    "A",
                    Alphabet.NT_STRICT,
                    id="self",
                    type="seqtype_2",
                    parent=Parent(
                        id="parent",
                        sequence_type="seqtype_1",
                        parent=Parent(id="grandparent", sequence_type="seqtype_1"),
                    ),
                ),
                "seqtype_3",
                True,
                False,
            ),
        ],
    )
    def test_has_ancestor_of_type(self, sequence, sequence_type, include_self, expected):
        assert sequence.has_ancestor_of_type(sequence_type, include_self=include_self) is expected

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            (Sequence("ATGCATATTTGGAAACCAA", Alphabet.NT_STRICT, id="test"), ">test\nATGCATATTT\nGGAAACCAA"),
            (Sequence("ATGCATATTTGGAAACCAA", Alphabet.NT_STRICT), ">None\nATGCATATTT\nGGAAACCAA"),
            (Sequence("GGAAACCAA", Alphabet.NT_STRICT, id="test"), ">test\nGGAAACCAA"),
            (
                Sequence("ATGCATATTTGGAAACCAAGGAAACCAA", Alphabet.NT_STRICT, id="test"),
                ">test\nATGCATATTT\nGGAAACCAAG\nGAAACCAA",
            ),
            (
                Sequence(
                    data="AAAAAAA",
                    alphabet=Alphabet.NT_STRICT,
                    id="test",
                    parent=Parent(location=SingleInterval(33, 40, Strand.MINUS)),
                ),
                ">test\nAAAAAAA",
            ),
        ],
    )
    def test_to_fasta(self, sequence, expected):
        s_fa = sequence.to_fasta(num_chars=10)
        assert s_fa == expected

    def test_empty_to_fasta(self):
        s = Sequence("", Alphabet.NT_STRICT)
        with pytest.raises(EmptySequenceFastaError):
            s.to_fasta()
