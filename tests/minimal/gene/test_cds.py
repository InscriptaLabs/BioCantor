import pytest

from inscripta.biocantor.exc import NoSuchAncestorException, MismatchedFrameException
from inscripta.biocantor.gene.cds import CDSInterval, TranslationTable
from inscripta.biocantor.gene.cds_frame import CDSPhase, CDSFrame
from inscripta.biocantor.gene.codon import Codon
from inscripta.biocantor.location.location_impl import CompoundInterval, SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import SequenceType
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence import Sequence
from inscripta.biocantor.sequence.alphabet import Alphabet


class TestCDSPhase:
    def test_from_int(self):
        assert CDSPhase.from_int(2) == CDSPhase.TWO
        with pytest.raises(ValueError):
            CDSPhase.from_int(3)

    @pytest.mark.parametrize(
        "phase,expected", [(CDSPhase.ONE, "1"), (CDSPhase.TWO, "2"), (CDSPhase.ZERO, "0"), (CDSPhase.NONE, ".")]
    )
    def test_to_gff(self, phase, expected):
        assert phase.to_gff() == expected


class TestCDSFrame:
    def test_from_int(self):
        assert CDSFrame.from_int(2) == CDSFrame.TWO
        with pytest.raises(ValueError):
            CDSFrame.from_int(3)

    @pytest.mark.parametrize(
        "frame,shift,expected",
        [
            (CDSFrame.ZERO, 0, CDSFrame.ZERO),
            (CDSFrame.ZERO, -1, CDSFrame.TWO),
            (CDSFrame.ZERO, 3, CDSFrame.ZERO),
            (CDSFrame.ONE, 1, CDSFrame.TWO),
            (CDSFrame.TWO, 5, CDSFrame.ONE),
            (CDSFrame.NONE, 1, CDSFrame.NONE),
        ],
    )
    def test_shift(self, frame, shift, expected):
        assert frame.shift(shift) == expected


class TestCDSInterval:

    alphabet = Alphabet.NT_STRICT

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # 2bp CDS
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 2, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                0,
            ),
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                3,
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                3,
            ),
            # Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                5,
            ),
        ],
    )
    def test_num_codons(self, cds, expected):
        assert cds.num_codons == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                [Codon.ATA, Codon.CGA, Codon.TCA],
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
            # Discontiguous CDS, plus strand, frame=1, 1bp deletion at start of exon 2
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon.GGA, Codon.CCC],  # GP
            ),
            # Discontiguous CDS, plus strand, frame=1,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO, CDSFrame.TWO],  # QGP
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
            # Discontiguous CDS, minus strand, frame=2
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
            # Discontiguous CDS, minus strand, frame=2, with frameshift that leads to truncation
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA],  # QG
            ),
            # Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),  # G gets repeated here
                [Codon.AGG, Codon.AAA, Codon.GGT, Codon.CCC, Codon.TGA],
            ),
        ],
    )
    def test_scan_codons(self, cds, expected):
        assert list(cds.scan_codons()) == expected
        # run again to ensure caching is not a problem
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("atacgatca", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                [Codon.ATA, Codon.CGA, Codon.TCA],
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("aaacaaaagggacccaaaaaa", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon.CAG, Codon.GGA, Codon.CCC],  # QGP
            ),
        ],
    )
    def test_scan_codons_lower_case(self, cds, expected):
        assert list(cds.scan_codons()) == expected
        # run again to ensure caching is not a problem
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds",
        [
            CDSInterval.from_location(
                SingleInterval(
                    0,
                    9,
                    Strand.PLUS,
                    parent=Sequence("ANACGATCA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                ),
                [CDSFrame.ZERO],
            ),
            CDSInterval.from_location(
                CompoundInterval(
                    [2, 8],
                    [5, 17],
                    Strand.PLUS,
                    parent=Sequence(
                        "AANNNNAAGGGTACCCAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                    ),
                ),
                [CDSFrame.ONE, CDSFrame.TWO],
            ),
        ],
    )
    def test_scan_codons_exception(self, cds):
        with pytest.raises(ValueError):
            _ = list(cds.scan_codons())

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # 2bp CDS
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 2, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                [],
            ),
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                [
                    SingleInterval(
                        0, 3, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),  # ATA
                    SingleInterval(
                        3, 6, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),  # CGA
                    SingleInterval(
                        6, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),  # TCA
                ],
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [
                    CompoundInterval(
                        [3, 8],
                        [5, 9],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CAG
                    SingleInterval(
                        9,
                        12,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                    SingleInterval(
                        12,
                        15,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CCC
                ],
            ),
            # Discontiguous CDS, plus strand, frame=1, 1bp deletion at start of exon 2
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [
                    SingleInterval(
                        8,
                        11,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                    SingleInterval(
                        11,
                        14,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CCC
                ],
            ),
            # Discontiguous CDS, plus strand, frame=0,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.ZERO],
                ),
                [
                    SingleInterval(
                        2,
                        5,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    SingleInterval(
                        8,
                        11,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    SingleInterval(
                        12,
                        15,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    SingleInterval(
                        15,
                        18,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                ],
            ),
            # Discontiguous CDS, plus strand, frame=1,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO, CDSFrame.TWO],
                ),
                [
                    CompoundInterval(
                        [3, 8],
                        [5, 9],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CAG
                    CompoundInterval(
                        [9, 12],
                        [11, 13],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                    SingleInterval(
                        13,
                        16,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CCC
                ],
            ),
            # Discontiguous CDS, plus strand, frame=2,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.TWO, CDSFrame.ONE, CDSFrame.ONE],
                ),
                [
                    CompoundInterval(
                        [4, 8],
                        [5, 10],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    CompoundInterval(
                        [10, 12],
                        [11, 14],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    SingleInterval(
                        14,
                        17,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                ],
            ),
            # Discontiguous CDS, minus strand, frame=2
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [
                    SingleInterval(
                        12,
                        15,
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CAG
                    SingleInterval(
                        9,
                        12,
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                    CompoundInterval(
                        [3, 8],
                        [5, 9],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CCC
                ],
            ),
            # Discontiguous CDS, minus strand, frame=2, with frameshift that leads to truncation
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                [
                    SingleInterval(
                        12,
                        15,
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CAG
                    SingleInterval(
                        9,
                        12,
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                ],
            ),
            # Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                [
                    SingleInterval(
                        2,
                        5,
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # AGG
                    SingleInterval(
                        5,
                        8,
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # AAA
                    CompoundInterval(
                        [8, 8],
                        [9, 10],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGT, G gets repeated
                    SingleInterval(
                        10,
                        13,
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CCC
                    SingleInterval(
                        13,
                        16,
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # TGA
                ],
            ),
            # Discontiguous CDS, plus strand, with +1 programmed frameshift that skips over a 1nt exon
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 6, 8],
                        [5, 7, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.ZERO],
                ),
                [
                    SingleInterval(
                        2,
                        5,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                    SingleInterval(
                        8,
                        11,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # GGA
                    SingleInterval(
                        11,
                        14,
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),  # CCC
                ],
            ),
            # 1bp exon on the negative strand gets removed due to being a partial codon
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [0, 7],
                        [5, 8],
                        Strand.MINUS,
                        parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.TWO],
                ),
                [
                    SingleInterval(
                        2, 5, Strand.MINUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                ],  # CGT
            ),
        ],
    )
    def test_scan_codon_locations(self, cds, expected):
        assert list(cds.scan_codon_locations()) == expected
        # run again to prove caching doesn't kill the iterator
        assert list(cds.scan_codon_locations()) == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # Contiguous CDS, plus strand, frame=0
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0,
                        9,
                        Strand.PLUS,
                        parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO],
                ),
                "IRS",
            ),
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                "QGP",
            ),
            # Discontiguous CDS, minus strand, frame=2, codons don't reach end of CDS
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                "QG",
            ),
        ],
    )
    def test_translate(self, cds, expected):
        assert str(cds.translate()) == expected
        # run again to ensure caching is not a problem
        assert str(cds.translate()) == expected

    def test_accessors(self):
        cds = CDSInterval.from_location(SingleInterval(0, 10, Strand.PLUS), [CDSFrame.ZERO])
        assert cds.start == cds._location.start
        assert cds.end == cds._location.end
        assert cds.strand == cds._location.strand

    @pytest.mark.parametrize(
        "cds,expected",
        [
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATGA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("atacgatca", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
        ],
    )
    def test_has_valid_stop(self, cds, expected):
        assert cds.has_valid_stop == expected
        # run again to ensure caching is not a problem
        assert cds.has_valid_stop == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATGTGAAAACCC", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("atacgatca", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
        ],
    )
    def test_in_frame_stop(self, cds, expected):
        assert cds.has_in_frame_stop == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATGTGAAAACCC", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.MINUS, parent=Sequence("ATGTGCCATCC", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.MINUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
        ],
    )
    def test_has_canonical_start_codon(self, cds, expected):
        assert cds.has_canonical_start_codon == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("TTGTGAAAACCC", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.MINUS, parent=Sequence("ATGTGCCAG", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                True,
            ),
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.MINUS, parent=Sequence("ATACGAGAT", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                False,
            ),
        ],
    )
    def test_has_start_codon_in_specific_translation_table(self, cds, expected):
        assert cds.has_start_codon_in_specific_translation_table(TranslationTable.STANDARD) == expected

    @pytest.mark.parametrize(
        "location,starting_frame,expected",
        [
            # 1bp exon on the negative strand gets removed due to being a partial codon
            (
                CompoundInterval(
                    [0, 7], [5, 8], Strand.MINUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                ),
                CDSFrame.TWO,
                [CDSFrame.TWO, CDSFrame.TWO],
            ),
            # 3 exons on plus strand with starting frame of 2
            (
                CompoundInterval(
                    [2, 8, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.TWO,
                [CDSFrame.TWO, CDSFrame.ONE, CDSFrame.ONE],
            ),
            # 3 exons on plus strand with starting frame of 0
            (
                CompoundInterval(
                    [2, 8, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ZERO,
                [CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.ZERO],
            ),
            # 3 exons on plus strand with starting frame of 0
            (
                CompoundInterval(
                    [0, 8, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ZERO,
                [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.TWO],
            ),
            # 3 exons on plus strand with starting frame of 2
            (
                CompoundInterval(
                    [0, 8, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.TWO,
                [CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
            ),
            # 3 exons on minus strand with starting frame of 0
            (
                CompoundInterval(
                    [2, 8, 12],
                    [5, 11, 18],
                    Strand.MINUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ZERO,
                [CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.ZERO],
            ),
            # 3 exons on minus strand with starting frame of 1
            (
                CompoundInterval(
                    [2, 8, 12],
                    [5, 11, 18],
                    Strand.MINUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ONE,
                [CDSFrame.TWO, CDSFrame.TWO, CDSFrame.ONE],
            ),
            # 3 exons on minus strand with starting frame of 2
            (
                CompoundInterval(
                    [2, 8, 12],
                    [5, 11, 18],
                    Strand.MINUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.TWO,
                [CDSFrame.ONE, CDSFrame.ONE, CDSFrame.TWO],
            ),
            (
                CompoundInterval(
                    [0, 7, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ZERO,
                [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO],
            ),
            (
                CompoundInterval(
                    [0, 7, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ONE,
                [CDSFrame.ONE, CDSFrame.ONE, CDSFrame.TWO],
            ),
            (
                CompoundInterval(
                    [0, 7, 12],
                    [5, 11, 18],
                    Strand.PLUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.TWO,
                [CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ONE],
            ),
            (
                CompoundInterval(
                    [0, 7, 12],
                    [5, 11, 18],
                    Strand.MINUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ZERO,
                [CDSFrame.ONE, CDSFrame.ZERO, CDSFrame.ZERO],
            ),
            (
                CompoundInterval(
                    [0, 7, 12],
                    [5, 11, 18],
                    Strand.MINUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.ONE,
                [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ONE],
            ),
            (
                CompoundInterval(
                    [0, 7, 12],
                    [5, 11, 18],
                    Strand.MINUS,
                    parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                ),
                CDSFrame.TWO,
                [CDSFrame.TWO, CDSFrame.ONE, CDSFrame.TWO],
            ),
        ],
    )
    def test_construct_frames_from_location(self, location, starting_frame, expected):
        frames = CDSInterval.construct_frames_from_location(location, starting_frame)
        assert frames == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # 1bp exon on the negative strand gets removed due to being a partial codon
            # this shifts the frame all the
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [0, 7],
                        [5, 8],
                        Strand.MINUS,
                        parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.TWO],
                ),
                [Codon.TAT],
            ),
            # Discontiguous CDS, plus strand, with +1 programmed frameshift that skips over a 1nt exon
            # last frame gets shifted as a result
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 6, 8],
                        [5, 7, 16],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.ZERO],
                ),
                [Codon.ACA, Codon.AGG, Codon.ACC, Codon.CAA],
            ),
            # Discontiguous CDS, plus strand, frame=2,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want(ed) to maintain original frame
            # frame is not changed
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 11, 18],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGTACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.TWO, CDSFrame.ONE, CDSFrame.ONE],
                ),
                [Codon.AGG, Codon.GAC, Codon.CCA],
            ),
        ],
    )
    def test_optimize_blocks(self, cds, expected):
        assert list(cds.optimize_blocks().scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # this overlapping interval gets merged into [0, 17], with frame 1, so we end up with 5 codons
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                [Codon.AGG, Codon.AAA, Codon.GTC, Codon.CCT, Codon.GAA],
            ),
            # this overlapping interval gets merged into [2,5], [8,18] with 0bp offset
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 13, 18],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon.AGG, Codon.GTC, Codon.CCT, Codon.GAA],
            ),
            # this overlapping interval gets merged into [2,5], [8,18] with 1bp offset
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 13, 18],
                        Strand.PLUS,
                        parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon.GGG, Codon.TCC, Codon.CTG, Codon.AAA],
            ),
        ],
    )
    def test_optimize_and_combine_blocks(self, cds, expected):
        assert list(cds.optimize_and_combine_blocks().scan_codons()) == expected

    def test_frame_exception(self):
        with pytest.raises(MismatchedFrameException):
            _ = CDSInterval.from_location(
                CompoundInterval(
                    [2, 8, 12],
                    [5, 13, 18],
                    Strand.PLUS,
                    parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                ),
                [CDSFrame.ONE, CDSFrame.ONE],
            )

    def test_frame_to_phase(self):
        cds = CDSInterval.from_location(
            CompoundInterval(
                [2, 8, 12],
                [5, 13, 18],
                Strand.PLUS,
                parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
            ),
            [CDSFrame.ONE.to_phase(), CDSFrame.ONE.to_phase(), CDSFrame.ZERO.to_phase()],
        )
        assert list(cds.scan_codons()) == [Codon.TCC, Codon.CTG, Codon.AAA]

    def test_frame_to_phase_mixed_exception(self):
        with pytest.raises(MismatchedFrameException):
            _ = CDSInterval.from_location(
                CompoundInterval(
                    [2, 8, 12],
                    [5, 13, 18],
                    Strand.PLUS,
                    parent=Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                ),
                [CDSFrame.ONE.to_phase(), CDSFrame.ONE, CDSFrame.ZERO.to_phase()],
            )

    @pytest.mark.parametrize(
        "cds",
        [
            dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ONE],
            ),
            dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ONE],
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        ],
    )
    def test_constructor(self, cds):
        cds = CDSInterval(**cds)
        cds2 = CDSInterval.from_location(cds._location, cds.frames)
        assert cds == cds2

    @pytest.mark.parametrize(
        "location,cds_frames,expected",
        [
            (
                SingleInterval(
                    5,
                    20,
                    Strand.PLUS,
                    parent=Parent(
                        id="NC_000913.3:222213-222241",
                        sequence=Sequence(
                            "AANAAATGGCGAGCACCTAACCCCCNCC",
                            Alphabet.NT_EXTENDED,
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                id="NC_000913.3",
                                location=SingleInterval(
                                    222213,
                                    222241,
                                    Strand.PLUS,
                                    parent=Parent(
                                        id="NC_000913.3",
                                        sequence_type=SequenceType.CHROMOSOME,
                                        location=SingleInterval(222213, 222241, Strand.PLUS),
                                    ),
                                ),
                                sequence_type=SequenceType.CHROMOSOME,
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ZERO],
                SingleInterval(222218, 222233, Strand.PLUS),
            )
        ],
    )
    def test_chunk_relative_constructor(self, location, cds_frames, expected):
        cds = CDSInterval.from_chunk_relative_location(location, cds_frames)
        assert cds.chromosome_location.reset_parent(None) == expected

    @pytest.mark.parametrize(
        "cds",
        [
            dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ONE],
            ),
            dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ONE],
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        ],
    )
    def test_dict(self, cds):
        cds = CDSInterval(**cds)
        cds2 = CDSInterval.from_dict(cds.to_dict())
        assert cds.to_dict() == cds2.to_dict()

    @pytest.mark.parametrize(
        "parent,expected",
        [
            # cuts 1bp off the 1st exon, shifting everything to zero
            (
                Parent(
                    id="test:3-21",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[3:],
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="test:3-21",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                3,
                                21,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                {
                    "cds_starts": (0, 5, 9),
                    "cds_ends": (2, 10, 15),
                    "strand": "PLUS",
                    "cds_frames": ["ZERO", "ONE", "ONE"],
                    "qualifiers": None,
                    "sequence_name": None,
                    "sequence_guid": None,
                    "protein_id": None,
                    "product": None,
                },
            ),
            # cuts off the 1st exon entirely
            (
                Parent(
                    id="test:3-21",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[6:],
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="test:3-21",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                6,
                                21,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                {
                    "cds_starts": (2, 6),
                    "cds_ends": (7, 12),
                    "strand": "PLUS",
                    "cds_frames": ["ONE", "ONE"],
                    "qualifiers": None,
                    "sequence_name": None,
                    "sequence_guid": None,
                    "protein_id": None,
                    "product": None,
                },
            ),
        ],
    )
    def test_dict_chunk_relative(self, parent, expected):
        cds = CDSInterval(
            **dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ONE],
            ),
            parent_or_seq_chunk_parent=parent,
        )
        assert cds.to_dict(chromosome_relative_coordinates=False) == expected

    @pytest.mark.parametrize(
        "parent",
        [
            # standard chromosome
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME)),
            # no type
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED)),
            # non-standard type
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED, type="nonstandard")),
        ],
    )
    def test_from_location(self, parent):
        _ = CDSInterval.from_location(SingleInterval(0, 9, Strand.PLUS, parent), [CDSFrame.ZERO])

    @pytest.mark.parametrize(
        "parent",
        [
            # flat chunk
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.SEQUENCE_CHUNK)),
            # proper chunk hierarchy
            Parent(
                id="test:0-9",
                sequence=Sequence(
                    "ATACGATCA",
                    alphabet,
                    id="test:0-9",
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            0,
                            9,
                            Strand.PLUS,
                            parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
        ],
    )
    def test_from_location_exception(self, parent):
        with pytest.raises(NoSuchAncestorException):
            _ = CDSInterval.from_location(SingleInterval(0, 9, Strand.PLUS, parent), [CDSFrame.ZERO])

    @pytest.mark.parametrize(
        "parent",
        [
            # proper chunk hierarchy
            Parent(
                id="test:0-9",
                sequence=Sequence(
                    "ATACGATCA",
                    alphabet,
                    id="test:0-9",
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            0,
                            9,
                            Strand.PLUS,
                            parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
            # negative strand chunk hierarchy works just fine
            Parent(
                id="test:0-9",
                location=SingleInterval(0, 9, Strand.MINUS),
                sequence=Sequence(
                    "ATACGATCA",
                    alphabet,
                    id="test:0-9",
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            0,
                            9,
                            Strand.PLUS,
                            parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
        ],
    )
    def test_from_chunk_relative_location(self, parent):
        _ = CDSInterval.from_chunk_relative_location(SingleInterval(0, 9, Strand.PLUS, parent), [CDSFrame.ZERO])

    @pytest.mark.parametrize(
        "parent",
        [
            # flat chunk
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.SEQUENCE_CHUNK)),
            # standard chromosome
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME)),
            # no type
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED)),
            # non-standard type
            Parent(sequence=Sequence("ATACGATCA", Alphabet.NT_EXTENDED_GAPPED, type="nonstandard")),
        ],
    )
    def test_from_chunk_relative_location_exception(self, parent):
        with pytest.raises(NoSuchAncestorException):
            _ = CDSInterval.from_chunk_relative_location(SingleInterval(0, 9, Strand.PLUS, parent), [CDSFrame.ZERO])

    @pytest.mark.parametrize(
        "parent,exp_frames,exp_codons",
        [
            # normal flat genome
            (
                Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                [CDSFrame.ZERO, CDSFrame.TWO],
                [Codon.AGG, Codon.AAA, Codon.GTT, Codon.CCC, Codon.TGA],
            ),
            # sequence subset that removes 1bp of exon
            (
                Parent(
                    id="test:3-21",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[3:],
                        alphabet,
                        id="test:3-21",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                3,
                                21,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                [CDSFrame.TWO, CDSFrame.TWO],
                [Codon.AAA, Codon.GTT, Codon.CCC, Codon.TGA],
            ),
            # sequence subset that removes 1st exon
            (
                Parent(
                    id="test:10-21",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[10:],
                        alphabet,
                        id="test:10-21",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                10,
                                21,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                [CDSFrame.ONE],
                # translation doesn't change because the original 1st exon is sliced off in codon iteration
                [Codon.CCT, Codon.GAA],
            ),
            # sequence subset that removes last exon and 1bp of 2nd exon
            (
                Parent(
                    id="test:0-9",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[:9],
                        alphabet,
                        id="test:0-9",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                0,
                                9,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                [CDSFrame.ZERO],
                [Codon.AGG, Codon.AAA],
            ),
        ],
    )
    def test_cds_frame_slicing_overlapping_cds(self, parent, exp_frames, exp_codons):
        """This CDS has a 1bp programmed frameshift, and so this must be maintained
        when it is sliced down.
        """
        cds = CDSInterval(
            [2, 9],
            [10, 18],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO],
            parent_or_seq_chunk_parent=parent,
        )
        assert cds.chunk_relative_frames == exp_frames
        assert list(cds.scan_codons()) == exp_codons

    @pytest.mark.parametrize(
        "parent,exp_frames,exp_codons",
        [
            # normal flat genome
            (
                Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.TWO],
                [Codon.AGA, Codon.CCC, Codon.GAA],
            ),
            # sequence subset that removes 1st exon
            (
                Parent(
                    id="test:5-21",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[5:],
                        alphabet,
                        id="test:5-21",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                5,
                                21,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                [CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.TWO],
                [Codon.CCC, Codon.GAA],
            ),
            # sequence subset that removes last block
            (
                Parent(
                    id="test:0-13",
                    sequence=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA"[:13],
                        alphabet,
                        id="test:0-13",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                0,
                                13,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO],
                [Codon.AGA, Codon.CCC],
            ),
        ],
    )
    def test_cds_frame_slicing_small_exons(self, parent, exp_frames, exp_codons):
        cds = CDSInterval(
            [2, 6, 10, 14, 16],
            [4, 7, 13, 16, 17],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.TWO],
            parent_or_seq_chunk_parent=parent,
        )
        assert cds.chunk_relative_frames == exp_frames
        assert list(cds.scan_codons()) == exp_codons

    def test_frameshift_overlapping(self):
        """This tests an edge case of overlapping CDS intervals and their codon scanning.

        This CDS should have every base of the transcript used in the translation, including the
        12th base (a C), which should be read twice.
        """
        cds = CDSInterval(
            **dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.TWO],
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        )
        assert str(cds.get_spliced_sequence()) == "AGGGTCCCCTGAAA"
        assert str(cds.extract_sequence()) == "AGGGTCCCCTGA"

    def test_equality_different_parents(self):
        cds1 = CDSInterval(
            **dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.TWO],
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        )
        cds2 = CDSInterval(
            **dict(
                cds_starts=[2, 8, 12],
                cds_ends=[5, 13, 18],
                strand=Strand.PLUS,
                frames_or_phases=[CDSFrame.ZERO, CDSFrame.ZERO, CDSFrame.TWO],
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        )
        assert cds1 != cds2  # Location equality comparison does compare sequence
        assert hash(cds1) == hash(cds2)  # Location hash comparison does not compare sequence
