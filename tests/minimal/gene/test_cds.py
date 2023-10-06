import pytest
from biocantor.exc import NoSuchAncestorException, MismatchedFrameException, InvalidPositionException
from biocantor.gene.cds import CDSInterval, TranslationTable
from biocantor.gene.cds_frame import CDSPhase, CDSFrame
from biocantor.gene.codon import Codon
from biocantor.location.location_impl import CompoundInterval, SingleInterval
from biocantor.location.strand import Strand
from biocantor.parent import SequenceType
from biocantor.parent.parent import Parent
from biocantor.sequence import Sequence
from biocantor.sequence.alphabet import Alphabet


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

    seq = "AAAGGAAAGTCCCTGAAAAAA"
    chrom_parent = Parent(
        sequence=Sequence(seq, alphabet, type=SequenceType.CHROMOSOME, id="genome2"),
        location=SingleInterval(0, len(seq), Strand.PLUS),
    )
    sequenceless_chrom_parent = Parent(
        sequence_type=SequenceType.CHROMOSOME, id="genome2", location=SingleInterval(0, len(seq), Strand.PLUS)
    )

    chunk_parent2_12 = Parent(
        id="genome2_2_12",
        sequence=Sequence(
            seq[2:12],
            Alphabet.NT_EXTENDED_GAPPED,
            type=SequenceType.SEQUENCE_CHUNK,
            parent=Parent(
                location=SingleInterval(
                    2,
                    12,
                    Strand.PLUS,
                    parent=Parent(id="genome2", sequence_type=SequenceType.CHROMOSOME),
                )
            ),
        ),
    )

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
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                5,
            ),
            # Discontiguous CDS with -1bp programmed frameshift on chunk-relative coordinates
            # chunk-relative coordinates slice off the last 2 codons
            # however, there should still be 5 codons despite the chunk slicing
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
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
            # Discontiguous CDS with -1bp programmed frameshift on chunk-relative coordinates
            # chunk-relative coordinates slice off the last 2 codons
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                3,
            ),
        ],
    )
    def test_num_chunk_relative_codons(self, cds, expected):
        assert cds.num_chunk_relative_codons == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # 0. Contiguous CDS, plus strand, frame=0
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0, 9, Strand.PLUS, parent=Sequence("ATACGATCA", alphabet, type=SequenceType.CHROMOSOME)
                    ),
                    [CDSFrame.ZERO],
                ),
                [Codon("ATA"), Codon("CGA"), Codon("TCA")],
            ),
            # 1 .Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            # Codon interval cleaning therefore removes the 1st codon entirely, because it is incomplete
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A       G G G  A  C  C  C  A  A  A  A
            # Zero Frame:     0 1 2       0 1 2  0  1  2  0  1  2  0  1
            # One Frame:      - 0 1       2 0 1  2  0  1  2  0  1  2  0  <- correct frame
            # Two Frame:      - - 0       1 2 0  1  2  0  1  2  0  1  2
            #
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
                [Codon("CAG"), Codon("GGA"), Codon("CCC")],
            ),
            # 2. Discontiguous CDS, plus strand, frame=1, 1bp deletion at start of exon 2
            # Codon interval cleaning therefore removes the 1st codon entirely, because it is incomplete
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A       G G G  A  C  C  C  A  A  A  A
            # Zero Frame:     0 1 2       0 1 2  0  1  2  0  1  2  0  1
            # One Frame:      - 0 1       2 0 1  2  0  1  2  0  1  2  0  <- correct frame if there was no deletion
            # Two Frame:      - - 0       1 2 0  1  2  0  1  2  0  1  2
            #
            # The codon cleaning algorithm therefore sees:
            # [A??] [GGG] [ACC] [CAA]
            # Because it detects that the Frame of the 2nd exon is 0, but it was expecting 2.
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence("AAACAAAAGGGACCCAAAAAA", alphabet, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon("GGG"), Codon("ACC"), Codon("CAA")],
            ),
            # 3. Discontiguous CDS, plus strand, frame=1,
            # 1bp insertion inside exon 2 relative to some canonical genome and we want to maintain original frame
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  T  C  C  C  A  A  A  A  A
            # Exons:          A C A       G G G     T  C  C  C  A  A  A  A
            # Zero Frame:     0 1 2       0 1 2     0  1  2  0  1  2  0  1
            # One Frame:      - 0 1       2 0 1     2  0  1  2  0  1  2  0  <- original frame
            # Two Frame:      - - 0       1 2 0     1  2  0  1  2  0  1  2
            #
            # The codon algorithm sees [CAG] [GGA] [CCC]
            # There is no frame trickery going on here, the insertion is modeled by the extra 'splice' at 11-12
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
                [Codon("CAG"), Codon("GGA"), Codon("CCC")],  # QGP
            ),
            # 4. Discontiguous CDS, minus strand, frame=2
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [Codon("CAG"), Codon("GGA"), Codon("CCC")],  # QGP
            ),
            # 5. Discontiguous CDS, minus strand, frame=2, with frameshift that leads to truncation
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.MINUS,
                        parent=chrom_parent,
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                [Codon("CAG"), Codon("GGA")],  # QG
            ),
            # 6. Discontiguous CDS, plus strand, with -1 bp programmed frameshift (overlapping interval)
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [9, 17],
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),  # G gets repeated here
                [Codon("AGG"), Codon("AAA"), Codon("GGT"), Codon("CCC"), Codon("TGA")],
            ),
            # N containing sequences now work
            (
                CDSInterval.from_location(
                    SingleInterval(
                        0,
                        9,
                        Strand.PLUS,
                        parent=Sequence("ANACGATCA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                    ),
                    [CDSFrame.ZERO],
                ),
                [Codon("ANA"), Codon("CGA"), Codon("TCA")],
            ),
            (
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
                [Codon("NNG"), Codon("GGT"), Codon("ACC")],
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
                [Codon("ATA"), Codon("CGA"), Codon("TCA")],
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
                [Codon("CAG"), Codon("GGA"), Codon("CCC")],  # QGP
            ),
        ],
    )
    def test_scan_codons_lower_case(self, cds, expected):
        assert list(cds.scan_codons()) == expected
        # run again to ensure caching is not a problem
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # chunk slices off only intergenic bases
            # frame for chunk should be same as frame for full
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     G G G  A  C  C  C  A  A  A  A
            # Zero Frame:     0 1 2 0     1 2 0  1  2  0  1  2  0  1  2
            # One Frame:      - 0 1 2     0 1 2  0  1  2  0  1  2  0  1
            # Two Frame:      - - 0 1     2 0 1  2  0  1  2  0  1  2  0
            #
            # Frame = ZERO, ONE
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[2:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    2,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("ACA"), Codon("AGG"), Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # chunk slices off first base of exon1
            # this removes the first codon
            # analysis of frame for resulting chunk:
            # Index:      0 1 2 | 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A | C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A | C A A     G G G  A  C  C  C  A  A  A  A
            # Zero Frame:       | 0 1 2     0 1 2  0  1  2  0  1  2  0  1
            # One Frame:        | - 0 1     2 0 1  2  0  1  2  0  1  2  0
            # Two Frame:        | - - 0     1 2 0  1  2  0  1  2  0  1  2      <- correct frame
            #
            # distance_from_tx_start = 0   input_frame = 0  chunk_shift = 1
            # Frame = TWO
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[3:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    3,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("AGG"), Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # chunk slices off 1st exon
            # Index:      0 1 2 3 4 5 6 7 | 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A | G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     | G G G  A  C  C  C  A  A  A  A
            # Zero Frame:                 | 0 1 2  0  1  2  0  1  2  0  1
            # One Frame:                  | - 0 1  2  0  1  2  0  1  2  0
            # Two Frame:                  | - - 0  1  2  0  1  2  0  1  2          <- correct frame
            #
            # frame for chunk:
            # distance_from_tx_start = 4   input_frame = 1  chunk_shift = 0
            # Frame = TWO
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[7:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    7,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # chunk slices of 1st exon and 1bp of exon 2
            # Index:      0 1 2 3 4 5 6 7 8 | 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G | G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     G | G G  A  C  C  C  A  A  A  A
            # Zero Frame:                   | 0 1  2  0  1  2  0  1  2  0
            # One Frame:                    | - 0  1  2  0  1  2  0  1  2          <- correct frame
            # Two Frame:                    | - -  0  1  2  0  1  2  0  1
            #
            # frame for chunk:
            # distance_from_tx_start = 4   input_frame = 1   chunk_shift = 1
            # Frame = ONE
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[9:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    9,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
        ],
    )
    def test_scan_codons_chunk_relative(self, cds, expected):
        assert list(cds.scan_codons()) == expected
        # run again to ensure caching is not a problem
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # chunk slices off only intergenic bases
            # frame for chunk should be same as frame for full
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T T     C C C  T  G  G  G  T  T  T  T
            # Zero Frame:     2 1 0 2     1 0 2  1  0  2  1  0  2  1  0
            # One Frame:      1 0 2 1     0 2 1  0  2  1  0  2  1  0  -
            # Two Frame:      0 2 1 0     2 1 0  2  1  0  2  1  0  -  -
            #
            # Frame = ZERO, TWO in transcription direction => [CDSFrame.TWO, CDSFrame.ZERO]
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.MINUS,
                    [CDSFrame.TWO, CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[2:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    2,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("TTT"), Codon("TGG"), Codon("GTC"), Codon("CCT"), Codon("TGT")],
            ),
            # chunk slices off first base of exon1
            # this removes the first codon
            # analysis of frame for resulting chunk:
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 | 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  | A  A  A
            # Exons:          T G T T     C C C  T  G  G  G  T  T  T  | T
            # Zero Frame:     2 1 0 2     1 0 2  1  0  2  1  0  2  1  | 0
            # One Frame:      1 0 2 1     0 2 1  0  2  1  0  2  1  0  | -
            # Two Frame:      0 2 1 0     2 1 0  2  1  0  2  1  0  -  | -
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.MINUS,
                    [CDSFrame.TWO, CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[1:18],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    1,
                                    18,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("TGG"), Codon("GTC"), Codon("CCT"), Codon("TGT")],
            ),
            # chunk slices off 1st exon
            # Index:      0 1 2 3 4 5 6 | 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A | A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T T   |   C C C  T  G  G  G  T  T  T  T
            # Zero Frame:     2 1 0 2   |   1 0 2  1  0  2  1  0  2  1  0
            # One Frame:      1 0 2 1   |   0 2 1  0  2  1  0  2  1  0  -
            # Two Frame:      0 2 1 0   |   2 1 0  2  1  0  2  1  0  -  -
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.MINUS,
                    [CDSFrame.TWO, CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[1:7],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    1,
                                    7,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("TGT")],
            ),
            # chunk slices of 1st exon and 1bp of exon 2
            # Index:      0 1 2 3 4 | 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A | A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T | T     C C C  T  G  G  G  T  T  T  T
            # Zero Frame:     2 1 0 | 2     1 0 2  1  0  2  1  0  2  1  0
            # One Frame:      1 0 2 | 1     0 2 1  0  2  1  0  2  1  0  -
            # Two Frame:      0 2 1 | 0     2 1 0  2  1  0  2  1  0  -  -
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.MINUS,
                    [CDSFrame.TWO, CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[1:5],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    1,
                                    5,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [Codon("TGT")],
            ),
        ],
    )
    def test_scan_codons_chunk_relative_negative_strand(self, cds, expected):
        assert list(cds.scan_codons()) == expected
        # run again to ensure caching is not a problem
        assert list(cds.scan_codons()) == expected

    @pytest.mark.parametrize(
        "cds,expected_frames,exp_codons",
        [
            # NOT CHUNK RELATIVE:
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A       G G G  A  C  C  C  A  A  A
            # Zero Frame:     0 1 2       0 1 2  0  1  2  0  1  2  0
            # One Frame:      - 0 1       2 0 1  2  0  1  2  0  1  2         <- correct frame
            # Two Frame:      - - 0       1 2 0  1  2  0  1  2  0  1
            #
            # frames are trivial; return input frames
            # constructing a new CDS with these frames returns the same translation
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8],
                        [5, 17],
                        Strand.PLUS,
                        parent=Sequence(
                            "AAACAAAAGGGACCCAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                        ),
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [CDSFrame.ONE, CDSFrame.TWO],
                [Codon("CAG"), Codon("GGA"), Codon("CCC")],
            ),
            # chunk slices off only intergenic bases
            # frame for chunk should be same as frame for full
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     G G G  A  C  C  C  A  A  A  A
            # Zero Frame:     0 1 2 0     1 2 0  1  2  0  1  2  0  1  2
            # One Frame:      - 0 1 2     0 1 2  0  1  2  0  1  2  0  1
            # Two Frame:      - - 0 1     2 0 1  2  0  1  2  0  1  2  0
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[2:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    2,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ZERO, CDSFrame.ONE],
                [Codon("ACA"), Codon("AGG"), Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # chunk slices off first base of exon1
            # this removes the first codon
            # analysis of frame for resulting chunk:
            # Index:      0 1 2 | 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A | C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A | C A A     G G G  A  C  C  C  A  A  A  A
            # Zero Frame:       | 0 1 2     0 1 2  0  1  2  0  1  2  0  1
            # One Frame:        | - 0 1     2 0 1  2  0  1  2  0  1  2  0
            # Two Frame:        | - - 0     1 2 0  1  2  0  1  2  0  1  2      <- correct frame
            #
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[3:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    3,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.TWO, CDSFrame.ONE],
                [Codon("AGG"), Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # chunk slices off 1st exon
            # Index:      0 1 2 3 4 5 6 7 | 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A | G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     | G G G  A  C  C  C  A  A  A  A
            # Zero Frame:                 | 0 1 2  0  1  2  0  1  2  0  1
            # One Frame:                  | - 0 1  2  0  1  2  0  1  2  0
            # Two Frame:                  | - - 0  1  2  0  1  2  0  1  2          <- correct frame
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[7:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    7,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.TWO],
                [Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # chunk slices of 1st exon and 1bp of exon 2; this does not change the translation from above
            # Index:      0 1 2 3 4 5 6 7 8 | 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G | G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     G | G G  A  C  C  C  A  A  A  A
            # Zero Frame:                   | 0 1  2  0  1  2  0  1  2  0
            # One Frame:                    | - 0  1  2  0  1  2  0  1  2          <- correct frame
            # Two Frame:                    | - -  0  1  2  0  1  2  0  1
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[9:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    9,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ONE],
                [Codon("GAC"), Codon("CCA"), Codon("AAA")],
            ),
            # The frame produced by chunk_relative_frames will be incorrect when the transcript has overlaps and
            # is chunk relative
            # THIS EXAMPLE IS NOT CHUNK RELATIVE, SO THE FRAMES ARE CORRECT
            # T at position 9 is repeated twice in the transcription
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A G G A A A G T C  C  C  T  G  A  A  A  A  A  A
            # Exons:          A G G A A A G T
            #                               T C  C  C  T  G  A  A  A
            # Zero Frame:     0 1 2 0 1 2 0 1
            #                               2 0  1  2  0  1  2  0  1
            # One Frame:      - 0 1 2 0 1 2 0
            #                               1 2  0  1  2  0  1  2  0
            # Two Frame:      - - 0 1 2 0 1 2
            #                               0 1  2  0  1  2  0  1  2
            (
                CDSInterval(
                    [2, 9],
                    [10, 18],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.TWO],
                    parent_or_seq_chunk_parent=Sequence(
                        "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                    ),
                ),
                [CDSFrame.ZERO, CDSFrame.TWO],
                [Codon("AGG"), Codon("AAA"), Codon("GTT"), Codon("CCC"), Codon("TGA")],
            ),
            # Now take the same interval and make it chunk relative, where the chunk only removes intergenic bases
            # this does not change the frames
            (
                CDSInterval(
                    [2, 9],
                    [10, 18],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.TWO],
                    parent_or_seq_chunk_parent=Parent(
                        id="test:1-19",
                        sequence=Sequence(
                            "AAAGGAAAGTCCCTGAAAAAA"[1:19],
                            alphabet,
                            id="test:1-19",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    1,
                                    19,
                                    Strand.PLUS,
                                    parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ZERO, CDSFrame.TWO],
                [Codon("AGG"), Codon("AAA"), Codon("GTT"), Codon("CCC"), Codon("TGA")],
            ),
            # slicing off the last exon and 1bp of the 1st exon (removing the double T)
            # Index:      0 1 2 3 4 5 6 7 8 | 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A G G A A A G | T C  C  C  T  G  A  A  A  A  A  A
            # Exons:          A G G A A A G | T
            #                               | T C  C  C  T  G  A  A  A
            # Zero Frame:     0 1 2 0 1 2 0 |
            #                               |
            # One Frame:      - 0 1 2 0 1 2 |
            #                               |
            # Two Frame:      - - 0 1 2 0 1 |
            #                               |
            (
                CDSInterval(
                    [2, 9],
                    [10, 18],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.TWO],
                    parent_or_seq_chunk_parent=Parent(
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
                ),
                [CDSFrame.ZERO],
                [Codon("AGG"), Codon("AAA")],
            ),
            # transcript starts in 1 frame
            # chunk slices off only intergenic bases
            # frame for chunk should be same as frame for full
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          A C A A     G G G  A  C  C  C  A  A  A  A
            # Zero Frame:     0 1 2 0     1 2 0  1  2  0  1  2  0  1  2
            # One Frame:      - 0 1 2     0 1 2  0  1  2  0  1  2  0  1
            # Two Frame:      - - 0 1     2 0 1  2  0  1  2  0  1  2  0
            (
                CDSInterval(
                    [2, 8],
                    [6, 19],
                    Strand.PLUS,
                    [CDSFrame.ONE, CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[2:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    2,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ONE, CDSFrame.ZERO],
                [Codon("CAA"), Codon("GGG"), Codon("ACC"), Codon("CAA")],
            ),
        ],
    )
    def test_chunk_relative_frames(self, cds, expected_frames, exp_codons):
        assert cds.chunk_relative_frames == expected_frames
        new_cds = CDSInterval.from_dict(
            cds.to_dict(chromosome_relative_coordinates=False),
            parent_or_seq_chunk_parent=Parent(
                sequence=Sequence(str(cds.chunk_relative_location.parent.sequence), Alphabet.NT_EXTENDED_GAPPED),
                sequence_type=SequenceType.CHROMOSOME,
            ),
        )
        assert list(new_cds.scan_codons()) == exp_codons

    @pytest.mark.parametrize(
        "cds,expected_frames,exp_codons",
        [
            # NOT CHUNK RELATIVE:
            # Discontiguous CDS, plus strand, frame=1, codons don't reach end of CDS
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T     T C C C  T  G  G  G  T  T  T
            # Zero Frame:     1 0 2     1 0 2 1  0  2  1  0  2  1  0
            # One Frame:      0 2 1     0 2 1 0  2  1  0  2  1  0  -         <- correct frame
            # Two Frame:      2 1 0     2 1 0 2  1  0  2  1  0  -  -
            #
            # transcription direction frame -> ZERO, ONE
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 7],
                        [5, 18],
                        Strand.MINUS,
                        parent=Sequence(
                            "AAACAAAAGGGACCCAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                        ),
                    ),
                    [CDSFrame.ONE, CDSFrame.ONE],
                ),
                [CDSFrame.ONE, CDSFrame.ONE],
                [Codon("TTG"), Codon("GGT"), Codon("CCC"), Codon("TTG")],
            ),
            # chunk slices off only intergenic bases
            # frame for chunk should be same as frame for full
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T     T C C C  T  G  G  G  T  T  T
            # Zero Frame:     1 0 2     1 0 2 1  0  2  1  0  2  1  0
            # One Frame:      0 2 1     0 2 1 0  2  1  0  2  1  0  -         <- correct frame
            # Two Frame:      2 1 0     2 1 0 2  1  0  2  1  0  -  -
            (
                CDSInterval(
                    [2, 7],
                    [5, 18],
                    Strand.MINUS,
                    [CDSFrame.ONE, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[2:20],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    2,
                                    20,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ONE, CDSFrame.ONE],
                [Codon("TTG"), Codon("GGT"), Codon("CCC"), Codon("TTG")],
            ),
            # chunk slices off 3 last bases of exon2
            # this removes the first codon
            # analysis of frame for resulting chunk:
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 | 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G G  A  C  C  C  | A  A  A  A  A  A
            # Exons:          T G T     T C C C  T  G  G  G  | T  T  T
            # Zero Frame:     1 0 2     1 0 2 1  0  2  1  0  | 2  1  0
            # One Frame:      0 2 1     0 2 1 0  2  1  0  2  | 1  0  -         <- correct frame
            # Two Frame:      2 1 0     2 1 0 2  1  0  2  1  | 0  -  -
            #
            (
                CDSInterval(
                    [2, 7],
                    [5, 18],
                    Strand.MINUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[1:15],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    1,
                                    15,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ONE, CDSFrame.ONE],
                [Codon("GGT"), Codon("CCC"), Codon("TTG")],
            ),
            # chunk slices off 1st exon (in coding orientation)
            # no translation at all because the last codon is only 2 bases
            # Index:      0 1 2 3 4 5 | 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A | A A G G G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T   |   T C C C  T  G  G  G  T  T  T
            # Zero Frame:     1 0 2   |   1 0 2 1  0  2  1  0  2  1  0
            # One Frame:      0 2 1   |   0 2 1 0  2  1  0  2  1  0  -         <- correct frame
            # Two Frame:      2 1 0   |   2 1 0 2  1  0  2  1  0  -  -
            (
                CDSInterval(
                    [2, 7],
                    [5, 18],
                    Strand.MINUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[:6],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    0,
                                    6,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.TWO],
                [],
            ),
            # chunk slices off most of the 1st coding direction exon
            # Index:      0 1 2 3 4 5 6 7 8 9 | 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A C A A A A G G | G  A  C  C  C  A  A  A  A  A  A
            # Exons:          T G T     T C C | C  T  G  G  G  T  T  T
            # Zero Frame:     1 0 2     1 0 2 | 1  0  2  1  0  2  1  0
            # One Frame:      0 2 1     0 2 1 | 0  2  1  0  2  1  0  -         <- correct frame
            # Two Frame:      2 1 0     2 1 0 | 2  1  0  2  1  0  -  -
            (
                CDSInterval(
                    [2, 7],
                    [5, 18],
                    Strand.MINUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=Parent(
                        id="chunk1",
                        sequence=Sequence(
                            "AAACAAAAGGGACCCAAAAAA"[0:10],
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="chunk1",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    0,
                                    10,
                                    Strand.PLUS,
                                    parent=Parent(id="seq", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
                [CDSFrame.ONE, CDSFrame.TWO],
                [Codon("TTG")],
            ),
        ],
    )
    def test_chunk_relative_frames_negative_strand(self, cds, expected_frames, exp_codons):
        assert cds.chunk_relative_frames == expected_frames
        new_cds = CDSInterval.from_dict(
            cds.to_dict(chromosome_relative_coordinates=False),
            parent_or_seq_chunk_parent=Parent(
                sequence=Sequence(str(cds.chunk_relative_location.parent.sequence), Alphabet.NT_EXTENDED_GAPPED),
                sequence_type=SequenceType.CHROMOSOME,
            ),
        )
        assert list(new_cds.scan_codons()) == exp_codons

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
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ONE, CDSFrame.TWO],
                ),
                [
                    SingleInterval(
                        12,
                        15,
                        Strand.MINUS,
                        parent=chrom_parent,
                    ),  # CAG
                    SingleInterval(
                        9,
                        12,
                        Strand.MINUS,
                        parent=chrom_parent,
                    ),  # GGA
                    CompoundInterval(
                        [3, 8],
                        [5, 9],
                        Strand.MINUS,
                        parent=chrom_parent,
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
                        parent=chrom_parent,
                    ),
                    [CDSFrame.TWO, CDSFrame.TWO],
                ),
                [
                    SingleInterval(
                        12,
                        15,
                        Strand.MINUS,
                        parent=chrom_parent,
                    ),  # CAG
                    SingleInterval(
                        9,
                        12,
                        Strand.MINUS,
                        parent=chrom_parent,
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
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                [
                    SingleInterval(
                        2,
                        5,
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),  # AGG
                    SingleInterval(
                        5,
                        8,
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),  # AAA
                    CompoundInterval(
                        [8, 8],
                        [9, 10],
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),  # GGT, G gets repeated
                    SingleInterval(
                        10,
                        13,
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),  # CCC
                    SingleInterval(
                        13,
                        16,
                        Strand.PLUS,
                        parent=chrom_parent,
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
    def test_scan_chunk_relative_codon_locations(self, cds, expected):
        assert list(cds.scan_chunk_relative_codon_locations()) == expected
        # run again to prove caching doesn't kill the iterator
        assert list(cds.scan_chunk_relative_codon_locations()) == expected
        # show that the tuple accessors works too
        assert cds.chunk_relative_codon_locations == tuple(expected)

    @pytest.mark.parametrize(
        "cds,expected",
        [
            # single exon, even though the CDS exists on a chunk because these are chromosome codons
            # we get the full set of codons, sequenceless
            (
                CDSInterval(
                    [2],
                    [17],
                    Strand.PLUS,
                    [CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                [
                    SingleInterval(2, 5, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # chunk relative, so the returned codons are sequenceless but still the whole thing
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                [
                    SingleInterval(2, 5, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, sequenceless_chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # genome relative, so returned codons are identical but have sequence
            (
                CDSInterval(
                    [2],
                    [17],
                    Strand.PLUS,
                    [CDSFrame.ZERO],
                    parent_or_seq_chunk_parent=chrom_parent,
                ),
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ),
            # genome relative, so returned codons are identical but have sequence
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chrom_parent,
                ),
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, chrom_parent),
                ],
            ),
        ],
    )
    def test_scan_chromosome_codon_locations(self, cds, expected):
        assert list(cds.scan_chromosome_codon_locations()) == expected
        # test that caching doesn't break the iterator
        assert list(cds.scan_chromosome_codon_locations()) == expected
        assert cds.chromosome_codon_locations == tuple(expected)

    @pytest.mark.parametrize(
        "cds,start,end,expected",
        [
            # single-exon, zero frame, no slicing
            (
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                None,
                None,
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ),
            # single-exon, zero frame, slice off first codon
            (
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                4,
                None,
                [
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ),
            # single-exon, zero frame, slice off last codon
            (
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                None,
                16,
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                ],
            ),
            # single-exon, one frame, no slicing
            (
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ONE], parent_or_seq_chunk_parent=chrom_parent),
                None,
                None,
                [
                    SingleInterval(3, 6, Strand.PLUS, chrom_parent),
                    SingleInterval(6, 9, Strand.PLUS, chrom_parent),
                    SingleInterval(9, 12, Strand.PLUS, chrom_parent),
                    SingleInterval(12, 15, Strand.PLUS, chrom_parent),
                ],
            ),
            # single-exon, one frame, slice off first codon
            (
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ONE], parent_or_seq_chunk_parent=chrom_parent),
                4,
                None,
                [
                    SingleInterval(6, 9, Strand.PLUS, chrom_parent),
                    SingleInterval(9, 12, Strand.PLUS, chrom_parent),
                    SingleInterval(12, 15, Strand.PLUS, chrom_parent),
                ],
            ),
            # single-exon, one frame, slice off last codon
            (
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ONE], parent_or_seq_chunk_parent=chrom_parent),
                None,
                14,
                [
                    SingleInterval(3, 6, Strand.PLUS, chrom_parent),
                    SingleInterval(6, 9, Strand.PLUS, chrom_parent),
                    SingleInterval(9, 12, Strand.PLUS, chrom_parent),
                ],
            ),
            # single-exon, zero frame, no slicing, minus strand
            (
                CDSInterval([2], [17], Strand.MINUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                None,
                None,
                [
                    SingleInterval(14, 17, Strand.MINUS, chrom_parent),
                    SingleInterval(11, 14, Strand.MINUS, chrom_parent),
                    SingleInterval(8, 11, Strand.MINUS, chrom_parent),
                    SingleInterval(5, 8, Strand.MINUS, chrom_parent),
                    SingleInterval(2, 5, Strand.MINUS, chrom_parent),
                ],
            ),
            # single-exon, zero frame, slice off last codon, minus strand
            (
                CDSInterval([2], [17], Strand.MINUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                4,
                None,
                [
                    SingleInterval(14, 17, Strand.MINUS, chrom_parent),
                    SingleInterval(11, 14, Strand.MINUS, chrom_parent),
                    SingleInterval(8, 11, Strand.MINUS, chrom_parent),
                    SingleInterval(5, 8, Strand.MINUS, chrom_parent),
                ],
            ),
            # single-exon, zero frame, slice off first codon, minus strand
            (
                CDSInterval([2], [17], Strand.MINUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                None,
                16,
                [
                    SingleInterval(11, 14, Strand.MINUS, chrom_parent),
                    SingleInterval(8, 11, Strand.MINUS, chrom_parent),
                    SingleInterval(5, 8, Strand.MINUS, chrom_parent),
                    SingleInterval(2, 5, Strand.MINUS, chrom_parent),
                ],
            ),
            # single-exon, one frame, no slicing, minus strand
            (
                CDSInterval([2], [17], Strand.MINUS, [CDSFrame.ONE], parent_or_seq_chunk_parent=chrom_parent),
                None,
                None,
                [
                    SingleInterval(13, 16, Strand.MINUS, chrom_parent),
                    SingleInterval(10, 13, Strand.MINUS, chrom_parent),
                    SingleInterval(7, 10, Strand.MINUS, chrom_parent),
                    SingleInterval(4, 7, Strand.MINUS, chrom_parent),
                ],
            ),
            # single-exon, one frame, slice off last codon, minus strand
            (
                CDSInterval([2], [17], Strand.MINUS, [CDSFrame.ONE], parent_or_seq_chunk_parent=chrom_parent),
                5,
                None,
                [
                    SingleInterval(13, 16, Strand.MINUS, chrom_parent),
                    SingleInterval(10, 13, Strand.MINUS, chrom_parent),
                    SingleInterval(7, 10, Strand.MINUS, chrom_parent),
                ],
            ),
            # single-exon, one frame, slice off first codon, minus strand
            (
                CDSInterval([2], [17], Strand.MINUS, [CDSFrame.ONE], parent_or_seq_chunk_parent=chrom_parent),
                None,
                14,
                [
                    SingleInterval(10, 13, Strand.MINUS, chrom_parent),
                    SingleInterval(7, 10, Strand.MINUS, chrom_parent),
                    SingleInterval(4, 7, Strand.MINUS, chrom_parent),
                ],
            ),
            # slice off no codons with all null
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                None,
                None,
                [
                    SingleInterval(2, 5, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, sequenceless_chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # slice off no codons with a large number
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                None,
                30,
                [
                    SingleInterval(2, 5, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, sequenceless_chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # slice off the first codon
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                3,
                None,
                [
                    SingleInterval(5, 8, Strand.PLUS, sequenceless_chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # slice off the first and last codon
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                3,
                15,
                [
                    SingleInterval(5, 8, Strand.PLUS, sequenceless_chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # slice off the first two codons
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chunk_parent2_12,
                ),
                8,
                17,
                [
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, sequenceless_chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, sequenceless_chrom_parent),
                ],
            ),
            # slice off the first codon in relative coordinates
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chrom_parent,
                ),
                3,
                None,
                [
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, chrom_parent),
                ],
            ),
            # slice off the first and last codons
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chrom_parent,
                ),
                3,
                15,
                [
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, chrom_parent),
                ],
            ),
            # slice off the first two codons
            (
                CDSInterval(
                    [2, 8],
                    [9, 17],
                    Strand.PLUS,
                    [CDSFrame.ZERO, CDSFrame.ONE],
                    parent_or_seq_chunk_parent=chrom_parent,
                ),
                8,
                17,
                [
                    CompoundInterval([8, 8], [9, 10], Strand.PLUS, chrom_parent),
                    SingleInterval(10, 13, Strand.PLUS, chrom_parent),
                    SingleInterval(13, 16, Strand.PLUS, chrom_parent),
                ],
            ),
        ],
    )
    def test_scan_chromosome_codon_locations_sliced(self, cds, start, end, expected):
        assert list(cds.scan_chromosome_codon_locations(start, end)) == expected

    @pytest.mark.parametrize(
        "cds,start,end,expected",
        (
            # no slicing, so comes back as full length
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                None,
                None,
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ],
            # 1bp slicing from front, so comes back as full length
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                3,
                None,
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ],
            # 1bp slicing from end, so comes back as full length
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                None,
                16,
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ],
            # 1bp slicing from both ends, so comes back as full length
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                3,
                16,
                [
                    SingleInterval(2, 5, Strand.PLUS, chrom_parent),
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                    SingleInterval(14, 17, Strand.PLUS, chrom_parent),
                ],
            ],
            # 3bp slicing from both ends, so comes back missing 1st and last codon
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                5,
                13,
                [
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                ],
            ],
            # 4bp slicing from both ends, so comes back missing 1st and last codon
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                6,
                12,
                [
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                    SingleInterval(11, 14, Strand.PLUS, chrom_parent),
                ],
            ],
            # 5bp slicing from both ends, so comes back missing 1st and last 2 codons
            # this is because these coordinates are half-open
            [
                CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=chrom_parent),
                7,
                11,
                [
                    SingleInterval(5, 8, Strand.PLUS, chrom_parent),
                    SingleInterval(8, 11, Strand.PLUS, chrom_parent),
                ],
            ],
        ),
    )
    def test_scan_chromosome_codon_locations_sliced_expand_coordinates(self, cds, start, end, expected):
        assert list(cds.scan_chromosome_codon_locations(start, end, expand_window_to_partial_codons=True)) == expected

    @pytest.mark.parametrize(
        "start,end,expected",
        (
            # no slicing
            [2, 17, (2, 17)],
            # 1bp off 1st codon -> full length
            [3, 17, (2, 17)],
            # 1bp off last codon -> full length
            [3, 16, (2, 17)],
            # 1bp off first codon and 2bp off last codon -> full length
            [3, 15, (2, 17)],
            # 4bp off 5' codon and 1bp off last codon -> lose 1st codon
            [6, 16, (5, 17)],
        ),
    )
    def test__expand_coordinates_to_codons(self, start, end, expected):
        cds = CDSInterval([2], [17], Strand.PLUS, [CDSFrame.ZERO], parent_or_seq_chunk_parent=self.chrom_parent)
        assert cds._expand_coordinates_to_codons(start, end) == expected

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
                        parent=chrom_parent,
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

    def test_translate_unknown_errors(self):
        # Makes sure ValueError raised when not allowing ambiguous translation (default action)
        cds = CDSInterval.from_location(
            SingleInterval(
                0,
                12,
                Strand.PLUS,
                parent=Sequence("ATGCTCGTNACG", Alphabet.NT_STRICT_UNKNOWN, type=SequenceType.CHROMOSOME),
            ),
            [CDSFrame.ZERO],
        )
        with pytest.raises(ValueError):
            cds.translate()

    @pytest.mark.parametrize(
        "seq, exp",
        [
            ("ATGCTCGCNACG", "MLAT"),
            ("ATGCTCGNNACG", "MLXT"),
        ],
    )
    def test_translate_ambiguous(self, seq, exp):
        # Makes sure translation works when allowing ambiguous translation
        cds = CDSInterval.from_location(
            SingleInterval(
                0,
                12,
                Strand.PLUS,
                parent=Sequence(seq, Alphabet.NT_STRICT_UNKNOWN, type=SequenceType.CHROMOSOME),
            ),
            [CDSFrame.ZERO],
        )
        obs = cds.translate(strict=False)
        assert obs == Sequence(exp, alphabet=Alphabet.AA_STRICT_UNKNOWN)

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
            # Index:      0 1 2 3 4 5 6 7 8
            # Sequence:   A T A C G A T C A
            # Exons:      T A T G C     G
            # Zero Frame: 2 1 0 2 1     0
            # One Frame:  1 0 2 1 0     -
            # Two Frame:  0 2 1 0 -     -
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
                [Codon("TAT")],
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
                [Codon("ACA"), Codon("AGG"), Codon("ACC"), Codon("CAA")],
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
                [Codon("AGG"), Codon("GAC"), Codon("CCA")],
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
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE],
                ),
                [Codon("AGG"), Codon("AAA"), Codon("GTC"), Codon("CCT"), Codon("GAA")],
            ),
            # this overlapping interval gets merged into [2,5], [8,18] with 0bp offset
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 13, 18],
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ZERO, CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon("AGG"), Codon("GTC"), Codon("CCT"), Codon("GAA")],
            ),
            # this overlapping interval gets merged into [2,5], [8,18] with 1bp offset
            (
                CDSInterval.from_location(
                    CompoundInterval(
                        [2, 8, 12],
                        [5, 13, 18],
                        Strand.PLUS,
                        parent=chrom_parent,
                    ),
                    [CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ZERO],
                ),
                [Codon("GGG"), Codon("TCC"), Codon("CTG"), Codon("AAA")],
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
        assert list(cds.scan_codons()) == [Codon("TCC"), Codon("CTG"), Codon("AAA")]

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
                CompoundInterval([222218], [222233], Strand.PLUS),
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
        "pos,exp_cds_pos",
        [(2, 0), (3, 1), (6, 2), (10, 3), (12, 5)],
    )
    def test_sequence_pos_to_cds_plus_strand(self, pos, exp_cds_pos):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
        )
        assert cds.sequence_pos_to_cds(pos) == exp_cds_pos

    @pytest.mark.parametrize(
        "pos",
        [0, 1, 5, 13, 20],
    )
    def test_sequence_pos_to_cds_plus_strand_exception(self, pos):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
        )
        with pytest.raises(InvalidPositionException):
            _ = cds.sequence_pos_to_cds(pos)

    @pytest.mark.parametrize(
        "pos,exp_aa_pos",
        [(2, 0), (3, 0), (6, 0), (10, 1), (12, 1)],
    )
    def test_sequence_pos_to_aa_pos_plus_strand(self, pos, exp_aa_pos):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
        )
        assert cds.sequence_pos_to_amino_acid(pos) == exp_aa_pos

    @pytest.mark.parametrize(
        "pos,exp_cds_pos",
        [(16, 0), (14, 2), (12, 3), (10, 5), (6, 6), (2, 8)],
    )
    def test_sequence_pos_to_cds_minus_strand(self, pos, exp_cds_pos):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.MINUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
        )
        assert cds.sequence_pos_to_cds(pos) == exp_cds_pos

    @pytest.mark.parametrize(
        "pos",
        [0, 1, 5, 13, 20],
    )
    def test_sequence_pos_to_cds_minus_strand_exception(self, pos):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.MINUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
        )
        with pytest.raises(InvalidPositionException):
            _ = cds.sequence_pos_to_cds(pos)

    @pytest.mark.parametrize(
        "pos,exp_aa_pos",
        [(16, 0), (14, 0), (12, 1), (10, 1), (6, 2), (2, 2)],
    )
    def test_sequence_pos_to_aa_pos_minus_strand(self, pos, exp_aa_pos):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.MINUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
        )
        assert cds.sequence_pos_to_amino_acid(pos) == exp_aa_pos

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
                    "cds_frames": ["ZERO", "TWO", "ONE"],
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
        "parent,exp_codons",
        [
            # normal flat genome
            # T at position 9 is repeated twice in the transcription
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A G G A A A G T C  C  C  T  G  A  A  A  A  A  A
            # Exons:          A G G A A A G T
            #                               T C  C  C  T  G  A  A  A
            # Zero Frame:     0 1 2 0 1 2 0 1
            #                               2 0  1  2  0  1  2  0  1
            # One Frame:      - 0 1 2 0 1 2 0
            #                               1 2  0  1  2  0  1  2  0
            # Two Frame:      - - 0 1 2 0 1 2
            #                               0 1  2  0  1  2  0  1  2
            (
                Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                [Codon("AGG"), Codon("AAA"), Codon("GTT"), Codon("CCC"), Codon("TGA")],
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
                [Codon("AAA"), Codon("GTT"), Codon("CCC"), Codon("TGA")],
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
                # translation doesn't change because the original 1st exon is sliced off in codon iteration
                [Codon("CCC"), Codon("TGA")],
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
                [Codon("AGG"), Codon("AAA")],
            ),
        ],
    )
    def test_cds_frame_slicing_overlapping_cds(self, parent, exp_codons):
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
        assert list(cds.scan_codons()) == exp_codons

    @pytest.mark.parametrize(
        "parent,exp_codons",
        [
            # normal flat genome
            # Index:      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A G G A A A G T C  C  C  T  G  A  A  A  A  A  A
            # Exons:          A G     A       C  C  C     G  A  A
            # Zero Frame:     0 1     2       0  1  2     0  1  2
            # One Frame:      - 0     1       2  0  1     2  0  1
            # Two Frame:      - -     0       1  2  0     1  2  0
            (
                Sequence("AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME),
                [Codon("AGA"), Codon("CCC"), Codon("GAA")],
            ),
            # sequence subset that removes 1st exon
            # Index:      0 1 2 3 4 | 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
            # Sequence:   A A A G G | A A A G T C  C  C  T  G  A  A  A  A  A  A
            # Exons:                |   A       C  C  C     G  A  A
            # Zero Frame:           |   0       1  2  0     1  2  0
            # One Frame:            |   -       0  1  2     0  1  2         <- correct frame
            # Two Frame:            |   -       -  0  1     2  0  1
            #
            # distance_from_tx_start = 2   input_frame = 2  chunk_shift = 0
            #
            # FRAME = ONE/TWO (does not matter), ZERO, ZERO
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
                [Codon("CCC"), Codon("GAA")],
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
                [Codon("AGA"), Codon("CCC")],
            ),
        ],
    )
    def test_cds_frame_slicing_small_exons(self, parent, exp_codons):
        cds = CDSInterval(
            [2, 6, 10, 14],
            [4, 7, 13, 17],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.ZERO, CDSFrame.ZERO],
            parent_or_seq_chunk_parent=parent,
        )
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

    def test_negative_strand_overlapping_codon_iterator(self):
        """
        This CDS has a 2bp overlap. In chunk relative coordinates, it looks like this:

        A T C G A T C G A T C  G
        0 1 2 3 4 5 6 7 8 9 10 11
                8 7 6 5 4 3 2  1
             1110 9

        which reverse complements into:
        CGATCGATATCGAT
                ^^        << duplicated "AT"
         GATCGATATCGA     << ORF is taken into account, these are actual codons in the chunk frame
        """
        cds = CDSInterval(
            **dict(
                cds_starts=[9, 202],
                cds_ends=[204, 587],
                strand=Strand.MINUS,
                frames_or_phases=[CDSFrame.ONE, CDSFrame.ZERO],
                parent_or_seq_chunk_parent=Parent(
                    id="chrVI.0:198-210",
                    sequence=Sequence(
                        "ATCGATCGATCG",
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="chrVI.0:198-210",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                198,
                                210,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
            )
        )
        orf = "GATCGATATCGA"
        assert str(cds.chunk_relative_location.extract_sequence()) == "CGATCGATATCGAT"
        assert str(cds.extract_sequence()) == orf
        # using codon iterator still produces the same result
        codons = (str(codon_location.extract_sequence()) for codon_location in cds.chunk_relative_codon_locations)
        old_seq = "".join(codons)
        assert old_seq == orf

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

    def test_codon_iteration_sequence_chunk(self):
        """Show that the heuristic that checks for codon iteration on spliced transcripts produces
        a valid translation by using the correct algorithm"""
        cds = CDSInterval(
            [143044, 143490],
            [143124, 144031],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO],
            parent_or_seq_chunk_parent=Parent(
                id="chrI:143480-144050",
                sequence=Sequence(
                    "TTTAAAATAGTACTGCTGTTTCTCAAGCTGACGTCACTGTCTTCAAGGCTTTCCAATCTGCTTACCCAGAATTCTCCAGATGGTTCAACCACATCGCTTCCAAGGCCGATGAATTCGACTCTTTCCCAGCTGCCTCTGCTGCCGCTGCCGAAGAAGAAGAAGATGACGATGTCGATTTATTCGGTTCCGACGATGAAGAAGCTGACGCTGAAGCTGAAAAGTTGAAGGCTGAAAGAATTGCCGCATACAACGCTAAGAAGGCTGCTAAGCCAGCTAAGCCAGCTGCTAAGTCCATTGTCACTCTAGATGTCAAGCCATGGGATGATGAAACCAATTTGGAAGAAATGGTTGCTAACGTCAAGGCCATCGAAATGGAAGGTTTGACCTGGGGTGCTCACCAATTTATCCCAATTGGTTTCGGTATCAAGAAGTTGCAAATTAACTGTGTTGTCGAAGATGACAAGGTTTCCTTGGATGACTTGCAACAAAGCATTGAAGAAGACGAAGACCACGTCCAATCTACCGATATTGCTGCTATGCAAAAATTATAAAAGGCTTTTTTATAAACTT",  # noqa: E501
                    Alphabet.NT_EXTENDED_GAPPED,
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            143480,
                            144050,
                            Strand.PLUS,
                            parent=Parent(id="chrI", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
        )
        # show that the function returns a sensible translation that starts from the 1st in-frame codon
        # need to re-build the CDS because the cache now contains the empty version caused by the mock
        cds = CDSInterval.from_dict(cds.to_dict(), parent_or_seq_chunk_parent=cds.chunk_relative_location.parent)
        translation = cds.translate()
        assert (
            str(translation)
            == "TAVSQADVTVFKAFQSAYPEFSRWFNHIASKADEFDSFPAASAAAAEEEEDDDVDLFGSDDEEADAEAEKLKAERIAAYNAKKAAKPAKPAAKSIVTLDVKPWDDETNLEEMVANVKAIEMEGLTWGAHQFIPIGFGIKKLQINCVVEDDKVSLDDLQQSIEEDEDHVQSTDIAAMQKL*"  # noqa: E501
        )

    def test_codon_iteration_sequence_chunk_negative_strand(self):
        """Show that the heuristic that checks for codon iteration on negative strand spliced transcripts produces
        a valid translation by using the correct algorithm"""
        cds = CDSInterval(
            [148464, 151967],
            [151877, 152037],
            Strand.MINUS,
            [CDSFrame.ONE, CDSFrame.ZERO],
            parent_or_seq_chunk_parent=Parent(
                id="chrI:148460-151890",
                sequence=Sequence(
                    "TTCATTATGTAGATTCATATATAGAATACCAATTATGATTGACCCAATAGCCATCAAAATCAGTAGTTATTAATACTTGTCTTTCTAGGAGCCATTTGCATATTTCTGATATTTCATGAAGCGAAAGTACTTCACGACACCTAGATTGCAATCTACTCAATGTTATCCCTGGATGAAATATTATTTCGTTAACGACCATAGTAACTACCTGCTTCCATATGTTTGGCCTAATGGAACCAGATCCATTCACCCATAAACGAGAAAATGGTTTGCCCAGTGGAACTTTGACAGCAGACTTCCTTGCTGTATTCAATTTTGTCTGAGAATTGGCATATATAATCAGAGGGGGAGTTAATGTTCGTATTTCAAATCTCCTTGAAGTATACGTTAAAGGTCGAACATTTCTCACCATTGGAATTACATCCATATTCAATAGCTCTCCCGAAATCAAATCAATTAAAACCCAAGAGGATATATCGGACGGCTCTTGATTGATAACAATAGCGTTTCCGGCCTCCAATAATTCATTAACCTTACATCTATACTGAAAAGCTACACCAAAATCTTTATAATTTCCTCTATTTTCCAAAATGTCTGGTAAAGTATCAGTACATTCAAGTTTTGAGCCATGGAGATAAATTTGCTTTTCCTTAGCCATATCCATGATGACGTTATCTATTGATTCGTTTCCAACGTTCTTCAACGCCTCTATTTCATTTCTAGTGGTCGAAGGACTTTCTATTAATATGGACCGGATCACTGTGCGAATATAATCGTCGCTTTGACTCTTCGATAAGTCCTTAGTAGAAGCGGAAATCTTTCTAGTGTAAGTTTTTTTTAAAGAAGAGATCTCTCTTTGAATCATAGAAGACATGGCCAGATCGTTGCCAGAATGTGTAAGTGTGTCATCACGTACCAGAGTAAATTTTTTTCTATTCTCTTCGTAGTTCTTATAAAGAAAGAGCGGCCTTTTTATTTCCTTTTCATCAAAAGAGGTCCAAATATCAAGCAATTTGATAAGATCTAGTTCTTCAACATCCCTCAGTGAAATCTTTTCACTTTTAATGGCTAGAACGAGCATTTTTTTCCACTTATCGACATATGCCCTCCAACCACTATGACCCATTCTTACTCGCCGTGCCGTCCATTTCTTTTTTAGGTTATCTAAAGAATTATTAGGAAATAATTTTGTTATTTTGTCCCACATTATTTCATTTTTAATACTTTTAGTAACTACAACAGCTCTGATTAAAGCCTGAACACCATCTTTAGTGCCTGCATGATAGACAGTTTTGTCTTCTTTAGTATTTTCCACGACCACCGTTGTCCTACCAGCAGACACTTTTTTGTCTCTCCTTTTGATCTTTCCATCTGATACGTTGACCGACGTACTCTTCTTAGAGATAGCGTCATCTGAAGCTTTGATCTTAGCATTCTTTTGGCGGTTCTGAAACTTACGAATTCTTTCAACTGATTTCTTTCCTCTGTGAAACCTATTTTTTTCCGTTTGGTCAAAGAAGTATATTTCCGTATCATGTATAACATCAGTAAAGGTTGCTTTTTTACTATCTTTTTCTTTCAGGATGTACCTTTGGATAGCGTCCTCTCCAACAGTGGGCAAAAAAATAATTTTTCTTCCTGATACAGGCTCTGTTCTGGCTCCTAATTTTTCGCTTTCTACCATCAAATCAACATCACCACGGACAGTCTTTTTATCTAATGTCGTTGTGGAGCCCATATATTTAGAAACGCTTTCGTAAAATTGTTCTCTCAGGTATGCTACCCCACCAATCGTATTCATAACTTTCAAAATGGCTCTCTGTCTCTGTAGTGAACGCAAAGAGCGGGCAGAAAAGCCGCCGAAGTTAACCACCTTGCCTTTGACAACAGTGCCTCCGTTTGAACTTGGACTATCTTCAGCAGATTTCGGCTCCTGTGCAGTGCTGACATGCTGCTCTAGCTTAATCCTTTTGGGATTCGAAATGTTTCCTGCAACAGAAGCATTAGTACTGTTTTTAACCTGCCTCTTCCGTTTGTTTTTATTCGGAGTTTTTTTTGAGTTTGGGGGAATTTTTAATTCACCGTGCCAGAAGAATATATCCTGTCCATCGCTGTCCGTTGTAAATCTAACAGTGTTGTTGAGTGCGACGAAATTATCCTCGTTGAGAGTTTTCAAATCGGTACGAGATTTGCCTAGCTCATCAAACCCTTTTGGAACGGATATTTCGTCTTCCGCATTTGTTAACTTTTGAAAGTTCTGAGCTGTGAACAGCCTAAAAAACTTCTTCTTTCCCTCAAAATCGTATATGCGAAAAAGCCTATACCCCCCTGTATTTTCTTTTTGCTTATCCACACTTTCTAAATAATATTCGCTTGATTTGGTAAAAGCTCGCTGAAATTCTTTTCCGGTAATTCGATTTACAACATCCATAGTTGAAATTCCTTTAAGGCCAGACTTATCTGCAATGTCATAAGTCTGATTTTGAAGTGGATAAAATCGATTAAGAAGAACTTCATTCTTTACAGCATCCTCTTTCTCTTCCATAACAAGGCCTTGATTTTGTAATAAATCAGTCGCATTGAAATTATCTAAACCTTCGACTAAGTCTTCATCTTCGAAAGCTGCCTTGCTATCTGATACAGAATCTTCATCCGCGCTATTGCTATCATACTCAAATGAAGGCGAGCCTTTAGAGTCTGGAATATCTTTCACGTATTTTACACATCTGATTTTAATGGCAGGATTCTTGGGTGATACTACAAGCACTTTCTTTAAGTACTCCTTTTCATCTAACCATGCAATAGCTGCAATAAAAGCTTTAGAAAGTCTTTTCTCTTTGTCAAATTTCAATTCACGCTTTAAATCAATTATCTGGCGAATACCATTTTTTGATCGTTTTACCACCTCAACTATTGTTGCTAAATGATCCCTAATATTAATATAGGGATTACTATCCACCCCGTCATGGCTGAATTTTTTTAGCTTCAATTGCTTCACGACGTGTCCCTTATAAATCAGTTGTGAACTTGTTAACAGGTGGTTTATTTTCTTGATACGTCCAGTCACACTTCTAGGATCTTGCCCAGTTACCTGCGCCAAATCCATAGTATTGATCCCTTTTTCTCCTGATTTGGCAACTTCGAGAAGTAGTTCAAATGCAGAATTTCCAATAGTTGACTCCTTTTTTGTGTATCCCGTTAATAATGTCCATAGGCTGTCCTCAGTAATCCCAACCGAGTATGAATGATTAGCGTCGCCTATAATATCAGTCACATTTTTAGTTGTTATAGCACCATCACAATACACCTCAATGTCCTTTTTCAATATCACGCATGAAAGCACGAACTGTTTAACTTTTTTATCAGACAAATCAAAATATTTACCAGATATATCCCACAGCTGATTCAAAGTGATTTCTTCAATGTGTCG",  # noqa: E501
                    Alphabet.NT_EXTENDED_GAPPED,
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            148460,
                            151890,
                            Strand.PLUS,
                            parent=Parent(id="chrI", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
        )
        # show that the function returns a sensible translation that starts from the 1st in-frame codon
        # need to re-build the CDS because the cache now contains the empty version caused by the mock
        cds = CDSInterval.from_dict(cds.to_dict(), parent_or_seq_chunk_parent=cds.chunk_relative_location.parent)
        translation = cds.translate()
        assert (
            str(translation)
            == "ITLNQLWDISGKYFDLSDKKVKQFVLSCVILKKDIEVYCDGAITTKNVTDIIGDANHSYSVGITEDSLWTLLTGYTKKESTIGNSAFELLLEVAKSGEKGINTMDLAQVTGQDPRSVTGRIKKINHLLTSSQLIYKGHVVKQLKLKKFSHDGVDSNPYINIRDHLATIVEVVKRSKNGIRQIIDLKRELKFDKEKRLSKAFIAAIAWLDEKEYLKKVLVVSPKNPAIKIRCVKYVKDIPDSKGSPSFEYDSNSADEDSVSDSKAAFEDEDLVEGLDNFNATDLLQNQGLVMEEKEDAVKNEVLLNRFYPLQNQTYDIADKSGLKGISTMDVVNRITGKEFQRAFTKSSEYYLESVDKQKENTGGYRLFRIYDFEGKKKFFRLFTAQNFQKLTNAEDEISVPKGFDELGKSRTDLKTLNEDNFVALNNTVRFTTDSDGQDIFFWHGELKIPPNSKKTPNKNKRKRQVKNSTNASVAGNISNPKRIKLEQHVSTAQEPKSAEDSPSSNGGTVVKGKVVNFGGFSARSLRSLQRQRAILKVMNTIGGVAYLREQFYESVSKYMGSTTTLDKKTVRGDVDLMVESEKLGARTEPVSGRKIIFLPTVGEDAIQRYILKEKDSKKATFTDVIHDTEIYFFDQTEKNRFHRGKKSVERIRKFQNRQKNAKIKASDDAISKKSTSVNVSDGKIKRRDKKVSAGRTTVVVENTKEDKTVYHAGTKDGVQALIRAVVVTKSIKNEIMWDKITKLFPNNSLDNLKKKWTARRVRMGHSGWRAYVDKWKKMLVLAIKSEKISLRDVEELDLIKLLDIWTSFDEKEIKRPLFLYKNYEENRKKFTLVRDDTLTHSGNDLAMSSMIQREISSLKKTYTRKISASTKDLSKSQSDDYIRTVIRSILIESPSTTRNEIEALKNVGNESIDNVIMDMAKEKQIYLHGSKLECTDTLPDILENRGNYKDFGVAFQYRCKVNELLEAGNAIVINQEPSDISSWVLIDLISGELLNMDVIPMVRNVRPLTYTSRRFEIRTLTPPLIIYANSQTKLNTARKSAVKVPLGKPFSRLWVNGSGSIRPNIWKQVVTMVVNEIIFHPGITLSRLQSRCREVLSLHEISEICKWLLERQVLITTDFDGYWVNHNWYSIYEST*"  # noqa: E501
        )

    def test_codon_iteration_negative_sequence_chunk_negative_gene(self):
        """
        The gene gatR identified issues with codon iteration. This test shows that this gene
        can be properly iterated over on the 2nd exon when sliced into a chunk that contains only the 2nd
        exon.

        The issue was identified as a problem with negative strand sequence chunks. On negative strand chunks,
        the strand of the chunk-relative Location will be resolved to be relative to the sequence, which means
        that it is a PLUS strand for a MINUS gene.

        The updated frame iterator takes this into account.

        This test is a sequence chunk that fully contains the first exon of the E. coli NC_000913.3 gene gatR.
        """
        cds = CDSInterval(
            [2169692, 2171428],
            [2170167, 2171727],
            Strand.MINUS,
            [CDSFrame.TWO, CDSFrame.ZERO],
            parent_or_seq_chunk_parent=Parent(
                id="NC_000913.3:2171278-2171877",
                sequence=Sequence(
                    "ATCGCTCGTAATGCTATGCCGGGCAAAGTGTTGCTCATTCCCTGAAACCGCGGGCCAGCGTGATGCTGGCCCGGTATTGTGCAAAACAGATCATTCACCAATGGTCCCCCTTCGTTTACACTAGCCACAATTGAATATGGGTAAATGACGATGAATTCATTCGAGCGAAGGAATAAGATCATCCAATTAGTGAATGAACAGGGAACCGTGCTTGTTCAGGATCTGGCGGGAGTATTTGCTGCCTCGGAAGCGACAATCCGTGCCGATTTGCGCTTTCTCGAACAAAAAGGCGTGGTTACGCGCTTTCATGGCGGTGCGGCGAAAATAATGTCTGGTAATAGTGAAACCGAGACCCAGGAAGTCGGGTTTAAAGAGCGATTTCAGCTCGCCAGCGCGCCAAAAAACAGAATAGCGCAGGCGGCAGTCAAAATGATCCACGAAGGGATGACTGATCCTACCCACGTAATATGGACACAGGCCTAAGCGAGGTTCTTGTTTTCAAATTGTTCCGGACTGAGGCCGCCACACCAACTGTGCCGCCGCCACCGATTGTAATCACATTCGATATAATTAAACACCGTTGCCCGCATTATTTCCCG",  # noqa: E501
                    Alphabet.NT_EXTENDED_GAPPED,
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            2171278,
                            2171877,
                            Strand.MINUS,
                            parent=Parent(id="chrI", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
        )
        cds = CDSInterval.from_dict(cds.to_dict(), parent_or_seq_chunk_parent=cds.chunk_relative_location.parent)
        # Genomic relative strand is MINUS
        assert cds.strand == Strand.MINUS
        # chunk relative strand is PLUS
        assert cds.chunk_relative_location.strand == Strand.PLUS
        translation = cds.translate()
        assert (
            str(translation)
            == "MNSFERRNKIIQLVNEQGTVLVQDLAGVFAASEATIRADLRFLEQKGVVTRFHGGAAKIMSGNSETETQEVGFKERFQLASAPKNRIAQAAVKMIHEGM"
        )

    def test_codon_iteration_negative_sequence_chunk_positive_gene(self):
        """
        This test is a sequence chunk that contains the first exon of the E. coli NC_000913.3 gene crl.
        """
        cds = CDSInterval(
            [257828, 258675],
            [257899, 259006],
            Strand.PLUS,
            [CDSFrame.ZERO, CDSFrame.TWO],
            parent_or_seq_chunk_parent=Parent(
                id="NC_000913.3:257820-257905",
                sequence=Sequence(
                    "CCTTCACGAATATACGGGCCTAGTGCGGTAAATTTTTTGATCAATCTGCTCTTCGGGTGTCCACTCGGTAACGTCATTGCTATCT",
                    Alphabet.NT_EXTENDED_GAPPED,
                    type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            257820,
                            257905,
                            Strand.MINUS,
                            parent=Parent(id="chrI", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
            ),
        )
        cds = CDSInterval.from_dict(cds.to_dict(), parent_or_seq_chunk_parent=cds.chunk_relative_location.parent)
        # Genomic relative strand is PLUS
        assert cds.strand == Strand.PLUS
        # chunk relative strand is MINUS
        assert cds.chunk_relative_location.strand == Strand.MINUS
        translation = cds.translate()
        assert str(translation) == "MTLPSGHPKSRLIKKFTALGPYI"

    @pytest.mark.parametrize(
        "start,end,cds_start,cds_end,frame,strand,expected",
        [
            # this chunk is beyond the CDS on both sides, so the translation is full length
            (0, 27, 3, 24, CDSFrame.ZERO, Strand.PLUS, "MPGFHP*"),
            # this chunk exactly bounds the CDS
            (3, 27, 3, 24, CDSFrame.ZERO, Strand.PLUS, "MPGFHP*"),
            # this chunk cuts the first off codon exactly
            (6, 27, 3, 24, CDSFrame.ZERO, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the 1st base of the 1st codon, so the translation starts from the 2nd codon
            (4, 27, 3, 24, CDSFrame.ZERO, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the 2nd base of the 1st codon, so the translation starts from the 2nd codon
            (5, 27, 3, 24, CDSFrame.ZERO, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the last base of the stop codon as well as the 2nd base of the start
            (5, 23, 3, 24, CDSFrame.ZERO, Strand.PLUS, "PGFHP"),
            # this chunk is beyond the CDS on both sides, so the translation is full length
            (0, 27, 3, 18, CDSFrame.ZERO, Strand.MINUS, "MKPGH"),
            # this chunk exactly bounds the CDS
            (3, 18, 3, 18, CDSFrame.ZERO, Strand.MINUS, "MKPGH"),
            # this chunk cuts the first off codon exactly
            (3, 15, 3, 18, CDSFrame.ZERO, Strand.MINUS, "KPGH"),
            # this chunk cuts off the 1st base of the 1st codon, so the translation starts from the 2nd codon
            (3, 17, 3, 18, CDSFrame.ZERO, Strand.MINUS, "KPGH"),
            # this chunk cuts off the 2nd base of the 1st codon, so the translation starts from the 2nd codon
            (3, 16, 3, 18, CDSFrame.ZERO, Strand.MINUS, "KPGH"),
            # this chunk cuts off the last base of the last codon as well as the 2nd base of the start
            (4, 16, 3, 18, CDSFrame.ZERO, Strand.MINUS, "KPG"),
            # this chunk is beyond the CDS on both sides, so the translation is full length
            (0, 27, 2, 24, CDSFrame.ONE, Strand.PLUS, "MPGFHP*"),
            # this chunk exactly bounds the CDS
            (2, 27, 2, 24, CDSFrame.ONE, Strand.PLUS, "MPGFHP*"),
            # this chunk cuts the first off codon exactly
            (5, 27, 2, 24, CDSFrame.ONE, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the 1st base of the 1st codon, so the translation starts from the 2nd codon
            (4, 27, 2, 24, CDSFrame.ONE, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the 2nd base of the 1st codon, so the translation starts from the 2nd codon
            (5, 27, 2, 24, CDSFrame.ONE, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the last base of the stop codon as well as the 2nd base of the start
            (5, 23, 2, 24, CDSFrame.ONE, Strand.PLUS, "PGFHP"),
            # this chunk is beyond the CDS on both sides, so the translation is full length
            (0, 27, 3, 19, CDSFrame.ONE, Strand.MINUS, "MKPGH"),
            # this chunk exactly bounds the CDS
            (3, 19, 3, 19, CDSFrame.ONE, Strand.MINUS, "MKPGH"),
            # this chunk cuts the first off codon exactly
            (3, 16, 3, 19, CDSFrame.ONE, Strand.MINUS, "KPGH"),
            # this chunk cuts off the 1st base of the 1st codon, so the translation starts from the 2nd codon
            (3, 17, 3, 19, CDSFrame.ONE, Strand.MINUS, "KPGH"),
            # this chunk cuts off the 2nd base of the 1st codon, so the translation starts from the 2nd codon
            (3, 16, 3, 19, CDSFrame.ONE, Strand.MINUS, "KPGH"),
            # this chunk cuts off the last base of the last codon as well as the 2nd base of the start
            (4, 16, 3, 19, CDSFrame.ONE, Strand.MINUS, "KPG"),
            # this chunk is beyond the CDS on both sides, so the translation is full length
            (0, 27, 1, 24, CDSFrame.TWO, Strand.PLUS, "MPGFHP*"),
            # this chunk exactly bounds the CDS
            (1, 27, 1, 24, CDSFrame.TWO, Strand.PLUS, "MPGFHP*"),
            # this chunk cuts the first off codon exactly
            (4, 27, 1, 24, CDSFrame.TWO, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the 1st base of the 1st codon, so the translation starts from the 2nd codon
            (3, 27, 1, 24, CDSFrame.TWO, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the 2nd base of the 1st codon, so the translation starts from the 2nd codon
            (4, 27, 1, 24, CDSFrame.TWO, Strand.PLUS, "PGFHP*"),
            # this chunk cuts off the last base of the stop codon as well as the 2nd base of the start
            (4, 23, 1, 24, CDSFrame.TWO, Strand.PLUS, "PGFHP"),
            # this chunk is beyond the CDS on both sides, so the translation is full length
            (0, 27, 3, 20, CDSFrame.TWO, Strand.MINUS, "MKPGH"),
            # this chunk exactly bounds the CDS
            (3, 20, 3, 20, CDSFrame.TWO, Strand.MINUS, "MKPGH"),
            # this chunk cuts the first off codon exactly
            (3, 17, 3, 20, CDSFrame.TWO, Strand.MINUS, "KPGH"),
            # this chunk cuts off the 1st base of the 1st codon, so the translation starts from the 2nd codon
            (3, 19, 3, 20, CDSFrame.TWO, Strand.MINUS, "KPGH"),
            # this chunk cuts off the 2nd base of the 1st codon, so the translation starts from the 2nd codon
            (3, 18, 3, 20, CDSFrame.TWO, Strand.MINUS, "KPGH"),
            # this chunk cuts off the last base of the last codon as well as the 2nd base of the start
            (4, 18, 3, 20, CDSFrame.TWO, Strand.MINUS, "KPG"),
        ],
    )
    def test_single_exon_chunk_relative_translation(self, start, end, cds_start, cds_end, frame, strand, expected):
        full_chunk = "AAAATGCCCGGGTTTCATCCCTGACCC"
        chunk_offset = 100
        parent = Parent(
            sequence=Sequence(
                full_chunk[start:end],
                Alphabet.NT_EXTENDED_GAPPED,
                type=SequenceType.SEQUENCE_CHUNK,
                parent=Parent(
                    location=SingleInterval(
                        chunk_offset + start,
                        chunk_offset + end,
                        Strand.PLUS,
                        parent=Parent(sequence_type=SequenceType.CHROMOSOME),
                    )
                ),
            ),
        )
        cds = CDSInterval(
            [cds_start + chunk_offset],
            [cds_end + chunk_offset],
            strand,
            [frame],
            parent_or_seq_chunk_parent=parent,
        )
        assert str(cds.translate()) == expected

    @pytest.mark.parametrize(
        "cds,cleaned_location,loc_on_chrom,expected_offset",
        [
            # + strand with full span chunk -> 0 offset
            (
                CDSInterval([0], [10], Strand.PLUS, [CDSFrame.ZERO]),
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(0, 10, Strand.PLUS),
                0,
            ),
            # + strand with 1bp cutoff from 5' end -> 2 offset
            (
                CDSInterval([0], [10], Strand.PLUS, [CDSFrame.ZERO]),
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(1, 10, Strand.PLUS),
                2,
            ),
            # + strand with 2bp cutoff from 5' end -> 1 offset
            (
                CDSInterval([0], [10], Strand.PLUS, [CDSFrame.ZERO]),
                SingleInterval(0, 10, Strand.PLUS),
                SingleInterval(2, 10, Strand.PLUS),
                1,
            ),
            # - strand with full span chunk -> 0 offset
            (
                CDSInterval([0], [10], Strand.MINUS, [CDSFrame.ZERO]),
                SingleInterval(0, 10, Strand.MINUS),
                SingleInterval(0, 10, Strand.MINUS),
                0,
            ),
            # - strand with 1bp cutoff from 5' end -> 2 offset
            (
                CDSInterval([0], [10], Strand.MINUS, [CDSFrame.ZERO]),
                SingleInterval(0, 10, Strand.MINUS),
                SingleInterval(0, 9, Strand.MINUS),
                2,
            ),
            # - strand with 2bp cutoff from 5' end -> 1 offset
            (
                CDSInterval([0], [10], Strand.MINUS, [CDSFrame.ZERO]),
                SingleInterval(0, 10, Strand.MINUS),
                SingleInterval(0, 8, Strand.MINUS),
                1,
            ),
        ],
    )
    def test__calculate_frame_offset(self, cds, cleaned_location, loc_on_chrom, expected_offset):
        assert cds._calculate_frame_offset(cleaned_location, loc_on_chrom) == expected_offset

    @pytest.mark.parametrize(
        "cds,expected",
        [
            (
                CDSInterval([0], [10], Strand.MINUS, [CDSFrame.ZERO]),
                ["None\tBioCantor\tCDS\t1\t10\t.\t-\t0\tID=072cb87f-e347-8702-ada3-20b519aa31e0-1"],
            )
        ],
    )
    def test_to_gff(self, cds, expected):
        assert [str(x) for x in cds.to_gff()] == expected


@pytest.mark.parametrize(
    "sequence,translation_table,expected",
    [
        # normal translation
        ("ATGCATATTTGGAAACCAA", TranslationTable.DEFAULT, "MHIWKP"),
        # CTG start codon
        ("CTGCATATTTGGAAACCAA", TranslationTable.STANDARD, "MHIWKP"),
        # GTG start codon (note that in-frame ATT does not change to M)
        ("CTGCATATTTGGAAACCAA", TranslationTable.STANDARD, "MHIWKP"),
        # invalid start codon is still translated
        ("AAGCATATTTGGAAACCAA", TranslationTable.PROKARYOTE, "KHIWKP"),
    ],
)
def test_translation_tables(sequence, translation_table, expected):
    s = Sequence(sequence, Alphabet.NT_EXTENDED, type=SequenceType.CHROMOSOME)
    i = SingleInterval(0, len(sequence), Strand.PLUS, parent=s)
    cds = CDSInterval.from_location(i, cds_frames=[CDSFrame.ZERO])
    assert str(cds.translate(translation_table=translation_table)) == expected
