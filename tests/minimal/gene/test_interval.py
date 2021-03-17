"""
Test static methods on the base classes in biocantor.gene.interval
"""
import pytest
from inscripta.biocantor.exc import (
    NoSuchAncestorException,
    NullSequenceException,
    LocationOverlapException,
    MismatchedParentException,
)
from inscripta.biocantor.gene.interval import AbstractInterval
from inscripta.biocantor.location.location_impl import SingleInterval, Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.sequence import SequenceType, Sequence, Alphabet


class TestAbstractInterval:
    @pytest.mark.parametrize(
        "location,parent,exp",
        [
            # no-op
            (
                SingleInterval(0, 10, Strand.PLUS),
                None,
                SingleInterval(0, 10, Strand.PLUS),
            ),
            # non-specified parent
            (SingleInterval(0, 10, Strand.PLUS), Parent(), SingleInterval(0, 10, Strand.PLUS, parent=Parent())),
            # specified parent that is non-standard
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(sequence_type="nonstandard"),
                SingleInterval(0, 10, Strand.PLUS, parent=Parent(sequence_type="nonstandard")),
            ),
            # chromosome parent
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(sequence_type=SequenceType.CHROMOSOME),
                SingleInterval(0, 10, Strand.PLUS, parent=Parent(sequence_type=SequenceType.CHROMOSOME)),
            ),
            # chromosome parent with real sequence
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(
                    sequence_type=SequenceType.CHROMOSOME, sequence=Sequence("ATACGATCAAA", Alphabet.NT_EXTENDED_GAPPED)
                ),
                SingleInterval(
                    0,
                    10,
                    Strand.PLUS,
                    parent=Parent(
                        sequence_type=SequenceType.CHROMOSOME,
                        sequence=Sequence("ATACGATCAAA", Alphabet.NT_EXTENDED_GAPPED),
                    ),
                ),
            ),
            # chromosome parent as would be built by seq_to_parent
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(
                    sequence=Sequence(
                        "ATACGATCAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME, id="test"
                    ),
                    location=SingleInterval(0, 11, Strand.PLUS),
                ),
                SingleInterval(
                    0,
                    10,
                    Strand.PLUS,
                    Parent(
                        sequence=Sequence(
                            "ATACGATCAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME, id="test"
                        ),
                        location=SingleInterval(0, 11, Strand.PLUS),
                    ),
                ),
            ),
            # chunk parent as would be built by seq_chunk_to_parent
            (
                SingleInterval(2, 10, Strand.PLUS),
                Parent(
                    id="test:1-11",
                    sequence=Sequence(
                        "TACGATCAAA",
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="test:1-11",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                1,
                                11,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                # chunk-relative location is 1-9
                SingleInterval(
                    1,
                    9,
                    Strand.PLUS,
                    Parent(
                        id="test:1-11",
                        sequence=Sequence(
                            "TACGATCAAA",
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="test:1-11",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                location=SingleInterval(
                                    1,
                                    11,
                                    Strand.PLUS,
                                    parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                                )
                            ),
                        ),
                    ),
                ),
            ),
        ],
    )
    def test_liftover_location_to_seq_chunk_parent(self, location, parent, exp):
        obs = AbstractInterval.liftover_location_to_seq_chunk_parent(location, parent)
        assert obs == exp

    @pytest.mark.parametrize(
        "location,parent,exception",
        [
            # cannot have a naked chunk parent
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(sequence_type=SequenceType.SEQUENCE_CHUNK),
                NoSuchAncestorException,
            ),
            # cannot have a sequence-less chunk parent
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(
                    id="test:1-11",
                    sequence_type=SequenceType.SEQUENCE_CHUNK,
                    parent=Parent(
                        location=SingleInterval(
                            1,
                            11,
                            Strand.PLUS,
                            parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                        )
                    ),
                ),
                NullSequenceException,
            ),
            # chunk location must overlap genomic location
            (
                SingleInterval(0, 10, Strand.PLUS),
                Parent(
                    id="test:11-20",
                    sequence=Sequence(
                        "ATACGATCA",
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="test:11-20",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                11,
                                20,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                LocationOverlapException,
            ),
            # if location is on a chunk, must be a proper chunk with a chromosome
            (
                SingleInterval(0, 10, Strand.PLUS, parent=Parent(sequence_type=SequenceType.SEQUENCE_CHUNK)),
                Parent(),
                NoSuchAncestorException,
            ),
            # if location is on a chunk, that chunk's chromosome must match the new parent
            (
                SingleInterval(
                    5,
                    8,
                    Strand.PLUS,
                    parent=Parent(
                        id="test:0-9",
                        sequence=Sequence(
                            "ATACGATCA",
                            Alphabet.NT_EXTENDED_GAPPED,
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
                Parent(
                    id="wrongchrom:0-9",
                    sequence=Sequence(
                        "ATACGATCA",
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="wrongchrom:0-9",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            location=SingleInterval(
                                0,
                                9,
                                Strand.PLUS,
                                parent=Parent(id="wrongchrom", sequence_type=SequenceType.CHROMOSOME),
                            )
                        ),
                    ),
                ),
                MismatchedParentException,
            ),
            # if location is on a chunk, and the chromosome parent has sequence, sequence must match
            (
                SingleInterval(
                    5,
                    8,
                    Strand.PLUS,
                    parent=Parent(
                        id="test:0-9",
                        sequence=Sequence(
                            "ATACGATCA",
                            Alphabet.NT_EXTENDED_GAPPED,
                            id="test:0-9",
                            type=SequenceType.SEQUENCE_CHUNK,
                            parent=Parent(
                                sequence=Sequence(
                                    "ATACGATCAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                                ),
                                location=SingleInterval(
                                    0,
                                    9,
                                    Strand.PLUS,
                                    parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                                ),
                            ),
                        ),
                    ),
                ),
                Parent(
                    id="test:0-9",
                    sequence=Sequence(
                        "ATACGATCA",
                        Alphabet.NT_EXTENDED_GAPPED,
                        id="test:0-9",
                        type=SequenceType.SEQUENCE_CHUNK,
                        parent=Parent(
                            sequence=Sequence(
                                "ATACGATCAATA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                            ),
                            location=SingleInterval(
                                0,
                                9,
                                Strand.PLUS,
                                parent=Parent(id="test", sequence_type=SequenceType.CHROMOSOME),
                            ),
                        ),
                    ),
                ),
                MismatchedParentException,
            ),
        ],
    )
    def test_liftover_location_to_seq_chunk_parent_exceptions(self, location, parent, exception):
        with pytest.raises(exception):
            _ = AbstractInterval.liftover_location_to_seq_chunk_parent(location, parent)
