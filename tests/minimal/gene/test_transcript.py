import pytest

from inscripta.biocantor.exc import (
    InvalidCDSIntervalError,
    EmptyLocationException,
    NullParentException,
    NoncodingTranscriptError,
    ValidationException,
    NullSequenceException,
)
from inscripta.biocantor.gene.cds_frame import CDSFrame
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.io.models import TranscriptIntervalModel
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval, EmptyLocation
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent, SequenceType
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence

# these features will be shared across all tests
genome = "GTATTCTTGGACCTAATT"
parent = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type=SequenceType.CHROMOSOME)
# offset the genome to show parent
genome2 = "AAGTATTCTTGGACCTAATT"
parent_genome2 = Parent(sequence=Sequence(genome2, Alphabet.NT_STRICT), sequence_type=SequenceType.CHROMOSOME)

# slice the genome down to contain some of the transcripts
parent_genome2_1_15 = Parent(
    sequence=Sequence(
        genome2[1:15],
        Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            location=SingleInterval(
                1, 15, Strand.PLUS, parent=Parent(id="genome_1_15", sequence_type=SequenceType.CHROMOSOME)
            )
        ),
    )
)

parent_no_seq = Parent(sequence_type=SequenceType.CHROMOSOME)
parent_nonstandard_type = Parent(sequence_type="SomeOtherType")
parent_nonstandard_type_with_sequence = Parent(
    sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type="SomeOtherType"
)

# Integer transcript definitions
# a single exon transcript that is this entire genome
se_unspliced = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[0],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[0],
        cds_ends=[18],
        cds_frames=[CDSFrame.ZERO.name],
    )
)
se_unspliced_repr = "TranscriptInterval((0-18:+), cds=[CDS((0-18:+), (CDSFrame.ZERO)], symbol=None)"
# a single exon transcript that is out of frame
se_unspliced_oof = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[0],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[0],
        cds_ends=[18],
        cds_frames=[CDSFrame.ONE.name],
    )
)
se_unspliced_oof_repr = "TranscriptInterval((0-18:+), cds=[CDS((0-18:+), (CDSFrame.ONE)], symbol=None)"
# a single exon noncoding transcript
se_noncoding = TranscriptIntervalModel.Schema().load(dict(exon_starts=[0], exon_ends=[18], strand=Strand.PLUS.name))
se_noncoding_repr = "TranscriptInterval((0-18:+), cds=[None], symbol=None)"
# a three exon transcript with a 2bp 5' UTR and a 2bp 3' UTR
e3_spliced = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
    )
)
e3_spliced_repr = (
    "TranscriptInterval((2-6:+, 7-10:+, 12-15:+), "
    "cds=[CDS((4-6:+, 7-10:+, 12-13:+), "
    "(CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.TWO)], symbol=None)"
)
# a three exon transcript with an entirely non-coding 1st and last exon
e3_spliced_utr = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[7],
        cds_ends=[10],
        cds_frames=[CDSFrame.ZERO.name],
    )
)
e3_spliced_utr_repr = "TranscriptInterval((2-6:+, 7-10:+, 12-15:+), cds=[CDS((7-10:+), (CDSFrame.ZERO)], symbol=None)"
# a three exon transcript with an entirely non-coding 1st and last exon; out of frame
# this means it has no translation
e3_spliced_notrans = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[7],
        cds_ends=[10],
        cds_frames=[CDSFrame.ONE.name],
    )
)
e3_spliced_notrans_repr = (
    "TranscriptInterval((2-6:+, 7-10:+, 12-15:+), cds=[CDS((7-10:+), (CDSFrame.ONE)], symbol=None)"
)

# Negative strand transcript definitions
# a single exon transcript that is this entire genome
se_unspliced_minus = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[0],
        exon_ends=[18],
        strand=Strand.MINUS.name,
        cds_starts=[0],
        cds_ends=[18],
        cds_frames=[CDSFrame.ZERO.name],
    )
)
se_unspliced_repr_minus = "TranscriptInterval((0-18:-), cds=[CDS((0-18:-), (CDSFrame.ZERO)], symbol=None)"
# a single exon transcript that is out of frame
se_unspliced_oof_minus = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[0],
        exon_ends=[18],
        strand=Strand.MINUS.name,
        cds_starts=[0],
        cds_ends=[18],
        cds_frames=[CDSFrame.ONE.name],
    )
)
se_unspliced_oof_repr_minus = "TranscriptInterval((0-18:-), cds=[CDS((0-18:-), (CDSFrame.ONE)], symbol=None)"
# a single exon noncoding transcript (explicitly defined)
se_noncoding_minus = TranscriptIntervalModel.Schema().load(
    dict(exon_starts=[0], exon_ends=[18], strand=Strand.MINUS.name)
)
se_noncoding_repr_minus = "TranscriptInterval((0-18:-), cds=[None], symbol=None)"
# a three exon transcript with a 2bp 5' UTR and a 2bp 3' UTR
e3_spliced_minus = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.MINUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ONE.name, CDSFrame.ONE.name, CDSFrame.ZERO.name],
    )
)
e3_spliced_repr_minus = (
    "TranscriptInterval((2-6:-, 7-10:-, 12-15:-), "
    "cds=[CDS((4-6:-, 7-10:-, 12-13:-), "
    "(CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ZERO)], symbol=None)"
)
# a three exon transcript with an entirely non-coding 1st and last exon
e3_spliced_utr_minus = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.MINUS.name,
        cds_starts=[7],
        cds_ends=[10],
        cds_frames=[CDSFrame.ZERO.name],
    )
)
e3_spliced_utr_repr_minus = (
    "TranscriptInterval((2-6:-, 7-10:-, 12-15:-), cds=[CDS((7-10:-), (CDSFrame.ZERO)], symbol=None)"
)
# a three exon transcript with an entirely non-coding 1st and last exon; out of frame
# this means it has no translation
e3_spliced_notrans_minus = TranscriptIntervalModel.Schema().load(
    dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.MINUS.name,
        cds_starts=[7],
        cds_ends=[10],
        cds_frames=[CDSFrame.ONE.name],
    )
)
e3_spliced_notrans_repr_minus = (
    "TranscriptInterval((2-6:-, 7-10:-, 12-15:-), cds=[CDS((7-10:-), (CDSFrame.ONE)], symbol=None)"
)


class TestTranscript:
    """Test constructing transcripts of various types"""

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (se_unspliced, se_unspliced_repr),
            (se_unspliced_oof, se_unspliced_oof_repr),
            (se_noncoding, se_noncoding_repr),
            (e3_spliced, e3_spliced_repr),
            (e3_spliced_utr, e3_spliced_utr_repr),
            (e3_spliced_notrans, e3_spliced_notrans_repr),
            (se_unspliced_minus, se_unspliced_repr_minus),
            (se_unspliced_oof_minus, se_unspliced_oof_repr_minus),
            (se_noncoding_minus, se_noncoding_repr_minus),
            (e3_spliced_minus, e3_spliced_repr_minus),
            (e3_spliced_utr_minus, e3_spliced_utr_repr_minus),
            (e3_spliced_notrans_minus, e3_spliced_notrans_repr_minus),
        ],
    )
    def test_transcript_constructor(self, schema, expected):
        assert str(schema.to_transcript_interval()) == expected
        assert str(schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)) == expected

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (se_unspliced, "VFLDLI"),
            (se_unspliced_oof, "YSWT*"),
            (e3_spliced, "SG"),
            (e3_spliced_utr, "W"),
            (e3_spliced_notrans, ""),
            (se_unspliced_minus, "N*VQEY"),
            (se_unspliced_oof_minus, "IRSKN"),
            (e3_spliced_minus, "AR"),
            (e3_spliced_utr_minus, "P"),
            (e3_spliced_notrans_minus, ""),
        ],
    )
    def test_translations(self, schema, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_protein_sequence()) == expected

    @pytest.mark.parametrize("schema,expected", [(se_unspliced_minus, "N*")])
    def test_truncations(self, schema, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_protein_sequence(truncate_at_in_frame_stop=True)) == expected

    @pytest.mark.parametrize("schema,expected_exception", [(se_unspliced_minus, NullParentException)])
    def test_translation_exceptions(self, schema, expected_exception):
        tx = schema.to_transcript_interval()
        with pytest.raises(expected_exception):
            assert str(tx.get_protein_sequence(truncate_at_in_frame_stop=True)) == expected_exception

    @pytest.mark.parametrize(
        "schema,expected_exception",
        [
            (
                TranscriptIntervalModel.Schema().load(
                    dict(exon_starts=[0], exon_ends=[18], strand=Strand.MINUS.name, cds_starts=[10], cds_ends=[0])
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(exon_starts=[0], exon_ends=[18], strand=Strand.MINUS.name, cds_starts=[0], cds_ends=[10])
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(exon_starts=[10], exon_ends=[18], strand=Strand.MINUS.name, cds_starts=[0], cds_ends=[10])
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(exon_starts=[10], exon_ends=[18], strand=Strand.MINUS.name, cds_starts=[0], cds_ends=[20])
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(exon_starts=[10], exon_ends=[18], strand=Strand.MINUS.name, cds_starts=[10], cds_ends=[20])
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.MINUS.name,
                        cds_starts=[0],
                        cds_ends=[10],
                        cds_frames=[CDSFrame.ZERO.name, CDSFrame.ZERO.name],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.MINUS.name,
                        cds_starts=[0],
                        cds_ends=[10],
                        cds_frames=None,
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.MINUS.name,
                        cds_starts=[0],
                        cds_ends=[10],
                        cds_frames=None,
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.MINUS.name,
                        cds_starts=[0],
                        cds_ends=[10],
                        cds_frames=None,
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=[0],
                        cds_frames=[CDSFrame.ZERO.name],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=[5],
                        cds_frames=[],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=None,
                        cds_ends=[0],
                        cds_frames=[CDSFrame.ZERO.name],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=None,
                        cds_frames=[CDSFrame.ZERO.name],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=[0, 5],
                        cds_frames=[CDSFrame.ZERO.name],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=[0],
                        cds_frames=[CDSFrame.ZERO.name, CDSFrame.ZERO.name],
                    )
                ),
                InvalidCDSIntervalError,
            ),
            (
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[5, 18],
                        strand=Strand.MINUS.name,
                        cds_ends=[0],
                    )
                ),
                ValidationException,
            ),
        ],
    )
    def test_transcript_from_transcript_exceptions(self, schema, expected_exception):
        with pytest.raises(expected_exception):
            _ = schema.to_transcript_interval()
            _ = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, 2, 0),
            (e3_spliced_utr, 7, 4),
            (e3_spliced_utr, 14, 9),
            (e3_spliced_utr_minus, 2, 9),
            (e3_spliced_utr_minus, 7, 5),
            (e3_spliced_utr_minus, 14, 0),
        ],
    )
    def test_sequence_pos_to_transcript(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.sequence_pos_to_transcript(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, (7, 13, Strand.PLUS), SingleInterval(4, 8, Strand.PLUS)),
            (e3_spliced_utr_minus, (7, 13, Strand.PLUS), (SingleInterval(2, 6, Strand.MINUS))),
        ],
    )
    def test_sequence_interval_to_transcript(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.sequence_interval_to_transcript(*value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, 0, 2),
            (e3_spliced_utr, 9, 14),
            (e3_spliced_utr, 4, 7),
            (e3_spliced_utr_minus, 0, 14),
            (e3_spliced_utr_minus, 9, 2),
            (e3_spliced_utr_minus, 5, 7),
        ],
    )
    def test_transcript_pos_to_sequence(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.transcript_pos_to_sequence(value) == expected
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.transcript_pos_to_sequence(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, (0, 5, Strand.PLUS), CompoundInterval([2, 7], [6, 8], Strand.PLUS)),
            (e3_spliced_utr_minus, (0, 5, Strand.PLUS), CompoundInterval([8, 12], [10, 15], Strand.MINUS)),
        ],
    )
    def test_transcript_interval_to_sequence(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.transcript_interval_to_sequence(*value).reset_parent(None) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, (0, 5, Strand.PLUS), CompoundInterval([2, 7], [6, 8], Strand.PLUS, parent=parent)),
            (
                e3_spliced_utr_minus,
                (0, 5, Strand.PLUS),
                CompoundInterval([8, 12], [10, 15], Strand.MINUS, parent=parent),
            ),
        ],
    )
    def test_transcript_interval_to_sequence_parent(self, schema, value, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.transcript_interval_to_sequence(*value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, 0, 7),
            (e3_spliced_utr, 1, 8),
            (e3_spliced_utr, 2, 9),
            (e3_spliced_utr_minus, 0, 9),
            (e3_spliced_utr_minus, 1, 8),
            (e3_spliced_utr_minus, 2, 7),
        ],
    )
    def test_cds_pos_to_sequence(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.cds_pos_to_sequence(value) == expected
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.cds_pos_to_sequence(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (
                e3_spliced_utr,
                (0, 3, Strand.PLUS),
                SingleInterval(7, 10, Strand.PLUS),
            ),
            (
                e3_spliced_utr_minus,
                (0, 3, Strand.PLUS),
                SingleInterval(7, 10, Strand.MINUS),
            ),
        ],
    )
    def test_cds_interval_to_sequence(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.cds_interval_to_sequence(*value).reset_parent(None) == expected
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.cds_interval_to_sequence(*value).reset_parent(None) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, 7, 0),
            (e3_spliced_utr, 8, 1),
            (e3_spliced_utr, 9, 2),
            (e3_spliced_utr_minus, 7, 2),
            (e3_spliced_utr_minus, 8, 1),
            (e3_spliced_utr_minus, 9, 0),
        ],
    )
    def test_sequence_pos_to_cds(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.sequence_pos_to_cds(value) == expected
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.sequence_pos_to_cds(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, (7, 10, Strand.PLUS), SingleInterval(0, 3, Strand.PLUS)),
            (
                e3_spliced_utr_minus,
                (7, 10, Strand.MINUS),
                SingleInterval(0, 3, Strand.PLUS),
            ),
        ],
    )
    def test_sequence_interval_to_cds(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.sequence_interval_to_cds(*value).reset_parent(None) == expected
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.sequence_interval_to_cds(*value).reset_parent(None) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, 0, 4),
            (e3_spliced_utr, 1, 5),
            (e3_spliced_utr, 2, 6),
            (e3_spliced_utr_minus, 0, 3),
            (e3_spliced_utr_minus, 1, 4),
            (e3_spliced_utr_minus, 2, 5),
        ],
    )
    def test_cds_pos_to_transcript(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.cds_pos_to_transcript(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced_utr, 4, 0),
            (e3_spliced_utr, 5, 1),
            (e3_spliced_utr, 6, 2),
            (e3_spliced_utr_minus, 3, 0),
            (e3_spliced_utr_minus, 4, 1),
            (e3_spliced_utr_minus, 5, 2),
        ],
    )
    def test_transcript_pos_to_cds(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        assert tx.transcript_pos_to_cds(value) == expected
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert tx.transcript_pos_to_cds(value) == expected

    def test_position_exceptions(self):
        tx = se_noncoding.to_transcript_interval()
        with pytest.raises(NoncodingTranscriptError):
            _ = tx.cds_pos_to_sequence(0)

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.cds_interval_to_sequence(0, 10, Strand.PLUS)

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.sequence_pos_to_cds(0)

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.sequence_interval_to_cds(0, 10, Strand.PLUS)

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.cds_pos_to_transcript(0)

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.transcript_pos_to_cds(0)

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.cds_location

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.get_5p_interval()

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.get_3p_interval()

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.get_cds_sequence()

        with pytest.raises(NoncodingTranscriptError):
            _ = tx.get_protein_sequence()

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (e3_spliced, SingleInterval(2, 4, Strand.PLUS)),
            (e3_spliced_utr, SingleInterval(2, 6, Strand.PLUS)),
            (se_unspliced, EmptyLocation()),
        ],
    )
    def test_get_5p_interval(self, schema, expected):
        tx = schema.to_transcript_interval()
        assert tx.get_5p_interval() == expected

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (e3_spliced, SingleInterval(13, 15, Strand.PLUS)),
            (e3_spliced_utr, SingleInterval(12, 15, Strand.PLUS)),
            (se_unspliced, EmptyLocation()),
        ],
    )
    def test_get_3p_interval(self, schema, expected):
        tx = schema.to_transcript_interval()
        assert tx.get_3p_interval() == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (
                e3_spliced_utr,
                SingleInterval(0, 10, Strand.PLUS, parent=parent),
                dict(
                    exon_starts=(2, 7),
                    exon_ends=(6, 10),
                    strand=Strand.PLUS.name,
                ),
            ),
            (
                e3_spliced_utr_minus,
                SingleInterval(0, 10, Strand.PLUS, parent=parent),
                dict(
                    exon_starts=(2, 7),
                    exon_ends=(6, 10),
                    strand=Strand.MINUS.name,
                ),
            ),
            (
                e3_spliced_utr_minus,
                SingleInterval(0, 5, Strand.PLUS, parent=parent),
                dict(exon_starts=(2,), exon_ends=(5,), strand=Strand.MINUS.name),
            ),
        ],
    )
    def test_intersection(self, schema, value, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        val = tx.intersect(value)
        expected = (
            TranscriptIntervalModel.Schema().load(expected).to_transcript_interval(parent_or_seq_chunk_parent=parent)
        )
        assert str(expected) == str(val)

    def test_intersection_exception(self):
        schema = TranscriptIntervalModel.Schema().load(dict(exon_starts=[10], exon_ends=[15], strand=Strand.MINUS.name))
        tx = schema.to_transcript_interval()
        s = SingleInterval(0, 5, Strand.PLUS)
        with pytest.raises(EmptyLocationException):
            _ = tx.intersect(s)

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (
                e3_spliced_utr,
                SingleInterval(0, 10, Strand.PLUS),
                dict(
                    exon_starts=(2, 7),
                    exon_ends=(6, 10),
                    strand=Strand.PLUS.name,
                ),
            ),
            (
                e3_spliced_utr_minus,
                SingleInterval(0, 10, Strand.PLUS),
                dict(
                    exon_starts=(2, 7),
                    exon_ends=(6, 10),
                    strand=Strand.MINUS.name,
                ),
            ),
        ],
    )
    def test_intersection_no_parent(self, schema, value, expected):
        tx = schema.to_transcript_interval()
        val = tx.intersect(value)
        expected = TranscriptIntervalModel.Schema().load(expected).to_transcript_interval()
        assert str(expected) == str(val)

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (
                e3_spliced,
                "IL",
            ),
            (
                e3_spliced_minus,
                "*D",
            ),
        ],
    )
    def test_reset_parent(self, schema, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        tx.reset_parent(parent_genome2)
        assert str(tx.get_protein_sequence()) == str(expected)

    def test_object_conversion(self):
        for model in [
            e3_spliced,
            e3_spliced_utr_minus,
            e3_spliced_minus,
            e3_spliced_notrans,
        ]:
            obj = model.to_transcript_interval()
            new_model = TranscriptIntervalModel.from_transcript_interval(obj)
            new_model.transcript_interval_guid = None
            assert model == new_model
            obj = model.to_transcript_interval(parent_or_seq_chunk_parent=parent)
            new_model = TranscriptIntervalModel.from_transcript_interval(obj)
            new_model.transcript_interval_guid = None
            assert model == new_model

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (e3_spliced, "ATTCTGGCTA"),
            (e3_spliced_minus, "TAGCCAGAAT"),
            (e3_spliced_notrans, "ATTCTGGCTA"),
            (e3_spliced_notrans_minus, "TAGCCAGAAT"),
        ],
    )
    def test_transcript_sequence(self, schema, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_transcript_sequence()) == str(tx.get_spliced_sequence()) == expected

    @pytest.mark.parametrize(
        "schema,expected_genomic,expected_stranded_genomic",
        [
            (e3_spliced, "ATTCTTGGACCTA", "ATTCTTGGACCTA"),
            (e3_spliced_minus, "ATTCTTGGACCTA", "TAGGTCCAAGAAT"),
            (e3_spliced_notrans, "ATTCTTGGACCTA", "ATTCTTGGACCTA"),
            (e3_spliced_notrans_minus, "ATTCTTGGACCTA", "TAGGTCCAAGAAT"),
        ],
    )
    def test_genomic_sequence(self, schema, expected_genomic, expected_stranded_genomic):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_reference_sequence()) == expected_genomic
        assert str(tx.get_genomic_sequence()) == expected_stranded_genomic

    @pytest.mark.parametrize(
        "schema,expected_cds",
        [
            (e3_spliced, "TCTGGC"),
            (e3_spliced_minus, "GCCAGA"),
            (e3_spliced_notrans, ""),
            (e3_spliced_notrans_minus, ""),
            (e3_spliced_utr, "TGG"),
            (e3_spliced_utr_minus, "CCA"),
        ],
    )
    def test_cds_sequence(self, schema, expected_cds):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_cds_sequence()) == expected_cds

    @pytest.mark.parametrize(
        "schema1,schema2,is_eq",
        [
            (e3_spliced, e3_spliced, True),
            (e3_spliced, e3_spliced_utr_minus, False),
            (
                se_unspliced,
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=[18],
                        cds_frames=[CDSFrame.ZERO.name],
                        qualifiers={"test": [1]},
                    )
                ),
                False,
            ),
            (
                se_unspliced,
                TranscriptIntervalModel.Schema().load(
                    dict(
                        exon_starts=[0],
                        exon_ends=[18],
                        strand=Strand.PLUS.name,
                        cds_starts=[0],
                        cds_ends=[18],
                        cds_frames=[CDSFrame.ZERO.name],
                    )
                ),
                True,
            ),
        ],
    )
    def test_equality(self, schema1, schema2, is_eq):
        tx1 = schema1.to_transcript_interval()
        tx2 = schema2.to_transcript_interval()
        assert (tx1 == tx2) == is_eq

    def test_cds_start_end(self):
        tx = e3_spliced.to_transcript_interval()
        assert tx.is_coding
        assert tx.cds_start == tx.chunk_relative_cds_start == 4
        assert tx.cds_end == tx.chunk_relative_cds_end == 13
        tx2 = se_noncoding.to_transcript_interval()
        assert not tx2.is_coding

        with pytest.raises(NoncodingTranscriptError):
            _ = tx2.cds_start
        with pytest.raises(NoncodingTranscriptError):
            _ = tx2.cds_end
        with pytest.raises(NoncodingTranscriptError):
            _ = tx2.chunk_relative_cds_start
        with pytest.raises(NoncodingTranscriptError):
            _ = tx2.chunk_relative_cds_end

        # empty translation but still coding
        tx3 = e3_spliced_notrans_minus.to_transcript_interval()
        assert tx3.is_coding
        assert tx3.cds_start

    def test_frameshifted(self):
        genome = "ATGGGGTGATGA"
        parent_genome = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type=SequenceType.CHROMOSOME)
        tx = dict(
            exon_starts=[0],
            exon_ends=[12],
            strand=Strand.PLUS.name,
            cds_starts=[0],
            cds_ends=[12],
            cds_frames=[CDSFrame.ZERO.name],
        )
        schema = TranscriptIntervalModel.Schema().load(tx)
        obj = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome)
        assert obj.has_in_frame_stop
        tx2 = dict(
            exon_starts=[0],
            exon_ends=[12],
            strand=Strand.PLUS.name,
            cds_starts=[0],
            cds_ends=[9],
            cds_frames=[CDSFrame.ZERO.name],
        )
        schema2 = TranscriptIntervalModel.Schema().load(tx2)
        obj2 = schema2.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome)
        assert not obj2.has_in_frame_stop

    def test_noncoding_frameshift(self):
        tx = se_noncoding.to_transcript_interval()
        with pytest.raises(NoncodingTranscriptError):
            _ = tx.has_in_frame_stop


class TestTranscriptWithoutModel:
    @pytest.mark.parametrize(
        "tx",
        [
            dict(
                exon_starts=[2, 8, 12],
                exon_ends=[5, 13, 18],
                strand=Strand.PLUS,
            ),
            dict(
                exon_starts=[2, 8, 12],
                exon_ends=[5, 13, 18],
                strand=Strand.PLUS,
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        ],
    )
    def test_constructor(self, tx):
        tx = TranscriptInterval(**tx)
        tx2 = TranscriptInterval.from_location(tx._location)
        assert tx == tx2

    @pytest.mark.parametrize(
        "location,expected",
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
                SingleInterval(222218, 222233, Strand.PLUS),
            )
        ],
    )
    def test_chunk_relative_constructor(self, location, expected):
        cds = TranscriptInterval.from_chunk_relative_location(location)
        assert cds.chromosome_location.reset_parent(None) == expected

    @pytest.mark.parametrize(
        "tx",
        [
            dict(
                exon_starts=[2, 8, 12],
                exon_ends=[5, 13, 18],
                strand=Strand.PLUS,
            ),
            dict(
                exon_starts=[2, 8, 12],
                exon_ends=[5, 13, 18],
                strand=Strand.PLUS,
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        ],
    )
    def test_dict(self, tx):
        tx = TranscriptInterval(**tx)
        tx2 = TranscriptInterval.from_dict(tx.to_dict())
        assert tx == tx2


class TestQualifiers:
    """Ensure that qualifier encoding is stable and digestable"""

    qualifiers = {"key1": ["a", "b", 2]}

    se_unspliced = TranscriptIntervalModel.Schema().load(
        dict(
            exon_starts=[0],
            exon_ends=[18],
            strand=Strand.PLUS.name,
            cds_starts=[0],
            cds_ends=[18],
            cds_frames=[CDSFrame.ZERO.name],
            qualifiers=qualifiers,
        )
    )

    def test_qualifier_import(self):
        """Importing qualifiers converts to a set of strings"""
        tx = self.se_unspliced.to_transcript_interval()
        assert tx.qualifiers == {"key1": {"2", "a", "b"}}

    def test_qualifier_export_to_list(self):
        """Exporting to a list of strings"""
        tx = self.se_unspliced.to_transcript_interval()
        # not the same as the input because the input was unsorted and not a string
        assert tx._export_qualifiers_to_list() != self.qualifiers
        assert tx._export_qualifiers_to_list() == {"key1": ["2", "a", "b"]}


class TestTranscriptIntervalSequenceSubset:
    """Test the ability to slice the genomic sequence of a feature interval and still get useful results."""

    @pytest.mark.parametrize(
        "schema,parent,expected",
        [
            (e3_spliced, parent_genome2_1_15, "GTATCTTACC"),
            (e3_spliced_minus, parent_genome2_1_15, "GGTAAGATAC"),
            (e3_spliced_notrans, parent_genome2_1_15, "GTATCTTACC"),
            (e3_spliced_notrans_minus, parent_genome2_1_15, "GGTAAGATAC"),
        ],
    )
    def test_transcript_sequence(self, schema, parent, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_transcript_sequence()) == str(tx.get_spliced_sequence()) == expected

    @pytest.mark.parametrize(
        "schema,parent,expected_genomic,expected_stranded_genomic",
        [
            (e3_spliced, parent_genome2_1_15, "GTATTCTTGGACC", "GTATTCTTGGACC"),
            (e3_spliced_minus, parent_genome2_1_15, "GTATTCTTGGACC", "GGTCCAAGAATAC"),
            (e3_spliced_notrans, parent_genome2_1_15, "GTATTCTTGGACC", "GTATTCTTGGACC"),
            (e3_spliced_notrans_minus, parent_genome2_1_15, "GTATTCTTGGACC", "GGTCCAAGAATAC"),
        ],
    )
    def test_genomic_sequence(self, schema, parent, expected_genomic, expected_stranded_genomic):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent)
        assert str(tx.get_reference_sequence()) == expected_genomic
        assert str(tx.get_genomic_sequence()) == expected_stranded_genomic

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (e3_spliced, "IL"),
            (e3_spliced_utr, "L"),
            (e3_spliced_notrans, ""),
            (e3_spliced_minus, "*D"),
            (e3_spliced_utr_minus, "K"),
            (e3_spliced_notrans_minus, ""),
        ],
    )
    def test_translations(self, schema, expected):
        tx = schema.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert str(tx.get_protein_sequence()) == expected

    def test_start_end(self):
        tx = e3_spliced.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.chunk_relative_start + 1 == tx.start
        assert tx.chunk_relative_end + 1 == tx.end
        assert tx.chunk_relative_cds_start + 1 == tx.cds_start
        assert tx.chunk_relative_cds_end + 1 == tx.cds_end

    @pytest.mark.parametrize("pos,expected", [(2, 0), (14, 9), (7, 4), (8, 5)])
    def test_sequence_pos_to_transcript(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.sequence_pos_to_transcript(pos) == expected

    @pytest.mark.parametrize("pos,expected", [(1, 0), (13, 9), (6, 4), (7, 5)])
    def test_chunk_relative_pos_to_transcript(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.chunk_relative_pos_to_transcript(pos) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (0, 10, Strand.PLUS, SingleInterval(0, 7, Strand.PLUS)),
            (0, 3, Strand.PLUS, SingleInterval(0, 1, Strand.PLUS)),
            (0, 13, Strand.PLUS, SingleInterval(0, 8, Strand.PLUS)),
        ],
    )
    def test_sequence_interval_to_transcript(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.sequence_interval_to_transcript(start, end, strand) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (0, 10, Strand.PLUS, SingleInterval(0, 7, Strand.PLUS)),
            (0, 3, Strand.PLUS, SingleInterval(0, 2, Strand.PLUS)),
            (0, 13, Strand.PLUS, SingleInterval(0, 9, Strand.PLUS)),
        ],
    )
    def test_chunk_relative_interval_to_transcript(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.chunk_relative_interval_to_transcript(start, end, strand) == expected

    @pytest.mark.parametrize("pos,expected", [(0, 2), (9, 14), (4, 7), (5, 8)])
    def test_transcript_pos_to_sequence(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.transcript_pos_to_sequence(pos) == expected

    @pytest.mark.parametrize("pos,expected", [(0, 1), (9, 13), (4, 6), (5, 7)])
    def test_transcript_pos_to_chunk_relative(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.transcript_pos_to_chunk_relative(pos) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (0, 7, Strand.PLUS, CompoundInterval([2, 7], [6, 10], Strand.PLUS)),
            (0, 1, Strand.PLUS, SingleInterval(2, 3, Strand.PLUS)),
            (0, 8, Strand.PLUS, CompoundInterval([2, 7, 12], [6, 10, 13], Strand.PLUS)),
        ],
    )
    def test_transcript_interval_to_sequence(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.transcript_interval_to_sequence(start, end, strand).reset_parent(new_parent=None) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (0, 7, Strand.PLUS, CompoundInterval([1, 6], [5, 9], Strand.PLUS)),
            (0, 1, Strand.PLUS, SingleInterval(1, 2, Strand.PLUS)),
            (0, 8, Strand.PLUS, CompoundInterval([1, 6, 11], [5, 9, 12], Strand.PLUS)),
        ],
    )
    def test_transcript_interval_to_chunk_relative(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.transcript_interval_to_chunk_relative(start, end, strand).reset_parent(new_parent=None) == expected

    @pytest.mark.parametrize("pos,expected", [(0, 7), (1, 8)])
    def test_cds_pos_to_sequence(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.cds_pos_to_sequence(pos) == expected

    @pytest.mark.parametrize("pos,expected", [(0, 6), (1, 7)])
    def test_cds_pos_to_chunk_relative(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.cds_pos_to_chunk_relative(pos) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (0, 1, Strand.PLUS, SingleInterval(7, 8, Strand.PLUS)),
            (0, 3, Strand.PLUS, SingleInterval(7, 10, Strand.PLUS)),
        ],
    )
    def test_cds_interval_to_sequence(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.cds_interval_to_sequence(start, end, strand).reset_parent(new_parent=None) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (0, 1, Strand.PLUS, SingleInterval(6, 7, Strand.PLUS)),
            (0, 3, Strand.PLUS, SingleInterval(6, 9, Strand.PLUS)),
        ],
    )
    def test_cds_interval_to_chunk_relative(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.cds_interval_to_chunk_relative(start, end, strand).reset_parent(new_parent=None) == expected

    @pytest.mark.parametrize("pos,expected", [(7, 0), (9, 2)])
    def test_sequence_pos_to_cds(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.sequence_pos_to_cds(pos) == expected

    @pytest.mark.parametrize("pos,expected", [(6, 0), (8, 2)])
    def chunk_relative_pos_to_cds(self, pos, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.chunk_relative_pos_to_cds(pos) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (7, 10, Strand.PLUS, SingleInterval(0, 3, Strand.PLUS)),
            (7, 8, Strand.PLUS, SingleInterval(0, 1, Strand.PLUS)),
        ],
    )
    def test_sequence_interval_to_cds(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.sequence_interval_to_cds(start, end, strand).reset_parent(new_parent=None) == expected

    @pytest.mark.parametrize(
        "start,end,strand,expected",
        [
            (6, 9, Strand.PLUS, SingleInterval(0, 3, Strand.PLUS)),
            (6, 7, Strand.PLUS, SingleInterval(0, 1, Strand.PLUS)),
        ],
    )
    def test_chunk_relative_interval_to_cds(self, start, end, strand, expected):
        tx = e3_spliced_utr.to_transcript_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert tx.chunk_relative_interval_to_cds(start, end, strand).reset_parent(new_parent=None) == expected

    def test_sequence_exceptions(self):
        """All sequence accessors should raise good errors when attempted without sequence info"""
        tx = e3_spliced.to_transcript_interval(parent_or_seq_chunk_parent=parent_no_seq)
        with pytest.raises(NullSequenceException):
            _ = tx.get_reference_sequence()
        with pytest.raises(NullSequenceException):
            _ = tx.get_spliced_sequence()
        with pytest.raises(NullSequenceException):
            _ = tx.get_genomic_sequence()

    def test_nonstandard_parents(self):
        tx0 = e3_spliced.to_transcript_interval(parent)
        seq0 = tx0.get_spliced_sequence()
        tx1 = e3_spliced.to_transcript_interval(parent_nonstandard_type)
        with pytest.raises(NullSequenceException):
            _ = tx1.get_spliced_sequence()
        tx2 = e3_spliced.to_transcript_interval(parent_no_seq)
        with pytest.raises(NullSequenceException):
            _ = tx2.get_spliced_sequence()
        tx3 = e3_spliced.to_transcript_interval(parent_nonstandard_type_with_sequence)
        seq = tx3.get_spliced_sequence()
        assert seq == seq0

        assert tx0.chromosome_location == tx0.chunk_relative_location
        assert tx1.chromosome_location == tx1.chunk_relative_location
        assert tx2.chromosome_location == tx2.chunk_relative_location
        assert tx3.chromosome_location == tx3.chunk_relative_location
        # OTOH, this is not the same
        tx4 = e3_spliced.to_transcript_interval(parent_genome2_1_15)
        assert tx4.chromosome_location != tx4.chunk_relative_location
