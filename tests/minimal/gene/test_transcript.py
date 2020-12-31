import pytest

from inscripta.biocantor.exc import (
    InvalidCDSIntervalError,
    LocationException,
    EmptyLocationException,
    NullParentException,
    NoncodingTranscriptError,
    NullSequenceException,
)
from inscripta.biocantor.gene.cds import CDSInterval, CDSFrame
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval, EmptyLocation
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.io.models import TranscriptIntervalModel
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence

# these features will be shared across all tests
genome = "GTATTCTTGGACCTAATT"
parent = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT))
# offset the genome to show parent
genome2 = "AAGTATTCTTGGACCTAATT"
parent_genome2 = Parent(sequence=Sequence(genome2, Alphabet.NT_STRICT))


class TestTranscript:
    """Test constructing transcripts of various types"""

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
    se_unspliced_repr = "TranscriptInterval((0-18:+), cds=[CDS((0-18:+), (CDSFrame.ZERO)], qualifiers=None)"
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
    se_unspliced_oof_repr = "TranscriptInterval((0-18:+), cds=[CDS((0-18:+), (CDSFrame.ONE)], qualifiers=None)"
    # a single exon noncoding transcript
    se_noncoding = TranscriptIntervalModel.Schema().load(dict(exon_starts=[0], exon_ends=[18], strand=Strand.PLUS.name))
    se_noncoding_repr = "TranscriptInterval((0-18:+), cds=[None], qualifiers=None)"
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
        "(CDSFrame.ZERO, CDSFrame.TWO, CDSFrame.TWO)], qualifiers=None)"
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
    e3_spliced_utr_repr = (
        "TranscriptInterval((2-6:+, 7-10:+, 12-15:+), cds=[CDS((7-10:+), (CDSFrame.ZERO)], qualifiers=None)"
    )
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
        "TranscriptInterval((2-6:+, 7-10:+, 12-15:+), cds=[CDS((7-10:+), (CDSFrame.ONE)], qualifiers=None)"
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
    se_unspliced_repr_minus = "TranscriptInterval((0-18:-), cds=[CDS((0-18:-), (CDSFrame.ZERO)], qualifiers=None)"
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
    se_unspliced_oof_repr_minus = "TranscriptInterval((0-18:-), cds=[CDS((0-18:-), (CDSFrame.ONE)], qualifiers=None)"
    # a single exon noncoding transcript (explicitly defined)
    se_noncoding_minus = TranscriptIntervalModel.Schema().load(
        dict(exon_starts=[0], exon_ends=[18], strand=Strand.MINUS.name)
    )
    se_noncoding_repr_minus = "TranscriptInterval((0-18:-), cds=[None], qualifiers=None)"
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
        "(CDSFrame.ONE, CDSFrame.ONE, CDSFrame.ZERO)], qualifiers=None)"
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
        "TranscriptInterval((2-6:-, 7-10:-, 12-15:-), cds=[CDS((7-10:-), (CDSFrame.ZERO)], qualifiers=None)"
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
        "TranscriptInterval((2-6:-, 7-10:-, 12-15:-), cds=[CDS((7-10:-), (CDSFrame.ONE)], qualifiers=None)"
    )

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
        assert str(schema.to_transcript_interval(parent=parent)) == expected

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
        tx = schema.to_transcript_interval(parent=parent)
        assert str(tx.get_protein_sequence()) == expected

    @pytest.mark.parametrize("schema,expected", [(se_unspliced_minus, "N*")])
    def test_truncations(self, schema, expected):
        tx = schema.to_transcript_interval(parent=parent)
        assert str(tx.get_protein_sequence(truncate_at_in_frame_stop=True)) == expected

    @pytest.mark.parametrize("schema,expected_exception", [(se_unspliced_minus, NullParentException)])
    def test_translation_exceptions(self, schema, expected_exception):
        tx = schema.to_transcript_interval()
        with pytest.raises(expected_exception):
            assert str(tx.get_protein_sequence(truncate_at_in_frame_stop=True)) == expected_exception

    @pytest.mark.parametrize(
        "location,cds,name,qualifiers,expected_exception",
        [
            (
                SingleInterval(0, 10, Strand.PLUS, parent),
                CDSInterval(SingleInterval(2, 8, Strand.PLUS), [CDSFrame.ZERO]),
                None,
                None,
                NullParentException,
            )
        ],
    )
    def test_transcript_constructor_exceptions(self, location, cds, name, qualifiers, expected_exception):
        with pytest.raises(expected_exception):
            TranscriptInterval(location, cds, name, qualifiers)

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
                LocationException,
            ),
        ],
    )
    def test_transcript_from_transcript_exceptions(self, schema, expected_exception):
        with pytest.raises(expected_exception):
            _ = schema.to_transcript_interval()
            _ = schema.to_transcript_interval(parent=parent)

    def test_no_such_ancestor(self):
        with pytest.raises(NullSequenceException):
            _ = self.se_unspliced.to_transcript_interval(parent=Parent())

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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
        assert tx.transcript_pos_to_cds(value) == expected

    def test_position_exceptions(self):
        tx = self.se_noncoding.to_transcript_interval()
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
            _ = tx.get_coding_interval()

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
        tx = schema.to_transcript_interval(parent=parent)
        val = tx.intersect(value)
        expected = TranscriptIntervalModel.Schema().load(expected).to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
        tx.update_parent(parent_genome2)
        assert str(tx.get_protein_sequence()) == str(expected)

    def test_object_conversion(self):
        for model in [
            self.e3_spliced,
            self.e3_spliced_utr_minus,
            self.e3_spliced_minus,
            self.e3_spliced_notrans,
        ]:
            obj = model.to_transcript_interval()
            new_model = TranscriptIntervalModel.from_transcript_interval(obj)
            new_model.transcript_guid = None
            assert model == new_model
            obj = model.to_transcript_interval(parent=parent)
            new_model = TranscriptIntervalModel.from_transcript_interval(obj)
            new_model.transcript_guid = None
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = schema.to_transcript_interval(parent=parent)
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
        tx = self.e3_spliced.to_transcript_interval()
        assert tx.is_coding
        assert tx.cds_start == 4
        assert tx.cds_end == 13
        tx2 = self.se_noncoding.to_transcript_interval()
        assert not tx2.is_coding

        with pytest.raises(NoncodingTranscriptError):
            _ = tx2.cds_start
        with pytest.raises(NoncodingTranscriptError):
            _ = tx2.cds_end

        # empty translation but still coding
        tx3 = self.e3_spliced_notrans_minus.to_transcript_interval()
        assert tx3.is_coding
        assert tx3.cds_start

    def test_overlapping(self):
        schema = TranscriptIntervalModel.Schema().load(
            dict(
                exon_starts=[0, 3],
                exon_ends=[5, 7],
                strand=Strand.PLUS.name,
                cds_starts=[1, 3],
                cds_ends=[5, 6],
                cds_frames=[CDSFrame.ZERO.name, CDSFrame.ZERO.name],
            )
        )
        tx = schema.to_transcript_interval()
        merged = tx.merge_overlapping()
        assert merged.cds is None
        assert merged.location == SingleInterval(0, 7, Strand.PLUS)

    def test_frameshifted(self):
        genome = "ATGGGGTGATGA"
        parent_genome = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT))
        tx = dict(
            exon_starts=[0],
            exon_ends=[12],
            strand=Strand.PLUS.name,
            cds_starts=[0],
            cds_ends=[12],
            cds_frames=[CDSFrame.ZERO.name],
        )
        schema = TranscriptIntervalModel.Schema().load(tx)
        obj = schema.to_transcript_interval(parent=parent_genome)
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
        obj2 = schema2.to_transcript_interval(parent=parent_genome)
        assert not obj2.has_in_frame_stop
