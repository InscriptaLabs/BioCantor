import pytest
from inscripta.biocantor.exc import ValidationException, EmptyLocationException, NullSequenceException
from inscripta.biocantor.io.models import FeatureIntervalModel
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence

# these features will be shared across all tests
genome = "GTATTCTTGGACCTAATT"
parent = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type="chromosome")
# offset the genome to show parent
genome2 = "AAGTATTCTTGGACCTAATT"
parent_genome2 = Parent(sequence=Sequence(genome2, Alphabet.NT_STRICT), sequence_type="chromosome")
parent_genome2_minus = Parent(
    sequence=Sequence(genome2, Alphabet.NT_STRICT), strand=Strand.MINUS, sequence_type="chromosome"
)

# slice the genome down to contain some of the transcripts
parent_genome2_1_15 = Parent(
    sequence=Sequence(
        genome2[1:15],
        Alphabet.NT_EXTENDED_GAPPED,
        type="sequence_chunk",
        parent=Parent(
            location=SingleInterval(1, 15, Strand.PLUS, parent=Parent(id="genome_1_15", sequence_type="chromosome"))
        ),
    )
)

# slice the genome down to contain none the transcripts
parent_genome2_2_8 = Parent(
    sequence=Sequence(
        genome2[2:8],
        Alphabet.NT_EXTENDED_GAPPED,
        type="sequence_chunk",
        parent=Parent(
            location=SingleInterval(2, 8, Strand.PLUS, parent=Parent(id="genome_2_8", sequence_type="chromosome"))
        ),
    )
)


# Integer feature definitions
# a single exon feature that is this entire genome
se_unspliced = FeatureIntervalModel.Schema().load(
    dict(interval_starts=[0], interval_ends=[18], strand=Strand.PLUS.name)
)
se_unspliced_repr = "FeatureInterval((0-18:+), name=None)"
# a three exon feature
e3_spliced = FeatureIntervalModel.Schema().load(
    dict(interval_starts=[2, 7, 12], interval_ends=[6, 10, 15], strand=Strand.PLUS.name)
)
e3_spliced_repr = "FeatureInterval((2-6:+, 7-10:+, 12-15:+), name=None)"
# Negative strand feature definitions
# a single exon feature that is this entire genome
se_unspliced_minus = FeatureIntervalModel.Schema().load(
    dict(interval_starts=[0], interval_ends=[18], strand=Strand.MINUS.name)
)
se_unspliced_repr_minus = "FeatureInterval((0-18:-), name=None)"
# a three exon negative strand feature
e3_spliced_minus = FeatureIntervalModel.Schema().load(
    dict(interval_starts=[2, 7, 12], interval_ends=[6, 10, 15], strand=Strand.MINUS.name)
)
e3_spliced_repr_minus = "FeatureInterval((2-6:-, 7-10:-, 12-15:-), name=None)"


class TestFeatureInterval:
    """Test constructing features of various types"""

    @pytest.mark.parametrize(
        "schema,expected",
        [
            (se_unspliced, se_unspliced_repr),
            (e3_spliced, e3_spliced_repr),
            (e3_spliced, e3_spliced_repr),
            (se_unspliced_minus, se_unspliced_repr_minus),
            (e3_spliced_minus, e3_spliced_repr_minus),
        ],
    )
    def test_feature_constructor(self, schema, expected):
        assert str(schema.to_feature_interval()) == expected
        assert str(schema.to_feature_interval(parent_or_seq_chunk_parent=parent)) == expected

    @pytest.mark.parametrize(
        "schema,expected_exception",
        [
            (
                FeatureIntervalModel.Schema().load(
                    dict(interval_starts=[0], interval_ends=[], strand=Strand.PLUS.name)
                ),
                ValidationException,
            ),
            (
                FeatureIntervalModel.Schema().load(
                    dict(interval_starts=[0], interval_ends=[2, 5], strand=Strand.PLUS.name)
                ),
                ValidationException,
            ),
            (
                FeatureIntervalModel.Schema().load(
                    dict(interval_starts=[], interval_ends=[2, 5], strand=Strand.PLUS.name)
                ),
                ValidationException,
            ),
        ],
    )
    def test_feature_excecptions(self, schema, expected_exception):
        with pytest.raises(expected_exception):
            _ = schema.to_feature_interval()

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced, 2, 0),
            (e3_spliced, 7, 4),
            (e3_spliced, 14, 9),
            (e3_spliced_minus, 2, 9),
            (e3_spliced_minus, 7, 5),
            (e3_spliced_minus, 14, 0),
        ],
    )
    def test_sequence_pos_to_feature(self, schema, value, expected):
        feat = schema.to_feature_interval()
        assert feat.sequence_pos_to_feature(value) == expected
        assert feat.chunk_relative_sequence_pos_to_feature(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced, (7, 13, Strand.PLUS), SingleInterval(4, 8, Strand.PLUS)),
            (e3_spliced_minus, (7, 13, Strand.PLUS), (SingleInterval(2, 6, Strand.MINUS))),
        ],
    )
    def test_sequence_interval_to_feature(self, schema, value, expected):
        feat = schema.to_feature_interval()
        assert feat.sequence_interval_to_feature(*value) == expected
        assert feat.chunk_relative_sequence_interval_to_feature(*value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced, 0, 2),
            (e3_spliced, 9, 14),
            (e3_spliced, 4, 7),
            (e3_spliced_minus, 0, 14),
            (e3_spliced_minus, 9, 2),
            (e3_spliced_minus, 5, 7),
        ],
    )
    def test_feature_pos_to_sequence(self, schema, value, expected):
        feat = schema.to_feature_interval()
        assert feat.feature_pos_to_sequence(value) == expected
        assert feat.feature_pos_to_chunk_relative_sequence(value) == expected
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert feat.feature_pos_to_sequence(value) == expected
        assert feat.feature_pos_to_chunk_relative_sequence(value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced, (0, 5, Strand.PLUS), CompoundInterval([2, 7], [6, 8], Strand.PLUS)),
            (e3_spliced_minus, (0, 5, Strand.PLUS), CompoundInterval([8, 12], [10, 15], Strand.MINUS)),
        ],
    )
    def test_feature_interval_to_sequence(self, schema, value, expected):
        feat = schema.to_feature_interval()
        assert feat.feature_interval_to_sequence(*value).reset_parent(None) == expected
        assert feat.feature_interval_to_chunk_relative_sequence(*value).reset_parent(None) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced, (0, 5, Strand.PLUS), CompoundInterval([2, 7], [6, 8], Strand.PLUS, parent=parent)),
            (e3_spliced_minus, (0, 5, Strand.PLUS), CompoundInterval([8, 12], [10, 15], Strand.MINUS, parent=parent)),
        ],
    )
    def test_feature_interval_to_sequence_parent(self, schema, value, expected):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert feat.feature_interval_to_sequence(*value) == expected
        assert feat.feature_interval_to_chunk_relative_sequence(*value) == expected

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (
                e3_spliced,
                SingleInterval(0, 10, Strand.PLUS, parent=parent),
                dict(interval_starts=(2, 7), interval_ends=(6, 10), strand=Strand.PLUS.name),
            ),
            (
                e3_spliced_minus,
                SingleInterval(0, 10, Strand.PLUS, parent=parent),
                dict(interval_starts=(2, 7), interval_ends=(6, 10), strand=Strand.MINUS.name),
            ),
        ],
    )
    def test_intersection(self, schema, value, expected):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        val = feat.intersect(value)
        expected = FeatureIntervalModel.Schema().load(expected).to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert str(expected) == str(val)

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (
                e3_spliced,
                SingleInterval(0, 10, Strand.PLUS),
                dict(interval_starts=(2, 7), interval_ends=(6, 10), strand=Strand.PLUS.name),
            ),
            (
                e3_spliced_minus,
                SingleInterval(0, 10, Strand.PLUS),
                dict(interval_starts=(2, 7), interval_ends=(6, 10), strand=Strand.MINUS.name),
            ),
        ],
    )
    def test_intersection_no_parent(self, schema, value, expected):
        feat = schema.to_feature_interval()
        val = feat.intersect(value)
        expected = FeatureIntervalModel.Schema().load(expected).to_feature_interval()
        assert str(expected) == str(val)

    def test_no_such_ancestor(self):
        with pytest.raises(NullSequenceException):
            _ = se_unspliced.to_feature_interval(parent_or_seq_chunk_parent=Parent(sequence_type="chromosome"))

    @pytest.mark.parametrize(
        "schema,parent,expected_spliced",
        [
            (e3_spliced, parent, "ATTCTGGCTA"),
            (e3_spliced, parent_genome2, "GTATCTTACC"),
            (
                # explicit strand on Parent has no effect
                e3_spliced,
                parent_genome2_minus,
                "GTATCTTACC",
            ),
            (e3_spliced_minus, parent, "TAGCCAGAAT"),
            (e3_spliced_minus, parent_genome2, "GGTAAGATAC"),
            (
                # explicit strand on Parent has no effect
                e3_spliced_minus,
                parent_genome2_minus,
                "GGTAAGATAC",
            ),
        ],
    )
    def test_spliced_sequence(self, schema, parent, expected_spliced):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert str(feat.get_spliced_sequence()) == expected_spliced

    @pytest.mark.parametrize(
        "schema,parent,expected_genomic,expected_stranded_genomic",
        [
            # positive strand, so genomic == genomic_stranded
            (e3_spliced, parent, "ATTCTTGGACCTA", "ATTCTTGGACCTA"),
            (e3_spliced, parent_genome2, "GTATTCTTGGACC", "GTATTCTTGGACC"),
            (
                # explicit strand on Parent has no effect
                e3_spliced,
                parent_genome2_minus,
                "GTATTCTTGGACC",
                "GTATTCTTGGACC",
            ),
            (e3_spliced_minus, parent, "ATTCTTGGACCTA", "TAGGTCCAAGAAT"),
            (e3_spliced_minus, parent_genome2, "GTATTCTTGGACC", "GGTCCAAGAATAC"),
            (
                # explicit strand on Parent flips it
                e3_spliced_minus,
                parent_genome2_minus,
                "GTATTCTTGGACC",
                "GGTCCAAGAATAC",
            ),
        ],
    )
    def test_genomic_sequence(self, schema, parent, expected_genomic, expected_stranded_genomic):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert str(feat.get_reference_sequence()) == expected_genomic
        assert str(feat.get_genomic_sequence()) == expected_stranded_genomic

    @pytest.mark.parametrize(
        "schema,expected",
        [(e3_spliced, "GTATCTTACC"), (e3_spliced_minus, "GGTAAGATAC")],
    )
    def test_reset_parent(self, schema, expected):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        feat.reset_parent(parent_genome2)
        assert str(feat.get_spliced_sequence()) == str(expected)

    def test_object_conversion(self):
        for model in [se_unspliced, e3_spliced_minus, e3_spliced]:
            obj = model.to_feature_interval()
            new_model = FeatureIntervalModel.from_feature_interval(obj)
            new_model.feature_interval_guid = None
            assert model == new_model
            obj = model.to_feature_interval(parent_or_seq_chunk_parent=parent)
            new_model = FeatureIntervalModel.from_feature_interval(obj)
            new_model.feature_interval_guid = None
            assert model == new_model

    def test_identifiers(self):
        feat = se_unspliced.to_feature_interval()
        feat.feature_name = "test"
        feat.feature_id = "testid"
        assert feat.identifiers == {"test", "testid"}
        assert feat.identifiers_dict == {"feature_name": "test", "feature_id": "testid"}

    def test_intersection_exception(self):
        schema = FeatureIntervalModel.Schema().load(
            dict(interval_starts=[10], interval_ends=[15], strand=Strand.MINUS.name)
        )
        feat = schema.to_feature_interval()
        s = SingleInterval(0, 5, Strand.PLUS)
        with pytest.raises(EmptyLocationException):
            _ = feat.intersect(s)


class TestFeatureIntervalSequenceSubset:
    """Test the ability to slice the genomic sequence of a feature interval and still get useful results."""

    @pytest.mark.parametrize(
        "schema,absolute_value,relative_value,expected,parent",
        [
            (e3_spliced, 2, 1, 0, parent_genome2_1_15),
            (e3_spliced, 7, 6, 4, parent_genome2_1_15),
            (e3_spliced, 14, 13, 9, parent_genome2_1_15),
            (e3_spliced_minus, 2, 1, 9, parent_genome2_1_15),
            (e3_spliced_minus, 7, 6, 5, parent_genome2_1_15),
            (e3_spliced_minus, 14, 13, 0, parent_genome2_1_15),
        ],
    )
    def test_sequence_pos_to_feature(self, schema, absolute_value, relative_value, expected, parent):
        feat = schema.to_feature_interval(parent)
        assert feat.sequence_pos_to_feature(absolute_value) == expected
        assert feat.chunk_relative_sequence_pos_to_feature(relative_value) == expected

    @pytest.mark.parametrize(
        "schema,absolute_value,relative_value,expected,parent",
        [
            (
                e3_spliced,
                (7, 13, Strand.PLUS),
                (6, 12, Strand.PLUS),
                SingleInterval(4, 8, Strand.PLUS),
                parent_genome2_1_15,
            ),
            (
                e3_spliced_minus,
                (7, 13, Strand.PLUS),
                (6, 12, Strand.PLUS),
                (SingleInterval(2, 6, Strand.MINUS)),
                parent_genome2_1_15,
            ),
        ],
    )
    def test_sequence_interval_to_feature(self, schema, absolute_value, relative_value, expected, parent):
        feat = schema.to_feature_interval(parent)
        assert feat.sequence_interval_to_feature(*absolute_value) == expected
        assert feat.chunk_relative_sequence_interval_to_feature(*relative_value) == expected

    @pytest.mark.parametrize(
        "schema,value,absolute_expected,relative_expected,parent",
        [
            (e3_spliced, 0, 2, 1, parent_genome2_1_15),
            (e3_spliced, 9, 14, 13, parent_genome2_1_15),
            (e3_spliced, 4, 7, 6, parent_genome2_1_15),
            (e3_spliced_minus, 0, 14, 13, parent_genome2_1_15),
            (e3_spliced_minus, 9, 2, 1, parent_genome2_1_15),
            (e3_spliced_minus, 5, 7, 6, parent_genome2_1_15),
        ],
    )
    def test_feature_pos_to_sequence(self, schema, value, absolute_expected, relative_expected, parent):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert feat.feature_pos_to_sequence(value) == absolute_expected
        assert feat.feature_pos_to_chunk_relative_sequence(value) == relative_expected

    @pytest.mark.parametrize(
        "schema,value,absolute_expected,relative_expected,parent",
        [
            (
                e3_spliced,
                (0, 5, Strand.PLUS),
                CompoundInterval([2, 7], [6, 8], Strand.PLUS),
                CompoundInterval([1, 6], [5, 7], Strand.PLUS),
                parent_genome2_1_15,
            ),
            (
                e3_spliced_minus,
                (0, 5, Strand.PLUS),
                CompoundInterval([8, 12], [10, 15], Strand.MINUS),
                CompoundInterval([7, 11], [9, 14], Strand.MINUS),
                parent_genome2_1_15,
            ),
        ],
    )
    def test_feature_interval_to_sequence(self, schema, value, absolute_expected, relative_expected, parent):
        feat = schema.to_feature_interval(parent)
        assert feat.feature_interval_to_sequence(*value).reset_parent(None) == absolute_expected
        assert feat.feature_interval_to_chunk_relative_sequence(*value).reset_parent(None) == relative_expected

    @pytest.mark.parametrize(
        "schema,parent,expected_spliced",
        [
            (e3_spliced, parent_genome2_1_15, "GTATCTTACC"),
            (
                e3_spliced,
                parent_genome2_1_15,
                "GTATCTTACC",
            ),
            (e3_spliced_minus, parent_genome2_1_15, "GGTAAGATAC"),
            (
                e3_spliced_minus,
                parent_genome2_1_15,
                "GGTAAGATAC",
            ),
        ],
    )
    def test_spliced_sequence(self, schema, parent, expected_spliced):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert str(feat.get_spliced_sequence()) == expected_spliced

    @pytest.mark.parametrize(
        "schema,parent,expected_genomic,expected_stranded_genomic",
        [
            (e3_spliced, parent_genome2_1_15, "GTATTCTTGGACC", "GTATTCTTGGACC"),
            (
                e3_spliced,
                parent_genome2_1_15,
                "GTATTCTTGGACC",
                "GTATTCTTGGACC",
            ),
            (e3_spliced_minus, parent_genome2_1_15, "GTATTCTTGGACC", "GGTCCAAGAATAC"),
        ],
    )
    def test_genomic_sequence(self, schema, parent, expected_genomic, expected_stranded_genomic):
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert str(feat.get_reference_sequence()) == expected_genomic
        assert str(feat.get_genomic_sequence()) == expected_stranded_genomic

    def test_start_end(self):
        feat = e3_spliced.to_feature_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert feat.chunk_relative_start + 1 == feat.start
        assert feat.chunk_relative_end + 1 == feat.end

    @pytest.mark.parametrize(
        "schema,parent,expected_exception",
        [
            [e3_spliced, parent_genome2_2_8, ValidationException],
            [e3_spliced_minus, parent_genome2_2_8, ValidationException],
            [se_unspliced, parent_genome2_1_15, ValidationException],
        ],
    )
    def test_position_exceptions(self, schema, parent, expected_exception):
        with pytest.raises(expected_exception):
            _ = schema.to_feature_interval(parent)

    @pytest.mark.parametrize(
        "schema,parent,expected_exception",
        [
            [e3_spliced, Parent(), ValidationException],
            [e3_spliced_minus, Parent(sequence_type="chromosome"), NullSequenceException],
            [se_unspliced, Parent(sequence_type="sequence_chunk"), NullSequenceException],
        ],
    )
    def test_constructor_exceptions(self, schema, parent, expected_exception):
        with pytest.raises(expected_exception):
            _ = schema.to_feature_interval(parent)

    def test_dict(self):
        feat = e3_spliced.to_feature_interval(parent_genome2_1_15)
        feat2 = e3_spliced.to_feature_interval()
        d = feat.to_dict()
        d2 = feat2.to_dict()
        del d["feature_interval_guid"]
        del d2["feature_interval_guid"]
        assert d == d2

        rel_d = feat.to_dict(chromosome_relative_coordinates=False)
        del rel_d["feature_interval_guid"]
        assert rel_d != d
        assert [x + 1 for x in rel_d["interval_starts"]] == list(d["interval_starts"])
