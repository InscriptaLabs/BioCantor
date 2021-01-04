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
parent = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT))
# offset the genome to show parent
genome2 = "AAGTATTCTTGGACCTAATT"
parent_genome2 = Parent(sequence=Sequence(genome2, Alphabet.NT_STRICT))
parent_genome2_minus = Parent(
    sequence=Sequence(genome2, Alphabet.NT_STRICT),
    strand=Strand.MINUS,
)


class TestFeatureInterval:
    """Test constructing features of various types"""

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
        assert str(schema.to_feature_interval(parent=parent)) == expected

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
        feat = schema.to_feature_interval(parent=parent)
        assert feat.feature_pos_to_sequence(value) == expected

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

    @pytest.mark.parametrize(
        "schema,value,expected",
        [
            (e3_spliced, (0, 5, Strand.PLUS), CompoundInterval([2, 7], [6, 8], Strand.PLUS, parent=parent)),
            (e3_spliced_minus, (0, 5, Strand.PLUS), CompoundInterval([8, 12], [10, 15], Strand.MINUS, parent=parent)),
        ],
    )
    def test_feature_interval_to_sequence_parent(self, schema, value, expected):
        feat = schema.to_feature_interval(parent=parent)
        assert feat.feature_interval_to_sequence(*value) == expected

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
        feat = schema.to_feature_interval(parent=parent)
        val = feat.intersect(value)
        expected = FeatureIntervalModel.Schema().load(expected).to_feature_interval(parent=parent)
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
            _ = self.se_unspliced.to_feature_interval(parent=Parent())

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
        feat = schema.to_feature_interval(parent=parent)
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
        feat = schema.to_feature_interval(parent=parent)
        assert str(feat.get_reference_sequence()) == expected_genomic
        assert str(feat.get_genomic_sequence()) == expected_stranded_genomic

    @pytest.mark.parametrize(
        "schema,expected",
        [(e3_spliced, "GTATCTTACC"), (e3_spliced_minus, "GGTAAGATAC")],
    )
    def test_reset_parent(self, schema, expected):
        feat = schema.to_feature_interval(parent=parent)
        feat.update_parent(parent_genome2)
        assert str(feat.get_spliced_sequence()) == str(expected)

    def test_object_conversion(self):
        for model in [self.se_unspliced, self.e3_spliced_minus, self.e3_spliced]:
            obj = model.to_feature_interval()
            new_model = FeatureIntervalModel.from_feature_interval(obj)
            new_model.feature_interval_guid = None
            assert model == new_model
            obj = model.to_feature_interval(parent=parent)
            new_model = FeatureIntervalModel.from_feature_interval(obj)
            new_model.feature_interval_guid = None
            assert model == new_model

    def test_identifiers(self):
        feat = self.se_unspliced.to_feature_interval()
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
