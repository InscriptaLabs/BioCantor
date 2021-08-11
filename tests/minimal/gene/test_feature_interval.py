from uuid import UUID

import pytest

from inscripta.biocantor.exc import (
    ValidationException,
    EmptyLocationException,
    NoSuchAncestorException,
    NullSequenceException,
    MismatchedParentException,
    NoncodingTranscriptError,
)
from inscripta.biocantor.gene.feature import FeatureInterval
from inscripta.biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from inscripta.biocantor.io.models import FeatureIntervalModel
from inscripta.biocantor.location.location_impl import SingleInterval, CompoundInterval, EmptyLocation
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent, SequenceType
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.util.object_validation import ObjectValidation

# these features will be shared across all tests
genome = "GTATTCTTGGACCTAATT"
parent = Parent(id="genome", sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type=SequenceType.CHROMOSOME)
# offset the genome to show parent
genome2 = "AAGTATTCTTGGACCTAATT"
parent_genome2 = Parent(
    id="genome2", sequence=Sequence(genome2, Alphabet.NT_STRICT), sequence_type=SequenceType.CHROMOSOME
)
parent_genome2_minus = Parent(
    id="genome2_minus",
    sequence=Sequence(genome2, Alphabet.NT_STRICT),
    strand=Strand.MINUS,
    sequence_type=SequenceType.CHROMOSOME,
)

# slice the genome down to contain some of the transcripts
parent_genome2_1_15 = Parent(
    id="genome2_1_15",
    sequence=Sequence(
        genome2[1:15],
        Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            location=SingleInterval(
                1, 15, Strand.PLUS, parent=Parent(id="genome2", sequence_type=SequenceType.CHROMOSOME)
            )
        ),
    ),
)

parent_genome2_2_8 = Parent(
    id="genome2_2_8",
    sequence=Sequence(
        genome2[2:8],
        Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            location=SingleInterval(
                2, 8, Strand.PLUS, parent=Parent(id="genome2", sequence_type=SequenceType.CHROMOSOME)
            )
        ),
    ),
)


parent_no_seq = Parent(sequence_type=SequenceType.CHROMOSOME)
parent_no_seq_with_id = Parent(sequence_type=SequenceType.CHROMOSOME, id="genome2")
parent_chunk_only = Parent(sequence_type=SequenceType.SEQUENCE_CHUNK)
parent_nonstandard_type = Parent(sequence_type="SomeOtherType", id="genome2")
parent_nonstandard_type_with_sequence = Parent(
    sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type="SomeOtherType"
)
parent_named = Parent(sequence_type=SequenceType.CHROMOSOME, id="chrom")


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
        assert feat.chunk_relative_pos_to_feature(value) == expected

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
        assert feat.chunk_relative_interval_to_feature(*value) == expected

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
        assert feat.feature_pos_to_chunk_relative(value) == expected
        feat = schema.to_feature_interval(parent_or_seq_chunk_parent=parent)
        assert feat.feature_pos_to_sequence(value) == expected
        assert feat.feature_pos_to_chunk_relative(value) == expected

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
        assert feat.feature_interval_to_chunk_relative(*value).reset_parent(None) == expected

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
        assert feat.feature_interval_to_chunk_relative(*value) == expected

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
        feat._reset_parent(parent_genome2)
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

    def test_gff_export_exceptions(self):
        feat = se_unspliced.to_feature_interval(parent_or_seq_chunk_parent=parent)
        with pytest.raises(GFF3MissingSequenceNameError):
            _ = "\n".join(str(x) for x in feat.to_gff())
        feat.sequence_name = "myseq"
        with pytest.raises(NoSuchAncestorException):
            _ = "\n".join(str(x) for x in feat.to_gff(chromosome_relative_coordinates=False))

    def test_gff_export(self):
        feat = se_unspliced.to_feature_interval(parent_or_seq_chunk_parent=parent)
        feat.sequence_name = "myseq"
        assert (
            "\n".join(str(x) for x in feat.to_gff())
            == "myseq\tBioCantor\tfeature_interval\t1\t18\t.\t+\t.\tID=6940c467-070a-3524-2dcb-a478a6fa0388\n"
            "myseq\tBioCantor\tsubregion\t1\t18\t.\t+\t.\tID=feature-6940c467-070a-3524-2dcb-a478a6fa0388-1;"
            "Parent=6940c467-070a-3524-2dcb-a478a6fa0388"
        )

    def test_gff_export_subset(self):
        feat = e3_spliced.to_feature_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        feat.sequence_name = "myseq"
        assert (
            "\n".join(str(x) for x in feat.to_gff())
            == "myseq\tBioCantor\tfeature_interval\t3\t15\t.\t+\t.\tID=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80\n"
            "myseq\tBioCantor\tsubregion\t3\t6\t.\t+\t.\tID=feature-c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80-1;"
            "Parent=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80\n"
            "myseq\tBioCantor\tsubregion\t8\t10\t.\t+\t.\tID=feature-c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80-2;"
            "Parent=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80\n"
            "myseq\tBioCantor\tsubregion\t13\t15\t.\t+\t.\tID=feature-c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80-3;"
            "Parent=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80"
        )
        assert (
            "\n".join(str(x) for x in feat.to_gff(chromosome_relative_coordinates=False))
            == "myseq\tBioCantor\tfeature_interval\t2\t14\t.\t+\t.\tID=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80\n"
            "myseq\tBioCantor\tsubregion\t2\t5\t.\t+\t.\tID=feature-c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80-1;"
            "Parent=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80\n"
            "myseq\tBioCantor\tsubregion\t7\t9\t.\t+\t.\tID=feature-c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80-2;"
            "Parent=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80\n"
            "myseq\tBioCantor\tsubregion\t12\t14\t.\t+\t.\tID=feature-c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80-3;"
            "Parent=c45e8d7b-cbd6-43b2-bb08-429d9cb7fe80"
        )

    @pytest.mark.parametrize(
        "feature,parent,expected_span",
        [
            (e3_spliced, None, SingleInterval(2, 15, Strand.PLUS)),
            (e3_spliced_minus, None, SingleInterval(2, 15, Strand.MINUS)),
            (se_unspliced, None, SingleInterval(0, 18, Strand.PLUS)),
            (se_unspliced, None, SingleInterval(0, 18, Strand.PLUS)),
            (e3_spliced, parent_genome2, SingleInterval(2, 15, Strand.PLUS, parent=parent_genome2)),
            (e3_spliced_minus, parent_genome2, SingleInterval(2, 15, Strand.MINUS, parent=parent_genome2)),
            (se_unspliced, parent_genome2, SingleInterval(0, 18, Strand.PLUS, parent=parent_genome2)),
            (se_unspliced, parent_genome2, SingleInterval(0, 18, Strand.PLUS, parent=parent_genome2)),
            (
                e3_spliced,
                parent_genome2_1_15,
                SingleInterval(
                    2, 15, Strand.PLUS, parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME)
                ),
            ),
            (
                e3_spliced_minus,
                parent_genome2_1_15,
                SingleInterval(
                    2, 15, Strand.MINUS, parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME)
                ),
            ),
            (
                e3_spliced,
                parent_genome2_2_8,
                SingleInterval(
                    2, 15, Strand.PLUS, parent=parent_genome2_2_8.first_ancestor_of_type(SequenceType.CHROMOSOME)
                ),
            ),
            (
                e3_spliced_minus,
                parent_genome2_2_8,
                SingleInterval(
                    2, 15, Strand.MINUS, parent=parent_genome2_2_8.first_ancestor_of_type(SequenceType.CHROMOSOME)
                ),
            ),
        ],
    )
    def test_chromosome_span(self, feature, parent, expected_span):
        feat = feature.to_feature_interval(parent)
        assert feat.chromosome_span == expected_span

    @pytest.mark.parametrize(
        "feature,parent,expected_gaps",
        [
            (e3_spliced, None, CompoundInterval([6, 10], [7, 12], Strand.PLUS)),
            (e3_spliced_minus, None, CompoundInterval([6, 10], [7, 12], Strand.MINUS)),
            (se_unspliced, None, EmptyLocation()),
            (se_unspliced, None, EmptyLocation()),
            (e3_spliced, parent_genome2, CompoundInterval([6, 10], [7, 12], Strand.PLUS, parent=parent_genome2)),
            (e3_spliced_minus, parent_genome2, CompoundInterval([6, 10], [7, 12], Strand.MINUS, parent=parent_genome2)),
            (se_unspliced, parent_genome2, EmptyLocation()),
            (se_unspliced, parent_genome2, EmptyLocation()),
            (
                e3_spliced,
                parent_genome2_1_15,
                CompoundInterval(
                    [6, 10],
                    [7, 12],
                    Strand.PLUS,
                    parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
            (
                e3_spliced_minus,
                parent_genome2_1_15,
                CompoundInterval(
                    [6, 10],
                    [7, 12],
                    Strand.MINUS,
                    parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
            (se_unspliced, parent_genome2_1_15, EmptyLocation()),
            (se_unspliced, parent_genome2_1_15, EmptyLocation()),
            (
                e3_spliced,
                parent_genome2_2_8,
                CompoundInterval(
                    [6, 10],
                    [7, 12],
                    Strand.PLUS,
                    parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
            (
                e3_spliced_minus,
                parent_genome2_2_8,
                CompoundInterval(
                    [6, 10],
                    [7, 12],
                    Strand.MINUS,
                    parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
        ],
    )
    def test_chromosome_gaps_location(self, feature, parent, expected_gaps):
        feat = feature.to_feature_interval(parent)
        assert feat.chromosome_gaps_location == expected_gaps


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
        assert feat.chunk_relative_pos_to_feature(relative_value) == expected

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
        assert feat.chunk_relative_interval_to_feature(*relative_value) == expected

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
        assert feat.feature_pos_to_chunk_relative(value) == relative_expected

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
        assert feat.feature_interval_to_chunk_relative(*value).reset_parent(None) == relative_expected

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

    def test_sequence_exceptions(self):
        """All sequence accessors should raise good errors when attempted without sequence info"""
        feat = e3_spliced.to_feature_interval(parent_or_seq_chunk_parent=parent_no_seq)
        with pytest.raises(NullSequenceException):
            _ = feat.get_reference_sequence()
        with pytest.raises(NullSequenceException):
            _ = feat.get_spliced_sequence()
        with pytest.raises(NullSequenceException):
            _ = feat.get_genomic_sequence()

    def test_start_end(self):
        feat = e3_spliced.to_feature_interval(parent_or_seq_chunk_parent=parent_genome2_1_15)
        assert feat.chunk_relative_start + 1 == feat.start
        assert feat.chunk_relative_end + 1 == feat.end

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

    def test_nonstandard_parents(self):
        feat0 = e3_spliced.to_feature_interval(parent)
        seq0 = feat0.get_spliced_sequence()
        feat1 = e3_spliced.to_feature_interval(parent_nonstandard_type)
        with pytest.raises(NullSequenceException):
            _ = feat1.get_spliced_sequence()
        feat2 = e3_spliced.to_feature_interval(parent_no_seq)
        with pytest.raises(NullSequenceException):
            _ = feat2.get_spliced_sequence()
        feat3 = e3_spliced.to_feature_interval(parent_nonstandard_type_with_sequence)
        seq = feat3.get_spliced_sequence()
        assert seq == seq0

        assert feat0.chromosome_location == feat0.chunk_relative_location
        assert feat1.chromosome_location != feat1.chunk_relative_location
        assert feat1._chunk_relative_bounded_chromosome_location == feat1.chunk_relative_location
        assert feat2.chromosome_location == feat2.chunk_relative_location
        assert feat3.chromosome_location != feat3.chunk_relative_location
        assert feat3._chunk_relative_bounded_chromosome_location == feat3.chunk_relative_location
        # OTOH, this is not the same
        feat4 = e3_spliced.to_feature_interval(parent_genome2_1_15)
        assert feat4.chromosome_location != feat4.chunk_relative_location

    def test_liftover_to_parent_or_seq_chunk_parent(self):
        feat0 = e3_spliced.to_feature_interval(parent_genome2)
        feat1 = feat0.liftover_to_parent_or_seq_chunk_parent(parent_genome2_1_15)
        assert str(feat0.get_spliced_sequence()) == str(feat1.get_spliced_sequence())
        assert feat0.chromosome_location.reset_parent(None) == feat1.chromosome_location.reset_parent(None)
        # bringing in a null parent means no sequence anymore
        feat2 = feat0.liftover_to_parent_or_seq_chunk_parent(parent_no_seq_with_id)
        with pytest.raises(NullSequenceException):
            _ = feat2.get_spliced_sequence()

        # we can also start in chunk coordinates, then lift to different chunk coordinates
        feat_chunk = e3_spliced.to_feature_interval(parent_genome2_1_15)
        feat_subchunk = feat_chunk.liftover_to_parent_or_seq_chunk_parent(parent_genome2_2_8)
        # this is now a subset
        assert str(feat_subchunk.get_spliced_sequence()) in str(feat_chunk.get_spliced_sequence())

    def test_liftover_to_parent_or_seq_chunk_parent_exception(self):
        feat0 = e3_spliced.to_feature_interval(parent_genome2_1_15)
        with pytest.raises(MismatchedParentException):
            _ = feat0.liftover_to_parent_or_seq_chunk_parent(parent_named)

        feat1 = e3_spliced.to_feature_interval(parent_genome2)
        with pytest.raises(MismatchedParentException):
            _ = feat1.liftover_to_parent_or_seq_chunk_parent(parent_named)

        feat2 = e3_spliced.to_feature_interval(parent_genome2_1_15)
        with pytest.raises(MismatchedParentException):
            _ = feat2.liftover_to_parent_or_seq_chunk_parent(parent_no_seq)

        with pytest.raises(MismatchedParentException):
            _ = feat2.liftover_to_parent_or_seq_chunk_parent(parent)

    @pytest.mark.parametrize(
        "feature,parent,expected",
        [
            (
                e3_spliced,
                parent_genome2_2_8,
                CompoundInterval(
                    [2, 7, 12],
                    [6, 10, 15],
                    Strand.PLUS,
                    parent=parent_genome2_2_8.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
            (
                e3_spliced,
                parent_genome2_1_15,
                CompoundInterval(
                    [2, 7, 12],
                    [6, 10, 15],
                    Strand.PLUS,
                    parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
        ],
    )
    def test_chromosome_location(self, feature, parent, expected):
        feat = feature.to_feature_interval(parent)
        assert feat.chromosome_location == expected

    @pytest.mark.parametrize(
        "feature,parent,expected",
        [
            # parent_genome2_2_8 slices off the first exon
            (
                e3_spliced,
                parent_genome2_2_8,
                CompoundInterval(
                    [2, 7],
                    [6, 8],
                    Strand.PLUS,
                    parent=parent_genome2_2_8.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
            # parent_genome2_1_15 does not slice off anything
            (
                e3_spliced,
                parent_genome2_1_15,
                CompoundInterval(
                    [2, 7, 12],
                    [6, 10, 15],
                    Strand.PLUS,
                    parent=parent_genome2_1_15.first_ancestor_of_type(SequenceType.CHROMOSOME),
                ),
            ),
        ],
    )
    def test__chunk_relative_bounded_chromosome_location(self, feature, parent, expected):
        feat = feature.to_feature_interval(parent)
        assert feat._chunk_relative_bounded_chromosome_location == expected

    @pytest.mark.parametrize(
        "feature,parent,expected_gaps",
        [
            (
                e3_spliced,
                parent_genome2_1_15,
                CompoundInterval([5, 9], [6, 11], Strand.PLUS, parent=parent_genome2_1_15),
            ),
            (
                e3_spliced_minus,
                parent_genome2_1_15,
                CompoundInterval([5, 9], [6, 11], Strand.MINUS, parent=parent_genome2_1_15),
            ),
            (e3_spliced, parent_genome2_2_8, SingleInterval(4, 5, Strand.PLUS, parent=parent_genome2_2_8)),
            (e3_spliced_minus, parent_genome2_2_8, SingleInterval(4, 5, Strand.MINUS, parent=parent_genome2_2_8)),
        ],
    )
    def test_chunk_relative_gaps_location(self, feature, parent, expected_gaps):
        feat = feature.to_feature_interval(parent)
        ObjectValidation.require_locations_have_same_nonempty_parent(feat.chunk_relative_gaps_location, expected_gaps)

    @pytest.mark.parametrize(
        "feature,parent,expected_span",
        [
            (e3_spliced, parent_genome2_1_15, SingleInterval(2, 14, Strand.PLUS, parent=parent_genome2_1_15)),
            (e3_spliced_minus, parent_genome2_1_15, SingleInterval(2, 14, Strand.MINUS, parent=parent_genome2_1_15)),
            (e3_spliced, parent_genome2_2_8, SingleInterval(2, 6, Strand.PLUS, parent=parent_genome2_2_8)),
            (e3_spliced_minus, parent_genome2_2_8, SingleInterval(2, 6, Strand.MINUS, parent=parent_genome2_2_8)),
        ],
    )
    def test_chunk_relative_span(self, feature, parent, expected_span):
        feat = feature.to_feature_interval(parent)
        ObjectValidation.require_locations_have_same_nonempty_parent(expected_span, feat.chunk_relative_span)


class TestFeatureWithoutModel:
    @pytest.mark.parametrize(
        "feat",
        [
            dict(
                interval_starts=[2, 8, 12],
                interval_ends=[5, 13, 18],
                strand=Strand.PLUS,
            ),
            dict(
                interval_starts=[2, 8, 12],
                interval_ends=[5, 13, 18],
                strand=Strand.PLUS,
                parent_or_seq_chunk_parent=Sequence(
                    "AAAGGAAAGTCCCTGAAAAAA", Alphabet.NT_EXTENDED_GAPPED, type=SequenceType.CHROMOSOME
                ),
            ),
        ],
    )
    def test_constructor(self, feat):
        feat = FeatureInterval(**feat)
        feat2 = FeatureInterval.from_location(feat._location)
        assert feat == feat2

    def test_equality_different_parents(self):
        feat1 = e3_spliced.to_feature_interval(parent_genome2)
        feat2 = e3_spliced.to_feature_interval(parent_genome2_1_15)
        assert feat1 != feat2


class TestOverlappingInterval:
    def test_overlapping_interval(self):
        exon_starts, exon_ends = ([2, 5], [6, 10])
        strand = Strand.PLUS
        i = FeatureInterval.initialize_location(exon_starts, exon_ends, strand, parent_genome2_1_15)
        assert i.num_blocks == 2


@pytest.mark.parametrize(
    "data,expected_exception",
    [
        # must be a dictionary
        (
            {
                "interval_starts": [0],
                "interval_ends": [18],
                "strand": "PLUS",
                "qualifiers": "abc",
                "feature_id": None,
                "feature_name": None,
                "feature_types": None,
                "sequence_name": None,
                "sequence_guid": None,
                "feature_interval_guid": UUID("6940c467-070a-3524-2dcb-a478a6fa0388"),
                "feature_guid": None,
                "is_primary_feature": None,
            },
            ValidationException,
        ),
        # must be a dictionary of lists
        (
            {
                "interval_starts": [0],
                "interval_ends": [18],
                "strand": "PLUS",
                "qualifiers": {"abc": [123], "xyz": 123},
                "feature_id": None,
                "feature_name": None,
                "feature_types": None,
                "sequence_name": None,
                "sequence_guid": None,
                "feature_interval_guid": UUID("6940c467-070a-3524-2dcb-a478a6fa0388"),
                "feature_guid": None,
                "is_primary_feature": None,
            },
            ValidationException,
        ),
    ],
)
def test_qualifiers_exceptions(data, expected_exception):
    """Qualifiers must be dictionaries of lists"""
    with pytest.raises(expected_exception):
        _ = FeatureInterval.from_dict(data)


@pytest.mark.parametrize(
    "property,expected_exception",
    [
        ("cds_start", NoncodingTranscriptError),
        ("cds_end", NoncodingTranscriptError),
        ("chunk_relative_cds_start", NoncodingTranscriptError),
        ("chunk_relative_cds_end", NoncodingTranscriptError),
        ("cds_location", NoncodingTranscriptError),
        ("cds_chunk_relative_location", NoncodingTranscriptError),
        ("is_coding", NoncodingTranscriptError),
        ("has_in_frame_stop", NoncodingTranscriptError),
        ("cds_size", NoncodingTranscriptError),
        ("chunk_relative_cds_size", NoncodingTranscriptError),
    ],
)
def test_property_exceptions(property, expected_exception):
    obj = e3_spliced.to_feature_interval()
    with pytest.raises(expected_exception):
        getattr(obj, property)
