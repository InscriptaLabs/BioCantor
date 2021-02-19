import pytest
from inscripta.biocantor.exc import InvalidAnnotationError, NoncodingTranscriptError, InvalidQueryError
from inscripta.biocantor.exc import ValidationException
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame
from inscripta.biocantor.io.models import (
    GeneIntervalModel,
    AnnotationCollectionModel,
    FeatureIntervalCollectionModel,
    TranscriptIntervalModel,
)
from inscripta.biocantor.location.location_impl import SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence

genome = "AAGTATTCTTGGACCTAATTAAAAAAAAAAAAAAAAAAA"
parent_genome = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type="chromosome")


class TestGene:
    """Test basic gene construction from Transcripts"""

    tx1 = dict(
        exon_starts=[2],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[5],
        cds_ends=[9],
        cds_frames=[CDSFrame.ZERO.name],
    )
    tx2 = dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
    )
    tx1_primary = dict(
        exon_starts=[2],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[5],
        cds_ends=[9],
        cds_frames=[CDSFrame.ZERO.name],
        is_primary_tx=True,
    )
    tx2_primary = dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        is_primary_tx=True,
    )
    tx_noncoding = dict(exon_starts=[2, 7, 12], exon_ends=[6, 10, 30], strand=Strand.PLUS.name)
    tx_noncoding_short = dict(exon_starts=[2], exon_ends=[6], strand=Strand.PLUS.name)

    gene = GeneIntervalModel.Schema().load(
        dict(transcripts=[tx1, tx2], qualifiers={}, gene_type=Biotype.protein_coding.name, gene_id="gene1")
    )

    gene_noncoding = GeneIntervalModel.Schema().load(dict(transcripts=[tx_noncoding, tx_noncoding_short]))

    def test_primary_inference(self):
        obj = self.gene.to_gene_interval()
        assert obj.primary_transcript == TranscriptIntervalModel.Schema().load(self.tx2).to_transcript_interval()

    def test_merged_interval(self):
        obj = self.gene.to_gene_interval()
        assert str(obj.get_merged_transcript()) == "FeatureInterval((2-18:+), name=None)"
        assert str(obj.get_merged_cds()) == "FeatureInterval((4-10:+, 12-13:+), name=None)"

    def test_failed_merge_interval(self):
        obj = self.gene_noncoding.to_gene_interval()
        with pytest.raises(NoncodingTranscriptError):
            _ = obj.get_merged_cds()

    def test_empty(self):
        gene = GeneIntervalModel.Schema().load(
            dict(transcripts=[], qualifiers={}, gene_type=Biotype.protein_coding.name, gene_id="gene1")
        )
        with pytest.raises(InvalidAnnotationError):
            _ = gene.to_gene_interval()

    def test_iter(self):
        obj = self.gene.to_gene_interval()
        assert list(obj) == obj.transcripts

    def test_primary_exception(self):
        with pytest.raises(ValidationException):
            gene = GeneIntervalModel.Schema().load(
                dict(
                    transcripts=[self.tx1_primary, self.tx2_primary],
                    qualifiers={},
                    gene_type=Biotype.protein_coding.name,
                    gene_id="gene1",
                )
            )
            _ = gene.to_gene_interval()

    @pytest.mark.parametrize(
        "gene,expected",
        [  # explicit primary overrules hierarchy
            (dict(transcripts=[tx1_primary, tx2]), tx1_primary),
            # without explicit primary, pick longest CDS
            (dict(transcripts=[tx1, tx2, tx_noncoding]), tx2),
            # longest non-coding
            (dict(transcripts=[tx_noncoding, tx_noncoding_short]), tx_noncoding),
        ],
    )
    def test_get_primary(self, gene, expected):
        obj = GeneIntervalModel.Schema().load(gene).to_gene_interval()
        expected = TranscriptIntervalModel.Schema().load(expected).to_transcript_interval()
        assert obj.get_primary_transcript() == expected

    def test_get_sequence(self):
        obj = self.gene.to_gene_interval(parent_or_seq_chunk_parent=parent_genome)
        assert str(obj.get_primary_transcript_sequence()) == "GTATCTTACC"
        assert str(obj.get_primary_cds_sequence()) == "ATCTTA"
        assert str(obj.get_primary_protein()) == "IL"


class TestFeatureIntervalCollection:
    feat1 = dict(
        interval_starts=[2], interval_ends=[5], strand=Strand.PLUS.name, feature_types=["a", "b"], feature_name="feat1"
    )
    feat2 = dict(
        interval_starts=[2, 7, 12],
        interval_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        feature_types=["b"],
        feature_name="feat2",
    )
    feat3 = dict(
        interval_starts=[25], interval_ends=[30], strand=Strand.MINUS.name, feature_types=["a"], feature_name="feat3"
    )
    collection1 = FeatureIntervalCollectionModel.Schema().load(
        dict(feature_intervals=[feat1, feat2], feature_collection_name="featgrp1")
    )
    collection2 = FeatureIntervalCollectionModel.Schema().load(
        dict(feature_intervals=[feat3], feature_collection_name="featgrp2")
    )

    def test_feature_collection(self):
        obj = self.collection1.to_feature_collection()
        model = FeatureIntervalCollectionModel.from_feature_collection(obj)
        # remove guids to make comparison work
        for item in model.feature_intervals:
            item.feature_interval_guid = None
        model.feature_collection_guid = None
        assert model == self.collection1

    def test_empty(self):
        feat = FeatureIntervalCollectionModel.Schema().load(dict(feature_intervals=[]))
        with pytest.raises(InvalidAnnotationError):
            _ = feat.to_feature_collection()

    def test_iter(self):
        obj = self.collection1.to_feature_collection()
        assert list(obj) == obj.feature_intervals


class TestAnnotationCollection:
    annot = AnnotationCollectionModel.Schema().load(
        dict(
            feature_collections=[
                dict(
                    feature_intervals=[TestFeatureIntervalCollection.feat1, TestFeatureIntervalCollection.feat2],
                    feature_collection_id="featgrp1",
                ),
                dict(feature_intervals=[TestFeatureIntervalCollection.feat3], feature_collection_id="featgrp2"),
            ],
            genes=[dict(transcripts=[TestGene.tx1, TestGene.tx2], gene_id="gene1")],
            start=0,
            end=30,
        )
    )

    annot_no_range = AnnotationCollectionModel.Schema().load(
        dict(
            feature_collections=[
                dict(
                    feature_intervals=[TestFeatureIntervalCollection.feat1, TestFeatureIntervalCollection.feat2],
                    feature_collection_id="featgrp1",
                ),
                dict(feature_intervals=[TestFeatureIntervalCollection.feat3], feature_collection_id="featgrp2"),
            ],
            genes=[dict(transcripts=[TestGene.tx1, TestGene.tx2], gene_id="gene1")],
        )
    )

    annot_no_features = AnnotationCollectionModel.Schema().load(
        dict(genes=[dict(transcripts=[TestGene.tx1, TestGene.tx2], gene_id="gene1")])
    )

    annot_no_genes = AnnotationCollectionModel.Schema().load(
        dict(
            feature_collections=[
                dict(
                    feature_intervals=[TestFeatureIntervalCollection.feat1, TestFeatureIntervalCollection.feat2],
                    feature_collection_id="featgrp1",
                ),
                dict(feature_intervals=[TestFeatureIntervalCollection.feat3], feature_collection_id="featgrp2"),
            ]
        )
    )

    empty_annot = AnnotationCollectionModel.Schema().load(
        dict(
            feature_collections=[],
            genes=[],
        )
    )

    def test_annotation(self):
        obj = self.annot.to_annotation_collection()
        model = AnnotationCollectionModel.from_annotation_collection(obj)

        # remove guids to make comparison work
        for feat_grp in model.feature_collections:
            for item in feat_grp.feature_intervals:
                item.feature_interval_guid = None
            feat_grp.feature_collection_guid = None

        # remove guids to make comparison work
        for gene in model.genes:
            for item in gene.transcripts:
                item.transcript_interval_guid = None
            gene.gene_guid = None

        assert model == self.annot

    def test_annot_no_range(self):
        obj = self.annot_no_range.to_annotation_collection()
        assert obj.location == SingleInterval(2, 30, Strand.PLUS)

    def test_annot_no_features(self):
        obj = self.annot_no_features.to_annotation_collection()
        assert len(obj.feature_collections) == 0

    def test_annot_no_genes(self):
        obj = self.annot_no_genes.to_annotation_collection()
        assert len(obj.genes) == 0

    def test_empty_annot(self):
        obj = self.empty_annot.to_annotation_collection()
        assert obj.is_empty

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected",
        [
            (
                None,
                None,
                False,
                True,
                #
                {"featgrp1", "featgrp2", "gene1"},
            ),
            (
                0,
                None,
                False,
                True,
                {"featgrp1", "featgrp2", "gene1"},
            ),
            (
                None,
                30,
                False,
                True,
                {"featgrp1", "featgrp2", "gene1"},
            ),
            (0, 0, False, True, {}),
            (
                0,
                20,
                False,
                True,
                {"featgrp1", "gene1"},
            ),
            (
                None,
                20,
                False,
                True,
                {"featgrp1", "gene1"},
            ),
            (25, None, False, True, {"featgrp2"}),
            (26, None, False, True, {}),
            (26, None, False, False, {"featgrp2"}),
            (
                0,
                3,
                False,
                False,
                {"featgrp1", "gene1"},
            ),
            (0, None, True, False, {"gene1"}),
        ],
    )
    def test_position_queries(self, start, end, coding_only, completely_within, expected):
        obj = self.annot.to_annotation_collection(parent_genome)
        r = obj.query_by_position(start, end, coding_only, completely_within)
        if r.is_empty:
            assert len(expected) == 0
        else:
            assert set.union(*[x.identifiers for x in r]) == expected
        for gene in r:
            orig_gene = next(obj.query_by_feature_identifiers(gene.identifiers).iter_children())
            assert orig_gene.get_primary_feature_sequence() == gene.get_primary_feature_sequence()

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected",
        [
            (None, None, False, True, SingleInterval(0, 30, Strand.PLUS)),
            (0, None, False, True, SingleInterval(0, 30, Strand.PLUS)),
            (None, 30, False, True, SingleInterval(0, 30, Strand.PLUS)),
            (0, 0, False, True, SingleInterval(0, 0, Strand.PLUS)),
            (0, 20, False, True, SingleInterval(0, 20, Strand.PLUS)),
            (None, 20, False, True, SingleInterval(0, 20, Strand.PLUS)),
            (25, None, False, True, SingleInterval(25, 30, Strand.PLUS)),
            (26, None, False, True, SingleInterval(26, 30, Strand.PLUS)),
            (0, 1, False, False, SingleInterval(0, 1, Strand.PLUS)),
        ],
    )
    def test_position_queries_location(self, start, end, completely_within, coding_only, expected):
        obj = self.annot.to_annotation_collection()
        r = obj.query_by_position(start, end, coding_only, completely_within)
        assert r.location == expected

    def test_query_position_exceptions(self):
        obj = self.annot.to_annotation_collection()
        with pytest.raises(InvalidQueryError):
            _ = obj.query_by_position(-1, 10)
        with pytest.raises(InvalidQueryError):
            _ = obj.query_by_position(15, 10)

    @pytest.mark.parametrize(
        "ids",
        (
            {"gene1"},
            {"gene1", "featgrp1"},
            {},
        ),
    )
    def test_query_by_identifiers(self, ids):
        obj = self.annot.to_annotation_collection()
        r = obj.query_by_feature_identifiers(ids)
        if r.is_empty:
            assert len(ids) == 0
        else:
            assert {x.gene_id for x in r.genes} | {x.feature_collection_id for x in r.feature_collections} == ids

    @pytest.mark.parametrize(
        "i",
        (
            "gene1",
            "featgrp1",
        ),
    )
    def test_query_by_identifiers_str(self, i):
        obj = self.annot.to_annotation_collection()
        r = obj.query_by_feature_identifiers(i)
        assert {x.gene_id for x in r.genes} | {x.feature_collection_id for x in r.feature_collections} == {i}

    def test_query_by_identifiers_with_extraneous(self):
        obj = self.annot.to_annotation_collection()
        r = obj.query_by_feature_identifiers(["gene1", "abc"])
        assert len(r.genes) == 1 and r.genes[0].gene_id == "gene1"

    def test_extract_sequence(self):
        obj = self.annot.to_annotation_collection(parent_or_seq_chunk_parent=parent_genome)
        seq = obj.get_reference_sequence()
        assert str(seq) == genome[:30]

    def test_query_by_ids(self):
        obj = self.annot.to_annotation_collection()
        my_ids = list(obj.hierarchical_children_guids.keys())

        # query them all
        assert obj.query_by_guids(my_ids).children_guids == set(my_ids)
        # query one
        assert obj.query_by_guids([my_ids[0]]).children_guids == {my_ids[0]}
        # query none
        assert obj.query_by_guids([]).children_guids == set()
