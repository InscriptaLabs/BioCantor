from copy import deepcopy

import pytest
from inscripta.biocantor.exc import (
    InvalidAnnotationError,
    NoncodingTranscriptError,
    InvalidQueryError,
    NoSuchAncestorException,
    ValidationException,
    NullSequenceException,
)
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds_frame import CDSFrame
from inscripta.biocantor.io.models import (
    GeneIntervalModel,
    AnnotationCollectionModel,
    FeatureIntervalCollectionModel,
    TranscriptIntervalModel,
)
from inscripta.biocantor.location.location_impl import SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent, SequenceType
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence

genome = "TTTTTTTTTTAAGTATTCTTGGACCTAATTAAAAAAAAAAAAAAAAAAACCCCC"
parent_genome = Parent(
    id="genome", sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type=SequenceType.CHROMOSOME
)
parent_genome_with_location = Parent(
    id="genome",
    sequence=Sequence(genome, Alphabet.NT_STRICT),
    sequence_type=SequenceType.CHROMOSOME,
    location=SingleInterval(0, len(genome), Strand.PLUS),
)
parent_genome_10_49 = Parent(
    id="genome_10_49",
    sequence=Sequence(
        genome[10:49],
        Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            location=SingleInterval(
                10, 49, Strand.PLUS, parent=Parent(id="genome", sequence_type=SequenceType.CHROMOSOME)
            )
        ),
    ),
)

parent_no_seq = Parent(sequence_type=SequenceType.CHROMOSOME, id="genome")
parent_nonstandard_type = Parent(sequence_type="SomeOtherType", id="genome")
parent_nonstandard_type_with_sequence = Parent(
    sequence=Sequence(genome, Alphabet.NT_STRICT), sequence_type="SomeOtherType", id="genome"
)


class TestGene:
    """Test basic gene construction from Transcripts"""

    tx1 = dict(
        exon_starts=[12],
        exon_ends=[28],
        strand=Strand.PLUS.name,
        cds_starts=[15],
        cds_ends=[19],
        cds_frames=[CDSFrame.ZERO.name],
        transcript_symbol="tx1",
    )
    tx2 = dict(
        exon_starts=[12, 17, 22],
        exon_ends=[16, 20, 25],
        strand=Strand.PLUS.name,
        cds_starts=[14, 17, 22],
        cds_ends=[16, 20, 23],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        transcript_symbol="tx2",
    )
    tx1_primary = dict(
        exon_starts=[12],
        exon_ends=[28],
        strand=Strand.PLUS.name,
        cds_starts=[15],
        cds_ends=[19],
        cds_frames=[CDSFrame.ZERO.name],
        is_primary_tx=True,
    )
    tx2_primary = dict(
        exon_starts=[12, 17, 22],
        exon_ends=[16, 20, 25],
        strand=Strand.PLUS.name,
        cds_starts=[14, 17, 22],
        cds_ends=[16, 20, 23],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        is_primary_tx=True,
    )
    tx_noncoding = dict(exon_starts=[12, 17, 22], exon_ends=[16, 20, 40], strand=Strand.PLUS.name)
    tx_noncoding_short = dict(exon_starts=[12], exon_ends=[16], strand=Strand.PLUS.name)

    gene = GeneIntervalModel.Schema().load(
        dict(transcripts=[tx1, tx2], qualifiers={}, gene_type=Biotype.protein_coding.name, gene_id="gene1")
    )

    gene_noncoding = GeneIntervalModel.Schema().load(dict(transcripts=[tx_noncoding, tx_noncoding_short]))

    def test_primary_inference(self):
        obj = self.gene.to_gene_interval()
        assert obj.primary_transcript == TranscriptIntervalModel.Schema().load(self.tx2).to_transcript_interval()

    def test_merged_interval(self):
        obj = self.gene.to_gene_interval()
        assert str(obj.get_merged_transcript()) == "FeatureInterval((12-28:+), name=None)"
        assert str(obj.get_merged_cds()) == "FeatureInterval((14-20:+, 22-23:+), name=None)"

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

    def test_query_by_guid(self):
        # query by all
        obj = self.gene.to_gene_interval(parent_or_seq_chunk_parent=parent_genome)
        assert obj.query_by_guids(list(obj.guid_map.keys())) == obj
        # query by none produces invalid gene
        with pytest.raises(InvalidAnnotationError):
            _ = obj.query_by_guids([])
        # query one
        guids = list(obj.guid_map.keys())[:1]
        assert {x.guid for x in obj.query_by_guids(guids)} == set(guids)

    def test_iterator(self):
        obj = self.gene.to_gene_interval(parent_or_seq_chunk_parent=parent_genome)
        assert list(obj.iter_children()) == list(obj)
        assert len(list(obj)) == 2

    def test_nonstandard_parents(self):
        obj0 = self.gene.to_gene_interval(parent_genome)
        obj1 = self.gene.to_gene_interval(parent_no_seq)
        obj2 = self.gene.to_gene_interval(parent_nonstandard_type)
        obj3 = self.gene.to_gene_interval(parent_nonstandard_type_with_sequence)
        with pytest.raises(NullSequenceException):
            _ = obj1.get_reference_sequence()
        with pytest.raises(NullSequenceException):
            _ = obj2.get_reference_sequence()
        assert obj0.get_reference_sequence() == obj3.get_reference_sequence()


class TestFeatureIntervalCollection:
    feat1 = dict(
        interval_starts=[12],
        interval_ends=[15],
        strand=Strand.PLUS.name,
        feature_types=["a", "b"],
        feature_name="feat1",
    )
    feat2 = dict(
        interval_starts=[12, 17, 22],
        interval_ends=[16, 20, 25],
        strand=Strand.PLUS.name,
        feature_types=["b"],
        feature_name="feat2",
    )
    feat3 = dict(
        interval_starts=[35], interval_ends=[40], strand=Strand.MINUS.name, feature_types=["a"], feature_name="feat3"
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

    def test_query_by_guid(self):
        # query by all
        obj = self.collection1.to_feature_collection()
        assert obj.query_by_guids(list(obj.guid_map.keys())) == obj
        # query by none produces invalid gene
        with pytest.raises(InvalidAnnotationError):
            _ = obj.query_by_guids([])
        # query one
        guids = list(obj.guid_map.keys())[:1]
        assert {x.guid for x in obj.query_by_guids(guids)} == set(guids)

    def test_iterator(self):
        obj = self.collection1.to_feature_collection()
        assert list(obj.iter_children()) == list(obj)
        assert len(list(obj)) == 2

    def test_nonstandard_parents(self):
        obj0 = self.collection1.to_feature_collection(parent_genome)
        obj1 = self.collection1.to_feature_collection(parent_no_seq)
        obj2 = self.collection1.to_feature_collection(parent_nonstandard_type)
        obj3 = self.collection1.to_feature_collection(parent_nonstandard_type_with_sequence)
        with pytest.raises(NullSequenceException):
            _ = obj1.get_reference_sequence()
        with pytest.raises(NullSequenceException):
            _ = obj2.get_reference_sequence()
        assert obj0.get_reference_sequence() == obj3.get_reference_sequence()


class TestAnnotationCollection:
    annot_full_range = AnnotationCollectionModel.Schema().load(
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
            end=54,
        )
    )

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
            start=2,
            end=40,
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
        assert obj._location == SingleInterval(12, 40, Strand.PLUS)

    def test_annot_no_features(self):
        obj = self.annot_no_features.to_annotation_collection()
        assert len(obj.feature_collections) == 0

    def test_annot_no_genes(self):
        obj = self.annot_no_genes.to_annotation_collection()
        assert len(obj.genes) == 0

    def test_empty_annot(self):
        obj = self.empty_annot.to_annotation_collection()
        assert obj.is_empty

    def test_equivalent_sequences(self):
        """Regardless of which start/end position was used, sequences are invariant"""
        obj = self.annot.to_annotation_collection(parent_genome)
        obj2 = self.annot_full_range.to_annotation_collection(parent_genome)
        obj3 = self.annot_no_range.to_annotation_collection(parent_genome)
        for gene in obj:
            gene2 = next(obj2.query_by_feature_identifiers(gene.identifiers).iter_children())
            gene3 = next(obj3.query_by_feature_identifiers(gene.identifiers).iter_children())
            for tx1 in gene:
                tx2 = gene2.guid_map[tx1.guid]
                assert tx1.get_spliced_sequence() == tx2.get_spliced_sequence()
                tx3 = gene3.guid_map[tx1.guid]
                assert tx1.get_spliced_sequence() == tx3.get_spliced_sequence()

    def test_equivalent_sequences_subset(self):
        """Regardless of which start/end position was used, sequences are invariant"""
        obj = self.annot.to_annotation_collection(parent_genome_10_49)
        obj2 = self.annot_full_range.to_annotation_collection(parent_genome_10_49)
        obj3 = self.annot_no_range.to_annotation_collection(parent_genome_10_49)
        for gene in obj:
            gene2 = next(obj2.query_by_feature_identifiers(gene.identifiers).iter_children())
            gene3 = next(obj3.query_by_feature_identifiers(gene.identifiers).iter_children())
            for tx1 in gene:
                tx2 = gene2.guid_map[tx1.guid]
                assert tx1.get_spliced_sequence() == tx2.get_spliced_sequence()
                tx3 = gene3.guid_map[tx1.guid]
                assert tx1.get_spliced_sequence() == tx3.get_spliced_sequence()

    def test_equivalent_sequences_all(self):
        """Regardless of which start/end position was used, sequences are invariant"""
        obj = self.annot.to_annotation_collection(parent_genome_10_49)
        obj2 = self.annot_full_range.to_annotation_collection(parent_genome_10_49)
        obj3 = self.annot_no_range.to_annotation_collection(parent_genome_10_49)
        obj4 = self.annot.to_annotation_collection(parent_genome)
        obj5 = self.annot_full_range.to_annotation_collection(parent_genome)
        obj6 = self.annot_no_range.to_annotation_collection(parent_genome)
        obj7 = self.annot.to_annotation_collection(parent_genome_with_location)
        obj8 = self.annot_full_range.to_annotation_collection(parent_genome_with_location)
        obj9 = self.annot_no_range.to_annotation_collection(parent_genome_with_location)

        objects = [obj, obj2, obj3, obj4, obj5, obj6, obj7, obj8, obj9]

        assert len({len(x) for x in objects}) == 1
        for gene in obj:
            for other_obj in objects[1:]:
                ogene = next(other_obj.query_by_feature_identifiers(gene.identifiers).iter_children())
                for tx1 in gene:
                    tx2 = ogene.guid_map[tx1.guid]
                    # could be a subset now
                    if len(tx1) < len(tx2):
                        assert str(tx1.get_spliced_sequence()) in str(tx2.get_spliced_sequence())
                    elif len(tx1) == len(tx2):
                        assert str(tx1.get_spliced_sequence()) == str(tx2.get_spliced_sequence())
                    else:
                        # cannot be longer, since obj is the subset here
                        assert False

    def test_annotation_collection_bounds(self):
        """Bounds will be inferred from children, depending on if a parent is provided"""
        obj = self.annot_no_range.to_annotation_collection()
        assert obj.start == 12 and obj.end == 40
        # if we provide a parent with no location information, then we do not shift boundaries
        obj = self.annot_no_range.to_annotation_collection(parent_genome)
        assert obj.start == 12 and obj.end == 40
        # but if we provide a parent with a location, then we retain those boundaries
        obj = self.annot_no_range.to_annotation_collection(parent_genome_with_location)
        assert obj.start == 0 and obj.end == 54
        obj = self.annot_no_range.to_annotation_collection(parent_genome_10_49)
        assert obj.start == 10 and obj.end == 49

    def test_annotation_bounds_exceptions(self):
        # hack a copy to avoid inconsistent state
        annot_no_range = deepcopy(self.annot_no_range)
        annot_no_range.start = 0
        with pytest.raises(InvalidAnnotationError):
            _ = annot_no_range.to_annotation_collection()
        annot_no_range.start = None
        annot_no_range.end = 10
        with pytest.raises(InvalidAnnotationError):
            _ = annot_no_range.to_annotation_collection()

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
                40,
                False,
                True,
                {"featgrp1", "featgrp2", "gene1"},
            ),
            (
                0,
                30,
                False,
                True,
                {"featgrp1", "gene1"},
            ),
            (
                None,
                30,
                False,
                True,
                {"featgrp1", "gene1"},
            ),
            (35, None, False, True, {"featgrp2"}),
            (36, None, False, True, {}),
            (36, None, False, False, {"featgrp2"}),
            (
                10,
                13,
                False,
                False,
                {"featgrp1", "gene1"},
            ),
            (0, None, True, False, {"gene1"}),
            (15, 30, False, False, {"featgrp1", "gene1"}),
            (10, 35, False, False, {"featgrp1", "gene1"}),
        ],
    )
    def test_position_queries(self, start, end, coding_only, completely_within, expected):
        obj = self.annot_full_range.to_annotation_collection(parent_genome)
        r = obj.query_by_position(start, end, coding_only, completely_within)
        if r.is_empty:
            assert len(expected) == 0
        else:
            assert set.union(*[x.identifiers for x in r]) == expected
        for gene_or_feature in r:
            orig_gene_or_feature = obj.guid_map[gene_or_feature.guid]
            for new_tx_or_feature in gene_or_feature:
                orig_tx_or_feature = orig_gene_or_feature.guid_map[new_tx_or_feature.guid]
                if len(new_tx_or_feature) == new_tx_or_feature.chunk_relative_size:
                    assert new_tx_or_feature.get_spliced_sequence() == orig_tx_or_feature.get_spliced_sequence()
                else:
                    assert str(new_tx_or_feature.get_spliced_sequence()) in str(
                        orig_tx_or_feature.get_spliced_sequence()
                    )

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
                None,
                40,
                False,
                True,
                {"featgrp1", "featgrp2", "gene1"},
            ),
            (
                10,
                30,
                False,
                True,
                {"featgrp1", "gene1"},
            ),
            (
                None,
                30,
                False,
                True,
                {"featgrp1", "gene1"},
            ),
            (35, None, False, True, {"featgrp2"}),
            (
                12,
                13,
                False,
                False,
                {"featgrp1", "gene1"},
            ),
            (10, None, True, False, {"gene1"}),
            (15, 30, False, False, {"featgrp1", "gene1"}),
            (10, 35, False, False, {"featgrp1", "gene1"}),
        ],
    )
    def test_position_queries_chunk_parent(self, start, end, coding_only, completely_within, expected):
        obj = self.annot_no_range.to_annotation_collection(parent_genome_10_49)
        r = obj.query_by_position(start, end, coding_only, completely_within)
        if r.is_empty:
            assert len(expected) == 0
        else:
            assert set.union(*[x.identifiers for x in r]) == expected
        for gene_or_feature in r:
            orig_gene_or_feature = obj.guid_map[gene_or_feature.guid]
            for new_tx_or_feature in gene_or_feature:
                orig_tx_or_feature = orig_gene_or_feature.guid_map[new_tx_or_feature.guid]
                if len(new_tx_or_feature) == new_tx_or_feature.chunk_relative_size:
                    assert new_tx_or_feature.get_spliced_sequence() == orig_tx_or_feature.get_spliced_sequence()
                else:
                    assert str(new_tx_or_feature.get_spliced_sequence()) in str(
                        orig_tx_or_feature.get_spliced_sequence()
                    )

    @pytest.mark.parametrize(
        "annot,parent",
        [(annot, parent) for annot in [annot, annot_full_range] for parent in [parent_genome_10_49, parent_genome]],
    )
    def test_nested_position_queries(self, annot, parent):
        obj = annot.to_annotation_collection(parent)
        r = obj.query_by_position(10, 35, completely_within=False)
        assert len(r) == 2
        assert r._location.parent
        for gene in r:
            orig_gene = next(obj.query_by_feature_identifiers(gene.identifiers).iter_children())
            for tx1 in gene:
                tx2 = orig_gene.guid_map[tx1.guid]
                assert tx1.get_spliced_sequence() == tx2.get_spliced_sequence()
        r2 = r.query_by_position(10, 20, completely_within=False)
        assert len(r2) == 2
        assert r2._location.parent
        # this slice cut some of transcripts into chunks, so now the sequences are a subset
        for gene in r2:
            orig_gene = next(obj.query_by_feature_identifiers(gene.identifiers).iter_children())
            for tx1 in gene:
                tx2 = orig_gene.guid_map[tx1.guid]
                assert str(tx1.get_spliced_sequence()) in str(tx2.get_spliced_sequence())
        r3 = r.query_by_position(10, 18, completely_within=False)
        assert r3._location.parent
        assert len(r3) == 2
        for gene in r3:
            orig_gene = next(obj.query_by_feature_identifiers(gene.identifiers).iter_children())
            for tx1 in gene:
                tx2 = orig_gene.guid_map[tx1.guid]
                assert str(tx1.get_spliced_sequence()) in str(tx2.get_spliced_sequence())

    def test_nested_position_query_out_of_bounds(self):
        obj = self.annot.to_annotation_collection(parent_genome)
        r = obj.query_by_position(10, 35, completely_within=False)
        with pytest.raises(InvalidQueryError):
            _ = r.query_by_position(10, 40)

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected",
        [
            (None, None, False, True, SingleInterval(0, 54, Strand.PLUS)),
            (0, None, False, True, SingleInterval(0, 54, Strand.PLUS)),
            (None, 30, False, True, SingleInterval(0, 30, Strand.PLUS)),
            (0, 20, False, True, SingleInterval(0, 20, Strand.PLUS)),
            (None, 20, False, True, SingleInterval(0, 20, Strand.PLUS)),
            (25, None, False, True, SingleInterval(25, 54, Strand.PLUS)),
            (26, None, False, True, SingleInterval(26, 54, Strand.PLUS)),
            (0, 1, False, False, SingleInterval(0, 1, Strand.PLUS)),
        ],
    )
    def test_position_queries_location_no_bounds(self, start, end, completely_within, coding_only, expected):
        obj = self.annot_full_range.to_annotation_collection()
        r = obj.query_by_position(start, end, coding_only, completely_within)
        assert r._location == expected

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected",
        [
            (None, None, False, True, SingleInterval(2, 40, Strand.PLUS)),
            (2, None, False, True, SingleInterval(2, 40, Strand.PLUS)),
            (None, 30, False, True, SingleInterval(2, 30, Strand.PLUS)),
            (2, 20, False, True, SingleInterval(2, 20, Strand.PLUS)),
            (None, 20, False, True, SingleInterval(2, 20, Strand.PLUS)),
            (25, None, False, True, SingleInterval(25, 40, Strand.PLUS)),
            (26, None, False, True, SingleInterval(26, 40, Strand.PLUS)),
            (2, 3, False, False, SingleInterval(2, 3, Strand.PLUS)),
        ],
    )
    def test_position_queries_location(self, start, end, completely_within, coding_only, expected):
        obj = self.annot.to_annotation_collection()
        r = obj.query_by_position(start, end, coding_only, completely_within)
        assert r._location == expected

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected",
        [
            (None, None, False, True, SingleInterval(12, 40, Strand.PLUS)),
            (12, None, False, True, SingleInterval(12, 40, Strand.PLUS)),
            (None, 30, False, True, SingleInterval(12, 30, Strand.PLUS)),
            (12, 20, False, True, SingleInterval(12, 20, Strand.PLUS)),
            (None, 20, False, True, SingleInterval(12, 20, Strand.PLUS)),
            (25, None, False, True, SingleInterval(25, 40, Strand.PLUS)),
            (26, None, False, True, SingleInterval(26, 40, Strand.PLUS)),
            (12, 13, False, False, SingleInterval(12, 13, Strand.PLUS)),
        ],
    )
    def test_position_queries_location_inferred(self, start, end, completely_within, coding_only, expected):
        obj = self.annot_no_range.to_annotation_collection()
        r = obj.query_by_position(start, end, coding_only, completely_within)
        assert r._location == expected

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,min_start,max_end",
        [
            (12, 13, False, False, 12, 28),
            (12, 15, False, False, 12, 28),
        ],
    )
    def test_position_queries_expanded_range(
        self,
        start,
        end,
        completely_within,
        coding_only,
        min_start,
        max_end,
    ):
        """The span of the AnnotationCollection object is the original query, but the full range
        of the child objects are larger"""
        obj = self.annot_no_range.to_annotation_collection()
        r = obj.query_by_position(start, end, coding_only, completely_within)
        assert min(x.start for x in r) == min_start
        assert max(x.end for x in r) == max_end

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected_identifiers",
        [
            # query removes everything due to completely_within=True
            (21, 22, False, True, set()),
            # feat1 and feat3 lost due to no overlap with query
            (21, 22, False, False, {"feat2", "tx1", "tx2"}),
            # tx1 hits 28, feat3 starts at 35, so nothing here
            (28, 35, False, False, set()),
            # moving to 36 now retains feat3
            (28, 36, False, False, {"feat3"}),
            # moving to 27 now retains tx1
            (27, 36, False, False, {"tx1", "feat3"}),
            # moving to 24 now retains all but feat1
            (24, 36, False, False, {"tx1", "tx2", "feat2", "feat3"}),
        ],
    )
    def test_position_queries_lose_isoforms(
        self,
        start,
        end,
        completely_within,
        coding_only,
        expected_identifiers,
    ):
        obj = self.annot_no_range.to_annotation_collection()
        r = obj.query_by_position(start, end, coding_only, completely_within)
        if len(r) == 0:
            assert expected_identifiers == set()
        else:
            assert set.union(*(y.identifiers for x in r for y in x)) == expected_identifiers

    def test_query_position_exceptions(self):
        obj = self.annot.to_annotation_collection()
        with pytest.raises(InvalidQueryError):
            _ = obj.query_by_position(-1, 10)
        with pytest.raises(InvalidQueryError):
            _ = obj.query_by_position(15, 10)
        with pytest.raises(InvalidQueryError):
            _ = obj.query_by_position(0, 55)
        with pytest.raises(InvalidQueryError):
            _ = obj.query_by_position(10, 10)

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
        assert str(seq) == genome[2:40]
        # with inferred range, now sequence is cut down based on bounds of transcripts
        obj = self.annot_no_range.to_annotation_collection(parent_or_seq_chunk_parent=parent_genome)
        seq = obj.get_reference_sequence()
        assert str(seq) == genome[12:40]
        # however, if we use a parent whose location is provided, then we retain more information
        obj = self.annot_no_range.to_annotation_collection(parent_or_seq_chunk_parent=parent_genome_with_location)
        seq = obj.get_reference_sequence()
        assert str(seq) == genome

    def test_query_by_ids(self):
        obj = self.annot.to_annotation_collection()
        my_ids = list(obj.hierarchical_children_guids.keys())

        # query them all
        assert obj.query_by_guids(my_ids).children_guids == set(my_ids)
        # query one
        assert obj.query_by_guids([my_ids[0]]).children_guids == {my_ids[0]}
        # query none
        assert obj.query_by_guids([]).children_guids == set()

    def test_gff3_export(self, test_data_dir):
        obj = self.annot.to_annotation_collection()
        # populate sequence names; normally this is done via the model constructors
        obj.sequence_name = "chr1"
        for item in obj:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
        with open(test_data_dir / "collection_gff3_export_chromosome_coordinates.gff") as fh:
            assert fh.read() == "\n".join(str(x) for x in obj.to_gff())

    def test_gff3_export_chunk_relative(self, test_data_dir):
        obj = self.annot.to_annotation_collection(parent_genome_10_49)
        # populate sequence names; normally this is done via the model constructors
        obj.sequence_name = "chr1"
        for item in obj:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
        with open(test_data_dir / "collection_gff3_export_chunk_relative.gff") as fh:
            assert fh.read() == "\n".join(str(x) for x in obj.to_gff(chromosome_relative_coordinates=False))

    def test_gff3_export_exception(self, test_data_dir):
        """Cannot export to GFF3 in relative coordinates without having sequence."""
        obj = self.annot.to_annotation_collection()
        obj.sequence_name = "chr1"
        for item in obj:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
        with pytest.raises(NoSuchAncestorException):
            _ = "\n".join(str(x) for x in obj.to_gff(chromosome_relative_coordinates=False))

    def test_reset_parent_noop(self):
        obj = self.annot.to_annotation_collection()
        # no-op
        obj._reset_parent()
        # equivalent
        obj._reset_parent(None)

    def test_reset_parent_null(self):
        obj = self.annot.to_annotation_collection(parent_genome)
        for child in obj:
            assert child._location.parent
        obj._reset_parent()
        for child in obj:
            assert not child._location.parent

    def test_reset_parent(self):
        obj = self.annot.to_annotation_collection(parent_genome)
        obj2 = obj.query_by_position(20, 40)
        obj2._reset_parent(parent_genome)
        # the coordinates are now broken, so the sequences are wrong
        for rec in obj2:
            orig_rec = next(obj.query_by_guids([rec.guid]).__iter__()).feature_intervals[0]
            assert orig_rec.get_spliced_sequence() != rec.feature_intervals[0].get_spliced_sequence()

    def test_iterator(self):
        obj = self.annot.to_annotation_collection(parent_genome)
        assert list(obj.iter_children()) == list(obj)
        assert len(list(obj)) == 3

    def test_nonstandard_parents(self):
        obj0 = self.annot.to_annotation_collection(parent_genome)
        obj1 = self.annot.to_annotation_collection(parent_no_seq)
        obj2 = self.annot.to_annotation_collection(parent_nonstandard_type)
        obj3 = self.annot.to_annotation_collection(parent_nonstandard_type_with_sequence)
        with pytest.raises(NullSequenceException):
            _ = obj1.get_reference_sequence()
        with pytest.raises(NullSequenceException):
            _ = obj2.get_reference_sequence()
        assert obj0.get_reference_sequence() == obj3.get_reference_sequence()

        assert obj0.chromosome_location == obj0.chunk_relative_location
        assert obj1.chromosome_location == obj1.chunk_relative_location
        assert obj2.chromosome_location == obj2.chunk_relative_location
        assert obj3.chromosome_location == obj3.chunk_relative_location
        # OTOH, this is not the same
        obj4 = self.annot.to_annotation_collection(parent_genome_10_49)
        assert obj4.chromosome_location != obj4.chunk_relative_location

    def test_lift_to_new_coordinates(self):
        obj0 = self.annot.to_annotation_collection(parent_genome)
        obj1 = obj0.liftover_to_parent_or_seq_chunk_parent(parent_genome_10_49)
        assert str(obj1.get_reference_sequence()) in str(obj0.get_reference_sequence())
        assert obj1.start == 2 and obj1.end == 40
        for gene in obj1:
            orig_gene = next(obj0.query_by_feature_identifiers(gene.identifiers).iter_children())
            orig_tx_or_feat = next(orig_gene.iter_children())
            tx_or_feat = next(gene.iter_children())
            assert str(orig_tx_or_feat.get_spliced_sequence()) == str(tx_or_feat.get_spliced_sequence())
