import pickle
from copy import deepcopy
from uuid import UUID

import pytest

from biocantor.exc import (
    InvalidAnnotationError,
    NoncodingTranscriptError,
    InvalidQueryError,
    NoSuchAncestorException,
    ValidationException,
    NullSequenceException,
)
from biocantor.io.gff3.exc import GTFExportException
from biocantor.gene.biotype import Biotype
from biocantor.gene.cds_frame import CDSFrame
from biocantor.gene.collections import AnnotationCollection
from biocantor.io.models import (
    GeneIntervalModel,
    AnnotationCollectionModel,
    FeatureIntervalCollectionModel,
    TranscriptIntervalModel,
)
from biocantor.io.parser import seq_chunk_to_parent
from biocantor.location.location_impl import SingleInterval
from biocantor.location.strand import Strand
from biocantor.parent.parent import Parent, SequenceType
from biocantor.sequence.alphabet import Alphabet
from biocantor.sequence.sequence import Sequence

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

genome_rev = "GGGGGTTTTTTTTTTTTTTTTTTTAATTAGGTCCAAGAATACTTAAAAAAAAAA"
parent_genome_rev = Parent(
    id="genome_rev",
    sequence_type=SequenceType.SEQUENCE_CHUNK,
    strand=None,
    location=None,
    sequence=Sequence(
        data=genome_rev,
        id="genome_rev",
        alphabet=Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            id="genome_rev",
            sequence_type=SequenceType.CHROMOSOME,
            strand=Strand.MINUS,
            location=SingleInterval(
                0,
                54,
                Strand.MINUS,
                parent=Parent(
                    id="genome_rev",
                    sequence_type=SequenceType.CHROMOSOME,
                    strand=Strand.MINUS,
                    location=SingleInterval(0, 54, Strand.MINUS),
                    sequence=None,
                    parent=None,
                ),
            ),
            sequence=None,
            parent=Parent(
                id="TestSeq",
                sequence_type=SequenceType.CHROMOSOME,
                strand=Strand.MINUS,
                location=SingleInterval(0, 54, Strand.MINUS, parent=None),
                sequence=None,
                parent=None,
            ),
        ),
    ),
)

parent_genome_rev_5_44 = Parent(
    id="genome_rev_5_44",
    sequence_type=SequenceType.SEQUENCE_CHUNK,
    strand=None,
    location=None,
    sequence=Sequence(
        data=genome_rev[5:44],
        id="genome_rev_5_44",
        alphabet=Alphabet.NT_EXTENDED_GAPPED,
        type=SequenceType.SEQUENCE_CHUNK,
        parent=Parent(
            id="genome_rev",
            sequence_type=SequenceType.CHROMOSOME,
            strand=Strand.MINUS,
            location=SingleInterval(
                10,
                49,
                Strand.MINUS,
                parent=Parent(
                    id="genome_rev",
                    sequence_type=SequenceType.CHROMOSOME,
                    strand=Strand.MINUS,
                    location=SingleInterval(10, 49, Strand.MINUS),
                    sequence=None,
                    parent=None,
                ),
            ),
            sequence=None,
            parent=Parent(
                id="TestSeq",
                sequence_type=SequenceType.CHROMOSOME,
                strand=Strand.MINUS,
                location=SingleInterval(10, 49, Strand.MINUS, parent=None),
                sequence=None,
                parent=None,
            ),
        ),
    ),
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
        assert (
            str(obj.get_merged_transcript()) == str(obj.get_merged_feature()) == "FeatureInterval((12-28:+), name=None)"
        )
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
        # query by none produces None
        assert obj.query_by_guids([]) is None
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

    def test_equality_different_parents(self):
        obj1 = self.gene.to_gene_interval(parent_genome_10_49)
        obj2 = self.gene.to_gene_interval(parent_genome_rev_5_44)
        assert obj1 != obj2
        assert hash(obj1) != hash(obj2)


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

    def test_merged_interval(self):
        obj = self.collection1.to_feature_collection()
        assert str(obj.get_merged_feature()) == "FeatureInterval((12-16:+, 17-20:+, 22-25:+), name=featgrp1)"

    def test_query_by_guid(self):
        # query by all
        obj = self.collection1.to_feature_collection()
        assert obj.query_by_guids(list(obj.guid_map.keys())) == obj
        # query by none produces None
        assert obj.query_by_guids([]) is None
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

    def test_equality_different_parents(self):
        obj1 = self.collection1.to_feature_collection(parent_genome_10_49)
        obj2 = self.collection1.to_feature_collection(parent_genome_rev_5_44)
        assert obj1 != obj2
        assert hash(obj1) != hash(obj2)


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
                if new_tx_or_feature.chunk_relative_location.is_empty:
                    continue
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
                # if there was no overlap between the query coordinates and an isoform (it was entirely intronic)
                # then the sequence information is inherently lost
                if new_tx_or_feature.chunk_relative_location.is_empty:
                    continue
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
        "start,end,coding_only,completely_within,expanded_start,expanded_end",
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
        expanded_start,
        expanded_end,
    ):
        """The span of the AnnotationCollection object is the original query, but the full range
        of the child objects are larger"""
        obj = self.annot_no_range.to_annotation_collection(parent_genome)
        r = obj.query_by_position(start, end, coding_only, completely_within, expand_location_to_children=True)
        assert min(x.chromosome_location.start for x in r) == expanded_start
        assert max(x.chromosome_location.end for x in r) == expanded_end
        # without expansion the children are all sliced; but this can only be noticed
        # if you look at the hidden _chunk_relative_bounded_chromosome_location
        r2 = obj.query_by_position(start, end, coding_only, completely_within, expand_location_to_children=False)
        assert min(x.chromosome_location.start for x in r2) == expanded_start
        assert max(x.chromosome_location.end for x in r2) == expanded_end
        assert min(x._chunk_relative_bounded_chromosome_location.start for x in r2) == start
        assert max(x._chunk_relative_bounded_chromosome_location.end for x in r2) == end

    def test_position_queries_expanded_range_exception(self):
        """If expand range is true, then an exception is thrown on sequence chunks if they expand beyond the chunk"""
        obj = self.annot_no_range.to_annotation_collection(parent_genome)
        # first query by range to produce a chunk-relative collection
        obj2 = obj.query_by_position(12, 20, completely_within=False)
        with pytest.raises(InvalidQueryError):
            _ = obj2.query_by_position(15, 20, expand_location_to_children=True, completely_within=False)
        # this is fine without range expansion
        _ = obj2.query_by_position(15, 20, completely_within=False)

        # this is fine without sequence info as well
        obj = self.annot_no_range.to_annotation_collection()
        obj2 = obj.query_by_position(12, 20, completely_within=False)
        _ = obj2.query_by_position(15, 20, expand_location_to_children=True, completely_within=False)

    @pytest.mark.parametrize(
        "start,end,expected_translations,expected_mrna",
        [
            (12, 13, {"tx1": "F", "tx2": "IL"}, {"tx1": "GTATTCTTGGACCTAA", "tx2": "GTATCTTACC"}),
            (12, 20, {"tx1": "F", "tx2": "IL"}, {"tx1": "GTATTCTTGGACCTAA", "tx2": "GTATCTTACC"}),
            (12, 25, {"tx1": "F", "tx2": "IL"}, {"tx1": "GTATTCTTGGACCTAA", "tx2": "GTATCTTACC"}),
        ],
    )
    def test_position_queries_expanded_range_translation(
        self,
        start,
        end,
        expected_translations,
        expected_mrna,
    ):
        """The span of the AnnotationCollection object is the original query, but the full range
        of the child objects are larger, and so the translation does not change even when the query is supposed to
        slice down the gene"""
        obj = self.annot_no_range.to_annotation_collection(parent_genome)
        r = obj.query_by_position(
            start, end, coding_only=True, completely_within=False, expand_location_to_children=True
        )
        for gene in r.genes:
            for tx in gene.transcripts:
                assert expected_translations[tx.transcript_symbol] == str(tx.get_protein_sequence())
                assert expected_mrna[tx.transcript_symbol] == str(tx.get_spliced_sequence())

    @pytest.mark.parametrize(
        "start,end,expected_translations,expected_mrna",
        [
            (12, 20, {"tx1": "F", "tx2": "I"}, {"tx1": "GTATTCTT", "tx2": "GTATCTT"}),
            (12, 25, {"tx1": "F", "tx2": "IL"}, {"tx1": "GTATTCTTGGACC", "tx2": "GTATCTTACC"}),
        ],
    )
    def test_position_queries_sliced_range_translation(
        self,
        start,
        end,
        expected_translations,
        expected_mrna,
    ):
        """The span of the AnnotationCollection object is the original query, but the full range
        of the child objects are larger, and so the translation does not change even when the query is supposed to
        slice down the gene"""
        obj = self.annot_no_range.to_annotation_collection(parent_genome)
        r = obj.query_by_position(start, end, coding_only=True, completely_within=False)
        for gene in r.genes:
            for tx in gene.transcripts:
                assert expected_translations[tx.transcript_symbol] == str(tx.get_protein_sequence())
                assert expected_mrna[tx.transcript_symbol] == str(tx.get_spliced_sequence())

    @pytest.mark.parametrize(
        "start,end,coding_only,completely_within,expected_nonempty_identifiers",
        [
            # query removes everything due to completely_within=True
            (21, 22, False, True, set()),
            # entirely intronic query makes non-overlapping isoforms tx2/feat2/feat1 empty
            (21, 22, False, False, {"tx1"}),
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
    def test_position_queries_intronic_queries(
        self,
        start,
        end,
        completely_within,
        coding_only,
        expected_nonempty_identifiers,
    ):
        # test this with different parents to show that the parent doesn't change the outcome
        for parent in [parent_genome, parent_genome_10_49]:
            obj = self.annot_no_range.to_annotation_collection(parent)
            r = obj.query_by_position(start, end, coding_only, completely_within)
            if len(r) == 0:
                assert expected_nonempty_identifiers == set()
            else:
                nonempty_identifiers = set()
                for gene_or_feature in r:
                    for child in gene_or_feature:
                        if not child.chunk_relative_location.is_empty:
                            nonempty_identifiers.update(child.identifiers)
                assert nonempty_identifiers == expected_nonempty_identifiers

    def test_position_queries_intronic(self):
        obj = self.annot_no_range.to_annotation_collection(parent_genome)
        r = obj.query_by_position(21, 22, completely_within=False)
        for gene in r.genes:
            for tx in gene.transcripts:
                if tx.transcript_symbol in ["tx2", "feat2", "feat1"]:
                    assert tx.chunk_relative_location.is_empty
                    assert not tx.chromosome_location.is_empty

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

    def test_hierarchical_children_guids(self):
        obj = self.annot.to_annotation_collection()
        assert obj.hierarchical_children_guids == {
            UUID("639c6178-5f15-935b-5085-dee5ed6badd5"): {UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")},
            UUID("94e30bde-d622-3b98-1745-ab022b6ae6ab"): {
                UUID("1e03f51a-5f3f-601c-1a27-2835c346d2bc"),
                UUID("848cf6c7-6867-c46b-d60f-f5e248febba4"),
            },
            UUID("af85efdb-05fd-fbab-62c4-27e8d0c874e3"): {
                UUID("043d7309-9036-7b27-d841-b7d6a2f70712"),
                UUID("2370657b-19cf-625f-566f-e1486d5dd163"),
            },
        }

    def test_query_by_interval_guids(self):
        obj = self.annot.to_annotation_collection()
        # only one isoform of gene1
        a = obj.query_by_interval_guids(UUID("043d7309-9036-7b27-d841-b7d6a2f70712"))
        assert (
            len(a.genes) == 1
            and a.genes[0].identifiers == {"gene1"}
            and a.genes[0].transcripts[0].guid == UUID("043d7309-9036-7b27-d841-b7d6a2f70712")
        )

        # both isoforms of gene1
        b = obj.query_by_interval_guids(
            [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("2370657b-19cf-625f-566f-e1486d5dd163")]
        )
        assert (
            len(b.genes) == 1
            and b.genes[0].identifiers == {"gene1"}
            and [x.guid for x in b.genes[0].transcripts]
            == [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("2370657b-19cf-625f-566f-e1486d5dd163")]
        )

        # one isoform of gene1 and one isoform of featgrp2
        c = obj.query_by_interval_guids(
            [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")]
        )
        assert (
            len(c.genes) == 1
            and c.genes[0].identifiers == {"gene1"}
            and c.genes[0].transcripts[0].guid == UUID("043d7309-9036-7b27-d841-b7d6a2f70712")
        )
        assert (
            len(c.feature_collections) == 1
            and c.feature_collections[0].identifiers == {"featgrp2"}
            and c.feature_collections[0].feature_intervals[0].guid == UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")
        )

    def test_query_by_transcript_interval_guids(self):
        obj = self.annot.to_annotation_collection()
        # only one isoform of gene1
        a = obj.query_by_transcript_interval_guids(UUID("043d7309-9036-7b27-d841-b7d6a2f70712"))
        assert (
            len(a.genes) == 1
            and a.genes[0].identifiers == {"gene1"}
            and a.genes[0].transcripts[0].guid == UUID("043d7309-9036-7b27-d841-b7d6a2f70712")
        )

        # both isoforms of gene1
        b = obj.query_by_transcript_interval_guids(
            [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("2370657b-19cf-625f-566f-e1486d5dd163")]
        )
        assert (
            len(b.genes) == 1
            and b.genes[0].identifiers == {"gene1"}
            and [x.guid for x in b.genes[0].transcripts]
            == [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("2370657b-19cf-625f-566f-e1486d5dd163")]
        )

        # one isoform of gene1 and one isoform of featgrp2 -- feature gets ignored
        c = obj.query_by_transcript_interval_guids(
            [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")]
        )
        assert (
            len(c.genes) == 1
            and c.genes[0].identifiers == {"gene1"}
            and c.genes[0].transcripts[0].guid == UUID("043d7309-9036-7b27-d841-b7d6a2f70712")
        )
        assert len(c.feature_collections) == 0

    def test_query_by_feature_interval_guids(self):
        obj = self.annot.to_annotation_collection()
        # only one isoform of gene1 -- genes are ignored
        a = obj.query_by_feature_interval_guids(UUID("043d7309-9036-7b27-d841-b7d6a2f70712"))
        assert len(a.genes) == 0

        # one isoform of gene1 and one isoform of featgrp2 (genes are ignored)
        c = obj.query_by_feature_interval_guids(
            [UUID("043d7309-9036-7b27-d841-b7d6a2f70712"), UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")]
        )
        assert len(c.genes) == 0
        assert (
            len(c.feature_collections) == 1
            and c.feature_collections[0].identifiers == {"featgrp2"}
            and c.feature_collections[0].feature_intervals[0].guid == UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")
        )

        # one isoform of featgrp2
        c = obj.query_by_feature_interval_guids([UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")])
        assert len(c.genes) == 0
        assert (
            len(c.feature_collections) == 1
            and c.feature_collections[0].identifiers == {"featgrp2"}
            and c.feature_collections[0].feature_intervals[0].guid == UUID("079c8c04-e2bd-590b-87f7-cb792ba67064")
        )

        # one isoform of featgrp1
        c = obj.query_by_feature_interval_guids([UUID("848cf6c7-6867-c46b-d60f-f5e248febba4")])
        assert len(c.genes) == 0
        assert (
            len(c.feature_collections) == 1
            and c.feature_collections[0].identifiers == {"featgrp1"}
            and c.feature_collections[0].feature_intervals[0].guid == UUID("848cf6c7-6867-c46b-d60f-f5e248febba4")
        )

        # both isoforms of featgrp1
        c = obj.query_by_feature_interval_guids(
            [UUID("848cf6c7-6867-c46b-d60f-f5e248febba4"), UUID("1e03f51a-5f3f-601c-1a27-2835c346d2bc")]
        )
        assert len(c.genes) == 0
        assert (
            len(c.feature_collections) == 1
            and c.feature_collections[0].identifiers == {"featgrp1"}
            and len(c.feature_collections[0].feature_intervals) == 2
        )
        assert [x.guid for x in c.feature_collections[0].feature_intervals] == [
            UUID("848cf6c7-6867-c46b-d60f-f5e248febba4"),
            UUID("1e03f51a-5f3f-601c-1a27-2835c346d2bc"),
        ]

    def test_interval_guids_to_collections(self):
        obj = self.annot.to_annotation_collection()
        m = obj.interval_guids_to_collections
        m = {key: val.to_dict() for key, val in m.items()}
        assert m == {
            UUID("043d7309-9036-7b27-d841-b7d6a2f70712"): {
                "transcripts": [
                    {
                        "exon_starts": [12],
                        "exon_ends": [28],
                        "strand": "PLUS",
                        "cds_starts": [15],
                        "cds_ends": [19],
                        "cds_frames": ["ZERO"],
                        "qualifiers": None,
                        "is_primary_tx": None,
                        "transcript_id": None,
                        "transcript_symbol": "tx1",
                        "transcript_type": None,
                        "sequence_name": None,
                        "sequence_guid": None,
                        "protein_id": None,
                        "product": None,
                        "transcript_guid": None,
                        "transcript_interval_guid": UUID("043d7309-9036-7b27-d841-b7d6a2f70712"),
                    },
                    {
                        "exon_starts": [12, 17, 22],
                        "exon_ends": [16, 20, 25],
                        "strand": "PLUS",
                        "cds_starts": [14, 17, 22],
                        "cds_ends": [16, 20, 23],
                        "cds_frames": ["ZERO", "TWO", "TWO"],
                        "qualifiers": None,
                        "is_primary_tx": None,
                        "transcript_id": None,
                        "transcript_symbol": "tx2",
                        "transcript_type": None,
                        "sequence_name": None,
                        "sequence_guid": None,
                        "protein_id": None,
                        "product": None,
                        "transcript_guid": None,
                        "transcript_interval_guid": UUID("2370657b-19cf-625f-566f-e1486d5dd163"),
                    },
                ],
                "gene_id": "gene1",
                "gene_symbol": None,
                "gene_type": None,
                "locus_tag": None,
                "qualifiers": None,
                "sequence_name": None,
                "sequence_guid": None,
                "gene_guid": UUID("af85efdb-05fd-fbab-62c4-27e8d0c874e3"),
            },
            UUID("2370657b-19cf-625f-566f-e1486d5dd163"): {
                "transcripts": [
                    {
                        "exon_starts": [12],
                        "exon_ends": [28],
                        "strand": "PLUS",
                        "cds_starts": [15],
                        "cds_ends": [19],
                        "cds_frames": ["ZERO"],
                        "qualifiers": None,
                        "is_primary_tx": None,
                        "transcript_id": None,
                        "transcript_symbol": "tx1",
                        "transcript_type": None,
                        "sequence_name": None,
                        "sequence_guid": None,
                        "protein_id": None,
                        "product": None,
                        "transcript_guid": None,
                        "transcript_interval_guid": UUID("043d7309-9036-7b27-d841-b7d6a2f70712"),
                    },
                    {
                        "exon_starts": [12, 17, 22],
                        "exon_ends": [16, 20, 25],
                        "strand": "PLUS",
                        "cds_starts": [14, 17, 22],
                        "cds_ends": [16, 20, 23],
                        "cds_frames": ["ZERO", "TWO", "TWO"],
                        "qualifiers": None,
                        "is_primary_tx": None,
                        "transcript_id": None,
                        "transcript_symbol": "tx2",
                        "transcript_type": None,
                        "sequence_name": None,
                        "sequence_guid": None,
                        "protein_id": None,
                        "product": None,
                        "transcript_guid": None,
                        "transcript_interval_guid": UUID("2370657b-19cf-625f-566f-e1486d5dd163"),
                    },
                ],
                "gene_id": "gene1",
                "gene_symbol": None,
                "gene_type": None,
                "locus_tag": None,
                "qualifiers": None,
                "sequence_name": None,
                "sequence_guid": None,
                "gene_guid": UUID("af85efdb-05fd-fbab-62c4-27e8d0c874e3"),
            },
            UUID("848cf6c7-6867-c46b-d60f-f5e248febba4"): {
                "feature_intervals": [
                    {
                        "interval_starts": [12],
                        "interval_ends": [15],
                        "strand": "PLUS",
                        "qualifiers": None,
                        "feature_id": None,
                        "feature_name": "feat1",
                        "feature_types": ["a", "b"],
                        "sequence_name": None,
                        "sequence_guid": None,
                        "feature_interval_guid": UUID("848cf6c7-6867-c46b-d60f-f5e248febba4"),
                        "feature_guid": None,
                        "is_primary_feature": None,
                    },
                    {
                        "interval_starts": [12, 17, 22],
                        "interval_ends": [16, 20, 25],
                        "strand": "PLUS",
                        "qualifiers": None,
                        "feature_id": None,
                        "feature_name": "feat2",
                        "feature_types": ["b"],
                        "sequence_name": None,
                        "sequence_guid": None,
                        "feature_interval_guid": UUID("1e03f51a-5f3f-601c-1a27-2835c346d2bc"),
                        "feature_guid": None,
                        "is_primary_feature": None,
                    },
                ],
                "feature_collection_name": None,
                "feature_collection_id": "featgrp1",
                "feature_collection_type": None,
                "locus_tag": None,
                "qualifiers": None,
                "sequence_name": None,
                "sequence_guid": None,
                "feature_collection_guid": UUID("94e30bde-d622-3b98-1745-ab022b6ae6ab"),
            },
            UUID("1e03f51a-5f3f-601c-1a27-2835c346d2bc"): {
                "feature_intervals": [
                    {
                        "interval_starts": [12],
                        "interval_ends": [15],
                        "strand": "PLUS",
                        "qualifiers": None,
                        "feature_id": None,
                        "feature_name": "feat1",
                        "feature_types": ["a", "b"],
                        "sequence_name": None,
                        "sequence_guid": None,
                        "feature_interval_guid": UUID("848cf6c7-6867-c46b-d60f-f5e248febba4"),
                        "feature_guid": None,
                        "is_primary_feature": None,
                    },
                    {
                        "interval_starts": [12, 17, 22],
                        "interval_ends": [16, 20, 25],
                        "strand": "PLUS",
                        "qualifiers": None,
                        "feature_id": None,
                        "feature_name": "feat2",
                        "feature_types": ["b"],
                        "sequence_name": None,
                        "sequence_guid": None,
                        "feature_interval_guid": UUID("1e03f51a-5f3f-601c-1a27-2835c346d2bc"),
                        "feature_guid": None,
                        "is_primary_feature": None,
                    },
                ],
                "feature_collection_name": None,
                "feature_collection_id": "featgrp1",
                "feature_collection_type": None,
                "locus_tag": None,
                "qualifiers": None,
                "sequence_name": None,
                "sequence_guid": None,
                "feature_collection_guid": UUID("94e30bde-d622-3b98-1745-ab022b6ae6ab"),
            },
            UUID("079c8c04-e2bd-590b-87f7-cb792ba67064"): {
                "feature_intervals": [
                    {
                        "interval_starts": [35],
                        "interval_ends": [40],
                        "strand": "MINUS",
                        "qualifiers": None,
                        "feature_id": None,
                        "feature_name": "feat3",
                        "feature_types": ["a"],
                        "sequence_name": None,
                        "sequence_guid": None,
                        "feature_interval_guid": UUID("079c8c04-e2bd-590b-87f7-cb792ba67064"),
                        "feature_guid": None,
                        "is_primary_feature": None,
                    }
                ],
                "feature_collection_name": None,
                "feature_collection_id": "featgrp2",
                "feature_collection_type": None,
                "locus_tag": None,
                "qualifiers": None,
                "sequence_name": None,
                "sequence_guid": None,
                "feature_collection_guid": UUID("639c6178-5f15-935b-5085-dee5ed6badd5"),
            },
        }

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
                if hasattr(subitem, "cds"):
                    subitem.cds.sequence_name = "chr1"
        with open(test_data_dir / "collection_gff3_export_chromosome_coordinates.gff") as fh:
            assert fh.read() == "\n".join(str(x) for x in obj.to_gff())

    def test_gtf_export_with_feature(self, test_data_dir):
        obj = self.annot.to_annotation_collection()
        obj.sequence_name = "chr1"
        for item in obj:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
                if hasattr(subitem, "cds"):
                    subitem.cds.sequence_name = "chr1"
        with pytest.raises(NotImplementedError):
            _ = "\n".join(str(x) for x in obj.to_gtf())
        # populate sequence names; normally this is done via the model constructors

    def test_gtf_export_no_transcript_id(self, test_data_dir):
        obj = self.annot_no_features.to_annotation_collection()
        obj.sequence_name = "chr1"
        for item in obj:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
                if hasattr(subitem, "cds"):
                    subitem.cds.sequence_name = "chr1"
        with pytest.raises(GTFExportException):
            _ = "\n".join(str(x) for x in obj.to_gtf())
        # populate sequence names; normally this is done via the model constructors

    def test_gtf_export(self, test_data_dir):
        tx1 = dict(
            exon_starts=[12],
            exon_ends=[28],
            strand=Strand.PLUS.name,
            cds_starts=[15],
            cds_ends=[19],
            cds_frames=[CDSFrame.ZERO.name],
            transcript_symbol="tx1",
            transcript_id="id_1",
        )
        tx2 = dict(
            exon_starts=[12, 17, 22],
            exon_ends=[16, 20, 25],
            strand=Strand.PLUS.name,
            cds_starts=[14, 17, 22],
            cds_ends=[16, 20, 23],
            cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
            transcript_symbol="tx2",
            transcript_id="id_2",
        )

        annot = (
            AnnotationCollectionModel.Schema()
            .load(
                dict(
                    genes=[dict(transcripts=[tx1, tx2], gene_id="gene1")],
                    start=2,
                    end=40,
                )
            )
            .to_annotation_collection()
        )

        annot.sequence_name = "chr1"
        for item in annot:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
                if hasattr(subitem, "cds"):
                    subitem.cds.sequence_name = "chr1"
        with open(test_data_dir / "collection_gtf_export_chromosome_coordinates.gtf") as fh:
            assert fh.read() == "\n".join(str(x) for x in annot.to_gtf())

    def test_gff3_export_chunk_relative(self, test_data_dir):
        obj = self.annot.to_annotation_collection(parent_genome_10_49)
        # populate sequence names; normally this is done via the model constructors
        obj.sequence_name = "chr1"
        for item in obj:
            item.sequence_name = "chr1"
            for subitem in item:
                subitem.sequence_name = "chr1"
                if hasattr(subitem, "cds"):
                    subitem.cds.sequence_name = "chr1"
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
        # not the same because of the non-standard parent
        assert obj2.chromosome_location != obj2.chunk_relative_location
        assert obj2._chunk_relative_bounded_chromosome_location == obj2.chunk_relative_location
        assert obj3.chromosome_location != obj3.chunk_relative_location
        assert obj3._chunk_relative_bounded_chromosome_location == obj3.chunk_relative_location
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

    def test_query_by_identifiers_subset(self):
        model = {
            "feature_collections": [],
            "genes": [
                {
                    "transcripts": [
                        {
                            "exon_starts": [2971596],
                            "exon_ends": [2972637],
                            "strand": "PLUS",
                            "cds_starts": [2971596],
                            "cds_ends": [2972637],
                            "cds_frames": ["ZERO"],
                            "qualifiers": {
                                "gene": ["tas"],
                                "locus_tag": ["b2834"],
                                "gene_synonym": ["ECK2830; JW2802; ygdS"],
                                "function": ["putative enzyme; Not classified"],
                                "codon_start": ["1"],
                                "transl_table": ["11"],
                                "product": ["putative NADP(H)-dependent aldo-keto reductase"],
                                "protein_id": ["NP_417311.1"],
                                "db_xref": [
                                    "ASAP:ABE-0009298",
                                    "EcoGene:EG13093",
                                    "GeneID:947306",
                                    "UniProtKB/Swiss-Prot:P0A9T4",
                                ],
                                "translation": [
                                    "MQYHRIPHSSLEVSTLGLGTMTFGEQNSEADAHAQLDYAVAQGINLIDVAEMYPVPPRPETQGLTETYVGNWLAKHGSREKLIIASKVSGPSRNNDKGIRPDQALDRKNIREALHDSLKRLQTDYLDLYQVHWPQRPTNCFGKLGYSWTDSAPAVSLLDTLDALAEYQRAGKIRYIGVSNETAFGVMRYLHLADKHDLPRIVTIQNPYSLLNRSFEVGLAEVSQYEGVELLAYSCLGFGTLTGKYLNGAKPAGARNTLFSRFTRYSGEQTQKAVAAYVDIARRHGLDPAQMALAFVRRQPFVASTLLGATTMDQLKTNIESLHLELSEDVLAEIEAVHQVYTYPAP"  # noqa: E501
                                ],
                            },
                            "is_primary_tx": False,
                            "transcript_id": None,
                            "protein_id": "NP_417311.1",
                            "product": "putative NADP(H)-dependent aldo-keto reductase",
                            "transcript_symbol": "tas",
                            "transcript_type": "protein_coding",
                            "sequence_name": "NC_000913.3",
                            "sequence_guid": None,
                            "transcript_interval_guid": "cdbed83f-cc17-b945-53df-620890e7e867",
                            "transcript_guid": None,
                        }
                    ],
                    "gene_id": None,
                    "gene_symbol": "tas",
                    "gene_type": "protein_coding",
                    "locus_tag": "b2834",
                    "qualifiers": {
                        "gene": ["tas"],
                        "locus_tag": ["b2834"],
                        "gene_synonym": ["ECK2830; JW2802; ygdS"],
                        "db_xref": ["EcoGene:EG13093", "GeneID:947306"],
                    },
                    "sequence_name": "NC_000913.3",
                    "sequence_guid": None,
                    "gene_guid": "6b729a41-3316-6cc8-ad99-db1597e4c68a",
                },
                {
                    "transcripts": [
                        {
                            "exon_starts": [2972668],
                            "exon_ends": [2973862],
                            "strand": "MINUS",
                            "cds_starts": [2972668],
                            "cds_ends": [2973862],
                            "cds_frames": ["ZERO"],
                            "qualifiers": {
                                "gene": ["lplT"],
                                "locus_tag": ["b2835"],
                                "gene_synonym": ["ECK2831; JW2803; ygeD"],
                                "function": ["orf; Drug/analog sensitivity"],
                                "GO_process": ["GO:0042493 - response to drug"],
                                "note": ["putative resistance proteins"],
                                "codon_start": ["1"],
                                "transl_table": ["11"],
                                "product": ["lysophospholipid transporter"],
                                "protein_id": ["NP_417312.1"],
                                "db_xref": [
                                    "ASAP:ABE-0009300",
                                    "EcoGene:EG12455",
                                    "GeneID:947317",
                                    "UniProtKB/Swiss-Prot:P39196",
                                ],
                                "translation": [
                                    "MSESVHTNTSLWSKGMKAVIVAQFLSAFGDNALLFATLALLKAQFYPEWSQPILQMVFVGAYILFAPFVGQVADSFAKGRVMMFANGLKLLGAASICFGINPFLGYTLVGVGAAAYSPAKYGILGELTTGSKLVKANGLMEASTIAAILLGSVAGGVLADWHVLVALAACALAYGGAVVANIYIPKLAAARPGQSWNLINMTRSFLNACTSLWRNGETRFSLVGTSLFWGAGVTLRFLLVLWVPVALGITDNATPTYLNAMVAIGIVVGAGAAAKLVTLETVSRCMPAGILIGVVVLIFSLQHELLPAYALLMLIGVMGGFFVVPLNALLQERGKKSVGAGNAIAVQNLGENSAMLLMLGIYSLAVMIGIPVVPIGIGFGALFALAITALWIWQRRH"  # noqa: E501
                                ],
                            },
                            "is_primary_tx": False,
                            "transcript_id": None,
                            "protein_id": "NP_417312.1",
                            "product": "lysophospholipid transporter",
                            "transcript_symbol": "lplT",
                            "transcript_type": "protein_coding",
                            "sequence_name": "NC_000913.3",
                            "sequence_guid": None,
                            "transcript_interval_guid": "91e4286f-9757-2bb9-1b50-2a42a9c668f0",
                            "transcript_guid": None,
                        }
                    ],
                    "gene_id": None,
                    "gene_symbol": "lplT",
                    "gene_type": "protein_coding",
                    "locus_tag": "b2835",
                    "qualifiers": {
                        "gene": ["lplT"],
                        "locus_tag": ["b2835"],
                        "gene_synonym": ["ECK2831; JW2803; ygeD"],
                        "db_xref": ["EcoGene:EG12455", "GeneID:947317"],
                    },
                    "sequence_name": "NC_000913.3",
                    "sequence_guid": None,
                    "gene_guid": "52745f9e-9ee8-7a2b-78cf-aa60aec4312b",
                },
                {
                    "transcripts": [
                        {
                            "exon_starts": [2973854],
                            "exon_ends": [2976014],
                            "strand": "MINUS",
                            "cds_starts": [2973854],
                            "cds_ends": [2976014],
                            "cds_frames": ["ZERO"],
                            "qualifiers": {
                                "gene": ["aas"],
                                "locus_tag": ["b2836"],
                                "gene_synonym": ["ECK2832; JW2804"],
                                "EC_number": ["2.3.1.40", "6.2.1.20"],
                                "function": ["enzyme; Fatty acid and phosphatidic acid biosynthesis"],
                                "GO_component": [
                                    "GO:0009274 - peptidoglycan-based cell wall",
                                    "GO:0019866 - organelle inner membrane",
                                ],
                                "GO_process": ["GO:0006464 - protein modification process"],
                                "note": [
                                    "2-acyl-glycerophospho-ethanolamine acyltransferase; acyl-acyl-carrier "
                                    "protein synthetase"
                                ],
                                "codon_start": ["1"],
                                "transl_table": ["11"],
                                "product": [
                                    "fused 2-acylglycerophospho-ethanolamine acyl transferase/acyl-acyl carrier "
                                    "protein synthetase"
                                ],
                                "protein_id": ["NP_417313.1"],
                                "db_xref": [
                                    "ASAP:ABE-0009302",
                                    "EcoGene:EG11679",
                                    "GeneID:947315",
                                    "UniProtKB/Swiss-Prot:P31119",
                                ],
                                "translation": [
                                    "MLFSFFRNLCRVLYRVRVTGDTQALKGERVLITPNHVSFIDGILLGLFLPVRPVFAVYTSISQQWYMRWLKSFIDFVPLDPTQPMAIKHLVRLVEQGRPVVIFPEGRITTTGSLMKIYDGAGFVAAKSGATVIPVRIEGAELTHFSRLKGLVKRRLFPQITLHILPPTQVAMPDAPRARDRRKIAGEMLHQIMMEARMAVRPRETLYESLLSAMYRFGAGKKCVEDVNFTPDSYRKLLTKTLFVGRILEKYSVEGERIGLMLPNAGISAAVIFGAIARRRMPAMMNYTAGVKGLTSAITAAEIKTIFTSRQFLDKGKLWHLPEQLTQVRWVYLEDLKADVTTADKVWIFAHLLMPRLAQVKQQPEEEALILFTSGSEGHPKGVVHSHKSILANVEQIKTIADFTTNDRFMSALPLFHSFGLTVGLFTPLLTGAEVFLYPSPLHYRIVPELVYDRSCTVLFGTSTFLGHYARFANPYDFYRLRYVVAGAEKLQESTKQLWQDKFGLRILEGYGVTECAPVVSINVPMAAKPGTVGRILPGMDARLLSVPGIEEGGRLQLKGPNIMNGYLRVEKPGVLEVPTAENVRGEMERGWYDTGDIVRFDEQGFVQIQGRAKRFAKIAGEMVSLEMVEQLALGVSPDKVHATAIKSDASKGEALVLFTTDNELTRDKLQQYAREHGVPELAVPRDIRYLKQMPLLGSGKPDFVTLKSWVDEAEQHDE"  # noqa: E501
                                ],
                            },
                            "is_primary_tx": False,
                            "transcript_id": None,
                            "protein_id": "NP_417313.1",
                            "product": "fused 2-acylglycerophospho-ethanolamine acyl transferase/acyl-acyl carrier"
                            " protein synthetase",
                            "transcript_symbol": "aas",
                            "transcript_type": "protein_coding",
                            "sequence_name": "NC_000913.3",
                            "sequence_guid": None,
                            "transcript_interval_guid": "e3189127-52ef-21aa-44c2-25e3ebd7834b",
                            "transcript_guid": None,
                        }
                    ],
                    "gene_id": None,
                    "gene_symbol": "aas",
                    "gene_type": "protein_coding",
                    "locus_tag": "b2836",
                    "qualifiers": {
                        "gene": ["aas"],
                        "locus_tag": ["b2836"],
                        "gene_synonym": ["ECK2832; JW2804"],
                        "db_xref": ["EcoGene:EG11679", "GeneID:947315"],
                    },
                    "sequence_name": "NC_000913.3",
                    "sequence_guid": None,
                    "gene_guid": "32710663-a671-5b98-b2e5-4d476811a91f",
                },
            ],
            "name": "NC_000913.3",
            "id": None,
            "sequence_name": "NC_000913.3",
            "sequence_guid": None,
            "sequence_path": None,
            "qualifiers": {
                "organism": ["Escherichia coli str. K-12 substr. MG1655"],
                "mol_type": ["genomic DNA"],
                "strain": ["K-12"],
                "sub_strain": ["MG1655"],
                "db_xref": ["taxon:511145"],
            },
            "start": 2972468,
            "end": 2974062,
            "completely_within": False,
        }
        chunk = seq_chunk_to_parent(
            "CAAAGGCGAGGCACTGGTGCTTTTCACCACAGATAACGAACTGACGCGCGATAAGTTGCAACAGTATGCCCGCGAGCACGGCGTGCCGGAGCTTGCTGTACC"
            "GCGCGATATTCGCTATCTGAAACAGATGCCATTACTTGGCAGCGGCAAACCTGACTTTGTCACGTTGAAAAGCTGGGTAGACGAAGCGGAACAACACGATGA"
            "GTGAGTCAGTGCACACTAACACTTCGTTGTGGTCGAAGGGGATGAAAGCGGTTATCGTGGCGCAGTTTCTCTCTGCGTTTGGCGATAATGCCCTACTGTTTG"
            "CCACTCTGGCGTTACTGAAAGCGCAGTTCTATCCGGAGTGGAGCCAGCCCATCCTGCAAATGGTGTTTGTAGGTGCTTACATTCTTTTTGCGCCGTTTGTCG"
            "GGCAGGTGGCGGATAGCTTCGCCAAAGGCCGGGTGATGATGTTTGCCAACGGCCTGAAGCTGCTGGGCGCAGCCAGTATCTGCTTTGGTATCAATCCGTTTC"
            "TCGGCTATACGCTGGTGGGTGTTGGTGCTGCAGCCTATTCACCGGCGAAATACGGTATTCTCGGCGAATTAACCACGGGTAGTAAGTTAGTGAAAGCTAACG"
            "GTTTAATGGAAGCTTCTACCATAGCGGCGATTTTGCTCGGTTCCGTAGCCGGTGGTGTGCTGGCTGACTGGCATGTCCTCGTCGCCCTGGCCGCATGCGCAC"
            "TGGCCTACGGTGGTGCGGTCGTTGCCAATATCTACATTCCCAAACTGGCGGCGGCGCGTCCGGGGCAGTCCTGGAATCTCATCAACATGACCCGCAGTTTCC"
            "TGAATGCCTGCACCTCGCTATGGCGCAATGGTGAAACGCGTTTTTCGCTGGTGGGCACCAGTTTATTCTGGGGAGCGGGTGTCACGCTGCGTTTCCTGTTGG"
            "TGCTGTGGGTACCGGTGGCGCTGGGCATTACCGATAACGCTACGCCCACCTATCTCAACGCGATGGTAGCGATTGGTATCGTGGTTGGCGCAGGTGCGGCAG"
            "CGAAGTTAGTTACGCTGGAAACCGTGTCACGCTGTATGCCAGCCGGGATTTTGATTGGCGTGGTGGTACTGATTTTTTCCCTGCAACACGAGCTGCTGCCAG"
            "CCTATGCCTTGTTGATGCTGATTGGCGTGATGGGGGGCTTTTTTGTCGTTCCGCTCAATGCGTTGCTACAGGAGCGGGGTAAAAAAAGCGTCGGGGCGGGGA"
            "ATGCGATTGCAGTACAAAACCTTGGCGAAAACAGCGCCATGTTGTTGATGCTGGGCATTTACTCGCTGGCGGTAATGATAGGCATCCCGGTCGTGCCCATTG"
            "GCATTGGCTTCGGTGCGCTGTTTGCGCTGGCAATAACGGCGCTGTGGATCTGGCAGCGCCGTCATTAATATTTAACGCCGGTTTTAACCGGCGTTAATCTTA"
            "TGGTGCCGGATAAGTATAAACCTGATGCACCGCTTCAATTTCAGCTAATACGTCTTCGCTTAACTCCAGATGCAAACTTTCGATGTTAGTTTTCAGCTGATC"
            "CATCGTGGTTGCGCCCAGCAGAGTGCTGGCAACAAACGGTTGACGGCGTACAAACGCGAGCGCC",
            strand=Strand.MINUS,
            sequence_name="'NC_000913.3'",
            start=2972468,
            end=2974062,
            alphabet=Alphabet.NT_EXTENDED,
        )
        a = AnnotationCollectionModel.Schema().load(model)
        aa = a.to_annotation_collection(chunk)
        assert not aa.query_by_feature_identifiers("lplT").genes[0].transcripts[0].has_in_frame_stop

    def test_get_children_by_type(self):
        obj = self.annot.to_annotation_collection()
        assert obj.get_children_by_type("transcript") == obj.genes
        assert obj.get_children_by_type("feature") == obj.feature_collections
        assert obj.get_children_by_type("Transcript") == obj.genes
        assert obj.get_children_by_type("Feature") == obj.feature_collections
        assert obj.get_children_by_type("TRANSCRIPT") == obj.genes
        assert obj.get_children_by_type("FEATURE") == obj.feature_collections
        with pytest.raises(InvalidQueryError):
            _ = obj.get_children_by_type("gene")

    @pytest.mark.parametrize(
        "parent",
        [
            parent_genome,
            parent_genome_10_49,
            parent_genome_with_location,
            parent_genome_rev,
            parent_genome_rev_5_44,
            parent_nonstandard_type_with_sequence,
        ],
    )
    def test_parent_to_dict(self, parent):
        as_dict = self.annot.to_annotation_collection(parent).to_dict(export_parent=True)
        obj = AnnotationCollectionModel.Schema().load(as_dict).to_annotation_collection()
        as_dict_cpy = as_dict.copy()
        obj2 = AnnotationCollection.from_dict(as_dict_cpy)
        assert as_dict_cpy == as_dict  # did not modify the dictionary
        assert obj.get_reference_sequence() == obj2.get_reference_sequence()
        assert obj == obj2

    @pytest.mark.parametrize(
        "parent",
        [
            parent_no_seq,
            parent_nonstandard_type,
        ],
    )
    def test_parent_to_dict_no_sequence(self, parent):
        as_dict = self.annot.to_annotation_collection(parent).to_dict(export_parent=True)
        obj = AnnotationCollectionModel.Schema().load(as_dict).to_annotation_collection()
        as_dict_cpy = as_dict.copy()
        obj2 = AnnotationCollection.from_dict(as_dict_cpy)
        assert as_dict_cpy == as_dict  # did not modify the dictionary
        assert obj == obj2
        with pytest.raises(NullSequenceException):
            _ = obj.get_reference_sequence()
        with pytest.raises(NullSequenceException):
            _ = obj2.get_reference_sequence()

    def test_parent_to_dict_exception(self):
        with pytest.raises(NotImplementedError):
            _ = self.annot.to_annotation_collection(parent_genome).to_dict(
                export_parent=True, chromosome_relative_coordinates=False
            )


class TestNegative:
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

    annot = AnnotationCollectionModel.Schema().load(
        dict(
            genes=[dict(transcripts=[tx1, tx2], gene_id="gene1")],
        )
    )

    tx1_minus = dict(
        exon_starts=[12],
        exon_ends=[28],
        strand=Strand.MINUS.name,
        cds_starts=[15],
        cds_ends=[19],
        cds_frames=[CDSFrame.ZERO.name],
        transcript_symbol="tx1",
    )
    tx2_minus = dict(
        exon_starts=[12, 17, 22],
        exon_ends=[16, 20, 25],
        strand=Strand.MINUS.name,
        cds_starts=[14, 17, 22],
        cds_ends=[16, 20, 23],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        transcript_symbol="tx2",
    )

    annot_minus = AnnotationCollectionModel.Schema().load(
        dict(
            genes=[dict(transcripts=[tx1_minus, tx2_minus], gene_id="gene1")],
        )
    )

    def test_translate(self):
        a = self.annot.to_annotation_collection(parent_genome)
        a_neg = self.annot.to_annotation_collection(parent_genome_rev)
        assert str(a.genes[0].get_reference_sequence()) == str(a_neg.genes[0].get_reference_sequence())
        assert str(a.genes[0].get_primary_transcript().get_spliced_sequence()) == str(
            a_neg.genes[0].get_primary_transcript().get_spliced_sequence()
        )

    def test_translate_chunk(self):
        a = self.annot.to_annotation_collection(parent_genome_10_49)
        a_neg = self.annot.to_annotation_collection(parent_genome_rev_5_44)
        assert str(a.genes[0].get_reference_sequence()) == str(a_neg.genes[0].get_reference_sequence())
        assert str(a.genes[0].get_primary_transcript().get_spliced_sequence()) == str(
            a_neg.genes[0].get_primary_transcript().get_spliced_sequence()
        )

    def test_translate_neg(self):
        a = self.annot_minus.to_annotation_collection(parent_genome)
        a_neg = self.annot_minus.to_annotation_collection(parent_genome_rev)
        assert str(a.genes[0].get_reference_sequence()) == str(a_neg.genes[0].get_reference_sequence())
        assert str(a.genes[0].get_primary_transcript().get_spliced_sequence()) == str(
            a_neg.genes[0].get_primary_transcript().get_spliced_sequence()
        )

    def test_translate_neg_chunk(self):
        a = self.annot_minus.to_annotation_collection(parent_genome_10_49)
        a_neg = self.annot_minus.to_annotation_collection(parent_genome_rev_5_44)
        assert str(a.genes[0].get_reference_sequence()) == str(a_neg.genes[0].get_reference_sequence())
        assert str(a.genes[0].get_primary_transcript().get_spliced_sequence()) == str(
            a_neg.genes[0].get_primary_transcript().get_spliced_sequence()
        )

    def test_subquery_neg_chunk(self):
        """Querying by interval still produces valid transcripts"""
        a_neg = self.annot_minus.to_annotation_collection(parent_genome_rev_5_44)
        a_neg_subquery = a_neg.query_by_position(12, 28)
        assert str(a_neg_subquery.genes[0].get_reference_sequence()) == str(a_neg.genes[0].get_reference_sequence())
        assert str(a_neg.genes[0].get_primary_transcript().get_spliced_sequence()) == str(
            a_neg_subquery.genes[0].get_primary_transcript().get_spliced_sequence()
        )

    def test_nested_parents(self):
        """Querying by interval should be robust to arbitrarily nested parents."""
        nested_parent = Parent(
            id="testseq:10-49",
            sequence_type=SequenceType.SEQUENCE_CHUNK,
            strand=Strand.PLUS,
            location=SingleInterval(10, 49, Strand.PLUS),
            sequence=Sequence(
                genome,
                type=SequenceType.SEQUENCE_CHUNK,
                id="testseq:10-49",
                alphabet=Alphabet.NT_EXTENDED,
                parent=Parent(
                    id="testseq",
                    sequence_type=SequenceType.CHROMOSOME,
                    location=SingleInterval(
                        0,
                        54,
                        Strand.PLUS,
                        parent=Parent(
                            id="testseq",
                            sequence_type=SequenceType.CHROMOSOME,
                            location=SingleInterval(
                                0,
                                54,
                                Strand.PLUS,
                                parent=Parent(
                                    id="testseq",
                                    sequence_type=SequenceType.CHROMOSOME,
                                    location=SingleInterval(0, 54, Strand.PLUS),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        )

        obj = self.annot.to_annotation_collection(nested_parent)
        assert obj.chromosome_location.parent.id == "testseq"

        ref_obj = self.annot.to_annotation_collection(parent_genome)
        assert str(obj.genes[0].get_reference_sequence()) == str(ref_obj.genes[0].get_reference_sequence())
        assert str(obj.genes[0].get_primary_transcript().get_spliced_sequence()) == str(
            ref_obj.genes[0].get_primary_transcript().get_spliced_sequence()
        )

        # the below test used to fail before `AnnotationCollection._subset_parent()` understood to look
        # for the sequence ID on the first chromosome ancestor type.
        subquery = obj.query_by_position(10, 28, completely_within=False)
        assert subquery.chromosome_location.parent.id == "testseq"

    def test_equality_different_parents(self):
        obj1 = self.annot.to_annotation_collection(parent_genome_10_49)
        obj2 = self.annot.to_annotation_collection(parent_genome_rev_5_44)
        assert obj1 != obj2
        assert hash(obj1) != hash(obj2)

    @pytest.mark.parametrize("parent", [parent_genome, parent_genome_10_49, parent_genome_rev_5_44])
    def test_pickling(self, parent):
        obj = self.annot.to_annotation_collection(parent)
        obj_str = pickle.dumps(obj)
        new_obj = pickle.loads(obj_str)
        assert obj.to_dict() == new_obj.to_dict()
        assert obj.get_reference_sequence() == new_obj.get_reference_sequence()

    def test_query_by_identifier_gene_exceeds_chunk_bounds(self):
        """
        Prior to this test, if an AnnotationCollection was chunk relative, with genes that exceed the chunk boundaries,
        then querying by identifier would raise an exception because the collection was trying to exist outside
        of its boundaries.
        """
        d = {
            "genes": [
                {
                    "transcripts": [
                        {
                            "exon_starts": [
                                111134,
                                112069,
                                112738,
                                116223,
                                116461,
                                116893,
                                117127,
                                117803,
                                118782,
                                119023,
                                119719,
                                120486,
                                121143,
                            ],
                            "exon_ends": [
                                111252,
                                112239,
                                112855,
                                116287,
                                116623,
                                117018,
                                117251,
                                118027,
                                118908,
                                119240,
                                120001,
                                120718,
                                123752,
                            ],
                            "strand": "PLUS",
                            "cds_starts": [
                                112083,
                                112738,
                                116223,
                                116461,
                                116893,
                                117127,
                                117803,
                                118782,
                                119023,
                                119719,
                                120486,
                                121143,
                            ],
                            "cds_ends": [
                                112239,
                                112855,
                                116287,
                                116623,
                                117018,
                                117251,
                                118027,
                                118908,
                                119240,
                                120001,
                                120718,
                                121519,
                            ],
                            "cds_frames": [
                                "ZERO",
                                "ZERO",
                                "ZERO",
                                "ONE",
                                "ONE",
                                "ZERO",
                                "ONE",
                                "ZERO",
                                "ZERO",
                                "ONE",
                                "ONE",
                                "TWO",
                            ],
                            "qualifiers": None,
                            "is_primary_tx": False,
                            "transcript_id": None,
                            "transcript_symbol": "rna-XM_017028208.1",
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr21:6000000-7000000",
                            "sequence_guid": None,
                            "protein_id": None,
                            "product": "salt inducible kinase 1B (putative)",
                            "transcript_guid": None,
                            "transcript_interval_guid": UUID("f1027f34-b941-52de-5b7e-de1257e9235a"),
                        },
                        {
                            "exon_starts": [
                                111130,
                                112069,
                                112738,
                                116223,
                                116461,
                                116893,
                                117127,
                                117803,
                                118255,
                                118782,
                                119023,
                                119719,
                                120486,
                                121143,
                            ],
                            "exon_ends": [
                                111252,
                                112239,
                                112855,
                                116287,
                                116623,
                                117018,
                                117251,
                                118027,
                                118402,
                                118908,
                                119240,
                                120001,
                                120718,
                                123778,
                            ],
                            "strand": "PLUS",
                            "cds_starts": [
                                112083,
                                112738,
                                116223,
                                116461,
                                116893,
                                117127,
                                117803,
                                118255,
                                118782,
                                119023,
                                119719,
                                120486,
                                121143,
                            ],
                            "cds_ends": [
                                112239,
                                112855,
                                116287,
                                116623,
                                117018,
                                117251,
                                118027,
                                118402,
                                118908,
                                119240,
                                120001,
                                120718,
                                121519,
                            ],
                            "cds_frames": [
                                "ZERO",
                                "ZERO",
                                "ZERO",
                                "ONE",
                                "ONE",
                                "ZERO",
                                "ONE",
                                "ZERO",
                                "ZERO",
                                "ZERO",
                                "ONE",
                                "ONE",
                                "TWO",
                            ],
                            "qualifiers": None,
                            "is_primary_tx": False,
                            "transcript_id": None,
                            "transcript_symbol": "rna-NM_001320643.3",
                            "transcript_type": "protein_coding",
                            "sequence_name": "chr21:6000000-7000000",
                            "sequence_guid": None,
                            "protein_id": None,
                            "product": "salt inducible kinase 1B (putative)",
                            "transcript_guid": None,
                            "transcript_interval_guid": UUID("50a31cc8-1d0e-9bf8-02ca-75c023007578"),
                        },
                    ],
                    "gene_id": "1e783983-3f9f-dee9-c5d7-86a903c32045",
                    "gene_symbol": "SIK1B",
                    "gene_type": "protein_coding",
                    "locus_tag": None,
                    "qualifiers": None,
                    "sequence_name": "chr21:6000000-7000000",
                    "sequence_guid": None,
                    "gene_guid": UUID("ebc74c11-fcef-3753-e3fc-a6923929ae08"),
                }
            ],
            "feature_collections": [],
            "variant_collections": [],
            "name": "chr21:6000000-7000000",
            "id": None,
            "qualifiers": None,
            "sequence_name": "chr21:6000000-7000000",
            "sequence_guid": None,
            "sequence_path": None,
            "start": 111933,
            "end": 121669,
            "completely_within": False,
            "parent_or_seq_chunk_parent": {
                "seq": "GTGGATAGAGTGGGGGCGACAGGCGTGGGAGCAGGTCCCTGCGTCCCGTGGGGACTGGGGGCTCCGCGGGCACGGATGGAGCCGACCGCGGGCGGCGGGGGCGCTGGTGGGCTCTGAGCTCTGTGCGGCCCCGCAGGTGCGCGCGGAGCCATGGTTATCATGTCGGAGTTCAGCGCGGACCCCGCGGGCCAGGGTCAGGGCCAGCAGAAGCCCCTCCGGGTGGGTTTTTACGACATCGAGCGGACCCTGGGCAAAGGCAACTTCGCGGTGGTGAAGCTGGCGCGGCATCGAGTCACCAAAACGCAGGTGCGTGGGGGCGGTGGGACCCAGCCGGCGGGGCCTCCCCACTGGACCCGGGGCGACCCGACCTTTGGCGGGTGGACCCCATGCAGATTCACTGCTCCAGGCTTGATTGTCGTGGTGGGTAAATAGTAACCGTTTTTAGTTCGGTAAAGAAAAATAATGCTTTATATTATGTTGCTGCATAATTTGGACATATATTGGGAGAAAGCAGAGAGGATAATAGCAAAAAAAAGGACCAACCTAAGATTGCAGTGGTTCATGGACTCAGAGCTAAACCCTGTAAAGTGAGCCTGCAAATACTTAAGTCACTTACTTAACTCTGATTATTTTAAATAACATTAGTGGTTACTTTGGTTATGTTTTCCCAACTTCTTGTGACCTTCTGGAGACAGAGTGTTGAGAAATTAACTTGCAAAAAATGAAATGTGTAAATTAGGTTGAGGGTTTCTTTTCTTTTTTTAAAAAAACCCACAGAATATCCTTTTTCTTTTTAAATTTGTAGGTTGCAATAAAAATAATTGATAAAACACGATTAGATTCAAGCAATTTGGAGAAAATCTATCGTGAGGTTCAGCTGATGAAGCTTCTGAACCATCCACACATCATAAAGCTTTACCAGGTAAGGGGTCAGCTTGCCTTTCTCTGCTAATCCTGGCAAATGTGTTATATTTTATTTCCAACAGCATAAGTCTGCAATTTCTGAAGGTGGCTTCCCTTTTGGGGAGTGTGTGGAAGTTGGGCATTTGGTTGTATAAGGTAAATAGATTGATTTTTATAGGCTGTGCAGATATATTTACAATTATTTTATTGCTTTTGGATCATAATTGAGATTTGCTTCATTAGAATTTTATTTAGTATCCTAACCTGATCCAGTTTTAAATGTGTAGAAGTCTGTTGTAAATTTGCTTTACTATGAAAACAGACAAACTTATCTGTAGGTATTGCCACAGGTCACGCATTTAGTATAAAAGGACAAATTTGAGTGAAACTGGAAACTTCTGCCTTTGCAGATAATTATTAAACTTTAAAGCCAGCTGCTCTTACTGGTTTTCATTAGTTTATTTTTGCAAGGACAGTAGTGAGCATTTGTTTTCCCTATTTAAACTCTTCCATGAAAAACGGAGAGATGAAAATACTTTTCCATAGTAAGAGCAAGGAGCTTGTCATTTAATCCATGCTAGGAACATGGCTTGCCCTGTCTCTCTGTTGGTTTATTTTGTGACTTAATTTTTGTTCTCATGTTGAAATTAATAAATTATATCATCAGGAAATAAGTACTAATTTAATTAGAACAGCAAAGCAATCATTTACACCTCCCTAATGACAAATGACCACTTTTGCAAGGTGAAATCAAGTCTTCTTGCAGAAGAAGTTCTTTGCAGGTGCTATTTCATCCATTTGTGTTGGGTTGAAGTTCCTCGTAATGATTCTTATGAAATGGAATGTGAGATGACTCTGGGATGTTGTTGTGTGTGTTTGTGCAATTAGGACAGCTCGGTGTACAGTGTGGTAGGAGAGAGAATTATGATTCGAAATTATCCCCCAAAAGCCTCTGCACTAGATTGAAAATAAAAGTAGCAGTATTGGTGGAGCTTACGACAGTTTGCTCAGTGTCTTTAAGCACAGAATTCATGCCGTGTGCGCCTGGTTTGTGTTTTGTTTCTGTTGTCATTTTTAAGTTGTGTGCTCCATGGCCACTGGTGTACTGCCTGTGCAGCGTGCCCAGCCCATTTCCTCAGGAGTGGGCTGAGCTGGTTTTGTTTTTGACTTTGCTCAAGGGCTGAACCCTAATGGTAGTACCTTCCTTTTCTTCCTGTTCTTCTCCTTTGCCCATTGTGTGGAATTAGCTCCGTTAAAGAAACAAATTGCCCATGTGTTTGCTGCTGAGGTCTTTTCTGGAGGTGCATTTTCTTACTGCTGCTGTATGTGATTCTTAACTGTTGGGAGATGATAGAATTCAGTGGTGTAGCAATGACTTTTACTCCTCTGGTTTTGATGCTAATGTCTAATTTTCTTTCTTTTCTTTCTTTCTTCCTGTCTTTTTTTTTTTTTTCTTTTTTGAGACGGAGTCTCTGTCCGTCACCCAGGCTGGAGTGCAGTGGCGCGATTTCGGCTCACTGCAACTGCCGCCTCCTGGGTTCAAGCAATTCTTCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACAACACCCAGCTAATTTTTTGTATTTTTAGTAGAGACAGGGTTTCACCGTGTTAGCCAGGATAGTCTGGATCTCCTGACCTTGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGAACCACCGCGCCCGGCCACTAGTGTCTAATTTTCATTCCGCTGGCCGGGGCCCCAGAACCCTGCTGAGACACTTCCCTCTTTTTTTCAAGTTTTCAAAAAGTTCAAACACTTTGACTTTTAATTCTTTTCCTCCTCCTCCCTTTCTTCAGTTGGGGACGGGTCAGCAGTGGTCTTAGAGTTATTGTACCACCCACAGGGGCGGGGCAGTTAAAATTCTTGCCAAACTGCTAGAATTTATAGAACACTTTTCTTTTTAGATTTATAACTCAAAGCAGCAATCTTTGAAGCATCTTTTTTTTTTTTTTTTTTTTTGAGATGGAGTCTCACTGTCGCCCAGGCTGGAGTGCAGTGGCGCTATGGATCGTGACTCACTGCAACCTCCACCTCCTGGGTTCAAGCGATTCTCCTACCTCAGCCTCCCAAGTAGCTGGGACTACAGGCACGCGCCACCACACCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGACCAGGCTGGTCTCGAACTCCTAACCTCAAGTGATCAGCCTGGGTGATCCCAAAGTACTGGGATTACAAGCATGAGCCACCACCTCCAGCACAAAGCATCTTTTTAAAGGAATAACTATATTTTGAACCCAAAAAAATCTTTTTTCATAAGAATATATTGCCTCATCCTCTCCAAAATTTTTCAATACTCACAGACATAAGTTATATTTTTAAATCTCTTTTTAAGTAAAATAAAGTTATTGATTAGTGAATGAACTTCAAAGGGCTATGTTGGCATGGCTGTTGTGTGGCTACATGCTGAGTGAATGGCTACATGCTGAGTGAATGAGGCTGGATTCTCACAGCCTCATTCATCCTGTCAATAAAGCTGGGTTCCCACAGCGTCATTCATCCTGTCAATCAAGCGAGTGGGAATGGGGTGTGCATTCTCTATTACCCAGGTGGCTCTGACCGTTGGCTGGAATCATCTGTTGCCTACTAGGCTTCTATAAAGCCCAGCAGACTTCCAGCAATCTCCTCCTCTGTAGCCTTCAGAATACTGTCTCTCCTTTTTTGAGATGGAGTCTCGCTCTGTCGCCCAAGTTGGAGTGCAGTGGTACGATCTCGGCTCACTGCAACCTCCGCCTTGCGGGTTCGAATGATCCTCTTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGTGCCCACCACCTCACCTAGCTAGTTTTTGTATTTTTTGTAGAGATGGGGTCTTGCCGTGTTGACCAGGCTGGTCTCGAACTCAAGCGATTCTCCTGCCTTGGCCTCCCGAAGTGCTGGGATTACAGGGCGTGAGCCACCATGCCCCGGCCTAAGCATTCTCAAATTGTATGTGAGGTGTCTTGAATTAAATATGACCACAGAGGGGAATTCTAATATATTTAGACCTTTTGTGAGAGCAAAGGTGTATCAGGATAGGCAAGTCAGTGAGGAATTCACCTCTGATTTTTCAAAGCTGTCTTTGTGAGTAGGAATCTCCTATTTTTTTTTTATTTGAGAAATTTTTCTTTTTCTTTGCATTTTTAAATGATGAATGCCTTTGTCCTAGAGAAGTTAGAAATTTCCTGGGGTTCTGGTGGATGGTTTGTAGCGCTGCCTTGTGCCTTTTTTGTGCGGTTAAATGTCAGGTCTTTCTTCTTTAGGTTATGGAAACAAAGGACATGCTTTACATCGTCACTGAATTTGCTAAAAATGGAGAAATGTTTGGTAAGCCACACTCGTTTTCAGTATCTGTTGAAAACCAATGGCACGAACGCGTCTGGATTTGCGCGGCCAGTGTCTGGGCCACACATCCCAGGACCCCGCGGTGTTCCCCAGGGCTCCATGGGTGGTGGGCAGGCTTTGCCCTGGGGGCTGTAGCTTTGTTTTGTGGCTCCTCAGATTATTTGACTTCCAACGGGCACCTGAGTGAGAACGAGGCGCGGAAGAAGTTCTGGCAAATCCTGTCGGCCGTGGAGTACTGTCACGACCATCACATCGTCCACCGGGACCTCAAGACCGAGAACCTCCTGCTGGATGGCAACATGGACATCAAGCTGGCAGGCACGGAGGGTCCGGGGTGGGAGCAGGGACATCAAGCTGGCAGGCACGGAGGGTCTGGGGTGGGAGCAGGGACATCAAGCTGGCAGGCACGGAGGGTCCGGGGTGGGAGCAGGGACATCAAGCTGGCAGGCACGGAGGCCCCGGGGTGGGAGTGTGCCCGCAAGCAGCCCCAGCTTCCCGGCATGTGCCGAAACCACAGCCCCTGTGCCGAAACCGCAGCCCCTTGCCATTTGCCGTTTAACCTCATCCATCTGTTTCTTTGCCGCTCAGATTTTGGATTTGGGAATTTCTACAAGTCAGGAGAGCCTCTGTCCACGTGGTGTGGGAGCCCCCCGTATGCCGCCCCGGAAGTCTTTGAGGGGAAGGAGTATGAAGGCCCCCAGCTGGACATCTGGGTAGGAGCCCTGTGCGCGCAGACCCCCTTCCCGAGGCCGCGTTCCCCGAGGGCGCCGTGTTCCCGGGGGCACCGCCCTGGCGCTGATGTGGCTCTGTGGTCCTCAACAGAGCCTGGGCGTGGTGCTGTACGTCCTGGTCTGCGGTTCTCTCCCCTTCGATGGGCCTAACCTGCCGACGCTGAGACAGCGGGTGCTGGAGGGCCGCTTCCGCATCCCCTTCTTCATGTCTCAAGGTGAGTCGTGTGGTCTCGCCTGGGCAGGGGCCTGTGTCTCCTGCAGCCCCTTCGTGGCTGGCTTCTGAGTCCTTGTTGAGTCTGGGAATCAGGCCCAGGTTACTGTGGGCTGCGGGAGAAACCCAGCCCATGGATGCCTTCTCAGTGCCCTATTTCGAGGGTGGACGGGTGCCCTTGCTGCACAGGCTTGTGCCTTGCGCCTGCCGTCCCCGCCCAGTACCCGCACCCAGGGTGCTGGCATCATCCCTGAGTCCAACCCACTGGCTGTGGCTCACCCGGCTTTTCTGTGTCATAAAGATGCCTGTCCACCTGGTGCTCAGGTGTGGGGCACGGCAGAAGGGTGAACCTAAGTGGTTCTCTGTGTGGTGGACACAAAGTGTGGCAGCCCCCAGGGTCTGGATCCCTGGCTGATGGCTCCCCTTGGGTGGCAGGACTGAGCCGATGGCCCGGCCCCTGCTCCTGGGCCCTGGCGTGTCAAACCACGTGGGGCGGTGGGTGGGAGCGCAGCCGAGGTCCCGGCCGCCCTGGCTCAGCCTCCTGGCCCCTCCACAGACTGTGAGAGCCTGATCCGCCGCATGCTGGTGGTGGACCCCGCCAGGCGCATCACCATCGCCCAGATCCGGCAGCACCGGTGGATGCGGGCTGAGCCCTGCTTGCCGGGACCCGCCTGCCCCGCCTTCTCCGCACACAGCTACACCTCCAACCTGGGCGACTACGATGAGCAGGCGCTGGGTATCATGCAGACCCTGGGCGTGGACCGGCAGAGGACGGTGGAGGTGAGCCTGGCCACACTTGCCCTGGCCTCACCCGCAGGGGCTAGGCAGCTGTCTCAGGGGAAGCAGGCACTCACCAGCTGAGTTAAAACGGGAGCACAGGCTAATTTAGGGGCCGGCTTGTCCACCCCACTAAAGGATATTCTCAGTCACCCAAGAATAACCTGGGATCGGGTGGGCCCTGGGGCTCCCCGCCATCCCCATGCTGATCTCTGCTCCTCTTGTTCCCAGTCACTGCAAAACAGCAGCTATAACCACTTTGCTGCCATTTATTACCTCCTCCTTGAGCGGCTCAAGGAGTATCGGAATGCCCAGTGCGCCCGCCCCGGGCCTGCCAGGCAGCCGCGGCCTCGGAGCTCGGACCTCAGTGGTTTGGAGGTGAGGGGGAGGAGTCTCCTCCCAGGCCCCAGGCTCCCTCCCCTGTCAGGCACCGGCTTGGAGGGCGGTTCCTTGCGTGGGCAGCGGGTCCCAGGCCTCGTGGGGAAGGGGGTGCCAGCTGCTGGGGGCTGGACTCTGCCCAGAGGCCACTTGTCCCTGACTCATCCCTGGGGCCGGCCTGTCACGCCACCTTCTGGCAGCGCAGCACCAGCACCTGGTCCTGAGCCCAGCCTGGCCAGCCAGTGTCCCTCCGTCAGCTGTCATGCCCACCAGCGAGGCAGCTTGCACTCCAGACAGAAGCCAGCTTGTTTCTTCTCTTGATGGCGCTGGTGGTCCGGAGTCGCTCCCCTGACAGCGTCTTTCCCGTCTGCCGGCCCCAGGTGCCTCAGGAAGGTCTTTCCACCGACCCTTTCCGACCTGCCTTGCTGTGCCCGCAGCCGCAGACCTTGGTGCAGTCCGTCCTCCAGGCCGAGATGGACTGTGAGCTCCAGAGCTCGCTGCAGTGGGTGAGTGCCCACAGCGGGTGTGCAGAGGGCTCGCCTCAGCCCAGCCCTGGTGCCCCCGGTGTGCCCGGGTCACCTGGAGTCCGAGAAGTCACTGGCTTGTGTCTCTCCAACGCAGCCCTTGTTCTTCCCGGTGGATGCCAGCTGCAGCGGAGTGTTCCGGCCCCGGCCCGTGTCCCCAAGCAGCCTGCTGGACACAGCCATCAGTGAGGAGGCCAGGCAGGGGCCGGGCCTAGAGGAGGAGCAGGACACGCAGGAGTCCCTGCCCAGCAGCACGGGCCGGAGGCACACCCTGGCCGAGGTCTCCACCCGCCTCTCCCCACTCACCGCGCCATGTAAGTGTCCCCGGGGGCCCAGGAGGACACCGGTGGATAGGCTTACGGTCGACGTGAGGGTGGGCTAATTTAGAATGGACGTTTTGCCCGGCAGCCTCTCAGGTTGGACTTCTCAGGATTTGCCATTTGTTTTAATCCCTGAGACCACACAGTTGATGTTTAGAGCCTGCCCTGCATGTGGTCGTTCCAGTGGAGGATACAGCATGGGGTCTGGCCTCCAGCAGGGTCCTCCCCAGGCCGCCCCTGGGTGCCGGGAGGGCAGCCCCTTGGCCTGAGGCCCACTATGACCTGCCCCCTGCAGCTGCACCGTGATGGTGGCTTGCCTTTGTGGCTCCCTGGGCTCTGGTGGCCTCGAGCCCTCTTCCCACCAGGTTGATGGTGGGGATGGGGAGGCCAGCGCAGCCATGTGTGCCCAGCAGTGGCCGGGGGAGCCTATTCCTTTGCACTGCAGCATCAAAAGCGCTGTCCTCCCCCCACAGGTATAGTCGTCTCCCCCTCCACCACGGCAAGTCCTGCAGAGGGAACCAGCTCTGACAGTTGTCTGACCTTCTCTGCGAGCAAAAGCCCCGCGGGGCTCAGTGGCACCCCGGCCACTCAGGGGCTGCTGGGCGCCTGCTCCCCGGTCAGGCTGGCCTCGCCCTTCCTGGGGTCGCAGTCCGCCACCCCAGTGCTGCAGGCTCAGGGGGGCTTGGGAGGAGCTGTTCTGCTCCCTGTCAGCTTCCAGGAGGGACGGCGGGCGTCGGACACCTCACTGACTCAAGGTGAGCCACGCTCCTCCCACACTTACCTCCACCTTCCCCAGGGGACCATGTGTGTCTCTGGCAGTACGTGACTTTGTCCGTGATGGCAGATGGCACCCCCTGTTCTCCACCGGGCCTGGGTGGGACCCTCAGTGCTCTGGGCAGGCTGGGGTGCTCAGTGCTCTGGTGGCTCGGGGCGTCACGGCCTGCTGGGATAGACACACATGGGTCCCTGAGCACGCGGCCTCCATGGCTGGTTCTGAAAGCACAGGAGACGACTCTGTGCTGGGCAGCACCTCCCGACTTGGAGGAGGGAGGGCCCCTGGCTGCCAGCTGCCTCCACGCCATCCTGGGGCTTAGGTGCCAGACTCCTGTCCAGACGTGCTTGTGTCACCCGTCCTTTCCTTACCCCCAACCTGAGGCTTGGAAACCCCTTAAGCCAAGGGCCGTGGATGCTGGGCTGAGAGCCGGGTGGCCGTTGACCTCCTGATTCATTCTCCCTGCAGGGCTGAAGGCCTTTCGGCAGCAGCTGAGGAAGACCACGCGGACCAAAGGGTTTCTGGGACTGAACAAAATCAAGGGGCTGGCTCGCCAGGTGTGCCAGGTCCCTGCCAGCCGGGCCAGCAGGGGCGGCCTGAGCCCCTTCCACGCCCCTGCACAGAGCCCAGGCCTGCACGGCGGCGCAGCCGGCAGCCGGGAGGGCTGGAGCCTGCTGGAGGAGGTGCTAGAGCAGCAGAGGTAGGGCCTGCCCCCGCCCTGGGACCCCGGGTGGGCACACGGCAGGTTATCTCCTCGAGGAACCTCATCTGCTAAGTGGTTCCCTCCTCTCTGTAGCCCAGTGCACACCCCCGCTCCCAGCCAGGGAGATGTGTGGGGCGTAGGTCCTAGGTGCTGAGCCATGGGGGTGCAGCAGGCGGGCGTGTCCTTTAAAGTCCCTGGGTGGGTGAGGGTGGCGGGGAGCGAGGGCGCCTTGTGGCCGCATCTCTGAGCTGCTGAGAAACCGGGTGGAGAATGAAAGGTGGGGCGCGGTCAGGGATCAGCCACGCACCTGCCCTCGGCAGCCGCGGCTGGCAGCTCCACGGGCGGGCCCTGCCACACGGGCACTCGGAAACCCGAGAACCCTGCGAGCCGGCGCAGTGACCACCTGTCCTCTGTTCCCACAGGCTGCTCCAGTTACAGCACCACCCGGCCGCTGCACCCGGCTGCTCCCAGGCCCCCCAGCCGGCCCCTGCCCCGTTTGTGATCGCCCCCTGTGATGGCCCTGGGGCTGCCCCGCTCCCCAGCACCCTCCTCACGTCGGGGCTCCCGCTGCTGCCGCCCCCACTCCTGCAGACCGGCGCGTCCCCGGTGGCCTCAGCGGCGCAGCTCCTGGACACACACCTGCACATTGGCACCGGCCCCACCGCCCTCCCCGCTGTGCCCCCACCACGCCTGGCCAGGCTGGCCCCAGGTTGTGAGCCCCTGGGGCTGCTGCAGGGGGACTGTGAGATGGAGGACCTGATGCCCTGCTCCCTAGGCACGTTTGTCCTGGTGCAGTGAGGGCAGCCCTGCATCCTGGCACGGACACTGACTCTTACAGCAATAACTTCAGAGGAGGTGAAGACATCTGGCCTCAAAGCCAAGAACTTTCTAGAAGCGAAATAAGCAATACGTTAGGTGTTTTGGCTTTTTAGTTTATTTTTGTTTTAT",  # noqa: E501
                "sequence_name": "chr21:6000000-7000000",
                "start": 111933,
                "end": 121669,
                "strand": "PLUS",
                "alphabet": "NT_EXTENDED_GAPPED",
                "type": "SEQUENCE_CHUNK",
            },
        }
        ac = AnnotationCollection.from_dict(d)
        new_ac = ac.query_by_feature_identifiers("SIK1B")
        assert new_ac.start == 111130
        assert new_ac.end == 123778
