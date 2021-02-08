"""
Tests for writing and reading feature intervals.
"""
import pytest
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame
from inscripta.biocantor.io.gff3.exc import GFF3LocusTagError, GFF3ChildParentMismatchError
from inscripta.biocantor.io.gff3.parser import parse_standard_gff3
from inscripta.biocantor.io.gff3.writer import collection_to_gff3
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence

genome = "AAGTATTCTTGGACCTAATTAAAAAAAAAAAAAAAAAAA"
parent_genome = Parent(sequence=Sequence(genome, Alphabet.NT_STRICT))


class TestFeatures:
    feat1 = dict(
        interval_starts=[2],
        interval_ends=[5],
        strand=Strand.PLUS.name,
        feature_types=["a", "b"],
        feature_name="feat1",
        sequence_name="chr1",
    )
    feat2 = dict(
        interval_starts=[2, 7, 12],
        interval_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        feature_types=["b"],
        feature_name="feat2",
        sequence_name="chr1",
    )
    featcollection = dict(
        feature_intervals=[feat1, feat2],
        feature_collection_name="featgrp1",
        feature_collection_id="abc123",
        feature_collection_type="group",
        sequence_name="chr1",
    )

    tx1 = dict(
        exon_starts=[2],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[5],
        cds_ends=[9],
        cds_frames=[CDSFrame.ZERO.name],
        sequence_name="chr1",
    )
    tx2 = dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        sequence_name="chr1",
    )
    gene = dict(
        transcripts=[tx1, tx2],
        qualifiers={},
        gene_type=Biotype.protein_coding.name,
        gene_id="gene1",
        sequence_name="chr1",
    )

    annot = AnnotationCollectionModel.Schema().load(
        dict(feature_collections=[featcollection], genes=[gene], sequence_name="chr1")
    )

    def test_gff3_export(self, test_data_dir):
        """Test that the GFF3 export did not change over time"""
        gff = test_data_dir / "feature_test_2.gff"
        lines = [x.rstrip() for x in open(gff) if not x.startswith("#")]
        assert lines == [str(x) for x in self.annot.to_annotation_collection().to_gff()]

    def test_gff3_locus_tag_error(self, test_data_dir):
        with pytest.raises(GFF3LocusTagError):
            _ = list(parse_standard_gff3(test_data_dir / "feature_test_locus_tag_error.gff"))

    def test_gff3_strand_error(self, tmp_path):
        tmp_gff3 = tmp_path / "tmp.gff3"

        annot = self.annot.to_annotation_collection().to_dict()
        annot["feature_collections"][0]["feature_intervals"][0]["strand"] = "MINUS"
        annot = AnnotationCollectionModel.Schema().load(annot).to_annotation_collection()

        with open(tmp_gff3, "w") as fh:
            collection_to_gff3([annot], fh)
        with pytest.raises(GFF3ChildParentMismatchError):
            _ = list(parse_standard_gff3(tmp_gff3))

    def test_gff3_chromosome(self, tmp_path):
        tmp_gff3 = tmp_path / "tmp.gff3"

        annot = self.annot.to_annotation_collection().to_dict()
        annot["feature_collections"][0]["feature_intervals"][0]["sequence_name"] = "chr2"
        annot = AnnotationCollectionModel.Schema().load(annot).to_annotation_collection()

        with open(tmp_gff3, "w") as fh:
            collection_to_gff3([annot], fh)
        with pytest.raises(GFF3ChildParentMismatchError):
            _ = list(parse_standard_gff3(tmp_gff3))
