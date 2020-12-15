"""
FeatureIntervals have the ability to write to BED.

"""
import pytest
from inscripta.biocantor.gene.cds import CDSFrame
from inscripta.biocantor.location.strand import Strand

from inscripta.biocantor.models import (
    TranscriptIntervalModel,
    FeatureIntervalModel,
)
from inscripta.biocantor.util.bed import RGB
from inscripta.biocantor.util.hashing import digest_object


class TestBedWriter:
    tx1 = dict(
        exon_starts=[2],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[5],
        cds_ends=[9],
        cds_frames=[CDSFrame.ZERO.name],
        sequence_symbol="chr1",
    )
    tx2 = dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        guid=digest_object(123),
        transcript_symbol="name",
    )
    feat1 = dict(interval_starts=[2], interval_ends=[5], strand=Strand.PLUS.name, sequence_symbol="chr10")
    feat2 = dict(interval_starts=[2, 7, 12], interval_ends=[6, 10, 15], strand=Strand.PLUS.name)
    feat3 = dict(interval_starts=[25], interval_ends=[30], strand=Strand.MINUS.name)

    @pytest.mark.parametrize(
        "tx,expected",
        [
            (tx1, ["chr1", "2", "18", "None", "0", "+", "5", "9", "0,0,0", "1", "16", "0"]),
            (
                tx2,
                [
                    "None",
                    "2",
                    "15",
                    "202cb962-ac59-075b-964b-07152d234b70",
                    "0",
                    "+",
                    "4",
                    "13",
                    "0,0,0",
                    "3",
                    "4,3,3",
                    "0,5,10",
                ],
            ),
        ],
    )
    def test_tx(self, tx, expected):
        model = TranscriptIntervalModel.Schema().load(tx)
        obj = model.to_transcript_interval()
        assert str(obj.to_bed12()) == "\t".join(expected)

    @pytest.mark.parametrize(
        "feat,expected",
        [
            (feat1, ["chr10", "2", "5", "None", "0", "+", "0", "0", "0,0,0", "1", "3", "0"]),
            (feat2, ["None", "2", "15", "None", "0", "+", "0", "0", "0,0,0", "3", "4,3,3", "0,5,10"]),
            (feat3, ["None", "25", "30", "None", "0", "-", "0", "0", "0,0,0", "1", "5", "0"]),
        ],
    )
    def test_feat(self, feat, expected):
        model = FeatureIntervalModel.Schema().load(feat)
        obj = model.to_feature_interval()
        assert str(obj.to_bed12()) == "\t".join(expected)

    @pytest.mark.parametrize(
        "tx,score,rgb,name,expected",
        [
            (
                tx2,
                10,
                RGB(128, 128, 128),
                "transcript_symbol",
                ["None", "2", "15", "name", "10", "+", "4", "13", "128,128,128", "3", "4,3,3", "0,5,10"],
            ),
            (  # if name is not an attribute, just pass it along
                tx2,
                10,
                RGB(128, 128, 128),
                "test",
                ["None", "2", "15", "test", "10", "+", "4", "13", "128,128,128", "3", "4,3,3", "0,5,10"],
            ),
        ],
    )
    def test_changed_metadata(self, tx, score, rgb, name, expected):
        model = TranscriptIntervalModel.Schema().load(tx)
        obj = model.to_transcript_interval()
        assert str(obj.to_bed12(score, rgb, name)) == "\t".join(expected)
