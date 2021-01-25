"""
Object representation of Transcripts.

Each object is capable of exporting itself to BED and GFF3.
"""
from itertools import count
from typing import Optional, Any, Dict, Iterable, Hashable, Set
from uuid import UUID

from inscripta.biocantor.exc import (
    EmptyLocationException,
)
from inscripta.biocantor.exc import NoncodingTranscriptError
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSInterval, CDSPhase
from inscripta.biocantor.gene.feature import AbstractFeatureInterval, QualifierValue
from inscripta.biocantor.io.bed import BED12, RGB
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorQualifiers, BioCantorFeatureTypes
from inscripta.biocantor.io.gff3.rows import GFFAttributes, GFFRow
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.location_impl import SingleInterval, EmptyLocation
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object
from inscripta.biocantor.util.object_validation import ObjectValidation
from methodtools import lru_cache


class TranscriptInterval(AbstractFeatureInterval):
    """Transcripts are a collection of Intervals. The canonical transcript has a 5' UTR,
    a CDS, and a 3' UTR. However, a transcript may be lacking any of those three features --
    it could be noncoding (no CDS), truncated (no 3' UTR or 5' UTR). It could also be a single exon.

    This object represents all of these possibilities through one Location that represents the Exons,
    and an optional CDS that represents the :class:`~biocantor.gene.cds.CDSInterval`.

    In addition to the interval information, this object has optional metadata and an optional understanding
    of sequence.
    """

    _identifiers = ["transcript_id", "transcript_symbol", "protein_id"]

    def __init__(
        self,
        location: Location,  # exons
        cds: Optional[CDSInterval] = None,  # optional CDS with frame
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        is_primary_tx: Optional[bool] = None,
        transcript_id: Optional[str] = None,
        transcript_symbol: Optional[str] = None,
        transcript_type: Optional[Biotype] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        protein_id: Optional[str] = None,
        guid: Optional[UUID] = None,
        transcript_guid: Optional[UUID] = None,
    ):
        self.location = location  # genomic CompoundInterval
        self._is_primary_feature = is_primary_tx
        self.transcript_id = transcript_id
        self.transcript_symbol = transcript_symbol
        self.transcript_type = transcript_type
        self.protein_id = protein_id
        self.sequence_guid = sequence_guid
        self.sequence_name = sequence_name
        self.bin = bins(self.start, self.end, fmt="bed")
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.cds = cds

        if guid is None:
            self.guid = digest_object(
                self.location,
                self.cds,
                self.qualifiers,
                self.transcript_id,
                self.transcript_symbol,
                self.transcript_type,
                self.protein_id,
                self.sequence_name,
                self.is_primary_tx,
            )
        else:
            self.guid = guid
        self.transcript_guid = transcript_guid

        if self.location.parent:
            ObjectValidation.require_location_has_parent_with_sequence(self.location)
            if self.cds:
                ObjectValidation.require_locations_have_same_nonempty_parent(location, cds.location)

    def __str__(self):
        return f"TranscriptInterval(({self.location}), cds=[{self.cds}], symbol={self.transcript_symbol})"

    def __repr__(self):
        return "<{}>".format(str(self))

    @property
    def is_primary_tx(self) -> bool:
        """Is this the primary transcript?"""
        return self.is_primary_feature

    def __len__(self):
        return len(self.location)

    @property
    def is_coding(self) -> bool:
        return self.cds is not None

    @property
    def has_in_frame_stop(self) -> bool:
        if not self.is_coding:
            raise NoncodingTranscriptError("Cannot have frameshifts on non-coding transcripts")
        return self.cds.has_in_frame_stop

    @property
    def cds_size(self) -> int:
        if self.is_coding:
            return len(self.cds)
        return 0

    @property
    def cds_start(self) -> int:
        if self.is_coding:
            return self.cds.location.start
        else:
            raise NoncodingTranscriptError("No CDS start for non-coding transcript")

    @property
    def cds_end(self) -> int:
        if self.is_coding:
            return self.cds.location.end
        else:
            raise NoncodingTranscriptError("No CDS end for non-coding transcript")

    def update_parent(self, parent: Parent):
        """Change parent to this new parent."""
        self.location = self.location.reset_parent(parent)
        self.cds.location = self.cds.location.reset_parent(parent)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`biocantor.io.models.TranscriptIntervalModel`."""
        exon_starts, exon_ends = list(zip(*([x.start, x.end] for x in self.location.blocks)))
        if self.cds:
            cds_starts, cds_ends = list(zip(*([x.start, x.end] for x in self.cds.location.blocks)))
            cds_frames = [f.name for f in self.cds.frames]
        else:
            cds_starts = None
            cds_ends = None
            cds_frames = None
        return dict(
            exon_starts=exon_starts,
            exon_ends=exon_ends,
            strand=self.strand.name,
            cds_starts=cds_starts,
            cds_ends=cds_ends,
            cds_frames=cds_frames,
            qualifiers=self._export_qualifiers_to_list(),
            is_primary_tx=self._is_primary_feature,
            transcript_id=self.transcript_id,
            transcript_symbol=self.transcript_symbol,
            transcript_type=self.transcript_type.name if self.transcript_type else None,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            protein_id=self.protein_id,
            transcript_guid=self.transcript_guid,
            transcript_interval_guid=self.guid,
        )

    def merge_overlapping(self) -> "TranscriptInterval":
        """Some tools can't handle overlapping intervals. This function discards CDS information."""
        if not self.location.is_overlapping:
            return self
        new_loc = self.location.merge_overlapping()
        tx = TranscriptInterval(
            new_loc,
            cds=None,
            qualifiers=self._export_qualifiers_to_list(),
            is_primary_tx=self.is_primary_tx,
            transcript_id=self.transcript_id,
            transcript_symbol=self.transcript_symbol,
            transcript_type=self.transcript_type,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            guid=self.guid,
        )
        return tx

    def intersect(
        self,
        location: Location,
        new_guid: Optional[UUID] = None,
        new_qualifiers: Optional[dict] = None,
    ) -> "TranscriptInterval":
        """Returns a new TranscriptInterval representing the intersection of this Transcript's location with the
        other location.

        Strand of the other location is ignored; returned Transcript is on the same strand as this Transcript.

        If this Transcript has a CDS, it will be dropped because CDS intersections are not currently supported.
        """
        if not new_qualifiers:
            new_qualifiers = self.qualifiers

        location_same_strand = location.reset_strand(self.location.strand)
        intersection = self.location.intersection(location_same_strand)

        if intersection.is_empty:
            raise EmptyLocationException("Can't intersect disjoint intervals")

        return TranscriptInterval(
            location=intersection,
            guid=new_guid,
            qualifiers=new_qualifiers,
        )

    def sequence_pos_to_transcript(self, pos: int) -> int:
        """Converts sequence position to relative position along this transcript."""
        return self.sequence_pos_to_feature(pos)

    def sequence_interval_to_transcript(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within this transcript."""
        return self.sequence_interval_to_feature(chr_start, chr_end, chr_strand)

    def transcript_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along this transcript to sequence coordinate."""
        return self.feature_pos_to_sequence(pos)

    def transcript_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to this transcript to a spliced location on the sequence."""
        return self.feature_interval_to_sequence(rel_start, rel_end, rel_strand)

    def cds_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along the CDS to sequence coordinate."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.cds.location.relative_to_parent_pos(pos)

    def cds_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to the CDS to a spliced location on the sequence."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.cds.location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    def sequence_pos_to_cds(self, pos: int) -> int:
        """Converts sequence position to relative position along the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.cds.location.parent_to_relative_pos(pos)

    def sequence_interval_to_cds(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.cds.location.parent_to_relative_location(
            SingleInterval(chr_start, chr_end, chr_strand, parent=self.location.parent)
        )

    def cds_pos_to_transcript(self, pos: int) -> int:
        """Converts a relative position along the CDS to a relative position along this transcript."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        chr_pos = self.cds_pos_to_sequence(pos)
        return self.sequence_pos_to_transcript(chr_pos)

    def transcript_pos_to_cds(self, pos: int) -> int:
        """Converts a relative position along this transcript to a relative position along the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        chr_pos = self.transcript_pos_to_sequence(pos)
        return self.sequence_pos_to_cds(chr_pos)

    def get_5p_interval(self) -> Location:
        """Return the 5' UTR as a Location, if it exists."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No 5' UTR on a non-coding transcript")
        # handle the edge case where the CDS is full length
        if self.cds.location == self.location:
            return EmptyLocation()
        cds_start_on_transcript = self.cds_pos_to_transcript(0)
        return self.location.relative_interval_to_parent_location(0, cds_start_on_transcript, Strand.PLUS)

    def get_3p_interval(self) -> Location:
        """Returns the 3' UTR as a location, if it exists."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No 3' UTR on a non-coding transcript")
        # handle the edge case where the CDS is full length
        if self.cds.location == self.location:
            return EmptyLocation()
        cds_inclusive_end_on_transcript = self.cds_pos_to_transcript(len(self.cds.location) - 1)
        return self.location.relative_interval_to_parent_location(
            cds_inclusive_end_on_transcript + 1, len(self.location), Strand.PLUS
        )

    def get_coding_interval(self) -> Location:
        """Get coding interval."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No coding interval on non-coding transcript")
        return self.cds.location

    @lru_cache(maxsize=1)
    def get_transcript_sequence(self) -> Sequence:
        """Returns the mRNA sequence."""
        return self.get_spliced_sequence()

    @lru_cache(maxsize=1)
    def get_cds_sequence(self) -> Sequence:
        """Returns the in-frame CDS sequence (always multiple of 3)."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS sequence on non-coding transcript")
        return self.cds.extract_sequence()

    @lru_cache(maxsize=2)
    def get_protein_sequence(self, truncate_at_in_frame_stop: Optional[bool] = False) -> Sequence:
        """Return the translation of this transcript, if possible."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No translation on non-coding transcript")
        return self.cds.translate(truncate_at_in_frame_stop)

    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None
    ) -> Dict[Hashable, Set[Hashable]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self._merge_qualifiers(parent_qualifiers)
        for key, val in [
            [BioCantorQualifiers.TRANSCRIPT_ID.value, self.transcript_id],
            [BioCantorQualifiers.TRANSCRIPT_NAME.value, self.transcript_symbol],
            [BioCantorQualifiers.TRANSCRIPT_TYPE.value, self.transcript_type.name if self.transcript_type else None],
            [BioCantorQualifiers.PROTEIN_ID.value, self.protein_id],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        return qualifiers

    def to_gff(
        self, parent: Optional[str] = None, parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None
    ) -> Iterable[GFFRow]:
        """Writes a GFF format list of lists for this gene.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.

        Yields:
            :class:`~biocantor.util.gff3.rows.GFFRow`
        """
        qualifiers = self.export_qualifiers(parent_qualifiers)

        tx_guid = str(self.guid) if self.guid else str(digest_object(self))

        attributes = GFFAttributes(id=tx_guid, qualifiers=qualifiers, name=self.transcript_symbol, parent=parent)

        # transcript feature
        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.TRANSCRIPT,
            self.start + 1,
            self.end,
            NULL_COLUMN,
            self.strand,
            CDSPhase.NONE,
            attributes,
        )
        yield row

        # start adding exon features
        # re-use qualifiers, updating ID each time
        for i, block in enumerate(self.location.blocks, 1):
            attributes = GFFAttributes(
                id=f"exon-{tx_guid}-{i}", qualifiers=qualifiers, name=self.transcript_symbol, parent=tx_guid
            )
            row = GFFRow(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.EXON,
                block.start + 1,
                block.end,
                NULL_COLUMN,
                self.strand,
                CDSPhase.NONE,
                attributes,
            )
            yield row

        # add CDS features, if applicable
        if self.cds:
            for i, block, frame in zip(count(1), self.cds.location.blocks, self.cds.frames):
                attributes = GFFAttributes(
                    id=f"cds-{tx_guid}-{i}",
                    qualifiers=qualifiers,
                    name=self.transcript_symbol,
                    parent=tx_guid,
                )
                row = GFFRow(
                    self.sequence_name,
                    GFF_SOURCE,
                    BioCantorFeatureTypes.CDS,
                    block.start + 1,
                    block.end,
                    NULL_COLUMN,
                    self.strand,
                    frame.to_phase(),
                    attributes,
                )
                yield row

    def to_bed12(
        self, score: Optional[int] = 0, rgb: Optional[RGB] = RGB(0, 0, 0), name: Optional[str] = "transcript_symbol"
    ) -> BED12:
        block_sizes = [b.end - b.start for b in self.location.blocks]
        block_starts = [b.start - self.start for b in self.location.blocks]
        return BED12(
            self.sequence_name,
            self.start,
            self.end,
            getattr(self, name, name),
            score,
            self.strand,
            self.cds.start if self.cds is not None else 0,
            self.cds.end if self.cds is not None else 0,
            rgb,
            len(self.location.blocks),
            block_sizes,
            block_starts,
        )
