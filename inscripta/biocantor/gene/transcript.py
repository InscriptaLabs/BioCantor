"""
Object representation of Transcripts.

Each object is capable of exporting itself to BED and GFF3.
"""
from itertools import count
from typing import Optional, Any, Dict, Iterable, Hashable, Set, List
from uuid import UUID

from methodtools import lru_cache

from inscripta.biocantor.exc import (
    EmptyLocationException,
    LocationOverlapException,
    NoncodingTranscriptError,
    InvalidCDSIntervalError,
)
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSInterval, CDSPhase, CDSFrame
from inscripta.biocantor.gene.feature import AbstractFeatureInterval, QualifierValue
from inscripta.biocantor.io.bed import BED12, RGB
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorQualifiers, BioCantorFeatureTypes
from inscripta.biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from inscripta.biocantor.io.gff3.rows import GFFAttributes, GFFRow
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.location_impl import SingleInterval, EmptyLocation
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent
from inscripta.biocantor.sequence.sequence import Sequence
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object
from inscripta.biocantor.util.object_validation import ObjectValidation


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
        exon_starts: List[int],
        exon_ends: List[int],
        strand: Strand,
        cds_starts: Optional[List[int]] = None,
        cds_ends: Optional[List[int]] = None,
        cds_frames: Optional[List[CDSFrame]] = None,
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
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        self._location = TranscriptInterval.initialize_location(
            exon_starts,
            exon_ends,
            strand,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

        if cds_starts is not None and cds_ends is None:
            raise InvalidCDSIntervalError("If CDS start is defined, CDS end must be defined")
        elif cds_starts is None and cds_ends is not None:
            raise InvalidCDSIntervalError("If CDS end is defined, CDS start must be defined")
        elif cds_starts is not None and cds_ends is not None:  # must be coding
            if len(cds_starts) != len(cds_ends):
                raise InvalidCDSIntervalError("Number of CDS starts does not number of CDS ends")
            elif cds_starts[0] < exon_starts[0]:
                raise InvalidCDSIntervalError("CDS start must be greater than or equal to exon start")
            elif cds_ends[-1] > exon_ends[-1]:
                raise InvalidCDSIntervalError("CDS end must be less than or equal to than exon end")
            elif cds_frames is None:
                raise InvalidCDSIntervalError("If CDS interval is defined, CDS frames must be defined")
            elif len(cds_frames) != len(cds_starts):
                raise InvalidCDSIntervalError("Number of CDS frames must match number of CDS starts/ends")

            # as a result of a parent or seq chunk parent constructor, it may be the case that this CDS is entirely
            # sliced out. Check this case, and then void out the CDS.
            try:
                cds_interval = TranscriptInterval.initialize_location(
                    cds_starts, cds_ends, strand, parent_or_seq_chunk_parent=parent_or_seq_chunk_parent
                )

                if len(cds_interval) == 0:
                    raise InvalidCDSIntervalError("Cannot have an empty CDS interval")

                # we may end up having a reduced number of CDSFrames now due to this region being a subset of this gene
                self._cds = CDSInterval(cds_interval, cds_frames[: len(cds_interval.blocks)])

            except LocationOverlapException:
                self._cds = None

            self._cds_frames = cds_frames

        else:
            self._cds = self._cds_frames = self._cds_start = self._cds_end = None

        if self._location.parent:
            ObjectValidation.require_location_has_parent_with_sequence(self._location)
            if self._cds:
                ObjectValidation.require_locations_have_same_nonempty_parent(self._location, self._cds._location)

        self._genomic_starts = exon_starts
        self._genomic_ends = exon_ends
        self.start = self.genomic_start = exon_starts[0]
        self.end = self.genomic_end = exon_ends[-1]
        self._genomic_cds_starts = cds_starts
        self._genomic_cds_ends = cds_ends

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

        if guid is None:
            self.guid = digest_object(
                self._genomic_starts,
                self._genomic_ends,
                self._genomic_cds_starts,
                self._genomic_cds_ends,
                self._cds_frames,
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

    def __str__(self):
        return f"TranscriptInterval(({self._location}), cds=[{self._cds}], symbol={self.transcript_symbol})"

    def __repr__(self):
        return "<{}>".format(str(self))

    @property
    def is_primary_tx(self) -> bool:
        """Is this the primary transcript?"""
        return self.is_primary_feature

    def __len__(self):
        return sum((end - start) for end, start in zip(self._genomic_ends, self._genomic_starts))

    @property
    def cds_location(self) -> Location:
        """Returns the Location of the CDS in *chromosome coordinates*"""
        return self.lift_cds_over_to_first_ancestor_of_type("chromosome")

    @property
    def cds_chunk_relative_location(self) -> Location:
        """Returns the Location of the CDS in *chunk relative coordinates*"""
        if not self.is_coding:
            raise NoncodingTranscriptError("No location on a non-coding transcript")
        return self._cds._location

    @property
    def cds_blocks(self) -> List[Location]:
        """Returns the blocks of the CDS, if it exists."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No blocks on a non-coding transcript")
        return self.cds_location.blocks

    @property
    def cds_chunk_relative_blocks(self) -> List[Location]:
        """Returns the chunk relative blocks of the CDS, if it exists."""
        return self.cds_chunk_relative_location.blocks

    @property
    def is_coding(self) -> bool:
        return self._cds is not None

    @property
    def has_in_frame_stop(self) -> bool:
        if not self.is_coding:
            raise NoncodingTranscriptError("Cannot have frameshifts on non-coding transcripts")
        return self._cds.has_in_frame_stop

    @property
    def cds_size(self) -> int:
        if self.is_coding:
            return sum((end - start) for end, start in zip(self._genomic_cds_ends, self._genomic_cds_starts))
        return 0

    @property
    def chunk_relative_cds_size(self) -> int:
        if self.is_coding:
            return len(self._cds)
        return 0

    @property
    def cds_start(self) -> int:
        if self.is_coding:
            return self._genomic_cds_starts[0]
        else:
            raise NoncodingTranscriptError("No CDS start for non-coding transcript")

    @property
    def cds_end(self) -> int:
        if self.is_coding:
            return self._genomic_cds_ends[-1]
        else:
            raise NoncodingTranscriptError("No CDS end for non-coding transcript")

    @property
    def chunk_relative_cds_start(self) -> int:
        if self.is_coding:
            return self._cds._location.start
        else:
            raise NoncodingTranscriptError("No CDS start for non-coding transcript")

    @property
    def chunk_relative_cds_end(self) -> int:
        if self.is_coding:
            return self._cds._location.end
        else:
            raise NoncodingTranscriptError("No CDS end for non-coding transcript")

    @property
    def cds_blocks(self) -> Iterable[SingleInterval]:
        """Wrapper for blocks function that reports blocks in chromosome coordinates"""
        if self.is_coding:
            yield from self.lift_cds_over_to_first_ancestor_of_type("chromosome").blocks
        else:
            raise NoncodingTranscriptError("No CDS blocks for non-coding transcript")

    @property
    def chunk_relative_cds_blocks(self) -> Iterable[SingleInterval]:
        """Wrapper for blocks function that reports blocks in chunk-relative coordinates"""
        if self.is_coding:
            yield from self._cds._location.blocks
        else:
            raise NoncodingTranscriptError("No relative CDS blocks for non-coding transcript")

    @property
    def id(self) -> str:
        """Returns the ID of this transcript. Provides a shared API across genes/transcripts and features."""
        return self.transcript_id

    @property
    def name(self) -> str:
        """Returns the name of this transcript. Provides a shared API across genes/transcripts and features."""
        return self.transcript_symbol

    def reset_parent(self, parent: Optional[Parent] = None):
        """
        Convenience function that wraps location.reset_parent(). Overrides the parent function in
        :class:`~biocantor.gene.feature.AbstractInterval` in order to also update the CDS interval.
        """
        if self.is_coding:
            self._cds._location = self._cds._location.reset_parent(parent)
        super().reset_parent(parent)

    def _liftover_this_location_to_seq_chunk_parent(
        self,
        seq_chunk_parent: Parent,
    ):
        """Lift *this* transcript interval to a new subset. Overrides the parent method in order to perform
        this task on the CDS interval as well.

        This could happen as the result of a subsetting operation.

        This will introduce chunk-relative coordinates to this interval, or reduce the size of existing chunk-relative
        coordinates.

        This function calls the parent static method :meth:`AbstractInterval.liftover_location_to_seq_chunk_parent()`,
        but differs in two key ways:

        1. It acts on an instantiated subclass of this abstract class, modifying the location.
        2. It handles the case where a subclass is already a slice, by first lifting up to genomic coordinates.

        For these reasons, and particularly #1, this is a private method that is intended to be used during
        construction of a subclass. Modifying the locations in-place are generally a bad idea after initial
        construction of a interval class.
        """
        super()._liftover_this_location_to_seq_chunk_parent(seq_chunk_parent)
        self._cds._location = self.liftover_location_to_seq_chunk_parent(
            self._cds._location.lift_over_to_first_ancestor_of_type("chromosome").reset_parent(seq_chunk_parent.parent)
        )

    def lift_cds_over_to_first_ancestor_of_type(self, sequence_type: Optional[str] = "chromosome") -> Location:
        """
        Lifts the CDS location member to another coordinate system. Is a no-op if there is no parent assigned.

        Returns:
            The lifted Location.
        """
        if self.is_coding:
            if self._cds._location.parent is None:
                return self._cds._location
            return self._cds._location.lift_over_to_first_ancestor_of_type(sequence_type)
        else:
            raise NoncodingTranscriptError("No CDS location for non-coding transcript")

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`biocantor.io.models.TranscriptIntervalModel`."""
        if chromosome_relative_coordinates:
            exon_starts = self._genomic_starts
            exon_ends = self._genomic_ends
        else:
            exon_starts, exon_ends = list(zip(*((x.start, x.end) for x in self.relative_blocks)))
        if self._cds:
            if chromosome_relative_coordinates:
                cds_starts = self._genomic_cds_starts
                cds_ends = self._genomic_cds_ends
            else:
                cds_starts, cds_ends = list(zip(*([x.start, x.end] for x in self.chunk_relative_cds_blocks)))
            cds_frames = [f.name for f in self._cds.frames]
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

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "TranscriptInterval":
        """Build a :class:`TranscriptInterval` from a dictionary."""

        return TranscriptInterval(
            exon_starts=vals["exon_starts"],
            exon_ends=vals["exon_ends"],
            strand=Strand[vals["strand"]],
            cds_starts=vals["cds_starts"] if vals["cds_starts"] else None,
            cds_ends=vals["cds_ends"] if vals["cds_ends"] else None,
            cds_frames=[CDSFrame[x] for x in vals["cds_frames"]] if vals["cds_frames"] else None,
            guid=vals["transcript_interval_guid"],
            transcript_guid=vals["transcript_guid"],
            qualifiers=vals["qualifiers"],
            is_primary_tx=vals["is_primary_tx"],
            transcript_id=vals["transcript_id"],
            transcript_symbol=vals["transcript_symbol"],
            transcript_type=Biotype[vals["transcript_type"]] if vals["transcript_type"] else None,
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            protein_id=vals["protein_id"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

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

        location_same_strand = location.reset_strand(self._location.strand)
        intersection = self._location.intersection(location_same_strand)

        if intersection.is_empty:
            raise EmptyLocationException("Can't intersect disjoint intervals")

        starts = [x.start for x in intersection.blocks]
        ends = [x.end for x in intersection.blocks]
        return TranscriptInterval(
            starts,
            ends,
            strand=intersection.strand,
            guid=new_guid,
            qualifiers=new_qualifiers,
            parent_or_seq_chunk_parent=intersection.parent,
        )

    def sequence_pos_to_transcript(self, pos: int) -> int:
        """Converts sequence position to relative position along this transcript."""
        return self.sequence_pos_to_feature(pos)

    def chunk_relative_pos_to_transcript(self, pos: int) -> int:
        """Converts chunk-relative sequence position to relative position along this transcript."""
        return self.chunk_relative_pos_to_feature(pos)

    def sequence_interval_to_transcript(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within this transcript."""
        return self.sequence_interval_to_feature(chr_start, chr_end, chr_strand)

    def chunk_relative_interval_to_transcript(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """
        Converts a contiguous interval on the chunk-relative sequence to a relative location within this transcript.
        """
        return self.chunk_relative_interval_to_feature(chr_start, chr_end, chr_strand)

    def transcript_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along this transcript to sequence coordinate."""
        return self.feature_pos_to_sequence(pos)

    def transcript_pos_to_chunk_relative(self, pos: int) -> int:
        """Converts a relative position along this transcript to chunk-relative sequence coordinate."""
        return self.feature_pos_to_chunk_relative(pos)

    def transcript_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to this transcript to a spliced location on the sequence."""
        return self.feature_interval_to_sequence(rel_start, rel_end, rel_strand)

    def transcript_interval_to_chunk_relative(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """
        Converts a contiguous interval relative to this transcript to a spliced location on the chunk-relative sequence.
        """
        return self.feature_interval_to_chunk_relative(rel_start, rel_end, rel_strand)

    def cds_pos_to_sequence(self, pos: int) -> int:
        """Converts a relative position along the CDS to sequence coordinate."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.lift_cds_over_to_first_ancestor_of_type("chromosome").relative_to_parent_pos(pos)

    def cds_pos_to_chunk_relative(self, pos: int) -> int:
        """Converts a relative position along the CDS to chunk-relative sequence coordinate."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self._cds._location.relative_to_parent_pos(pos)

    def cds_interval_to_sequence(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to the CDS to a spliced location on the sequence."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.lift_cds_over_to_first_ancestor_of_type("chromosome").relative_interval_to_parent_location(
            rel_start, rel_end, rel_strand
        )

    def cds_interval_to_chunk_relative(self, rel_start: int, rel_end: int, rel_strand: Strand) -> Location:
        """Converts a contiguous interval relative to the CDS to a spliced location on the chunk-relative sequence."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self._cds._location.relative_interval_to_parent_location(rel_start, rel_end, rel_strand)

    def sequence_pos_to_cds(self, pos: int) -> int:
        """Converts sequence position to relative position along the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self.lift_cds_over_to_first_ancestor_of_type("chromosome").parent_to_relative_pos(pos)

    def chunk_relative_pos_to_cds(self, pos: int) -> int:
        """Converts chunk-relative sequence position to relative position along the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self._cds._location.parent_to_relative_pos(pos)

    def sequence_interval_to_cds(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the sequence to a relative location within the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        loc = self.lift_cds_over_to_first_ancestor_of_type("chromosome")
        i = SingleInterval(chr_start, chr_end, chr_strand, parent=loc.parent)
        return loc.parent_to_relative_location(i)

    def chunk_relative_interval_to_cds(self, chr_start: int, chr_end: int, chr_strand: Strand) -> Location:
        """Converts a contiguous interval on the chunk-relative sequence to a relative location within the CDS."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS positions on non-coding transcript")
        return self._cds._location.parent_to_relative_location(
            SingleInterval(chr_start, chr_end, chr_strand, parent=self._cds._location.parent)
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
        if self._cds._location == self._location:
            return EmptyLocation()
        cds_start_on_transcript = self.cds_pos_to_transcript(0)
        return self._location.relative_interval_to_parent_location(0, cds_start_on_transcript, Strand.PLUS)

    def get_3p_interval(self) -> Location:
        """Returns the 3' UTR as a location, if it exists."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No 3' UTR on a non-coding transcript")
        # handle the edge case where the CDS is full length
        if self._cds._location == self._location:
            return EmptyLocation()
        cds_inclusive_end_on_transcript = self.cds_pos_to_transcript(len(self._cds._location) - 1)
        return self._location.relative_interval_to_parent_location(
            cds_inclusive_end_on_transcript + 1, len(self._location), Strand.PLUS
        )

    def get_coding_interval(self) -> Location:
        """Get coding interval."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No coding interval on non-coding transcript")
        return self._cds._location

    @lru_cache(maxsize=1)
    def get_transcript_sequence(self) -> Sequence:
        """Returns the mRNA sequence."""
        return self.get_spliced_sequence()

    @lru_cache(maxsize=1)
    def get_cds_sequence(self) -> Sequence:
        """Returns the in-frame CDS sequence (always multiple of 3)."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No CDS sequence on non-coding transcript")
        return self._cds.extract_sequence()

    @lru_cache(maxsize=2)
    def get_protein_sequence(self, truncate_at_in_frame_stop: Optional[bool] = False) -> Sequence:
        """Return the translation of this transcript, if possible."""
        if not self.is_coding:
            raise NoncodingTranscriptError("No translation on non-coding transcript")
        return self._cds.translate(truncate_at_in_frame_stop)

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
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterable[GFFRow]:
        """Writes a GFF format list of lists for this transcript.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
            GFF3MissingSequenceNameError: If there are no sequence names associated with this transcript.
        """

        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        qualifiers = self.export_qualifiers(parent_qualifiers)

        tx_guid = str(self.guid)

        attributes = GFFAttributes(id=tx_guid, qualifiers=qualifiers, name=self.transcript_symbol, parent=parent)

        # transcript feature
        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.TRANSCRIPT,
            (self.start if chromosome_relative_coordinates else self.chunk_relative_start) + 1,
            self.end if chromosome_relative_coordinates else self.chunk_relative_end,
            NULL_COLUMN,
            self.strand,
            CDSPhase.NONE,
            attributes,
        )
        yield row

        # start adding exon features
        # re-use qualifiers, updating ID each time
        if chromosome_relative_coordinates:
            blocks = zip(self._genomic_starts, self._genomic_ends)
        else:
            blocks = [[x.start, x.end] for x in self.relative_blocks]

        for i, (start, end) in enumerate(blocks, 1):
            attributes = GFFAttributes(
                id=f"exon-{tx_guid}-{i}", qualifiers=qualifiers, name=self.transcript_symbol, parent=tx_guid
            )
            row = GFFRow(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.EXON,
                start + 1,
                end,
                NULL_COLUMN,
                self.strand,
                CDSPhase.NONE,
                attributes,
            )
            yield row

        # add CDS features, if applicable
        if self._cds:
            if chromosome_relative_coordinates:
                cds_blocks = zip(self._genomic_cds_starts, self._genomic_cds_ends)
            else:
                cds_blocks = [[x.start, x.end] for x in self.chunk_relative_cds_blocks]

            for i, block, frame in zip(count(1), cds_blocks, self._cds.frames):
                start, end = block
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
                    start + 1,
                    end,
                    NULL_COLUMN,
                    self.strand,
                    frame.to_phase(),
                    attributes,
                )
                yield row

    def to_bed12(
        self,
        score: Optional[int] = 0,
        rgb: Optional[RGB] = RGB(0, 0, 0),
        name: Optional[str] = "transcript_symbol",
        chromosome_relative_coordinates: bool = True,
    ) -> BED12:
        """Write a BED12 format representation of this :class:`TranscriptInterval`.

        Both of these optional arguments are specific to the BED12 format.

        Args:
            score: An optional score associated with a interval. UCSC requires an integer between 0 and 1000.
            rgb: An optional RGB string for visualization on a browser. This allows you to have multiple colors
                on a single UCSC track.
            name: Which identifier in this record to use as 'name'. feature_name to guid. If the supplied string
                is not a valid attribute, it is used directly.
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Return:
            A :class:`~biocantor.io.bed.BED12` object.

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        if chromosome_relative_coordinates:
            blocks = list(zip(self._genomic_starts, self._genomic_ends))
            num_blocks = len(self._genomic_starts)
        else:
            blocks = [[x.start, x.end] for x in self.relative_blocks]
            num_blocks = self._location.num_blocks
        block_sizes = [end - start for start, end in blocks]
        block_starts = [start - self.start for start, _ in blocks]

        if chromosome_relative_coordinates:
            start = self.start
            end = self.end
            if self._cds:
                cds_start = self.cds_start
                cds_end = self.cds_end
            else:
                cds_start = cds_end = 0
        else:
            start = self.chunk_relative_start
            end = self.chunk_relative_end
            if self._cds:
                cds_start = self.chunk_relative_cds_start
                cds_end = self.chunk_relative_cds_end
            else:
                cds_start = cds_end = 0

        return BED12(
            self.sequence_name,
            start,
            end,
            getattr(self, name, name),
            score,
            self.strand,
            cds_start,
            cds_end,
            rgb,
            num_blocks,
            block_sizes,
            block_starts,
        )
