"""
Object representation of features. Includes an abstract feature class that is also used by transcripts.

Each object is capable of exporting itself to BED and GFF3.
"""
from typing import Optional, Any, Dict, List, Set, Iterable, Hashable
from uuid import UUID

from inscripta.biocantor.exc import (
    EmptyLocationException,
    NoSuchAncestorException,
)
from inscripta.biocantor.gene.cds_frame import CDSPhase
from inscripta.biocantor.gene.interval import AbstractFeatureInterval, QualifierValue, IntervalType
from inscripta.biocantor.io.bed import BED12, RGB
from inscripta.biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorFeatureTypes, BioCantorQualifiers
from inscripta.biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from inscripta.biocantor.io.gff3.rows import GFFAttributes, GFFRow
from inscripta.biocantor.location.location import Location
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent.parent import Parent, SequenceType
from inscripta.biocantor.util.bins import bins
from inscripta.biocantor.util.hashing import digest_object


class FeatureInterval(AbstractFeatureInterval):
    """FeatureIntervals are generic intervals. These can be used to model genome promoters,
    open chromatin sites, etc.
    """

    interval_type = IntervalType.FEATURE
    _identifiers = ["feature_name", "feature_id"]

    def __init__(
        self,
        interval_starts: List[int],
        interval_ends: List[int],
        strand: Strand,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        feature_types: Optional[List[str]] = None,
        feature_name: Optional[str] = None,
        feature_id: Optional[str] = None,
        guid: Optional[UUID] = None,
        feature_guid: Optional[UUID] = None,
        is_primary_feature: Optional[bool] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        self._location = self.initialize_location(interval_starts, interval_ends, strand, parent_or_seq_chunk_parent)
        self._genomic_starts = interval_starts
        self._genomic_ends = interval_ends
        self.start = self.genomic_start = interval_starts[0]
        self.end = self.genomic_end = interval_ends[-1]
        self._strand = strand
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent
        self.sequence_guid = sequence_guid
        self.sequence_name = sequence_name
        self.feature_types = set(feature_types) if feature_types else set()  # stored as a set of types
        self.feature_name = feature_name
        self.feature_id = feature_id
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.bin = bins(self.start, self.end, fmt="bed")
        self._is_primary_feature = is_primary_feature

        if guid is None:
            self.guid = digest_object(
                self._genomic_starts,
                self._genomic_ends,
                self.qualifiers,
                self.sequence_name,
                self.feature_types,
                self.feature_name,
                self.feature_id,
                self.is_primary_feature,
            )
        else:
            self.guid = guid
        self.feature_guid = feature_guid

    def __str__(self):
        return f"FeatureInterval(({self.chromosome_location}), name={self.feature_name})"

    def __repr__(self):
        return "<{}>".format(str(self))

    @property
    def id(self) -> str:
        """Returns the ID of this feature. Provides a shared API across genes/transcripts and features."""
        return self.feature_id

    @property
    def name(self) -> str:
        """Returns the name of this feature. Provides a shared API across genes/transcripts and features."""
        return self.feature_name

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`biocantor.io.models.FeatureIntervalModel`."""
        if chromosome_relative_coordinates:
            interval_starts = self._genomic_starts
            interval_ends = self._genomic_ends
        else:
            interval_starts, interval_ends = list(zip(*((x.start, x.end) for x in self.relative_blocks)))

        return dict(
            interval_starts=interval_starts,
            interval_ends=interval_ends,
            strand=self.strand.name,
            qualifiers=self._export_qualifiers_to_list(),
            feature_id=self.feature_id,
            feature_name=self.feature_name,
            feature_types=sorted(self.feature_types) if self.feature_types else None,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            feature_interval_guid=self.guid,
            feature_guid=self.feature_guid,
            is_primary_feature=self._is_primary_feature,
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "FeatureInterval":
        """Build a :class:`FeatureInterval` from a dictionary."""
        return FeatureInterval(
            interval_starts=vals["interval_starts"],
            interval_ends=vals["interval_ends"],
            strand=Strand[vals["strand"]],
            qualifiers=vals["qualifiers"],
            sequence_guid=vals["sequence_guid"],
            sequence_name=vals["sequence_name"],
            feature_types=vals["feature_types"],
            feature_name=vals["feature_name"],
            feature_id=vals["feature_id"],
            guid=vals["feature_interval_guid"],
            feature_guid=vals["feature_guid"],
            is_primary_feature=vals["is_primary_feature"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_location(
        location: Location,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        guid: Optional[UUID] = None,
        feature_guid: Optional[UUID] = None,
        feature_types: Optional[List[str]] = None,
        feature_id: Optional[str] = None,
        feature_name: Optional[str] = None,
        is_primary_feature: Optional[str] = None,
    ) -> "FeatureInterval":

        return FeatureInterval(
            interval_starts=[x.start for x in location.blocks],
            interval_ends=[x.end for x in location.blocks],
            strand=location.strand,
            guid=guid,
            feature_guid=feature_guid,
            qualifiers=qualifiers,
            sequence_name=sequence_name,
            sequence_guid=sequence_guid,
            feature_types=feature_types,
            feature_id=feature_id,
            feature_name=feature_name,
            is_primary_feature=is_primary_feature,
            parent_or_seq_chunk_parent=location.parent,
        )

    def intersect(
        self,
        location: Location,
        new_guid: Optional[UUID] = None,
        new_qualifiers: Optional[dict] = None,
    ) -> "FeatureInterval":
        """Returns a new FeatureInterval representing the intersection of this FeatureInterval's location with the
        other location.

        Strand of the other location is ignored; returned FeatureInterval is on the same strand as this FeatureInterval.
        """
        if not new_qualifiers:
            new_qualifiers = self.qualifiers

        location_same_strand = location.reset_strand(self.chromosome_location.strand)
        intersection = self.chromosome_location.intersection(location_same_strand)

        if intersection.is_empty:
            raise EmptyLocationException("Can't intersect disjoint intervals")

        starts = [x.start for x in intersection.blocks]
        ends = [x.end for x in intersection.blocks]
        return FeatureInterval(
            starts,
            ends,
            strand=intersection.strand,
            guid=new_guid,
            qualifiers=new_qualifiers,
            parent_or_seq_chunk_parent=intersection.parent,
        )

    def export_qualifiers(
        self, parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None
    ) -> Dict[Hashable, Set[str]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self._merge_qualifiers(parent_qualifiers)
        for key, val in [
            [BioCantorQualifiers.FEATURE_SYMBOL.value, self.feature_name],
            [BioCantorQualifiers.FEATURE_ID.value, self.feature_id],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        if self.feature_types:
            qualifiers[BioCantorQualifiers.FEATURE_TYPE.value] = self.feature_types
        return qualifiers

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterable[GFFRow]:
        """Writes a GFF format list of lists for this feature.

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
            GFF3MissingSequenceNameError: If there are no sequence names associated with this feature.
        """

        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers(parent_qualifiers)

        feature_id = str(self.guid)

        attributes = GFFAttributes(id=feature_id, qualifiers=qualifiers, name=self.feature_name, parent=parent)

        # "transcript" (feature interval) feature
        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.FEATURE_INTERVAL,
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
                id=f"feature-{feature_id}-{i}", qualifiers=qualifiers, name=self.feature_name, parent=feature_id
            )
            row = GFFRow(
                self.sequence_name,
                GFF_SOURCE,
                BioCantorFeatureTypes.FEATURE_INTERVAL_REGION,
                start + 1,
                end,
                NULL_COLUMN,
                self.strand,
                CDSPhase.NONE,
                attributes,
            )
            yield row

    def to_bed12(
        self,
        score: Optional[int] = 0,
        rgb: Optional[RGB] = RGB(0, 0, 0),
        name: Optional[str] = "feature_name",
        chromosome_relative_coordinates: bool = True,
    ) -> BED12:
        """Write a BED12 format representation of this :class:`FeatureInterval`.

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
            num_blocks = self.chunk_relative_location.num_blocks
        block_sizes = [end - start for start, end in blocks]
        block_starts = [start - self.start for start, _ in blocks]

        if chromosome_relative_coordinates:
            start = self.start
            end = self.end
        else:
            start = self.chunk_relative_start
            end = self.chunk_relative_end

        return BED12(
            self.sequence_name,
            start,
            end,
            getattr(self, name, name),
            score,
            self.strand,
            0,  # thickStart always 0 for non-coding
            0,  # thickEnd always 0 for non-coding
            rgb,
            num_blocks,
            block_sizes,
            block_starts,
        )
