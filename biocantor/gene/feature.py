"""
Object representation of features. Includes an abstract feature class that is also used by transcripts.

Each object is capable of exporting itself to BED and GFF3.
"""
from functools import reduce
from typing import Optional, Any, Dict, List, Set, Iterable, Iterator, Hashable, Union, TYPE_CHECKING, Type
from uuid import UUID

from biocantor.exc import (
    EmptyLocationException,
    NoSuchAncestorException,
    NoncodingTranscriptError,
    InvalidAnnotationError,
    DuplicateFeatureError,
)
from biocantor.gene.cds_frame import CDSPhase
from biocantor.gene.interval import (
    AbstractFeatureInterval,
    QualifierValue,
    IntervalType,
    AbstractFeatureIntervalCollection,
)
from biocantor.io.bed import BED12, RGB
from biocantor.io.gff3.constants import GFF_SOURCE, NULL_COLUMN, BioCantorFeatureTypes, BioCantorQualifiers
from biocantor.io.gff3.exc import GFF3MissingSequenceNameError
from biocantor.io.gff3.rows import GFFAttributes, GFFRow, GTFAttributes, GTFRow
from biocantor.location.location import Location
from biocantor.location.strand import Strand
from biocantor.parent.parent import Parent, SequenceType
from biocantor.sequence import Sequence
from biocantor.util.bins import bins
from biocantor.util.hashing import digest_object

if TYPE_CHECKING:
    from biocantor.gene.variants import VariantIntervalCollection, VariantInterval


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
                self.strand,
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

    @property
    def cds_start(self) -> int:
        raise NoncodingTranscriptError("No CDS start for non-transcribed features")

    @property
    def cds_end(self) -> int:
        raise NoncodingTranscriptError("No CDS end for non-transcribed features")

    @property
    def chunk_relative_cds_start(self) -> int:
        raise NoncodingTranscriptError("No CDS start for non-transcribed features")

    @property
    def chunk_relative_cds_end(self) -> int:
        raise NoncodingTranscriptError("No CDS end for non-transcribed features")

    @property
    def cds_location(self) -> Location:
        """Returns the Location of the CDS in *chromosome coordinates*"""
        raise NoncodingTranscriptError("No location on a non-transcribed feature")

    @property
    def cds_chunk_relative_location(self) -> Location:
        """Returns the Location of the CDS in *chunk relative coordinates*"""
        raise NoncodingTranscriptError("No location on a non-transcribed feature")

    @property
    def is_coding(self) -> bool:
        raise NoncodingTranscriptError("Non-transcribed features cannot be coding")

    @property
    def has_in_frame_stop(self) -> bool:
        raise NoncodingTranscriptError("Cannot have frameshifts on non-transcribed features")

    @property
    def cds_size(self) -> int:
        """CDS size, regardless of chunk relativity (does not shrink)"""
        raise NoncodingTranscriptError("No cds size on a non-transcribed feature")

    @property
    def chunk_relative_cds_size(self) -> int:
        """Chunk relative CDS size (can shrink if the Location is a slice of the full transcript)"""
        raise NoncodingTranscriptError("No chunk-relative CDS size on a non-transcribed feature")

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
        is_primary_feature: Optional[bool] = None,
    ) -> "FeatureInterval":
        if location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot call from_location with a chunk-relative location. Use from_chunk_relative_location()."
            )

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

    @staticmethod
    def from_chunk_relative_location(
        location: Location,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_name: Optional[str] = None,
        guid: Optional[UUID] = None,
        feature_guid: Optional[UUID] = None,
        feature_types: Optional[List[str]] = None,
        feature_id: Optional[str] = None,
        feature_name: Optional[str] = None,
        is_primary_feature: Optional[bool] = None,
    ) -> "FeatureInterval":
        """
        Allows construction of a FeatureInterval from a chunk-relative location. This is a location
        present on a sequence chunk, which should be built by the convenience function seq_chunk_to_parent:

        .. code-block:: python

            from biocantor.io.parser import seq_chunk_to_parent
            parent = seq_chunk_to_parent('AANAAATGGCGAGCACCTAACCCCCNCC', "NC_000913.3", 222213, 222241)
            loc = SingleInterval(5, 20, Strand.PLUS, parent=parent)

        And then, this can be lifted back to chromosomal coordinates like such:

        .. code-block:: python

            loc.lift_over_to_first_ancestor_of_type("chromosome")

        """
        if not location.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException("Must have a sequence chunk in the parent hierarchy.")

        chromosome_location = location.lift_over_to_first_ancestor_of_type("chromosome")

        return FeatureInterval(
            interval_starts=[x.start for x in chromosome_location.blocks],
            interval_ends=[x.end for x in chromosome_location.blocks],
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

    def _to_gff_or_gtf(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
        row_type: Union[Type[GFFRow], Type[GTFRow]] = GFFRow,
        attribute_type: Union[Type[GFFAttributes], Type[GTFAttributes]] = GFFAttributes,
    ) -> Iterator[Union[GFFRow, GTFRow]]:
        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers(parent_qualifiers)

        feature_id = str(self.guid)

        attributes = attribute_type(
            id=feature_id,
            qualifiers=qualifiers,
            name=self.feature_name,
            parent=parent,
            raise_on_reserved_attributes=raise_on_reserved_attributes,
        )

        # "transcript" (feature interval) feature
        row = row_type(
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
            attributes = attribute_type(
                id=f"feature-{feature_id}-{i}",
                qualifiers=qualifiers,
                name=self.feature_name,
                parent=feature_id,
                raise_on_reserved_attributes=raise_on_reserved_attributes,
            )
            row = row_type(
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

    def to_gff(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
    ) -> Iterator[GFFRow]:
        """Writes a GFF format list of lists for this feature.

        The additional qualifiers are used when writing a hierarchical relationship back to files. GFF files
        are easier to work with if the children features have the qualifiers of their parents.

        Args:
            parent: ID of the Parent of this transcript.
            parent_qualifiers: Directly pull qualifiers in from this dictionary.
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.
            raise_on_reserved_attributes: If ``True``, then GFF3 reserved attributes such as ``ID`` and ``Name`` present
                in the qualifiers will lead to an exception and not a warning.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
            GFF3MissingSequenceNameError: If there are no sequence names associated with this feature.
        """

        yield from self._to_gff_or_gtf(
            parent,
            parent_qualifiers,
            chromosome_relative_coordinates,
            raise_on_reserved_attributes,
            GFFRow,
            GFFAttributes,
        )

    def to_gtf(
        self,
        parent: Optional[str] = None,
        parent_qualifiers: Optional[Dict[Hashable, Set[str]]] = None,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GFFRow]:
        raise NotImplementedError("Cannot export features to GTF")

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

    def incorporate_variants(
        self, variants: Union["VariantInterval", "VariantIntervalCollection"]
    ) -> "FeatureInterval":
        """
        Incorporate all of the variant(s) for an input VariantInterval or VariantIntervalCollection,
        producing a new FeatureInterval with those changes incorporated.
        """
        new_loc = variants.lift_over_location(self.chunk_relative_location)
        if new_loc.is_empty:
            raise EmptyLocationException("Variant incorporation led to an EmptyLocation")
        fn = FeatureInterval.from_chunk_relative_location if self.is_chunk_relative else FeatureInterval.from_location
        return fn(
            new_loc,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            guid=None,
            feature_guid=self.feature_guid,
            feature_types=sorted(self.feature_types) if self.feature_types else None,
            feature_id=self.feature_id,
            feature_name=self.feature_name,
            is_primary_feature=self.is_primary_feature,
        )


class FeatureIntervalCollection(AbstractFeatureIntervalCollection):
    """A FeatureIntervalCollection is arbitrary container of intervals.

    This can be thought of to be analogous to a :class:`GeneInterval`, but for non-transcribed
    features that are grouped in some fashion. An example is transcription factor
    binding sites for a specific transcription factor.

    The :class:`~biocantor.location.strand.Strand` of this feature interval collection is always the *plus* strand.

    This cannot be empty; it must have at least one feature interval.
    """

    interval_type = IntervalType.FEATURE
    _identifiers = ["feature_collection_id", "feature_collection_name", "locus_tag"]

    def __init__(
        self,
        feature_intervals: List[FeatureInterval],
        feature_collection_name: Optional[str] = None,
        feature_collection_id: Optional[str] = None,
        feature_collection_type: Optional[str] = None,
        locus_tag: Optional[str] = None,
        sequence_name: Optional[str] = None,
        sequence_guid: Optional[UUID] = None,
        guid: Optional[UUID] = None,
        qualifiers: Optional[Dict[Hashable, List[QualifierValue]]] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        if not feature_intervals:
            raise InvalidAnnotationError("Must have at least one feature interval.")

        self.feature_intervals = feature_intervals
        self.feature_collection_name = feature_collection_name
        self.feature_collection_id = feature_collection_id
        self.locus_tag = locus_tag
        self.sequence_name = sequence_name
        self.sequence_guid = sequence_guid
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)
        self.feature_collection_type = feature_collection_type
        self.feature_types = set.union(*[x.feature_types for x in feature_intervals])

        if not self.feature_intervals:
            raise InvalidAnnotationError("FeatureCollection must have features")

        # start/end are assumed to be in genomic coordinates, and then _initialize_location
        # will transform them into chunk-relative coordinates if necessary
        self.start = self.genomic_start = min(f.start for f in self.feature_intervals)
        self.end = self.genomic_end = max(f.end for f in self.feature_intervals)
        self._initialize_location(self.start, self.end, parent_or_seq_chunk_parent)
        self.bin = bins(self.start, self.end, fmt="bed")

        self.primary_feature = AbstractFeatureIntervalCollection._find_primary_feature(self.feature_intervals)

        if guid is None:
            self.guid = digest_object(
                self.chunk_relative_location,
                self.feature_collection_name,
                self.feature_collection_id,
                self.feature_collection_type,
                self.feature_types,
                self.locus_tag,
                self.sequence_name,
                self.qualifiers,
                self.children_guids,
            )
        else:
            self.guid = guid

        self.guid_map = {}
        for feat in self.feature_intervals:
            if feat.guid in self.guid_map:
                raise DuplicateFeatureError(f"Guid {feat.guid} found more than once in this FeatureIntervalCollection")
            self.guid_map[feat.guid] = feat

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(identifiers={self.identifiers}, "
            f"Intervals:{','.join(str(f) for f in self.feature_intervals)})"
        )

    def iter_children(self) -> Iterable[FeatureInterval]:
        """Iterate over all intervals in this collection."""
        yield from self.feature_intervals

    @property
    def is_coding(self) -> bool:
        """Never coding."""
        return False

    @property
    def children_guids(self) -> Set[UUID]:
        return {x.guid for x in self.feature_intervals}

    @property
    def id(self) -> str:
        """Returns the ID of this feature collection. Provides a shared API across genes/transcripts and features."""
        return self.feature_collection_id

    @property
    def name(self) -> str:
        """Returns the name of this feature collection. Provides a shared API across genes/transcripts and features."""
        return self.feature_collection_name

    def get_primary_feature(self) -> Union[FeatureInterval, None]:
        """Get the primary feature, if it exists."""

        return self.primary_feature

    def get_primary_feature_sequence(self) -> Union[Sequence, None]:
        """Convenience function that provides shared API between features and transcripts"""
        if self.get_primary_feature() is not None:
            return self.get_primary_feature().get_spliced_sequence()

    def get_merged_feature(self) -> FeatureInterval:
        """Generate a single :class:`~biocantor.gene.feature.FeatureInterval` that merges all intervals together."""
        intervals = []
        for tx in self.feature_intervals:
            for i in tx.chromosome_location.blocks:
                intervals.append(i)
        merged = reduce(lambda x, y: x.union(y), intervals)
        interval_starts = [x.start for x in merged.blocks]
        interval_ends = [x.end for x in merged.blocks]

        return FeatureInterval(
            interval_starts=interval_starts,
            interval_ends=interval_ends,
            strand=self.chunk_relative_location.strand,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            feature_types=list(self.feature_types),
            feature_name=self.feature_collection_name,
            feature_id=self.feature_collection_id,
            guid=self.guid,
            parent_or_seq_chunk_parent=self.chunk_relative_location.parent,
        )

    def to_dict(self, chromosome_relative_coordinates: bool = True) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`~biocantor.io.models.FeatureIntervalCollectionModel`."""
        return dict(
            feature_intervals=[feat.to_dict(chromosome_relative_coordinates) for feat in self.feature_intervals],
            feature_collection_name=self.feature_collection_name,
            feature_collection_id=self.feature_collection_id,
            feature_collection_type=self.feature_collection_type,
            locus_tag=self.locus_tag,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            feature_collection_guid=self.guid,
        )

    @staticmethod
    def from_dict(
        vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None
    ) -> "FeatureIntervalCollection":
        """Build a :class:`FeatureIntervalCollection` from a dictionary representation"""
        return FeatureIntervalCollection(
            feature_intervals=[
                FeatureInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["feature_intervals"]
            ],
            feature_collection_name=vals["feature_collection_name"],
            feature_collection_id=vals["feature_collection_id"],
            feature_collection_type=vals["feature_collection_type"],
            locus_tag=vals["locus_tag"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            guid=vals["feature_collection_guid"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def export_qualifiers(self) -> Dict[Hashable, Set[str]]:
        """Exports qualifiers for GFF3/GenBank export"""
        qualifiers = self.qualifiers.copy()
        for key, val in [
            [BioCantorQualifiers.FEATURE_COLLECTION_ID.value, self.feature_collection_id],
            [BioCantorQualifiers.FEATURE_COLLECTION_NAME.value, self.feature_collection_name],
            [BioCantorQualifiers.LOCUS_TAG.value, self.locus_tag],
            [BioCantorQualifiers.FEATURE_COLLETION_TYPE.value, self.feature_collection_type],
        ]:
            if not val:
                continue
            if key not in qualifiers:
                qualifiers[key] = set()
            qualifiers[key].add(val)
        if self.feature_types:
            qualifiers[BioCantorQualifiers.FEATURE_TYPE.value] = self.feature_types
        return qualifiers

    def query_by_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> Optional["FeatureIntervalCollection"]:
        """Filter this feature collection object by a list of unique IDs.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`FeatureIntervalCollection`, or None if there are no matching GUIDs.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        features = [self.guid_map[i] for i in ids if i in self.guid_map]
        if features:
            return FeatureIntervalCollection(
                feature_intervals=features,
                feature_collection_name=self.feature_collection_name,
                feature_collection_id=self.feature_collection_id,
                feature_collection_type=self.feature_collection_type,
                locus_tag=self.locus_tag,
                qualifiers=self._export_qualifiers_to_list(),
                sequence_name=self.sequence_name,
                sequence_guid=self.sequence_guid,
                guid=self.guid,
                parent_or_seq_chunk_parent=self.chunk_relative_location.parent,
            )

    def _to_gff_or_gtf(
        self,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
        row_type: Union[Type[GFFRow], Type[GTFRow]] = GFFRow,
        attribute_type: Union[Type[GFFAttributes], Type[GTFAttributes]] = GFFAttributes,
    ) -> Iterator[Union[GFFRow, GTFRow]]:
        if not self.sequence_name:
            raise GFF3MissingSequenceNameError("Must have sequence names to export to GFF3.")

        if not chromosome_relative_coordinates and not self.has_ancestor_of_type(SequenceType.SEQUENCE_CHUNK):
            raise NoSuchAncestorException(
                "Cannot export GFF in relative coordinates without a sequence_chunk ancestor."
            )

        qualifiers = self.export_qualifiers()

        feat_group_id = str(self.guid)

        attributes = GFFAttributes(
            id=feat_group_id,
            qualifiers=qualifiers,
            name=self.feature_collection_name,
            parent=None,
            raise_on_reserved_attributes=raise_on_reserved_attributes,
        )

        row = GFFRow(
            self.sequence_name,
            GFF_SOURCE,
            BioCantorFeatureTypes.FEATURE_COLLECTION,
            (self.start if chromosome_relative_coordinates else self.chunk_relative_start) + 1,
            self.end if chromosome_relative_coordinates else self.chunk_relative_end,
            NULL_COLUMN,
            self.chunk_relative_location.strand,
            CDSPhase.NONE,
            attributes,
        )
        yield row

        for feature in self.feature_intervals:
            if row_type == GFFRow:
                yield from feature.to_gff(
                    feat_group_id,
                    qualifiers,
                    chromosome_relative_coordinates=chromosome_relative_coordinates,
                    raise_on_reserved_attributes=raise_on_reserved_attributes,
                )
            else:
                yield from feature.to_gtf(
                    feat_group_id,
                    qualifiers,
                    chromosome_relative_coordinates=chromosome_relative_coordinates,
                )

    def to_gff(
        self,
        chromosome_relative_coordinates: bool = True,
        raise_on_reserved_attributes: Optional[bool] = True,
    ) -> Iterator[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this feature collection and its
        children.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.
            raise_on_reserved_attributes: If ``True``, then GFF3 reserved attributes such as ``ID`` and ``Name`` present
                in the qualifiers will lead to an exception and not a warning.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        yield from self._to_gff_or_gtf(
            chromosome_relative_coordinates,
            raise_on_reserved_attributes,
            GFFRow,
            GFFAttributes,
        )

    def to_gtf(
        self,
        chromosome_relative_coordinates: bool = True,
    ) -> Iterator[GTFRow]:
        raise NotImplementedError("Cannot export features to GTF")

    def incorporate_variants(
        self, variants: Union["VariantInterval", "VariantIntervalCollection"]
    ) -> "FeatureIntervalCollection":
        """
        Incorporate all of the variant(s) for an input VariantInterval or VariantIntervalCollection,
        producing a new FeatureIntervalCollection with those changes incorporated on every child.
        """
        new_features = [feature.incorporate_variants(variants) for feature in self.feature_intervals]
        if variants.has_sequence:
            new_parent = variants.parent_with_alternative_sequence
        else:
            new_parent = variants.chunk_relative_location.parent
        return FeatureIntervalCollection(
            new_features,
            feature_collection_name=self.feature_collection_name,
            feature_collection_id=self.feature_collection_id,
            feature_collection_type=self.feature_collection_type,
            locus_tag=self.locus_tag,
            sequence_guid=self.sequence_guid,
            sequence_name=self.sequence_name,
            guid=None,  # generate a new Interval GUID based on updated data
            qualifiers=self._export_qualifiers_to_list(),
            parent_or_seq_chunk_parent=new_parent,
        )
