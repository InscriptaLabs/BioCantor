"""
Collection classes. The data model is structured into two general categories,
transcripts and features. Each of those are wrapped into genes and feature collections,
respectively. These are then wrapped up into one :class:`AnnotationIntervalCollection`.

:class:`AnnotationIntervalCollections` are the topmost class and hold all possible annotations
for a given interval, as well as the place to find their sequence information.

It is useful to think of transcripts/genes as *transcriptional units*, which mean these data structures
model *transcribed sequence*. In contrast, features are *non-transcribed*, and are meant to model things
such as promoters or transcription factor binding sites.

Each object is capable of exporting itself to BED and GFF3.
"""
import itertools
from typing import List, Any, Dict, Set, Hashable, Optional, Union, Iterator, Tuple
from uuid import UUID

from methodtools import lru_cache

from biocantor.exc import (
    InvalidAnnotationError,
    InvalidQueryError,
)
from biocantor.gene.feature import FeatureIntervalCollection, FeatureInterval
from biocantor.gene.gene import GeneInterval
from biocantor.gene.interval import QualifierValue, IntervalType, AbstractFeatureIntervalCollection
from biocantor.gene.transcript import TranscriptInterval
from biocantor.gene.variants import VariantIntervalCollection, VariantInterval
from biocantor.io.gff3.rows import GFFRow, GTFRow
from biocantor.location import SingleInterval, EmptyLocation, Strand
from biocantor.parent import Parent, SequenceType
from biocantor.sequence import Alphabet
from biocantor.util.bins import bins
from biocantor.util.hashing import digest_object

try:
    import cgranges
except ImportError:
    HAS_CGRANGES = False
else:
    HAS_CGRANGES = True


class AnnotationCollection(AbstractFeatureIntervalCollection):
    """An AnnotationCollection is a container to contain :class:`GeneInterval`,
    :class:`FeatureIntervalCollection` and :class:`VariantIntervalCollection`.

    Encapsulates all possible annotations for a given interval on a specific source.

    If no start/end points are provided, the interval for this collection is the min/max of the data it contains. The
    interval for an AnnotationCollection is always on the plus strand.

    An AnnotationCollection can be empty (``feature_collections``, ``genes``, and ``variant_collections``
    can be ``None``).

    The object provided to ``parent_or_seq_chunk_parent`` must have a ``chromosome`` sequence-type in its ancestry,
    and there must be associated sequence. This object should look like the object produced by the function
    :meth:`biocantor.io.parser.seq_to_parent()`, and represent a full chromosome sequence. This will be automatically
    instantiated if you use the constructor method in :class:`biocantor.io.parser.ParsedAnnotationRecord`,
    which will import the sequence from a BioPython ``SeqRecord`` object.

    If you are using file parsers, then if the associated file types have sequence information (GenBank or GFF3+FASTA),
    then the sequences will also be automatically included when the :class:`~biocantor.io.parser.ParsedAnnotationRecord`
    is returned.

    *Object Bounds*: If `start` is provided, `end` must be provided, and vice versa. If neither are provided, and a
    `parent_or_seq_chunk_parent` is provided, then the bounds of this collection will be inferred from that object,
    if possible. If not possible, the bounds of the collection will be the bounds of the child objects associated.

    It is possible to instantiate a :class:`AnnotationCollection` with a ``sequence_chunk`` as well. A
    ``sequence_chunk`` is a slice of a chromosomal sequence that allows operations without loading an entire chromosome
    into memory. The easiest way to produce the parental relationship required for this object to operate on
    ``sequence_chunk`` is to instantiate via the constructor :meth:`biocantor.io.parser.seq_chunk_to_parent()`,
    to which you provide the slice of sequence, the chromosomal start/end positions of that slice, and a sequence name,
    and the returned Parent object will be suitable for passing to this class.
    """

    _identifiers = ["name"]

    def __init__(
        self,
        feature_collections: Optional[List[FeatureIntervalCollection]] = None,
        genes: Optional[List[GeneInterval]] = None,
        variant_collections: Optional[List[VariantIntervalCollection]] = None,
        name: Optional[str] = None,
        id: Optional[str] = None,
        sequence_name: Optional[str] = None,
        sequence_guid: Optional[UUID] = None,
        sequence_path: Optional[str] = None,
        qualifiers: Optional[Dict[Hashable, QualifierValue]] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        completely_within: Optional[bool] = None,
        parent_or_seq_chunk_parent: Optional[Parent] = None,
    ):
        self.feature_collections = feature_collections if feature_collections else []
        self.genes = genes if genes else []
        self.variant_collections = variant_collections if variant_collections else []
        self.sequence_name = sequence_name
        self.sequence_guid = sequence_guid
        self.sequence_path = sequence_path
        self._name = name
        self._id = id
        self._parent_or_seq_chunk_parent = parent_or_seq_chunk_parent
        # qualifiers come in as a List, convert to Set
        self._import_qualifiers_from_list(qualifiers)

        # we store the sequence explicitly, because this is how we can retain sequence information
        # for empty collections
        if parent_or_seq_chunk_parent and parent_or_seq_chunk_parent.sequence:
            self.sequence = parent_or_seq_chunk_parent.sequence
        else:
            self.sequence = None

        if start is None and end is not None:
            raise InvalidAnnotationError("If end is provided, start must also be provided.")
        elif end is None and start is not None:
            raise InvalidAnnotationError("If start is provided, end must also be provided.")
        # must both be unset
        elif start is None:
            # build the start/end coordinates of this collection. Start by looking at the provided parent
            # to use the coordinates provided there.
            if parent_or_seq_chunk_parent and parent_or_seq_chunk_parent.has_ancestor_of_type(SequenceType.CHROMOSOME):
                chrom_parent = parent_or_seq_chunk_parent.first_ancestor_of_type(SequenceType.CHROMOSOME)
                if chrom_parent.location:
                    start = chrom_parent.location.start
                    end = chrom_parent.location.end

            # if we have children, and the above did not work, then use the children
            # cannot infer a range for an empty collection
            if start is None and not self.is_empty:
                start = min(f.start for f in self.iter_children())
                end = max(f.end for f in self.iter_children())

        if start is None and end is None:
            # if we still have nothing, we are empty
            self._location = EmptyLocation()
        else:
            self._initialize_location(start, end, parent_or_seq_chunk_parent)
            self.start = start
            self.end = end
            self.bin = bins(self.start, self.end, fmt="bed")

        self.completely_within = completely_within

        self.guid_map: Dict[UUID, Union[GeneInterval, FeatureIntervalCollection, VariantIntervalCollection]] = {
            x.guid: x for x in self.iter_children()
        }

        self.guid: UUID = digest_object(
            self._location, self.name, self.sequence_name, self.qualifiers, self.completely_within, self.children_guids
        )

        self._associate_intervals_with_variant_intervals()

    def __repr__(self):
        return f"{self.__class__.__name__}({','.join(str(f) for f in self.iter_children())})"

    def __len__(self):
        return len(self.feature_collections) + len(self.genes)

    def __getstate__(self):
        return self.to_dict(export_parent=True)

    def __setstate__(self, state):
        ac = self.from_dict(state)
        self.__init__(
            ac.feature_collections,
            ac.genes,
            ac.variant_collections,
            ac.name,
            ac.id,
            ac.sequence_name,
            ac.sequence_guid,
            ac.sequence_path,
            ac._export_qualifiers_to_list(),
            ac.start,
            ac.end,
            ac.completely_within,
            ac._parent_or_seq_chunk_parent,
        )

    def _associate_intervals_with_variant_intervals(self):
        """
        If the constructor for this AnnotationCollection was passed one or more VariantIntervalCollections,
        then construct a mapping that associates them together. This produces new GeneInterval/FeatureIntervalCollection
        objects whose Parent are the alternative haplotype defined by the VariantIntervalCollection.

        If a GeneInterval or FeatureIntervalCollection overlap multiple VariantIntervalCollections, then they will
        exist on sequence chunks that define the sub-fraction of the interval that the haplotype represents.

        This is different from :meth:`incorporate_variants` because it is applying the variants within this collection,
        enabling comparison of haplotypes rather than generating an entirely new AnnotationCollection centered
        on the alternative haplotypes. :meth:`incorporate_variants` can only apply one VariantIntervalCollection
        at a time to an entire Interval, whereas this will instead group them by haplotype.
        """
        if not self.variant_collections:
            self.alternative_haplotype_mapping = None
            return

        self.alternative_haplotype_mapping = {}

        if HAS_CGRANGES is False:
            for variant_collection in self.variant_collections:
                for gene_or_feature in itertools.chain(self.genes, self.feature_collections):
                    if gene_or_feature.chunk_relative_location.has_overlap(variant_collection.chunk_relative_location):
                        new_gene_or_feature = gene_or_feature.incorporate_variants(variant_collection)
                        if variant_collection.guid not in self.alternative_haplotype_mapping:
                            self.alternative_haplotype_mapping[variant_collection.guid] = []
                        self.alternative_haplotype_mapping[variant_collection.guid].append(new_gene_or_feature)
        else:
            tree = cgranges.cgranges()
            for i, variant_collection in enumerate(self.variant_collections):
                tree.add(
                    "",
                    variant_collection.genomic_start,
                    variant_collection.genomic_end,
                    i,
                )
            tree.index()
            for gene_or_feature in itertools.chain(self.genes, self.feature_collections):
                for _, __, variant_collection_idx in tree.overlap(
                    "", gene_or_feature.genomic_start, gene_or_feature.genomic_end
                ):
                    variant_collection = self.variant_collections[variant_collection_idx]
                    new_gene_or_feature = gene_or_feature.incorporate_variants(variant_collection)
                    if variant_collection.guid not in self.alternative_haplotype_mapping:
                        self.alternative_haplotype_mapping[variant_collection.guid] = []
                    self.alternative_haplotype_mapping[variant_collection.guid].append(new_gene_or_feature)

    @property
    def is_empty(self) -> bool:
        """Is this an empty collection?"""
        return len(self) == 0

    @property
    def children_guids(self) -> set:
        return {x.guid for x in self.iter_children()}

    @lru_cache(maxsize=1)
    @property
    def hierarchical_children_guids(self) -> Dict[UUID, Set[UUID]]:
        """Returns children GUIDs in their hierarchical structure."""
        retval = {}
        for child in self.iter_children():
            if child.guid in retval:
                raise InvalidAnnotationError(f"Found multiple interval collections with the same GUID: {child.guid}")
            retval[child.guid] = child.children_guids
        return retval

    @lru_cache(maxsize=1)
    @property
    def interval_guids_to_collections(self) -> Dict[UUID, Union[GeneInterval, FeatureIntervalCollection]]:
        """
        For example, if this collection had a gene with two transcripts with GUID ABC and 123, and the gene
        had GUID XYZ, this would return:

        .. code-block::

            {
              "ABC": GeneInterval(guid=XYZ),
              "123": GeneInterval(guid=XYZ)
            }


        Returns:
            A map of sub-feature GUIDs to their containing elements.
        """
        retval = {}
        for child in self.iter_children():
            for interval in child.iter_children():
                if interval.guid in retval:
                    raise InvalidAnnotationError(f"Found multiple child intervals with the same GUID: {interval.guid}")
                retval[interval.guid] = child
        return retval

    @property
    def id(self) -> str:
        """Returns the ID of this collection. Provides a shared API across genes/transcripts and features."""
        return self._id

    @property
    def name(self) -> str:
        """Returns the name of this collection. Provides a shared API across genes/transcripts and features."""
        return self._name

    @lru_cache(maxsize=1)
    @property
    def children(self) -> List[Union[GeneInterval, FeatureIntervalCollection, VariantIntervalCollection]]:
        """
        Sorted list of all children. Cached.
        """
        chain_iter = itertools.chain(self.genes, self.feature_collections, self.variant_collections)
        return sorted(chain_iter, key=lambda x: x.start)

    @lru_cache(maxsize=1)
    @property
    def non_variant_children(self) -> List[Union[GeneInterval, FeatureIntervalCollection, VariantIntervalCollection]]:
        """
        Sorted list of all non-variant children. Cached.
        """
        chain_iter = itertools.chain(self.genes, self.feature_collections)
        return sorted(chain_iter, key=lambda x: x.start)

    def iter_children(self) -> Iterator[Union[GeneInterval, FeatureIntervalCollection, VariantIntervalCollection]]:
        """Iterate over all intervals in this collection, in sorted order."""
        yield from self.children

    def iter_non_variant_children(self) -> Iterator[Union[GeneInterval, FeatureIntervalCollection]]:
        """Iterate over all intervals in this collection, in sorted order."""
        yield from self.non_variant_children

    def to_dict(self, chromosome_relative_coordinates: bool = True, export_parent: bool = False) -> Dict[str, Any]:
        """Convert to a dict usable by :class:`~biocantor.io.models.AnnotationCollectionModel`.

        Allows export of the parent object as well, which allows for sequence information to be serialized to disk.

        It is not currently possible to export the parent in chunk-relative coordinates.

        Raises:
            NotImplementedError if chromosome_relative_coordinates is ``False`` and export_parent is ``True``.
        """

        if export_parent is True:
            parent_or_seq_chunk_parent = self._parent_to_dict(chromosome_relative_coordinates)
        else:
            parent_or_seq_chunk_parent = None

        return dict(
            genes=([gene.to_dict(chromosome_relative_coordinates) for gene in self.genes]),
            feature_collections=(
                [feature.to_dict(chromosome_relative_coordinates) for feature in self.feature_collections]
            ),
            variant_collections=(
                [variant.to_dict(chromosome_relative_coordinates) for variant in self.variant_collections]
            ),
            name=self.name,
            id=self.id,
            qualifiers=self._export_qualifiers_to_list(),
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            sequence_path=self.sequence_path,
            start=self.start if chromosome_relative_coordinates else self.chunk_relative_start,
            end=self.end if chromosome_relative_coordinates else self.chunk_relative_end,
            completely_within=self.completely_within,
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    @staticmethod
    def from_dict(vals: Dict[str, Any], parent_or_seq_chunk_parent: Optional[Parent] = None) -> "AnnotationCollection":
        """Build a :class:`AnnotationCollection` from a dictionary representation.

        Will use the ``parent_or_seq_chunk_parent`` value encoded in the dict if it exists,
        but this will be overridden by anything passed to the parameter.
        """

        def extract_parent_or_seq_chunk_parent_from_parent_dict(parent_dict: Dict[str, Any]) -> Optional[Parent]:
            """
            Extract a ``parent_or_seq_chunk_parent`` from a dictionary representation. This function is called
            if no ``parent_or_seq_chunk_parent`` is provided explicitly to ``from_dict``.

            NOTE: This function modifies the contents of ``parent_dict``. It is expected that it is only called
            from ``convert_parent_dict_to_parent``, which copies ``parent_dict`` out of the input ``vals``
            dictionary.

            When the original :class:`AnnotationCollection` was exported to dictionary, it may or may not have had
            a ``parent_or_seq_chunk_parent``, and that :class:`~biocantor.parent.Parent` may or may not have had
            any sequence. The sequence may or may not have been a sequence-chunk.

            This function first checks for the presence of any sequence information, then infers if the sequence
            was a chunk or not.
            """
            if parent_dict.get("seq"):
                # have to import here to avoid circular imports
                from biocantor.io.parser import seq_chunk_to_parent, seq_to_parent

                # use dictionary to prevent seq_to_parent or seq_chunk_to_parent from retaining their default parameters
                if parent_dict.get("seq_type") and parent_dict["seq_type"] == SequenceType.SEQUENCE_CHUNK:
                    fn = seq_chunk_to_parent
                    # not an argument to seq_chunk_to_parent because it is implicit
                    del parent_dict["seq_type"]
                else:
                    # seq_to_parent uses slightly different kwargs, unfortunately
                    fn = seq_to_parent
                    parent_dict["seq_id"] = parent_dict["sequence_name"]
                    # remove incorrectly named or invalid parameters
                    parent_dict = {
                        k: v
                        for k, v in parent_dict.items()
                        if k not in ["sequence_name", "type", "start", "end", "strand"]
                    }
                parent_or_seq_chunk_parent = fn(**parent_dict)
            elif "seq_type" in parent_dict or "sequence_name" in parent_dict:
                # the previous AnnotationCollection had a sequenceless Parent
                parent_or_seq_chunk_parent = Parent(
                    id=parent_dict.get("sequence_name"), sequence_type=parent_dict.get("seq_type")
                )
            else:
                parent_or_seq_chunk_parent = None
            return parent_or_seq_chunk_parent

        def convert_parent_dict_to_parent(vals: Dict[str, Any]) -> Optional[Parent]:
            """Converts the optional dictionary representation of a ``parent_or_seq_chunk_parent`` to
            a Parent object.

            NOTE: This function copies the input dictionary, stripping out null values.
            """
            # copy the input parent dictionary, stripping out null values
            if "parent_or_seq_chunk_parent" in vals and vals["parent_or_seq_chunk_parent"] is not None:
                parent_dict = {k: v for k, v in vals["parent_or_seq_chunk_parent"].items() if v is not None}
            else:
                parent_dict = {}

            if parent_dict.get("alphabet"):
                parent_dict["alphabet"] = Alphabet[parent_dict["alphabet"]]
            if parent_dict.get("strand"):
                parent_dict["strand"] = Strand[parent_dict["strand"]]
            if parent_dict.get("type"):
                parent_dict["seq_type"] = SequenceType.sequence_type_str_to_type(parent_dict["type"])
                del parent_dict["type"]

            return extract_parent_or_seq_chunk_parent_from_parent_dict(parent_dict)

        if not parent_or_seq_chunk_parent and "parent_or_seq_chunk_parent" in vals:
            parent_or_seq_chunk_parent = convert_parent_dict_to_parent(vals)

        return AnnotationCollection(
            genes=[GeneInterval.from_dict(x, parent_or_seq_chunk_parent) for x in vals["genes"]]
            if vals["genes"]
            else None,
            feature_collections=[
                FeatureIntervalCollection.from_dict(x, parent_or_seq_chunk_parent) for x in vals["feature_collections"]
            ]
            if vals["feature_collections"]
            else None,
            variant_collections=[
                VariantIntervalCollection.from_dict(x, parent_or_seq_chunk_parent) for x in vals["variant_collections"]
            ]
            if vals["variant_collections"]
            else None,
            name=vals["name"],
            id=vals["id"],
            qualifiers=vals["qualifiers"],
            sequence_name=vals["sequence_name"],
            sequence_guid=vals["sequence_guid"],
            sequence_path=vals["sequence_path"],
            start=vals["start"],
            end=vals["end"],
            completely_within=vals["completely_within"],
            parent_or_seq_chunk_parent=parent_or_seq_chunk_parent,
        )

    def _subset_parent(self, start: int, end: int) -> Optional[Parent]:
        """
        Subset the Parent of this collection to a new interval, building a chunk parent.

        Args:
            start: Genome relative start position.
            end: Genome relative end position.

        Returns:
            A parent, or ``None`` if this location has no parent, or if start == end (empty interval).
        """
        if not self.chunk_relative_location.parent:
            return None
        # edge case for a now null interval
        elif start == end:
            return None
        # edge case -- we are not actually subsetting at all
        if start == self.start and end == self.end:
            return self.chunk_relative_location.parent

        chrom_ancestor = self.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME)

        # if this subset operation is about to walk off the edge of the chunk this collection exists on,
        # don't allow this
        if self.is_chunk_relative and start < self.chromosome_location.start:
            start = self.chromosome_location.start
        chunk_relative_start = chrom_ancestor.parent_to_relative_pos(start)

        # handle the edge case where the end is the end of the current chunk
        if end == self.end:
            chunk_relative_end = (
                self.lift_over_to_first_ancestor_of_type(SequenceType.CHROMOSOME).parent_to_relative_pos(end - 1) + 1
            )
        else:
            # if this subset operation is about to walk off the edge of the chunk this collection exists on,
            # don't allow this
            if self.is_chunk_relative and end > self.chromosome_location.end:
                end = self.chromosome_location.end - 1
            chunk_relative_end = self.lift_over_to_first_ancestor_of_type(
                SequenceType.CHROMOSOME
            ).parent_to_relative_pos(end)

        seq_subset = self.chunk_relative_location.extract_sequence()[chunk_relative_start:chunk_relative_end]

        parent_id = chrom_ancestor.parent.id
        # TODO: FIXME: handle circular imports by doing this import within the function
        from biocantor.io.parser import seq_chunk_to_parent

        return seq_chunk_to_parent(
            str(seq_subset),
            parent_id,
            start,
            end,
            self.chromosome_location.strand,
            self.chunk_relative_location.parent.sequence.alphabet,
        )

    def _build_new_collection_from_query(
        self,
        genes_to_keep: List[GeneInterval],
        features_collections_to_keep: List[FeatureIntervalCollection],
        variants_to_keep: List[VariantIntervalCollection],
        start: Optional[int],
        end: Optional[int],
        completely_within: Optional[bool],
    ) -> "AnnotationCollection":
        """Convenience function that wraps functionality to build new collections"""
        seq_chunk_parent = self._subset_parent(start, end)
        return AnnotationCollection.from_dict(
            dict(
                feature_collections=[x.to_dict() for x in features_collections_to_keep],
                genes=[x.to_dict() for x in genes_to_keep],
                variant_collections=[x.to_dict() for x in variants_to_keep],
                name=self.name,
                id=self.id,
                sequence_name=self.sequence_name,
                sequence_guid=self.sequence_guid,
                sequence_path=self.sequence_path,
                qualifiers=self._export_qualifiers_to_list(),
                start=start,
                end=end,
                completely_within=completely_within,
            ),
            parent_or_seq_chunk_parent=seq_chunk_parent,
        )

    def query_by_position(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
        coding_only: Optional[bool] = False,
        completely_within: Optional[bool] = True,
        expand_location_to_children: Optional[bool] = False,
    ) -> "AnnotationCollection":
        """Filter this annotation collection object based on positions, sequence, and boolean flags.

        In all cases, the comparisons are made without considering strand.  Intronic queries are still valid.
        In other words, a query from ``[10,20]`` would still return a transcript whose intervals were
        ``[0,9], [21, 30]``.

        The resulting :class:`AnnotationCollection` returned will have a `._location` member whose bounds
        exactly match the query. If ``expand_location_to_children`` is ``True``, then the
        child genes/feature collections will potentially extend beyond this range,
        in order to encapsulate their full length. The resulting gene/feature collections will potentially have a
        reduced set of transcripts/features, if those transcripts/features are outside the query range.
        However, if ``expand_location_to_children`` is ``False``, then the child genes/feature collections
        will have location objects that represent the exact bounds of the query, which means that they
        may be sliced down. If the sliced down coordinates are entirely intronic for any isoform, then
        this isoform will have an EmptyLocation `chunk_relative_location` member, because it is no longer
        possible to have a relationship to the location object associated with this collection.

        Here is an example (equals are exons, dashes are introns):

        .. code-block::

                          10      15      20      25      30      35      40
            Gene1: Tx1:     12============20
                   Tx2:     12======16-17=20--22==25
            Fc1:    F1:     12====15
                    F2:     12======16-17=20--22==25
                    F3:                                           35======40


        Results:

        +------------+------------+--------------------+------------------+
        | start      | end        | completely_within  | result           |
        +============+============+====================+==================+
        | 21         | 22         | True               | EmptyCollection  |
        +------------+------------+--------------------+------------------+
        | 21         | 22         | False              | Tx1,Tx2,F2       |
        +------------+------------+--------------------+------------------+
        | 28         | 35         | False              | EmptyCollection  |
        +------------+------------+--------------------+------------------+
        | 28         | 36         | False              | F3               |
        +------------+------------+--------------------+------------------+
        | 27         | 36         | False              | Tx1,F3           |
        +------------+------------+--------------------+------------------+
        | 24         | 36         | False              | Tx1,Tx2,F2,F3    |
        +------------+------------+--------------------+------------------+

        Args:
            start: Genome relative start position. If not set, will be 0.
            end: Genome relative end position. If not set, will be unbounded.
            coding_only: Filter for coding genes only?
            completely_within: Strict *query* boundaries? If ``False``, features that partially overlap
                will be included in the output. Bins optimization cannot be used, so these queries are slower.
            expand_location_to_children: Should the underlying location objects be expanded so that no
                child gene/transcripts get sliced? If this is ``False``, then the constituent objects may not
                actually represent their full lengths, although the original position information is retained.

        Returns:
           :class:`AnnotationCollection` that may be empty, and otherwise will contain new copies of every
            constituent member.

        Raises:
            InvalidQueryError: If the start/end bounds are not valid. This could be because they exceed the
            bounds of the current interval. It could also happen if ``expand_location_to_children`` is ``True``
            and the new expanded range would exceed the range of an associated sequence chunk.
        """
        # after bins were decided, we can now force start/end to min/max values
        # for exact checking
        start = self.start if start is None else start
        end = self.end if end is None else end
        if start < 0:
            raise InvalidQueryError("Start must be positive")
        elif start > end:
            raise InvalidQueryError("Start must be less than or equal to end")
        elif start < self.start:
            raise InvalidQueryError(
                f"Start {start} must be within bounds of current interval [{self.start}-{self.end})"
            )
        elif end > self.end:
            raise InvalidQueryError(f"End {end} must be within bounds of current interval [{self.start}-{self.end})")
        elif start == end:
            raise InvalidQueryError("Cannot query a 0bp interval (start must not be the same as end).")

        if HAS_CGRANGES:
            (
                genes_to_keep,
                features_collections_to_keep,
                variant_collections_to_keep,
            ) = self._optimized_query_by_position(start, end, completely_within, coding_only)
        else:
            genes_to_keep, features_collections_to_keep, variant_collections_to_keep = self._query_by_position(
                start, end, completely_within, coding_only
            )

        # if completely within is False, expand the range of seq_chunk_parent to retain the full span
        # of all child intervals. This prevents features getting cut in half.
        if expand_location_to_children is True and completely_within is False:
            for g_or_fc in itertools.chain(features_collections_to_keep, genes_to_keep):
                if g_or_fc.start < start:
                    start = g_or_fc.start
                if g_or_fc.end > end:
                    end = g_or_fc.end

        # if there is a sequence chunk, then some validation checks must be performed
        if self.chunk_relative_location.parent and self.chunk_relative_location.parent.sequence:
            # not possible to expand range if it exceeds parent bounds
            if start < self.start or end > self.end:
                raise InvalidQueryError(
                    f"Cannot expand range of location to {start}-{end} because the associated sequence chunk "
                    f"lies from {self.start}-{self.end}"
                )

        return self._build_new_collection_from_query(
            genes_to_keep, features_collections_to_keep, variant_collections_to_keep, start, end, completely_within
        )

    @lru_cache(maxsize=1)
    def _build_position_interval_tree(self):
        """
        Build a position tree of every child interval.
        """
        tree = cgranges.cgranges()
        for i, child in enumerate(self.children):
            tree.add("", child.genomic_start, child.genomic_end, i)
        tree.index()
        return tree

    def _optimized_query_by_position(
        self, start: int, end: int, completely_within: bool, coding_only: bool
    ) -> Tuple[List[GeneInterval], List[FeatureIntervalCollection], List[VariantIntervalCollection]]:
        """
        Optimized implementation of position query. Used when `cgranges` is installed. The tree is cached
        if this is the first time it is being built.
        """
        if not HAS_CGRANGES:
            raise RuntimeError("Cannot use this query mode without cgranges")
        tree = self._build_position_interval_tree()

        # cgranges only does .overlap() so need to restrict search when completely_within is True
        if completely_within is True:
            query_loc = SingleInterval(start, end, Strand.PLUS, parent=self.chromosome_location.parent)
            coordinate_fn = query_loc.contains

        genes_to_keep = []
        features_collections_to_keep = []
        variant_collections_to_keep = []
        for _, __, result_idx in tree.overlap("", start, end):
            child = self.children[result_idx]
            # if completely_within is true, need to restrict further
            if completely_within is True and not coordinate_fn(
                child.chromosome_location, match_strand=False, full_span=True, strict_parent_compare=True
            ):
                continue
            elif coding_only is True and child.is_coding is False:
                continue
            if child.interval_type == IntervalType.FEATURE:
                features_collections_to_keep.append(child)
            elif child.interval_type == IntervalType.TRANSCRIPT:
                genes_to_keep.append(child)
            else:
                variant_collections_to_keep.append(child)
        return genes_to_keep, features_collections_to_keep, variant_collections_to_keep

    def _query_by_position(
        self, start: int, end: int, completely_within: bool, coding_only: bool
    ) -> Tuple[List[GeneInterval], List[FeatureIntervalCollection], List[VariantIntervalCollection]]:
        """
        Non-optimized implementation of position query. Used when `cgranges` is not installed.
        """
        # bins are only valid if we have start, end and completely_within
        if completely_within and start and end:
            my_bins = bins(start, end, fmt="bed", one=False)
        else:
            my_bins = None

        # coordinate_fn will be applied when filtering specific transcripts/features
        query_loc = SingleInterval(start, end, Strand.PLUS, parent=self.chromosome_location.parent)
        if completely_within:
            coordinate_fn = query_loc.contains
        else:
            coordinate_fn = query_loc.has_overlap

        genes_to_keep = []
        features_collections_to_keep = []
        variant_collections_to_keep = []
        for child in self.iter_children():
            if coding_only and not child.is_coding:
                continue

            # my_bins only exists if completely_within, start and end
            # if no children match these bins, skip
            elif my_bins and not any(grandchild.bin in my_bins for grandchild in child):
                continue

            # regardless of completely_within flag, first just look for overlaps on the gene/feature collection level
            elif coordinate_fn(
                child.chromosome_location, match_strand=False, full_span=True, strict_parent_compare=True
            ):
                if child.interval_type == IntervalType.FEATURE:
                    features_collections_to_keep.append(child)
                elif child.interval_type == IntervalType.TRANSCRIPT:
                    genes_to_keep.append(child)
                else:
                    variant_collections_to_keep.append(child)
        return genes_to_keep, features_collections_to_keep, variant_collections_to_keep

    def _return_collection_for_id_queries(
        self,
        genes_to_keep: List[GeneInterval],
        features_collections_to_keep: List[FeatureIntervalCollection],
        variant_collections_to_keep: List[VariantIntervalCollection],
    ) -> "AnnotationCollection":
        """Convenience function shared by all functions that query by identifiers or GUIDs."""

        if genes_to_keep or features_collections_to_keep or variant_collections_to_keep:
            start = min(
                self.start,
                min(
                    x.start
                    for x in itertools.chain(genes_to_keep, features_collections_to_keep, variant_collections_to_keep)
                ),
            )
            end = max(
                self.end,
                max(
                    x.end
                    for x in itertools.chain(genes_to_keep, features_collections_to_keep, variant_collections_to_keep)
                ),
            )
        else:
            start = self.start
            end = self.end

        return self._build_new_collection_from_query(
            genes_to_keep, features_collections_to_keep, variant_collections_to_keep, start, end, self.completely_within
        )

    def query_by_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique IDs.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        genes_to_keep = []
        features_collections_to_keep = []
        variant_collections_to_keep = []
        for i in ids:
            child = self.guid_map.get(i)
            if child is None:
                continue
            elif child.interval_type == IntervalType.FEATURE:
                features_collections_to_keep.append(child)
            elif child.interval_type == IntervalType.TRANSCRIPT:
                genes_to_keep.append(child)
            else:
                variant_collections_to_keep.append(child)

        return self._return_collection_for_id_queries(
            genes_to_keep, features_collections_to_keep, variant_collections_to_keep
        )

    @lru_cache(maxsize=1)
    @property
    def _child_interval_guid_map(
        self,
    ) -> Dict[
        UUID,
        Tuple[
            Union[GeneInterval, FeatureIntervalCollection, VariantIntervalCollection],
            Union[TranscriptInterval, FeatureInterval, VariantInterval],
        ],
    ]:
        """
        Construct a dictionary mapping grandchildren (interval GUIDs) to the children themselves.
        """
        guid_map = {}
        for child in self.iter_children():
            for grandchild in child.iter_children():
                guid_map[grandchild.guid] = (child, grandchild)
        return guid_map

    def query_by_interval_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique *interval* IDs.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        gene_guids_to_keep = set()
        features_collection_guids_to_keep = set()
        variant_collection_guids_to_keep = set()
        for guid in ids:
            if guid not in self._child_interval_guid_map:
                continue
            child, _ = self._child_interval_guid_map[guid]
            if child.interval_type == IntervalType.FEATURE:
                features_collection_guids_to_keep.add(child.guid)
            elif child.interval_type == IntervalType.TRANSCRIPT:
                gene_guids_to_keep.add(child.guid)
            else:
                variant_collection_guids_to_keep.add(child.guid)

        genes_to_keep = [self.guid_map[x].query_by_guids(ids) for x in gene_guids_to_keep]
        features_collections_to_keep = [self.guid_map[x].query_by_guids(ids) for x in features_collection_guids_to_keep]
        variant_collections_to_keep = [self.guid_map[x].query_by_guids(ids) for x in variant_collection_guids_to_keep]
        return self._return_collection_for_id_queries(
            genes_to_keep, features_collections_to_keep, variant_collections_to_keep
        )

    def query_by_transcript_interval_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique *TranscriptInterval* IDs.

        This function wraps the ``query_by_guid`` function of child GeneInterval objects.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        gene_guids_to_keep = set()
        for guid in ids:
            if guid not in self._child_interval_guid_map:
                continue
            child, _ = self._child_interval_guid_map[guid]
            if child.interval_type == IntervalType.TRANSCRIPT:
                gene_guids_to_keep.add(child.guid)

        genes_to_keep = [self.guid_map[x].query_by_guids(ids) for x in gene_guids_to_keep]
        return self._return_collection_for_id_queries(genes_to_keep, [], [])

    def query_by_feature_interval_guids(self, id_or_ids: Union[UUID, List[UUID]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of unique *interval* IDs.

        This function wraps the ``query_by_guid`` function of child FeatureIntervalCollection objects.

        NOTE: If the children of this collection have GUID collisions, either across genes or features or
        within genes and features, this function will return all members with the matching GUID.

        Args:
            id_or_ids: List of GUIDs, or unique IDs. Can also be a single ID.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, UUID):
            ids = [id_or_ids]
        else:
            ids = id_or_ids

        features_collection_guids_to_keep = set()
        for guid in ids:
            if guid not in self._child_interval_guid_map:
                continue
            child, _ = self._child_interval_guid_map[guid]
            if child.interval_type == IntervalType.FEATURE:
                features_collection_guids_to_keep.add(child.guid)

        features_collections_to_keep = [self.guid_map[x].query_by_guids(ids) for x in features_collection_guids_to_keep]
        return self._return_collection_for_id_queries([], features_collections_to_keep, [])

    def query_by_feature_identifiers(self, id_or_ids: Union[str, List[str]]) -> "AnnotationCollection":
        """Filter this annotation collection object by a list of identifiers, or a single identifier.

        Identifiers are not necessarily unique; if your identifier matches more than one interval,
        all matching intervals will be returned. These ambiguous results will be adjacent in the resulting collection,
        but are not grouped or signified in any way.

        This method is ``O(n_ids * m_identifiers)``.

        Args:
            id_or_ids: List of identifiers, or a single identifier.

        Returns:
           :class:`AnnotationCollection` that may be empty.
        """
        if isinstance(id_or_ids, str):
            ids = {id_or_ids}
        else:
            ids = set(id_or_ids)

        genes_to_keep = []
        features_collections_to_keep = []
        variant_collections_to_keep = []
        for child in self.iter_children():
            if ids & child.identifiers:
                if child.interval_type == IntervalType.FEATURE:
                    features_collections_to_keep.append(child)
                elif child.interval_type == IntervalType.TRANSCRIPT:
                    genes_to_keep.append(child)
                else:
                    variant_collections_to_keep.append(child)

        return self._return_collection_for_id_queries(
            genes_to_keep, features_collections_to_keep, variant_collections_to_keep
        )

    def get_children_by_type(
        self, child_type: str
    ) -> Union[List[GeneInterval], List[FeatureIntervalCollection], List[VariantIntervalCollection]]:
        if child_type.lower() == IntervalType.FEATURE:
            return self.feature_collections
        elif child_type.lower() == IntervalType.TRANSCRIPT:
            return self.genes
        elif child_type.lower() == IntervalType.VARIANT:
            return self.variant_collections
        else:
            raise InvalidQueryError(f"Cannot get children of type {child_type}")

    def _unsorted_gff_iter(
        self, chromosome_relative_coordinates: bool = True, raise_on_reserved_attributes: bool = True
    ) -> Iterator[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this annotation collection and its
        children.

        The positions of the genes will be ordered by genomic position, but may not be globally position sorted
        because it could be the case that children gene/features will overlap. This private function
        exists to provide an iterator to sort in the main ``to_gff()`` function.

        Args:
            chromosome_relative_coordinates: Output GFF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.
            raise_on_reserved_attributes: If ``True``, then GFF3 reserved attributes such as ``ID`` and ``Name`` present
                in the qualifiers will lead to an exception and not a warning.

        Yields:
            :class:`~biocantor.io.gff3.rows.GFFRow`
        """
        for item in self.iter_children():
            yield from item.to_gff(
                chromosome_relative_coordinates=chromosome_relative_coordinates,
                raise_on_reserved_attributes=raise_on_reserved_attributes,
            )

    def to_gff(
        self, chromosome_relative_coordinates: bool = True, raise_on_reserved_attributes: Optional[bool] = True
    ) -> Iterator[GFFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GFFRow` for this annotation collection and its
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
        yield from sorted(
            self._unsorted_gff_iter(chromosome_relative_coordinates, raise_on_reserved_attributes),
            key=lambda x: x.start,
        )

    def _unsorted_gtf_iter(self, chromosome_relative_coordinates: bool = True) -> Iterator[GTFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GTFRow` for this annotation collection and its
        children.

        The positions of the genes will be ordered by genomic position, but may not be globally position sorted
        because it could be the case that children gene/features will overlap. This private function
        exists to provide an iterator to sort in the main ``to_gtf()`` function.

        Args:
            chromosome_relative_coordinates: Output GTF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.
            raise_on_reserved_attributes: If ``True``, then GFF3 reserved attributes such as ``ID`` and ``Name`` present
                in the qualifiers will lead to an exception and not a warning.

        Yields:
            :class:`~biocantor.io.gff3.rows.GTFRow`
        """
        for item in self.iter_children():
            yield from item.to_gtf(
                chromosome_relative_coordinates=chromosome_relative_coordinates,
            )

    def to_gtf(self, chromosome_relative_coordinates: bool = True) -> Iterator[GTFRow]:
        """Produces iterable of :class:`~biocantor.io.gff3.rows.GTFRow` for this annotation collection and its
        children.

        Args:
            chromosome_relative_coordinates: Output GTF in chromosome-relative coordinates? Will raise an exception
                if there is not a ``sequence_chunk`` ancestor type.

        Yields:
            :class:`~biocantor.io.gff3.rows.GTFRow`

        Raises:
            NoSuchAncestorException: If ``chromosome_relative_coordinates`` is ``False`` but there is no
            ``sequence_chunk`` ancestor type.
        """
        yield from sorted(
            self._unsorted_gtf_iter(chromosome_relative_coordinates),
            key=lambda x: x.start,
        )

    def incorporate_variants(
        self, variants: Union[VariantInterval, VariantIntervalCollection]
    ) -> "AnnotationCollection":
        """
        Incorporate all of the variant(s) for an input VariantInterval or VariantIntervalCollection,
        producing a new AnnotationCollection with those changes incorporated on every child.
        """
        new_genes = [tx.incorporate_variants(variants) for tx in self.genes]
        new_features = [feature.incorporate_variants(variants) for feature in self.feature_collections]
        if variants.has_sequence:
            new_parent = variants.parent_with_alternative_sequence
        else:
            new_parent = variants.chunk_relative_location.parent
        return AnnotationCollection(
            new_features,
            new_genes,
            variant_collections=None,
            name=self.name,
            id=self.id,
            sequence_name=self.sequence_name,
            sequence_guid=self.sequence_guid,
            sequence_path=self.sequence_path,
            qualifiers=self._export_qualifiers_to_list(),
            start=None,  # recalculate start
            end=None,  # recalculate end
            completely_within=None,  # no longer can be sure
            parent_or_seq_chunk_parent=new_parent,
        )
