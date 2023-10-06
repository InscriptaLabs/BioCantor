"""
Parse GenBank files. Biopython provides the core parsing functionality, but is not capable of producing a
hierarchical model. Thus this module does this, by depending on the ordering of the GenBank file.

There are two ways to infer hierarchy from GenBank files, that are not always followed.

The first (Model 1A) is sort order: so that it always goes

gene -> {mRNA, tRNA, rRNA} -> CDS (for coding genes only)

Each transcript feature can repeat. Each mRNA feature must be followed by a CDS feature.
The presence of a new gene feature is the divider between genes.

In some genomes (often Prokaryotic), there is no transcript level feature for coding genes.
That is, it goes from gene -> CDS. This is Model 1B.

The second way that a GenBank file can be grouped is via the locus_tag qualifiers. This method is the
default for this parsing module.

The generic parsing function that interprets the BioPython results to BioCantor data models is implemented in
:meth:`GeneFeature.to_gene_model()`. This function can be over-ridden to provide custom parsing implementations.
"""
import itertools
import pathlib
import warnings
from abc import ABC, abstractmethod
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from typing import Optional, TextIO, Iterator, List, Dict, Callable, Any, Union

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from biocantor.gene import Biotype, CDSInterval, CDSFrame
from biocantor.io.exc import (
    DuplicateSequenceException,
    InvalidCDSIntervalWarning,
    InvalidIntervalWarning,
    StrandViolationWarning,
    DuplicateFeatureWarning,
    DuplicateTranscriptWarning,
    InvalidInputError,
)
from biocantor.io.features import extract_feature_types, extract_feature_name_id, merge_qualifiers
from biocantor.io.genbank.constants import (
    GeneFeatures,
    TranscriptFeatures,
    GeneIntervalFeatures,
    MetadataFeatures,
    GenBankParserType,
    KnownQualifiers,
    NonCodingTranscriptFeatures,
    GENBANK_GENE_FEATURES,
)
from biocantor.io.vcf.parser import parse_vcf_file, VariantIntervalCollectionModel
from biocantor.io.genbank.exc import (
    GenBankParserError,
    EmptyGenBankError,
    GenBankLocusTagError,
    GenBankLocationException,
    GenBankNullStrandException,
    GenBankUnknownFeatureWarning,
    GenBankEmptyGeneWarning,
    GenBankDuplicateLocusTagWarning,
    UnknownGenBankFeatureWarning,
)
from biocantor.io.models import (
    GeneIntervalModel,
    AnnotationCollectionModel,
    FeatureIntervalCollectionModel,
)
from biocantor.io.parser import ParsedAnnotationRecord
from biocantor.location import (
    Location,
    Strand,
    CompoundInterval,
    SingleInterval,
    EmptyLocation,
)


class Feature(ABC):
    """Generic feature."""

    types = set()

    def __init__(self, feature: SeqFeature, record: SeqRecord):
        if feature.type not in self.types:
            raise GenBankParserError(f"Invalid feature type {feature.type}")
        if feature.location is None:
            raise GenBankLocationException(f"Feature {feature} did not have parseable coordinates.")
        if not feature.strand:
            raise GenBankNullStrandException(f"Feature {feature} is unstranded or has multiple strands.")
        self._seq_feature = feature
        self.record = record

    @property
    def type(self) -> str:
        return self._seq_feature.type

    @property
    def strand(self) -> int:
        return Strand.from_int(self._seq_feature.strand)

    @property
    def start(self) -> int:
        return int(self._seq_feature.location.start)


class FeatureIntervalGenBankCollection:
    """A collection of generic (non-transcribed) feature intervals."""

    def __init__(self, features: List[SeqFeature], record: SeqRecord):
        """
        Build a generic feature from a grouping of features found in a genbank parsing event.

        Args:
            features: One or more ``SeqFeature``s found in a GenBank file that are associated together, but which
                could not be interpreted as a gene.
            record: The ``SeqRecord`` these features were found on.
        """
        for feature in features:
            if feature.location is None:
                raise GenBankLocationException(f"Feature {feature} did not have parseable coordinates.")

        self.types = {feature.type for feature in features}
        self.record = record
        self._seq_features = features

    @property
    def start(self) -> int:
        return min(int(x.location.start) for x in self._seq_features)

    @staticmethod
    def to_feature_model(cls: "FeatureIntervalGenBankCollection") -> Dict[str, Any]:
        """Convert to a Dict representation of a :class:`biocantor.gene.collections.FeatureIntervalCollection`
        that can be used for analyses.

        This is the default function, that can be over-ridden by specific implementations.

        Looks for identifiers in the hierarchy defined by the Enum
        :class:`biocantor.io.genbank.constants.FeatureIntervalIdentifierKeys`.

        The feature collection produced will be named either the locus tag if provided, and otherwise by definition of
        the parser we have only one feature, so the first name is chosen.
        """
        features = []
        feature_names = Counter()
        feature_ids = Counter()
        locus_tag = None
        merged_qualifiers = {}
        for feature in cls._seq_features:
            interval_starts = []
            interval_ends = []
            for loc in sorted(feature.location.parts, key=lambda p: p.start):
                interval_starts.append(int(loc.start))
                interval_ends.append(int(loc.end))
            strand = Strand.from_int(feature.location.strand)

            # extract feature types, including the base type
            feature_types = {feature.type}
            extract_feature_types(feature_types, feature.qualifiers)

            # extract primary identifier
            feature_name, feature_id = extract_feature_name_id(feature.qualifiers)
            # keep track of feature names seen to identify consensus feature name for collection
            if feature_name:
                feature_names[feature_name] += 1
            if feature_id:
                feature_ids[feature_id] += 1

            # try to find a locus tag
            if KnownQualifiers.LOCUS_TAG.value in feature.qualifiers:
                locus_tag = feature.qualifiers[KnownQualifiers.LOCUS_TAG.value][0]

            feature_model = dict(
                interval_starts=interval_starts,
                interval_ends=interval_ends,
                strand=strand.name,
                qualifiers=feature.qualifiers,
                feature_id=feature_id,
                feature_name=feature_name,
                feature_types=sorted(feature_types),
                sequence_name=cls.record.id,
                is_primary_feature=False,
            )
            if feature_model in features:
                warnings.warn(DuplicateFeatureWarning(f"Feature {feature} found twice within a feature collection"))
            else:
                features.append(feature_model)

            merged_qualifiers = merge_qualifiers(merged_qualifiers, feature.qualifiers)

        if len(feature_names) > 0:
            feature_name = feature_names.most_common(1)[0][0]
        else:
            feature_name = locus_tag

        if len(feature_ids) > 0:
            feature_id = feature_ids.most_common(1)[0][0]
        else:
            feature_id = locus_tag

        feature_collection = FeatureIntervalCollectionModel.Schema().load(
            dict(
                feature_intervals=features,
                feature_collection_name=feature_name,
                feature_collection_id=feature_id,
                locus_tag=locus_tag,
                qualifiers=merged_qualifiers,
                sequence_name=cls.record.id,
            )
        )
        # construct a FeatureIntervalCollection to run validations
        feature_collection = feature_collection.to_feature_collection()
        return feature_collection.to_dict()


class GeneFeature(Feature):
    """A gene."""

    types = {x.value for x in GeneFeatures}

    def __init__(self, feature: SeqFeature, record: SeqRecord):
        super().__init__(feature, record)
        self.children = []

    def __str__(self):
        return self._seq_feature.__repr__()

    def __repr__(self):
        return "\n".join((str(x) for x in itertools.chain((self,), self.children)))

    @property
    def type(self) -> str:
        return self._seq_feature.type

    @property
    def has_children(self) -> bool:
        return len(self.children) > 0

    @staticmethod
    def from_transcript_or_cds_feature(feature: SeqFeature, seqrecord: SeqRecord) -> "GeneFeature":
        """Some GenBank files lack a gene-level feature, but have transcript-level features or CDS-level features only.

        Construct a GeneFeature from such records."""
        old_type = feature.type
        feature.type = GeneFeatures.GENE.value
        gene = GeneFeature(feature, seqrecord)
        feature.type = old_type
        gene.add_child(feature)
        return gene

    def add_child(self, feature: SeqFeature, cds_feature: Optional[SeqFeature] = None):
        """Add a new feature as a child. Infer Transcripts if this child is a CDS feature."""
        if feature.type in TranscriptFeature.types:
            self.children.append(TranscriptFeature(feature, self.record, cds_feature=cds_feature))
        elif feature.type in CDSFeature.types:
            # infer a transcript
            tx_feature = deepcopy(feature)
            if tx_feature.type == GeneIntervalFeatures.CDS.value:
                tx_feature.type = TranscriptFeatures.CODING_TRANSCRIPT.value
            # this means we have an exon as a direct child of a gene
            elif tx_feature.type == GeneIntervalFeatures.EXON.value:
                tx_feature.type = TranscriptFeatures.CODING_TRANSCRIPT.value
            else:
                tx_feature.type = TranscriptFeatures[tx_feature.type].value
            tx = TranscriptFeature(tx_feature, self.record, cds_feature=feature)
            self.children.append(tx)
        else:
            raise GenBankParserError(f"Invalid feature type {feature.type}")

    def infer_child(self):
        """If this is an isolated gene feature, then construct a child transcript feature that is a copy"""
        tx_feature = deepcopy(self._seq_feature)
        tx_feature.type = TranscriptFeatures.NONCODING_TRANSCRIPT.value
        tx = TranscriptFeature(tx_feature, self.record)
        self.children.append(tx)

    @staticmethod
    def to_gene_model(cls: "GeneFeature") -> Dict[str, Any]:
        """Convert to a Dict representation of a :class:`biocantor.gene.collections.GeneInterval`
        that can be used for analyses.

        This is the default function, that can be over-ridden by specific implementations.

        Looks for /transcript_id, /protein_id, and /gene on the transcript level, and
        looks for /gene_id, /gene, and /locus_tag on the gene level.
        """
        transcripts = []
        tx_biotypes = Counter()
        for tx in cls.children:
            exon_interval: CompoundInterval = tx.find_exon_interval()
            exon_starts = exon_interval._starts
            exon_ends = exon_interval._ends
            strand = Strand.from_int(tx.strand)

            cds_interval = tx.find_cds_interval()
            if cds_interval.is_empty:
                cds_starts = None
                cds_ends = None
                cds_frames = []
            else:
                cds_starts = cds_interval._starts
                cds_ends = cds_interval._ends
                # TODO: this will produce invalid Frames if there is a programmed frameshift or other such problem.
                # This is a limitation of genbank files. HOWEVER, this code could detect the presence of a
                # overlapping interval and try all possible frames there to see if any produce a valid translation.
                cds_frames = tx.construct_frames(cds_interval)

            if "pseudo" in tx._seq_feature.qualifiers:
                transcript_biotype = Biotype.pseudogene
            elif tx._seq_feature.type == TranscriptFeatures.CODING_TRANSCRIPT.value:
                transcript_biotype = Biotype.protein_coding
            else:
                transcript_biotype = Biotype[tx._seq_feature.type]
            tx_biotypes[transcript_biotype] += 1

            tx_model = dict(
                exon_starts=exon_starts,
                exon_ends=exon_ends,
                strand=strand.name,
                cds_starts=cds_starts,
                cds_ends=cds_ends,
                cds_frames=cds_frames,
                qualifiers=tx.merge_cds_qualifiers_to_transcript(),
                is_primary_tx=False,
                transcript_id=tx.get_qualifier_from_tx_or_cds_features(KnownQualifiers.TRANSCRIPT_ID.value),
                protein_id=tx.get_qualifier_from_tx_or_cds_features(KnownQualifiers.PROTEIN_ID.value),
                product=tx.get_qualifier_from_tx_or_cds_features(KnownQualifiers.PRODUCT.value),
                transcript_symbol=tx.get_qualifier_from_tx_or_cds_features(KnownQualifiers.GENE.value),
                transcript_type=transcript_biotype.name,
                sequence_name=tx.record.id,
            )
            if tx_model in transcripts:
                warnings.warn(DuplicateTranscriptWarning(f"Transcript {tx_model} found twice within a gene"))
            else:
                transcripts.append(tx_model)

        # pick most common transcript type; hacky
        gene_biotype = tx_biotypes.most_common(1)[0][0]
        gene = GeneIntervalModel.Schema().load(
            dict(
                transcripts=transcripts,
                gene_id=cls._seq_feature.qualifiers.get(KnownQualifiers.GENE_ID.value, [None])[0],
                gene_symbol=cls._seq_feature.qualifiers.get(KnownQualifiers.GENE.value, [None])[0],
                locus_tag=cls._seq_feature.qualifiers.get(KnownQualifiers.LOCUS_TAG.value, [None])[0],
                gene_type=gene_biotype.name,
                qualifiers=cls._seq_feature.qualifiers,
                sequence_name=cls.record.id,
            )
        )
        # construct a GeneInterval to run validations
        gene = gene.to_gene_interval()
        return gene.to_dict()


class TranscriptFeature(Feature):
    """A transcript"""

    types = {x.value for x in TranscriptFeatures}
    _exon_interval = None
    _cds_interval = None

    def __init__(self, feature: SeqFeature, record: SeqRecord, cds_feature: Optional[SeqFeature] = None):
        super().__init__(feature, record)
        if cds_feature:
            self.cds_feature = CDSFeature(cds_feature, self.record)
        else:
            self.cds_feature = None

    def __str__(self):
        as_str = f"--> {self._seq_feature.__repr__()}"
        if self.cds_feature:
            as_str += f"\n---> {self.cds_feature._seq_feature.__repr__()}"
        return as_str

    def construct_frames(self, cds_interval: Location) -> List[str]:
        """We need to build frames. Since GenBank lacks this info, do our best"""
        # make 0 based offset, if possible, otherwise assume always in frame
        frame = int(self.cds_feature._seq_feature.qualifiers.get(KnownQualifiers.CODON_START.value, [1])[0]) - 1
        frame = CDSFrame.from_int(frame)
        frames = CDSInterval.construct_frames_from_location(cds_interval, frame)
        return [x.name for x in frames]

    def get_qualifier_from_tx_or_cds_features(self, qualifier: str) -> Optional[str]:
        """Get a specific qualifier, if it exists. Look at tx first, then children"""
        if qualifier in self._seq_feature.qualifiers:
            return self._seq_feature.qualifiers[qualifier][0]
        if self.cds_feature:
            if qualifier in self.cds_feature._seq_feature.qualifiers:
                return self.cds_feature._seq_feature.qualifiers[qualifier][0]

    def find_exon_interval(self) -> CompoundInterval:
        """Finds the Location of the Exons."""
        if self._exon_interval:
            return self._exon_interval
        exon_starts = []
        exon_ends = []
        for part in sorted(self._seq_feature.location.parts, key=lambda p: p.start):
            exon_starts.append(int(part.start))
            exon_ends.append(int(part.end))
        self._exon_interval = CompoundInterval(
            exon_starts,
            exon_ends,
            Strand.from_int(self.strand),
        )
        return self._exon_interval

    def find_transcript_interval(self) -> Location:
        """Finds the Location that spans the full length of the Transcript"""
        exon_interval = self.find_exon_interval()
        return SingleInterval(exon_interval.blocks[0].start, exon_interval.blocks[-1].end, Strand.from_int(self.strand))

    def find_cds_interval(self) -> Location:
        """Finds the Location of the CDS.

        Handle edge cases from tools like Geneious where the CDS Interval exceeds the bounds of the Exons.
        """
        if self._cds_interval:
            return self._cds_interval

        # this is a non-coding transcript
        if not self.cds_feature:
            self._cds_interval = EmptyLocation()
            return self._cds_interval

        cds_starts = []
        cds_ends = []
        for part in sorted(self.cds_feature._seq_feature.location.parts, key=lambda p: p.start):
            cds_starts.append(int(part.start))
            cds_ends.append(int(part.end))
            cds_interval = CompoundInterval(
                cds_starts,
                cds_ends,
                Strand.from_int(self.strand),
            )
        cds_interval_intersection_with_exons = cds_interval.intersection(
            self.find_transcript_interval(), optimize_blocks=False
        )
        if cds_interval_intersection_with_exons != cds_interval:
            warnings.warn(
                InvalidCDSIntervalWarning(
                    f"CDS Interval {cds_interval} was sliced down to {cds_interval_intersection_with_exons} "
                    f"because that is the exon boundaries found"
                )
            )
        # convert to CompoundInterval, if SingleInterval, for consistent ability to extract position information
        if isinstance(cds_interval_intersection_with_exons, SingleInterval):
            cds_interval_intersection_with_exons = CompoundInterval.from_single_intervals(
                [cds_interval_intersection_with_exons]
            )
        self._cds_interval = cds_interval_intersection_with_exons
        return self._cds_interval

    def merge_cds_qualifiers_to_transcript(self) -> Dict[str, List[str]]:
        """
        If there were distinct transcript-level features, the qualifiers on the CDS feature will be lost
        when converting to the BioCantor data model unless those qualifiers are rolled into the qualifiers on
        the transcript feature.
        """
        qualifiers = {key: set(vals) for key, vals in self._seq_feature.qualifiers.items()}
        if self.cds_feature:
            for key, vals in self.cds_feature._seq_feature.qualifiers.items():
                if key not in qualifiers:
                    qualifiers[key] = set(vals)
                else:
                    qualifiers[key].update(vals)
        return {key: list(vals) for key, vals in qualifiers.items()}


class CDSFeature(Feature):
    """A CDS interval"""

    types = {GeneIntervalFeatures.CDS.value}

    def __str__(self):
        return f"----> {self._seq_feature.__repr__()}"


@dataclass
class GroupedGeneFeatures:
    """
    Container class for a grouping of gene-like SeqFeatures on their associated SeqRrecord.

    This class is used by implementations of the :class:`BaseGenBankParser` to store groupings of features
    that are considered to be part of a gene unit together with the `SeqRecord` they came from.

    Due to the various flavors of GenBank files out there, any of the gene, transcript or CDS features might not
    be present. The downstream usage of this class will infer the missing feature types.
    """

    seqrecord: SeqRecord
    gene_feature: Optional[SeqFeature] = None
    transcript_features: Optional[List[SeqFeature]] = None
    cds_features: Optional[List[SeqFeature]] = None


class BaseGenBankParser(ABC):
    """
    Base class for GenBank parsing.
    """

    genbank_parser_type: GenBankParserType

    def __init__(
        self,
        seq_records: List[SeqRecord],
        parsed_variants: Dict[str, List[VariantIntervalCollectionModel]],
        gene_parse_func: Callable[[GeneFeature], Dict[str, Any]],
        feature_parse_func: Callable[[FeatureIntervalGenBankCollection], Dict[str, Any]],
    ):
        self.seq_records = seq_records
        self.parsed_variants = parsed_variants
        self.gene_parse_func = gene_parse_func
        self.feature_parse_func = feature_parse_func

        self.sources: List[Optional[SeqFeature]] = [None] * len(seq_records)
        # this is the filtered set of SeqFeatures in the same order as the SeqRecords
        self.gene_filtered_features: List[List[SeqFeature]] = [list() for _ in self.seq_records]
        self.feature_features: List[List[SeqFeature]] = [list() for _ in self.seq_records]

        self.grouped_gene_features: List[List[GroupedGeneFeatures]] = [list() for _ in self.seq_records]

        self.genes: List[List[GeneFeature]] = [list() for _ in self.seq_records]
        self.feature_collections: List[List[FeatureIntervalGenBankCollection]] = [list() for _ in self.seq_records]

    @property
    def num_genes(self) -> int:
        return sum(len(x) for x in self.genes)

    @property
    def num_feature_collections(self):
        return sum(len(x) for x in self.feature_collections)

    @abstractmethod
    def parse(self) -> Iterator[ParsedAnnotationRecord]:
        """Parse features"""

    @staticmethod
    def validate_seqfeature(feature: SeqFeature) -> bool:
        """
        Perform validation checks on a SeqFeature. If there are issues with the feature,
        warnings and raised and this function returns ``False``.
        """
        if not feature.location:
            warnings.warn(
                InvalidIntervalWarning(f"Invalid location for feature {feature}. This feature will be skipped.")
            )
            return False
        elif not feature.strand:
            warnings.warn(
                StrandViolationWarning(f"Found multiple strands for feature {feature}. This feature will be skipped.")
            )
            return False
        return True

    @staticmethod
    def _construct_gene_from_feature(feature: SeqFeature, seqrecord: SeqRecord) -> Optional[GeneFeature]:
        """Convenience function for deciding which function to use when converting a feature to a gene"""
        if feature.type in GeneFeature.types:
            return GeneFeature(feature, seqrecord)
        elif feature.type in TranscriptFeature.types or feature.type in CDSFeature.types:
            return GeneFeature.from_transcript_or_cds_feature(feature, seqrecord)
        else:
            warnings.warn(
                UnknownGenBankFeatureWarning(
                    f"Feature {feature} was associated with a gene but was not of a known type and will be ignored"
                )
            )

    @staticmethod
    def _sort_features_by_position_and_type(features: List[SeqFeature]) -> List[SeqFeature]:
        """
        sort features first by position then by the expected order of features for a coding gene
        all non-coding transcripts as a result end up at the end of the group
        """
        return sorted(
            features,
            key=lambda x: (
                int(x.location.start),
                x.type != GeneFeatures.GENE.value,
                x.type != TranscriptFeatures.CODING_TRANSCRIPT.value,
                x.type != GeneIntervalFeatures.CDS.value,
            ),
        )

    @staticmethod
    def _group_sorted_features_by_type(features: List[SeqFeature]) -> Iterator[List[SeqFeature]]:
        """
        Iterator for grouping sorted features by each time a ``gene`` feature appears.

        This function identifies groups of canonical genes, as well as interspersed non-coding genes,
        while also trying to identify isolated CDS records.

        As an example, take this set of ordered features:

        .. code-block::

            gene
            mRNA
            CDS
            tRNA
            rRNA
            gene
            CDS

        This would be grouped as:

        .. code-block::

            [gene, mRNA, CDS]
            [tRNA]
            [rRNA]
            [gene, CDS]

        However, sometimes non-coding genes still have the gene feature. Take this example:

        .. code-block::

            gene
            mRNA
            CDS
            gene
            ncRNA
            gene
            CDS

        This would be grouped as:

        .. code-block::

            [gene, mRNA, CDS]
            [gene, ncRNA]
            [gene, CDS]

        Now consider this problematic ordering, where a non-coding object interrupts a coding object:

        .. code-block::

            gene
            tRNA
            CDS

        This is ambiguous -- is the gene meant to go with the tRNA, or the CDS? This would be grouped as:

        .. code-block::

            [gene, tRNA]
            [CDS]

        And there will be a new ``gene`` object inferred for the ``CDS``, regardless of if the ``gene`` was supposed
        to go with it or not.

        This is an example of why the Sorted parser is inherently more fragile than the LocusTag parser.
        """
        group = []
        for feature in features:
            # new gene group
            if not group:
                if feature.type in GeneFeature.types:
                    group.append(feature)
                    # base case for starting with a isolated ncRNA or CDS feature; immediately add them
                    # and reset the gene to None
                elif feature.type in TranscriptFeature.types or feature.type in CDSFeature.types:
                    group.append(feature)
                    yield group
                    group = []
                else:
                    warnings.warn(
                        UnknownGenBankFeatureWarning(
                            f"Feature {feature} was associated with a gene "
                            f"but was not of a known type and will be ignored"
                        )
                    )
            # next gene; reset
            elif feature.type in GeneFeature.types:
                yield group
                group = [feature]
            elif feature.type in TranscriptFeature.types:
                # if the current gene is non-empty, and the feature is not a mRNA, then this is a isolated ncRNA
                # finish this gene and start a new one
                if group and feature.type != TranscriptFeatures.CODING_TRANSCRIPT:
                    if group[-1].type in GeneFeature.types:
                        group.append(feature)
                    else:
                        yield group
                        group = [feature]
                else:
                    group.append(feature)
            elif feature.type in CDSFeature.types:
                # are we about to associate this CDS with a non-coding transcript type? If so, don't do so
                if any(NonCodingTranscriptFeatures.has_value(x.type) for x in group):
                    yield group
                    group = [feature]
                else:
                    group.append(feature)
            else:
                warnings.warn(
                    UnknownGenBankFeatureWarning(
                        f"Feature {feature} was associated with a gene but was not of a known type and will be ignored"
                    )
                )
        if group:
            yield group

    def _group_features_by_position(self, features: List[SeqFeature], seqrecord: SeqRecord, idx: int):
        """
        Group features by position. Since this function is always called in `seqrecord` order,
        `self.grouped_gene_features` is incremented here.
        """
        grouped_features = []
        for feature_group in self._group_sorted_features_by_type(features):
            gene_features = []
            transcript_features = []
            cds_features = []
            for feature in feature_group:
                if feature.type in GeneFeature.types:
                    gene_features.append(feature)
                elif feature.type in TranscriptFeature.types:
                    transcript_features.append(feature)
                elif feature.type in CDSFeature.types:
                    cds_features.append(feature)
                else:
                    warnings.warn(
                        UnknownGenBankFeatureWarning(
                            f"Feature {feature} was associated with a gene but was not of a known type "
                            f"and will be ignored"
                        )
                    )
            # due to the grouping predicate in _group_position_sorted_by_gene_feature it is not possible
            # for there to be more than one gene feature
            if gene_features:
                gene_feature = gene_features[0]
            else:
                gene_feature = None
            if len(transcript_features) > 1 and len(cds_features) > 1:
                warnings.warn(
                    GenBankUnknownFeatureWarning(
                        f"Grouping by position found multiple transcript features associated with the same gene\n"
                        f"{gene_feature}: {transcript_features}\n"
                        f"Extra transcript will be ignored"
                    )
                )
                transcript_features = [transcript_features[0]]
            grouped_features.append(GroupedGeneFeatures(seqrecord, gene_feature, transcript_features, cds_features))
        self.grouped_gene_features[idx].extend(grouped_features)

    def _group_features_by_locus_tag(self, features: List[SeqFeature], seqrecord: SeqRecord, idx: int):
        """
        Group features by locus tag. Since this function is always called in `seqrecord` order,
        `self.grouped_gene_features` is incremented here.
        """
        grouped_features = []
        for locus_tag, gene_features in itertools.groupby(
            features,
            key=lambda f: f.qualifiers[KnownQualifiers.LOCUS_TAG.value][0],
        ):
            gene_feature = None
            transcript_features = []
            cds_features = []
            for feature in gene_features:
                if feature.type == GeneFeatures.GENE.value:
                    if gene_feature is not None:
                        raise GenBankLocusTagError(
                            f"Grouping by locus tag found multiple gene features on the same sequence "
                            f"with the same locus tag:"
                            f"\n{gene_feature}\n{feature}"
                        )
                    else:
                        gene_feature = feature
                elif feature.type in TranscriptFeature.types:
                    transcript_features.append(feature)
                elif feature.type in CDSFeature.types:
                    cds_features.append(feature)
                else:
                    warnings.warn(
                        UnknownGenBankFeatureWarning(
                            f"Feature {feature} was associated with a gene but was not of a known type "
                            f"and will be ignored"
                        )
                    )

            if len(transcript_features) > 1 and len(cds_features) > 1:
                warnings.warn(
                    DuplicateTranscriptWarning(
                        f"Grouping by locus tag found multiple coding transcript features on the same sequence "
                        f"with the same locus tag:"
                        f"\n{transcript_features[0]}\n{transcript_features[1]}\n"
                        f"Extra transcripts will be skipped"
                    )
                )
                transcript_features = [transcript_features[0]]

            grouped_features.append(GroupedGeneFeatures(seqrecord, gene_feature, transcript_features, cds_features))
        self.grouped_gene_features[idx].extend(grouped_features)

    def _parse_features(
        self,
    ):
        """
        Extract all generic features from a SeqRecord. These are anything that did not qualify as a gene, based
        on the feature type being one of the known members of :class:`biocantor.io.genbank.constants.GenBankFeatures`.

        Feature collections are inferred through the ``locus_tag`` field. Any items without such a tag are treated
        separately.

        """
        for idx, seqrecord in enumerate(self.seq_records):
            for locus_tag, features in itertools.groupby(
                self.feature_features[idx],
                key=lambda f: f.qualifiers.get(KnownQualifiers.LOCUS_TAG.value, [""])[0],
            ):
                if not locus_tag:
                    # we are in the null scenario,
                    # meaning that there are no locus tag information and thus no groupings.
                    for feature in features:
                        feature_collection = FeatureIntervalGenBankCollection([feature], seqrecord)
                        if feature_collection:
                            self.feature_collections[idx].append(feature_collection)
                else:
                    feature_collection = FeatureIntervalGenBankCollection(list(features), seqrecord)
                    if feature_collection:
                        self.feature_collections[idx].append(feature_collection)

    def _convert_seqfeature_to_gene(
        self, grouped_gene_features: GroupedGeneFeatures, seqrecord: SeqRecord
    ) -> GeneFeature:
        # if there is no gene level feature, then it must be inferred
        if not grouped_gene_features.gene_feature:
            if grouped_gene_features.transcript_features:
                gene = GeneFeature.from_transcript_or_cds_feature(
                    grouped_gene_features.transcript_features[0], seqrecord
                )
            # CDS feature(s) must exist
            else:
                gene = GeneFeature.from_transcript_or_cds_feature(grouped_gene_features.cds_features[0], seqrecord)
        else:
            gene = GeneFeature(grouped_gene_features.gene_feature, seqrecord)

        if grouped_gene_features.transcript_features and grouped_gene_features.cds_features:
            # there must be exactly one transcript_feature because cds_features exist
            for cds_feature in grouped_gene_features.cds_features:
                gene.add_child(grouped_gene_features.transcript_features[0], cds_feature)
        # this must be a non-coding gene
        elif grouped_gene_features.transcript_features:
            for transcript_feature in grouped_gene_features.transcript_features:
                gene.add_child(transcript_feature)
        # gene -> CDS without a transcript level feature
        elif grouped_gene_features.cds_features:
            for cds_feature in grouped_gene_features.cds_features:
                gene.add_child(cds_feature)
        if not gene.has_children:
            gene.infer_child()
            warnings.warn(GenBankEmptyGeneWarning(f"Gene {gene} has no valid children features"))
        return gene

    def _convert_seqfeatures_to_genes(self):
        """
        After the gene-like features have been grouped by either position, locus tag, or both, the groups
        are evaluated and converted into :class:`GeneFeature` objects.

        Gene-level objects with no children (transcripts or CDSes) are skipped.
        """
        for i, seqrecord_grouped_gene_features in enumerate(self.grouped_gene_features):
            genes = []
            for grouped_gene_features in seqrecord_grouped_gene_features:
                seqrecord = grouped_gene_features.seqrecord
                gene = self._convert_seqfeature_to_gene(grouped_gene_features, seqrecord)
                genes.append(gene)
            self.genes[i].extend(genes)

    def _export_annotation_collections(self) -> Iterator[ParsedAnnotationRecord]:
        """
        The final step of GenBank parsing is exporting the annotations as :class:`ParsedAnnotationRecord`.

        These objects contain both the annotations as a :class:`AnnotationCollectionModel`, as well as the associated
        :class:`SeqRecord`, which allows for construction of a `AnnotationCollection` with sequence information.
        """
        if self.num_genes + self.num_feature_collections == 0:
            raise EmptyGenBankError("GenBank parsing produced zero genes and zero features.")

        for i, seqrecord in enumerate(self.seq_records):
            genes = [GeneFeature.to_gene_model(x) for x in sorted(self.genes[i], key=lambda x: x.start)]
            feature_collections = [
                FeatureIntervalGenBankCollection.to_feature_model(x)
                for x in sorted(self.feature_collections[i], key=lambda x: x.start)
            ]

            if self.parsed_variants:
                variant_collections = VariantIntervalCollectionModel.Schema().dump(
                    self.parsed_variants.get(seqrecord.id), many=True
                )
            else:
                variant_collections = []

            annotation = AnnotationCollectionModel.Schema().load(
                dict(
                    genes=genes,
                    feature_collections=feature_collections,
                    variant_collections=variant_collections,
                    name=seqrecord.id,
                    sequence_name=seqrecord.id,
                    start=0,
                    end=len(seqrecord),
                    qualifiers=self.sources[i].qualifiers if self.sources[i] else None,
                )
            )
            yield ParsedAnnotationRecord(annotation=annotation, seqrecord=seqrecord)


class SortedGenBankParser(BaseGenBankParser):
    """
    The Sorted GenBank parser relies entirely on increasing genomic position to partition features into genes or
    feature groups. This is inherently challenging because of issues like overlapping genes or multiple isoforms.
    """

    genbank_parser_type = GenBankParserType.SORTED

    def parse(self) -> Iterator[ParsedAnnotationRecord]:
        self._extract_seqfeatures_from_seqrecords()
        self._group_gene_features_by_position()
        self._convert_seqfeatures_to_genes()
        self._parse_features()
        yield from self._export_annotation_collections()

    def _extract_seqfeatures_from_seqrecords(self):
        for i, seqrecord in enumerate(self.seq_records):
            gene_filtered_features = []
            remaining_features = []
            for f in seqrecord.features:
                if not self.validate_seqfeature(f):
                    continue
                elif f.type in GENBANK_GENE_FEATURES:
                    gene_filtered_features.append(f)
                elif f.type == MetadataFeatures.SOURCE.value:
                    self.sources[i] = f
                else:
                    remaining_features.append(f)

            self.gene_filtered_features[i].extend(self._sort_features_by_position_and_type(gene_filtered_features))
            self.feature_features[i].extend(remaining_features)

    def _group_gene_features_by_position(self):
        """
        Soerted parser groups features using :meth:`BaseGenBankParser._group_features_by_position()`.
        """
        for i, seqrecord in enumerate(self.seq_records):
            self._group_features_by_position(self.gene_filtered_features[i], seqrecord, i)


class LocusTagGenBankParser(BaseGenBankParser):
    """
    The LocusTag parser expects that every gene feature in the GenBank file contains a /locus_tag qualifier.

    Gene-type features without a locus tag qualifier are ignored, and an exception is raised if multiple
    gene features have the same locus tag.
    """

    genbank_parser_type = GenBankParserType.LOCUS_TAG

    def parse(self) -> Iterator[ParsedAnnotationRecord]:
        self._extract_seqfeatures_from_seqrecords()
        self._group_gene_features_by_locus_tag()
        self._convert_seqfeatures_to_genes()
        self._parse_features()
        yield from self._export_annotation_collections()

    def _extract_seqfeatures_from_seqrecords(self):
        for i, seqrecord in enumerate(self.seq_records):
            gene_filtered_features = []
            remaining_features = []
            for f in seqrecord.features:
                if not self.validate_seqfeature(f):
                    continue
                elif f.type in GENBANK_GENE_FEATURES and KnownQualifiers.LOCUS_TAG.value in f.qualifiers:
                    gene_filtered_features.append(f)
                elif f.type == MetadataFeatures.SOURCE.value:
                    self.sources[i] = f
                else:
                    remaining_features.append(f)

            locus_tag_sorted_gene_filtered_features = sorted(
                gene_filtered_features, key=lambda f: f.qualifiers[KnownQualifiers.LOCUS_TAG.value]
            )
            self.gene_filtered_features[i].extend(locus_tag_sorted_gene_filtered_features)
            self.feature_features[i].extend(remaining_features)

    def _group_gene_features_by_locus_tag(self):
        """
        Locus tag parser groups features using :meth:`BaseGenBankParser._group_features_by_locus_tag()`.
        """
        for i, seqrecord in enumerate(self.seq_records):
            self._group_features_by_locus_tag(self.gene_filtered_features[i], seqrecord, i)


class HybridGenBankParser(LocusTagGenBankParser, SortedGenBankParser):
    """
    The Hybrid parsing mode combines both LocusTag and Sorted parsing. LocusTag is preferentially used,
    with features that either lack a locus tag or with duplicate tags are be sent to the Sorted parser.
    """

    genbank_parser_type = GenBankParserType.HYBRID

    def __init__(
        self,
        seq_records: List[SeqRecord],
        parsed_variants: Dict[str, List[VariantIntervalCollectionModel]],
        gene_parse_func: Callable[[GeneFeature], Dict[str, Any]],
        feature_parse_func: Callable[[FeatureIntervalGenBankCollection], Dict[str, Any]],
    ):
        super().__init__(seq_records, parsed_variants, gene_parse_func, feature_parse_func)
        self.gene_filtered_features_without_locus_tag: List[List[SeqFeature]] = [list() for _ in self.seq_records]

    def parse(self) -> Iterator[ParsedAnnotationRecord]:
        self._extract_seqfeatures_from_seqrecords()
        self._identify_locus_tag_collisions()
        self._group_gene_features_by_locus_tag_and_position()
        self._convert_seqfeatures_to_genes()
        self._parse_features()
        yield from self._export_annotation_collections()

    def _extract_seqfeatures_from_seqrecords(self):
        """
        Hybrid parser partitions all features by locus tag, if they exist.
        Any features without a locus tag or with duplicate tags will be sent to the Sorted parser.
        """
        for i, seqrecord in enumerate(self.seq_records):
            gene_filtered_features = []
            gene_filtered_features_without_locus_tag = []
            remaining_features = []
            for f in seqrecord.features:
                if not self.validate_seqfeature(f):
                    continue
                elif f.type in GENBANK_GENE_FEATURES:
                    if KnownQualifiers.LOCUS_TAG.value in f.qualifiers:
                        gene_filtered_features.append(f)
                    else:
                        gene_filtered_features_without_locus_tag.append(f)
                elif f.type == MetadataFeatures.SOURCE.value:
                    self.sources[i] = f
                else:
                    remaining_features.append(f)

            sorted_gene_filtered_features = sorted(
                gene_filtered_features, key=lambda f: f.qualifiers[KnownQualifiers.LOCUS_TAG.value]
            )

            self.gene_filtered_features[i].extend(sorted_gene_filtered_features)
            self.gene_filtered_features_without_locus_tag[i].extend(gene_filtered_features_without_locus_tag)
            self.feature_features[i].extend(remaining_features)

    def _identify_locus_tag_collisions(self):
        """
        Identify duplicated locus tags on gene features across all of the SeqRecords. Raise a warning for each
        one, then reassign it to the Sorted parser.

        Also sorts the objects in ``gene_filtered_features_without_locus_tag`` to prepare for sorted parsing.
        """
        bad_locus_tags = set()
        for locus_tag, features in itertools.groupby(
            sorted(
                itertools.chain.from_iterable(self.gene_filtered_features),
                key=lambda x: x.qualifiers.get(KnownQualifiers.LOCUS_TAG.value, None)[0],
            ),
            key=lambda x: x.qualifiers[KnownQualifiers.LOCUS_TAG.value][0],
        ):
            if len([f for f in features if f.type == GeneFeatures.GENE.value]) > 1:
                warnings.warn(
                    GenBankDuplicateLocusTagWarning(
                        f"Locus tag {locus_tag} seen on two or more gene features. "
                        f"These features will be sent to the Sorted parser."
                    )
                )
                bad_locus_tags.add(locus_tag)

        for i, features in enumerate(self.gene_filtered_features):
            good_features = []
            bad_features = []
            for feature in features:
                if feature.qualifiers[KnownQualifiers.LOCUS_TAG.value][0] in bad_locus_tags:
                    bad_features.append(feature)
                else:
                    good_features.append(feature)
            self.gene_filtered_features[i] = good_features
            self.gene_filtered_features_without_locus_tag[i] = self._sort_features_by_position_and_type(
                list(itertools.chain(self.gene_filtered_features_without_locus_tag[i], bad_features))
            )

    def _group_gene_features_by_locus_tag_and_position(self):
        """
        Hybrid parser groups genes by both position and locus tag, after they were partitioned by
        :meth:`HybridGenBankParser._extract_seqfeatures_from_seqrecords()`.
        """
        for i, seqrecord in enumerate(self.seq_records):
            self._group_features_by_position(self.gene_filtered_features_without_locus_tag[i], seqrecord, i)
            self._group_features_by_locus_tag(self.gene_filtered_features[i], seqrecord, i)


def parse_genbank(
    genbank_handle_or_path: Union[TextIO, str, pathlib.Path],
    variant_handle_or_path: Optional[Union[TextIO, str, pathlib.Path]] = None,
    parsed_variants: Optional[Dict[str, List[VariantIntervalCollectionModel]]] = None,
    gene_parse_func: Callable[[GeneFeature], Dict[str, Any]] = GeneFeature.to_gene_model,
    feature_parse_func: Callable[
        [FeatureIntervalGenBankCollection], Dict[str, Any]
    ] = FeatureIntervalGenBankCollection.to_feature_model,
    gbk_type: GenBankParserType = GenBankParserType.HYBRID,
    allow_duplicate_sequence_identifiers: bool = False,
) -> Iterator[ParsedAnnotationRecord]:
    """This is the main GenBank parsing function. The parse function implemented in :class:`GeneFeature` can be
    over-ridden to provide a custom implementation.

    Args:
        genbank_handle_or_path: An open GenBank file or a path to a locally stored GenBank file.
        variant_handle_or_path: Optional open handle to a VCF file. Mutually exclusive with ``parsed_variants``.
        parsed_variants: Optional parsed variants. Mutually exclusive with ``variant_handle_or_path``.
        gene_parse_func: Optional gene parse function implementation. Defaults to :meth:`GeneFeature.to_gene_model()`
            implemented in this module.
        feature_parse_func: Optional feature interval parse function implementation.
            Defaults to :meth:`FeatureIntervalGenBankCollection.to_feature_model()` implemented in this module.
        gbk_type: Use Hybrid, Sorted or LocusTag based parsing? Defaults to Hybrid.
        allow_duplicate_sequence_identifiers: Should this parser raise an exception if the same identifier is seen
            twice? Defaults to `False`.

    Yields:
         :class:`ParsedAnnotationRecord`.
    """
    seq_records = list(SeqIO.parse(genbank_handle_or_path, format="genbank"))

    if variant_handle_or_path and parsed_variants:
        raise InvalidInputError("Cannot pass both a variant file and parsed variants")
    elif variant_handle_or_path:
        parsed_variants = parse_vcf_file(variant_handle_or_path)

    if allow_duplicate_sequence_identifiers is False:
        seqrecords_dict = {}
        for rec in seq_records:
            if rec.id in seqrecords_dict:
                raise DuplicateSequenceException(f"Sequence {rec.id} found twice in GenBank file.")
            seqrecords_dict[rec.id] = rec

    if parsed_variants:
        seqrecord_ids = {x.id for x in seq_records}
        for seq_id in parsed_variants:
            if seq_id not in seqrecord_ids:
                raise InvalidInputError(f"Sequence ID {seq_id} found in variant file but not found in GenBank file.")

    if gbk_type == GenBankParserType.SORTED:
        parser = SortedGenBankParser(seq_records, parsed_variants, gene_parse_func, feature_parse_func)
    elif gbk_type == GenBankParserType.HYBRID:
        parser = HybridGenBankParser(seq_records, parsed_variants, gene_parse_func, feature_parse_func)
    else:
        parser = LocusTagGenBankParser(seq_records, parsed_variants, gene_parse_func, feature_parse_func)
    yield from parser.parse()
