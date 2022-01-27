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
default for this parsing module. This does not work for genomes with alternative isoforms.

The generic parsing function that interprets the BioPython results to BioCantor data models is implemented in
:meth:`GeneFeature.to_gene_model()`. This function can be over-ridden to provide custom parsing implementations.
"""
import itertools
import pathlib
import warnings
from abc import ABC
from collections import Counter
from copy import deepcopy
from typing import Optional, TextIO, Iterator, List, Dict, Callable, Tuple, Any, Union

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from inscripta.biocantor.gene import Biotype, CDSInterval, CDSFrame
from inscripta.biocantor.io.exc import (
    DuplicateSequenceException,
    InvalidCDSIntervalWarning,
    StrandViolationWarning,
    DuplicateFeatureWarning,
    DuplicateTranscriptWarning,
)
from inscripta.biocantor.io.features import extract_feature_types, extract_feature_name_id, merge_qualifiers
from inscripta.biocantor.io.genbank.constants import (
    GeneFeatures,
    TranscriptFeatures,
    GeneIntervalFeatures,
    MetadataFeatures,
    GenBankParserType,
    KnownQualifiers,
    GENBANK_GENE_FEATURES,
)
from inscripta.biocantor.io.genbank.exc import (
    GenBankParserError,
    EmptyGenBankError,
    GenBankLocusTagError,
    GenBankLocationException,
    GenBankNullStrandException,
)
from inscripta.biocantor.io.models import (
    GeneIntervalModel,
    AnnotationCollectionModel,
    FeatureIntervalCollectionModel,
)
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.location import (
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
        self.feature = feature
        self.record = record
        self.children = []

    def __iter__(self) -> Iterator[SeqFeature]:
        """Pre-order depth first traversal of this feature."""
        yield self
        for child in self.children:
            yield from child

    def __repr__(self):
        return "\n".join((str(x) for x in self))

    @property
    def type(self) -> str:
        return str(self.feature.type)

    @property
    def has_children(self) -> bool:
        return len(self.children) > 0

    @property
    def strand(self) -> int:
        return self.feature.strand


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
        self.features = features

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
        for feature in cls.features:
            interval_starts = []
            interval_ends = []
            for loc in sorted(feature.location.parts, key=lambda p: p.start):
                interval_starts.append(loc.nofuzzy_start)
                interval_ends.append(loc.nofuzzy_end)
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

    def __str__(self):
        return self.feature.__repr__()

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

    def add_child(self, feature: SeqFeature):
        """Add a new feature as a child. Infer Transcripts if this child is a CDS or exon feature."""
        if feature.type in TranscriptFeature.types:
            self.children.append(TranscriptFeature(feature, self.record))
        elif feature.type in IntervalFeature.types:
            # infer a transcript
            tx_feature = deepcopy(feature)
            if tx_feature.type == GeneIntervalFeatures.CDS.value:
                tx_feature.type = TranscriptFeatures.CODING_TRANSCRIPT.value
            # this means we have an exon as a direct child of a gene
            elif tx_feature.type == GeneIntervalFeatures.EXON.value:
                tx_feature.type = TranscriptFeatures.CODING_TRANSCRIPT.value
            else:
                tx_feature.type = TranscriptFeatures[tx_feature.type].value
            tx = TranscriptFeature(tx_feature, self.record)
            self.children.append(tx)
            tx.add_child(feature)
        else:
            raise GenBankParserError(f"Invalid feature type {feature.type}")

    def finalize(self):
        """Make sure we have a full hierarchy; infer children if necessary.

        This is often needed for non-coding genes which lack an explicit exon.
        """
        for transcript in self.children:
            transcript.infer_exon_features()

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
            exon_interval = tx.find_exon_interval()
            exon_starts = []
            exon_ends = []
            for block in exon_interval.blocks:
                exon_starts.append(block.start)
                exon_ends.append(block.end)

            strand = Strand.from_int(tx.strand)

            cds_interval = tx.find_cds_interval()
            if cds_interval.is_empty:
                cds_starts = None
                cds_ends = None
                cds_frames = []
            else:
                cds_starts = []
                cds_ends = []
                for block in cds_interval.blocks:
                    cds_starts.append(block.start)
                    cds_ends.append(block.end)
                # TODO: this will produce invalid Frames if there is a programmed frameshift or other such problem.
                # This is a limitation of genbank files. HOWEVER, this code could detect the presence of a
                # overlapping interval and try all possible frames there to see if any produce a valid translation.
                cds_frames = tx.construct_frames(cds_interval)

            if "pseudo" in tx.feature.qualifiers:
                transcript_biotype = Biotype.pseudogene
            elif tx.feature.type == TranscriptFeatures.CODING_TRANSCRIPT.value:
                transcript_biotype = Biotype.protein_coding
            else:
                transcript_biotype = Biotype[tx.feature.type]
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
                warnings.warn(DuplicateTranscriptWarning("Transcript {tx_model} found twice within a gene"))
            else:
                transcripts.append(tx_model)

        # pick most common transcript type; hacky
        gene_biotype = tx_biotypes.most_common(1)[0][0]
        gene = GeneIntervalModel.Schema().load(
            dict(
                transcripts=transcripts,
                gene_id=cls.feature.qualifiers.get(KnownQualifiers.GENE_ID.value, [None])[0],
                gene_symbol=cls.feature.qualifiers.get(KnownQualifiers.GENE.value, [None])[0],
                locus_tag=cls.feature.qualifiers.get(KnownQualifiers.LOCUS_TAG.value, [None])[0],
                gene_type=gene_biotype.name,
                qualifiers=cls.feature.qualifiers,
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

    def __str__(self):
        return f"--> {self.feature.__repr__()}"

    def add_child(self, feature: SeqFeature):
        self.children.append(IntervalFeature(feature, self.record))

    def infer_exon_features(self):
        """Commonly, non-coding genes lack IntervalFeatures, coding genes only have CDS features"""
        # no children means this is likely a non-coding transcript with no exon features
        # no exon features means this is likely a coding transcript with no exon features
        if len(self.children) == 0 or len(list(self.exon_features)) == 0:
            # add an exon with the same interval as the transcript
            feature = deepcopy(self.feature)
            feature.type = GeneIntervalFeatures.EXON.value
            self.add_child(feature)

    def construct_frames(self, cds_interval: Location) -> List[str]:
        """We need to build frames. Since GenBank lacks this info, do our best"""
        # make 0 based offset, if possible, otherwise assume always in frame
        frame = int(self.children[0].feature.qualifiers.get(KnownQualifiers.CODON_START.value, [1])[0]) - 1
        frame = CDSFrame.from_int(frame)
        frames = CDSInterval.construct_frames_from_location(cds_interval, frame)
        return [x.name for x in frames]

    @property
    def exon_features(self) -> SeqFeature:
        for f in self.children:
            if f.type == GeneIntervalFeatures.EXON.value:
                yield f

    @property
    def cds_features(self) -> SeqFeature:
        for f in self.children:
            if f.type == GeneIntervalFeatures.CDS.value:
                yield f

    def get_qualifier_from_tx_or_cds_features(self, qualifier: str) -> Optional[str]:
        """Get a specific qualifier, if it exists. Look at tx first, then children"""
        if qualifier in self.feature.qualifiers:
            return self.feature.qualifiers[qualifier][0]
        for feature in self.cds_features:
            if qualifier in feature.feature.qualifiers:
                return feature.feature.qualifiers[qualifier][0]

    def iterate_exon_intervals(self) -> Iterator[Tuple[int, int]]:
        """Iterate over the location parts"""
        for exon in sorted(self.exon_features, key=lambda e: e.feature.location.nofuzzy_start):
            for part in sorted(exon.feature.location.parts, key=lambda p: p.start):
                yield int(part.nofuzzy_start), int(part.nofuzzy_end)

    def iterate_cds_intervals(self) -> Iterator[Tuple[int, int]]:
        """Iterate over the CDS location parts"""
        for cds in sorted(self.cds_features, key=lambda e: e.feature.location.nofuzzy_start):
            for part in sorted(cds.feature.location.parts, key=lambda p: p.start):
                yield int(part.nofuzzy_start), int(part.nofuzzy_end)

    def find_exon_interval(self) -> CompoundInterval:
        """Finds the Location of the Exons."""
        if self._exon_interval:
            return self._exon_interval
        exon_starts = []
        exon_ends = []
        for exon_start, exon_end in self.iterate_exon_intervals():
            exon_starts.append(exon_start)
            exon_ends.append(exon_end)
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
        cds_starts = []
        cds_ends = []
        for cds_start, cds_end in self.iterate_cds_intervals():
            cds_starts.append(cds_start)
            cds_ends.append(cds_end)
        if len(cds_starts) == 0:
            return EmptyLocation()
        # must use SingleInterval here because otherwise the optimization step of cds_interval.intersection below
        # will return a SingleInterval, and the equality comparison will raise a spurious InvalidCDSIntervalWarning
        elif len(cds_starts) == 1:
            cds_interval = SingleInterval(cds_starts[0], cds_ends[0], Strand.from_int(self.strand))
        else:
            cds_interval = CompoundInterval(
                cds_starts,
                cds_ends,
                Strand.from_int(self.strand),
            )
        cds_interval_intersection_with_exons = cds_interval.intersection(self.find_transcript_interval())
        if cds_interval_intersection_with_exons != cds_interval:
            warnings.warn(
                InvalidCDSIntervalWarning(
                    f"CDS Interval {cds_interval} was sliced down to {cds_interval_intersection_with_exons} "
                    f"because that is the exon boundaries found"
                )
            )
        return cds_interval_intersection_with_exons

    def merge_cds_qualifiers_to_transcript(self) -> Dict[str, List[str]]:
        """
        If there were distinct transcript-level features, the qualifiers on the CDS feature will be lost
        when converting to the BioCantor data model unless those qualifiers are rolled into the qualifiers on
        the transcript feature.
        """
        qualifiers = {key: set(vals) for key, vals in self.feature.qualifiers.items()}
        for cds_feature in self.cds_features:
            for key, vals in cds_feature.feature.qualifiers.items():
                if key not in qualifiers:
                    qualifiers[key] = set(vals)
                else:
                    qualifiers[key].update(vals)
        return {key: list(vals) for key, vals in qualifiers.items()}


class IntervalFeature(Feature):
    """A set of intervals"""

    types = {"CDS", "exon"}

    def __str__(self):
        return f"----> {self.feature.__repr__()}"


def _construct_gene_from_feature(
    feature: SeqFeature,
    seqrecord: SeqRecord,
    cls_or_fn: Callable[[SeqFeature, SeqRecord], GeneFeature],
) -> Optional[GeneFeature]:
    """
    Convenience function for producing :class:`GeneFeature` from a `SeqFeature`, handling both
    possible constructor routes (construction from a ``gene`` feature as found in the GenBank,
    or inference from a transcript/interval level feature in case no ``gene`` level feature was found).

    This wrapper function catches exceptions raised for common errors, and converts them to warnings
    as appropriate.
    """
    try:
        return cls_or_fn(feature, seqrecord)
    except GenBankNullStrandException:
        warnings.warn(
            StrandViolationWarning(f"Found multiple strands for feature {feature}. This feature will be skipped.")
        )


def _construct_feature_collection_from_features(
    features: List[SeqFeature],
    seqrecord: SeqRecord,
) -> Optional[FeatureIntervalGenBankCollection]:
    """
    Convenience function for producing :class:`FeatureIntervalGenBankCollection` from a `SeqFeature`.

    This wrapper function catches exceptions raised for common errors, and converts them to warnings
    as appropriate.
    """
    try:
        return FeatureIntervalGenBankCollection(features, seqrecord)
    except GenBankNullStrandException:
        warnings.warn(
            StrandViolationWarning(
                f"Found multiple strands for feature group {features}. " f"This feature collection will be skipped."
            )
        )


def parse_genbank(
    genbank_handle_or_path: Union[TextIO, str, pathlib.Path],
    parse_func: Optional[Callable[[GeneFeature], Dict[str, Any]]] = GeneFeature.to_gene_model,
    feature_parse_func: Optional[
        Callable[[FeatureIntervalGenBankCollection], Dict[str, Any]]
    ] = FeatureIntervalGenBankCollection.to_feature_model,
    gbk_type: Optional[GenBankParserType] = GenBankParserType.HYBRID,
) -> Iterator[ParsedAnnotationRecord]:
    """This is the main GenBank parsing function. The parse function implemented in :class:`GeneFeature` can be
    over-ridden to provide a custom implementation.

    Args:
        genbank_handle_or_path: An open GenBank file or a path to a locally stored GenBank file.
        parse_func: Optional parse function implementation.
        feature_parse_func: Optional feature interval parse function implementation.
        gbk_type: Do we want to use model 1 or model 2? Must be one of ``sorted``, ``locus_tag``.

    Yields:
         :class:`ParsedAnnotationRecord`.
    """
    seq_records = list(SeqIO.parse(genbank_handle_or_path, format="genbank"))

    seqrecords_dict = {}
    for rec in seq_records:
        if rec.id in seqrecords_dict:
            raise DuplicateSequenceException(f"Sequence {rec.id} found twice in GenBank file.")
        seqrecords_dict[rec.id] = rec

    if gbk_type == GenBankParserType.SORTED:
        gene_records = group_gene_records_from_sorted_genbank(seq_records, parse_func, feature_parse_func)
    else:
        gene_records = group_gene_records_by_locus_tag(seq_records, parse_func, feature_parse_func, gbk_type)
    yield from gene_records


def group_gene_records_from_sorted_genbank(
    record_iter: Iterator[SeqRecord],
    parse_func: Callable[[GeneFeature], Dict[str, Any]],
    feature_parse_func: Callable[[FeatureIntervalGenBankCollection], Dict[str, Any]],
) -> Iterator[ParsedAnnotationRecord]:
    """Model 1: position sorted GenBank.

    This function looks for canonical gene records:
        gene -> Optional(mRNA) -> CDS records
    It also looks for canonical non-coding records:
        gene -> {misc_RNA,tRNA,rRNA,etc)

    It also will infer non-canonical record types, including non-coding transcripts and coding genes
    from isolated CDS/non-coding features (those without a gene feature before them in the sort order).

    Any features that do not fit the above bins are interpreted as generic features.

    Some GenBank files are improperly ordered, and will have things like the CDS feature first, or the mRNA feature
    first. To try and capture this, the full set of records are sorted first by position, then in the order:

    gene
    mRNA
    CDS
    exon
    anything else

    Args:
        record_iter: Iterator of SeqRecord objects.
        parse_func: Optional parse function implementation.
        feature_parse_func: Optional feature interval parse function implementation.

    Yields:
        :class:`ParsedAnnotationRecord`.
    """
    tot_genes = 0
    tot_features = 0
    for seqrecord in record_iter:
        gene = None
        source = None
        genes = []
        # capture non-gene intervals downstream
        feature_features = []

        # sort features to try to capture weirdly ordered genbank files
        sorted_features = sorted(
            seqrecord.features,
            key=lambda x: (
                # features with invalid coordinates will have no location
                # this will raise GenBankLocationException once it is consumed by a subclass of Feature
                x.location.nofuzzy_start if x.location else -1,
                x.type != GeneFeatures.GENE.value,
                x.type != TranscriptFeatures.CODING_TRANSCRIPT.value,
                x.type != GeneIntervalFeatures.CDS.value,
                x.type != GeneIntervalFeatures.EXON.value,
            ),
        )
        for feature in sorted_features:
            # try to capture the Source field, if it exists
            if feature.type == MetadataFeatures.SOURCE.value:
                source = feature
            # base case for start; iterate until we find a gene
            elif gene is None:
                if feature.type in GeneFeature.types:
                    gene = _construct_gene_from_feature(feature, seqrecord, GeneFeature)
                    # gene is None if it was not parseable
                    if not gene:
                        continue
                # base case for starting with a isolated ncRNA or CDS feature; immediately add them
                # and reset the gene to None
                elif feature.type in TranscriptFeature.types or feature.type in IntervalFeature.types:
                    gene = _construct_gene_from_feature(feature, seqrecord, GeneFeature.from_transcript_or_cds_feature)
                    # gene is None if it was not parseable
                    if gene:
                        gene.finalize()
                        gene = parse_func(gene)
                        genes.append(gene)
                        gene = None
                # this must be a generic feature
                else:
                    feature_features.append(feature)
            # next gene; re-set the gene object and report out the collection
            elif feature.type in GeneFeature.types:
                if gene.has_children:
                    gene.finalize()
                    gene = parse_func(gene)
                    genes.append(gene)
                gene = _construct_gene_from_feature(feature, seqrecord, GeneFeature)
                if not gene:
                    continue
            elif feature.type in TranscriptFeature.types:
                # if the current gene is non-empty, and the feature is not a mRNA, then this is a isolated ncRNA
                # finish this gene and start a new one
                if feature.type != TranscriptFeatures.CODING_TRANSCRIPT and gene.has_children:
                    gene.finalize()
                    gene = parse_func(gene)
                    genes.append(gene)
                    gene = _construct_gene_from_feature(feature, seqrecord, GeneFeature.from_transcript_or_cds_feature)
                    # gene is None if it was not parseable
                    if not gene:
                        continue
                else:
                    gene.add_child(feature)
            elif feature.type in IntervalFeature.types:
                if not gene.has_children:
                    gene.add_child(feature)
                else:
                    gene.children[-1].add_child(feature)
            else:
                feature_features.append(feature)

        # gene could be None if this record has no annotations
        if gene is not None and gene.has_children:
            gene.finalize()
            gene = parse_func(gene)
            genes.append(gene)

        if source is not None:
            source_qualifiers = source.qualifiers
        else:
            source_qualifiers = None

        feature_collections = _extract_generic_features(seqrecord, feature_features, feature_parse_func)

        tot_features += len(feature_collections) if feature_collections else 0
        tot_genes += len(genes) if genes else 0

        annotation = AnnotationCollectionModel.Schema().load(
            dict(
                genes=genes,
                feature_collections=feature_collections,
                sequence_name=seqrecord.id,
                start=0,
                end=len(seqrecord),
                qualifiers=source_qualifiers,
            )
        )
        yield ParsedAnnotationRecord(annotation=annotation, seqrecord=seqrecord)

    if tot_genes + tot_features == 0:
        raise EmptyGenBankError("GenBank parsing produced zero genes and zero features.")


def group_gene_records_by_locus_tag(
    record_iter: Iterator[SeqRecord],
    parse_func: Callable[[GeneFeature], Dict[str, Any]],
    feature_parse_func: Callable[[FeatureIntervalGenBankCollection], Dict[str, Any]],
    genbank_parser_type: GenBankParserType = GenBankParserType.LOCUS_TAG,
) -> Iterator[ParsedAnnotationRecord]:
    """Model 2: ``locus_tag`` defined GenBank.

    All feature types that qualify within the hierarchical structure, possess a locus_tag, and whose feature type
    are valid for a known transcribed interval type, will be included in the gene parsing.

    All other feature types will become generic features (FeatureIntervals), unless we are in hybrid mode.

    In hybrid mode, locus_tag is used first, then all of the remaining features are sent to the
    sorted parser.

    Args:
        record_iter: Iterator of SeqRecord objects.
        parse_func: Optional parse function implementation.
        feature_parse_func: Optional feature interval parse function implementation.
        genbank_parser_type: Optional parser type. Changing this to GenBankParserType.HYBRID
            will enable hybrid parsing mode.

    Yields:
        :class:`ParsedAnnotationRecord`.
    """
    if genbank_parser_type not in [GenBankParserType.LOCUS_TAG, GenBankParserType.HYBRID]:
        raise GenBankParserError("Must use either locus_tag or hybrid")

    tot_genes = 0
    tot_features = 0
    for seqrecord in record_iter:
        gene_filtered_features = []
        remaining_features = []
        source = None
        for f in seqrecord.features:
            if f.type in GENBANK_GENE_FEATURES and KnownQualifiers.LOCUS_TAG.value in f.qualifiers:
                gene_filtered_features.append(f)
            elif f.type == MetadataFeatures.SOURCE.value:
                source = f
            else:
                remaining_features.append(f)

        sorted_gene_filtered_features = sorted(
            gene_filtered_features, key=lambda f: f.qualifiers[KnownQualifiers.LOCUS_TAG.value]
        )

        genes = []
        for locus_tag, gene_features in itertools.groupby(
            sorted_gene_filtered_features, key=lambda f: f.qualifiers[KnownQualifiers.LOCUS_TAG.value][0]
        ):
            # sort the features for this locus tag to bubble the "gene" feature to the top, if it exists
            gene_features = sorted(gene_features, key=lambda f: f.type != GeneFeatures.GENE.value)

            # do we have more than one gene with this locus_tag?
            if len(gene_features) > 1 and gene_features[1].type == GeneFeatures.GENE.value:
                raise GenBankLocusTagError(
                    f"Grouping by locus tag found multiple gene features with the same locus tag:"
                    f"\n{gene_features[0]}\n{gene_features[1]}"
                )

            gene_feature = gene_features[0]
            if gene_feature.type == GeneFeatures.GENE.value:
                gene = _construct_gene_from_feature(gene_feature, seqrecord, GeneFeature)
            else:
                gene = _construct_gene_from_feature(gene_feature, seqrecord, GeneFeature.from_transcript_or_cds_feature)
            # gene is None if it was not parseable
            if not gene:
                continue

            for feature in gene_features[1:]:
                if feature.type in TranscriptFeature.types:
                    gene.add_child(feature)
                elif feature.type in IntervalFeature.types:
                    if len(gene.children) == 0:
                        gene.add_child(feature)
                    else:
                        gene.children[-1].add_child(feature)

            if gene.has_children:
                gene.finalize()
                gene = parse_func(gene)
                genes.append(gene)

        if source is not None:
            source_qualifiers = source.qualifiers
        else:
            source_qualifiers = None

        if genbank_parser_type == GenBankParserType.LOCUS_TAG:
            feature_collections = _extract_generic_features(seqrecord, remaining_features, feature_parse_func)
        else:
            # hybrid parsing mode
            tmp_seqrecord = deepcopy(seqrecord)
            tmp_seqrecord.features = remaining_features
            tmp_annotation = next(
                group_gene_records_from_sorted_genbank((tmp_seqrecord,), parse_func, feature_parse_func)
            )
            if tmp_annotation.annotation.feature_collections:
                feature_collections = [
                    FeatureIntervalCollectionModel.Schema().dump(x)
                    for x in tmp_annotation.annotation.feature_collections
                ]
            else:
                feature_collections = None
            if tmp_annotation.annotation.genes:
                genes.extend([GeneIntervalModel.Schema().dump(x) for x in tmp_annotation.annotation.genes])

        tot_features += len(feature_collections) if feature_collections else 0
        tot_genes += len(genes) if genes else 0

        annotation = AnnotationCollectionModel.Schema().load(
            dict(
                genes=genes,
                feature_collections=feature_collections,
                name=seqrecord.id,
                sequence_name=seqrecord.id,
                start=0,
                end=len(seqrecord),
                qualifiers=source_qualifiers,
            )
        )
        yield ParsedAnnotationRecord(annotation=annotation, seqrecord=seqrecord)

    if tot_genes + tot_features == 0:
        raise EmptyGenBankError("GenBank parsing produced zero genes and zero features.")


def _extract_generic_features(
    seqrecord: SeqRecord,
    filtered_features: List[SeqFeature],
    feature_parse_func: Callable[[FeatureIntervalGenBankCollection], Dict[str, Any]],
) -> Optional[List[Dict[str, Any]]]:
    """
    Extract all generic features from a SeqRecord. These are anything that did not qualify as a gene, based
    on the feature type being one of the known members of :class:`biocantor.io.genbank.constants.GenBankFeatures`.

    Feature collections are inferred through the ``locus_tag`` field. Any items without such a tag are treated
    separately.

    Args:
        seqrecord: A SeqRecord object.
        filtered_features: List of SeqFeature objects associated with the SeqRecord that are not gene-like.
        feature_parse_func: Optional feature interval parse function implementation.

    Returns:
        A list of dictionary representations of feature interval collections, or ``None`` if no feature intervals were
        found.
    """

    # sort by locus tag, or null if no locus tag is provided.
    sorted_filtered_features = sorted(
        filtered_features, key=lambda f: f.qualifiers.get(KnownQualifiers.LOCUS_TAG.value, [""])[0]
    )

    feature_collections = []
    for locus_tag, features in itertools.groupby(
        sorted_filtered_features, key=lambda f: f.qualifiers.get(KnownQualifiers.LOCUS_TAG.value, [""])[0]
    ):
        if not locus_tag:
            # we are in the null scenario, meaning that there are no locus tag information and thus no groupings.
            for feature in features:
                feature_collection = _construct_feature_collection_from_features([feature], seqrecord)
                if feature_collection:
                    feature_collections.append(feature_collection)
        else:
            feature_collection = _construct_feature_collection_from_features(list(features), seqrecord)
            if feature_collection:
                feature_collections.append(feature_collection)

    return [feature_parse_func(fc) for fc in feature_collections] if feature_collections else None
