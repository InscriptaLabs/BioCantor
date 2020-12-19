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
import logging
from abc import ABC
from collections import Counter
from copy import deepcopy
from typing import Optional, TextIO, Iterable, List, Dict, Callable, Tuple, Any

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSFrame
from inscripta.biocantor.io.genbank.constants import (
    GeneFeatures,
    TranscriptFeatures,
    IntervalFeatures,
    GenBankFeatures,
    MetadataFeatures,
    GenBankParserType,
)
from inscripta.biocantor.io.genbank.exc import GenBankValidationError
from inscripta.biocantor.io.models import GeneIntervalModel, AnnotationCollectionModel
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.location import Location
from inscripta.biocantor.location.location_impl import (
    SingleInterval,
    CompoundInterval,
    EmptyLocation,
)
from inscripta.biocantor.location.strand import Strand

logger = logging.getLogger(__name__)


class Feature(ABC):
    """Generic feature."""

    types = set()

    def __init__(self, feature: SeqFeature, record: SeqRecord):
        if feature.type not in self.types:
            raise GenBankValidationError(f"Invalid feature type {feature.type}")
        self.feature = feature
        self.record = record
        self.children = []

    def __iter__(self) -> Iterable[SeqFeature]:
        """Pre-order depth first traversal of this feature."""
        yield self
        for child in self.children:
            yield from child

    def __repr__(self):
        r = []
        for x in self:
            r.append(str(x))
        return "\n".join(r)

    @property
    def type(self) -> str:
        return str(self.feature.type)

    @property
    def has_children(self) -> bool:
        return len(self.children) > 0

    @property
    def strand(self) -> int:
        return self.feature.strand


class GeneFeature(Feature):
    """A gene."""

    types = {x.value for x in GeneFeatures}

    def __str__(self):
        return self.feature.__repr__()

    def add_child(self, feature: SeqFeature):
        """Add a new feature as a child. Infer Transcripts if this child is a CDS or exon feature."""
        if feature.type in TranscriptFeature.types:
            self.children.append(TranscriptFeature(feature, self.record))
        elif feature.type in IntervalFeature.types:
            # infer a transcript
            logger.debug(f"Inferring a transcript for {feature}")
            tx_feature = deepcopy(feature)
            if tx_feature.type == GenBankFeatures.CDS.value:
                tx_feature.type = TranscriptFeatures.CODING_TRANSCRIPT.value
            else:
                tx_feature.type = TranscriptFeatures[tx_feature.type].value
            tx = TranscriptFeature(tx_feature, self.record)
            self.children.append(tx)
            tx.add_child(feature)
        else:
            raise GenBankValidationError(f"Invalid feature type {feature.type}")

    def finalize(self):
        """Make sure we have a full hierarchy; infer children if necessary.

        This is often needed for non-coding genes which lack an explicit exon.
        """
        for transcript in self.children:
            transcript.infer_exon_features()

    @staticmethod
    def to_gene_model(cls: "GeneFeature") -> Dict[str, Any]:
        """Convert to a Dict representation of a :class:`biocantor.gene.GeneInterval` that can be used for analyses.

        This is the default function, that can be over-ridden by specific implementations.

        Looks for /transcript_id, /protein_id, and /gene on the transcript level, and
        looks for /gene_id, /gene, and /locus_tag on the gene level.
        """
        transcripts = []
        tx_biotypes = Counter()
        for tx in cls.children:
            exon_starts = []
            exon_ends = []
            for start, end in tx.iterate_intervals():
                exon_starts.append(start)
                exon_ends.append(end)

            strand = Strand.from_int(tx.strand)
            exon_interval = CompoundInterval(exon_starts, exon_ends, strand)

            cds_interval = tx.find_cds_interval(exon_interval)
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
                qualifiers=tx.feature.qualifiers,
                is_primary_tx=False,
                transcript_id=tx.feature.qualifiers.get("transcript_id", [None])[0],
                protein_id=tx.feature.qualifiers.get("protein_id", [None])[0],
                transcript_symbol=tx.feature.qualifiers.get("gene", [None])[0],
                transcript_type=transcript_biotype.name,
                sequence_name=tx.record.id,
            )
            transcripts.append(tx_model)

        # pick most common transcript type; hacky
        gene_biotype = tx_biotypes.most_common(1)[0][0]
        gene = GeneIntervalModel.Schema().load(
            dict(
                transcripts=transcripts,
                gene_id=cls.feature.qualifiers.get("gene_id", [None])[0],
                gene_symbol=cls.feature.qualifiers.get("gene", [None])[0],
                locus_tag=cls.feature.qualifiers.get("locus_tag", [None])[0],
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
            feature.type = GenBankFeatures.EXON.value
            self.add_child(feature)

    def construct_frames(self, cds_interval: Location) -> List[str]:
        """We need to build frames. Since GenBank lacks this info, do our best"""
        if self.feature.type != GenBankFeatures.CODING_TRANSCRIPT.value:
            return []
        # make 0 based offset, if possible, otherwise assume always in frame
        frame = int(self.children[0].feature.qualifiers.get("codon_start", [1])[0]) - 1
        frame = CDSFrame.from_int(frame)
        frames = []
        for block in cds_interval.blocks:
            frames.append(frame.name)
            frame = frame.shift(len(block))
        return frames

    @property
    def exon_features(self) -> SeqFeature:
        for f in self.children:
            if f.type == GenBankFeatures.EXON.value:
                yield f

    @property
    def cds_features(self) -> SeqFeature:
        for f in self.children:
            if f.type == GenBankFeatures.CDS.value:
                yield f

    def iterate_intervals(self) -> Iterable[Tuple[int, int]]:
        """Iterate over the location parts"""
        for exon in sorted(self.exon_features, key=lambda e: e.feature.location.nofuzzy_start):
            for part in sorted(exon.feature.location.parts, key=lambda p: p.start):
                yield int(part.start), int(part.end)

    def find_cds_interval(self, exon_interval: Location) -> Location:
        cds = sorted(self.cds_features, key=lambda c: c.feature.location.nofuzzy_start)
        if len(cds) == 0:
            return EmptyLocation()
        cds_i = SingleInterval(
            cds[0].feature.location.nofuzzy_start,
            cds[-1].feature.location.nofuzzy_end,
            Strand.from_int(self.feature.strand),
        )
        return exon_interval.intersection(cds_i)


class IntervalFeature(Feature):
    """A set of intervals"""

    types = {x.value for x in IntervalFeatures}

    def __str__(self):
        return f"----> {self.feature.__repr__()}"


def parse_genbank(
    genbank_handle: TextIO,
    parse_func: Optional[Callable[[GeneFeature], Dict[str, Any]]] = GeneFeature.to_gene_model,
    gbk_type: Optional[GenBankParserType] = GenBankParserType.LOCUS_TAG,
) -> Iterable[ParsedAnnotationRecord]:
    """This is the main GenBank parsing function. The parse function implemented in :class:`GeneFeature` can be
    over-ridden to provide a custom implementation.

    Args:
        genbank_handle: An open GenBank file.
        parse_func: Optional parse function implementation.
        gbk_type: Do we want to use model 1 or model 2? Must be one of ``sorted``, ``locus_tag``.

    Yields:
         :class:`ParsedAnnotationRecord`.
    """
    seq_records = list(SeqIO.parse(genbank_handle, format="genbank"))
    if gbk_type == GenBankParserType.SORTED:
        gene_records = group_gene_records_from_sorted_genbank(seq_records, parse_func)
    else:
        gene_records = group_gene_records_by_locus_tag(seq_records, parse_func)
    yield from gene_records


def group_gene_records_from_sorted_genbank(
    record_iter: Iterable[SeqRecord],
    parse_func: Callable[[GeneFeature], Dict[str, Any]],
) -> Iterable[ParsedAnnotationRecord]:
    """Model 1: position sorted GenBank.

    Args:
        record_iter: Iterable of SeqRecord objects.
        parse_func: Optional parse function implementation.

    Yields:
        :class:`ParsedAnnotationRecord`.
    """
    for seqrecord in record_iter:
        gene = None
        source = None
        genes = []
        for feature in seqrecord.features:

            # try to capture the Source field, if it exists
            if feature.type == MetadataFeatures.SOURCE.value:
                source = feature

            # base case for start; iterate until we find a gene
            if gene is None:
                if feature.type in GeneFeature.types:
                    gene = GeneFeature(feature, seqrecord)
                else:
                    continue
            elif feature.type in GeneFeature.types:  # next gene
                if gene.has_children:
                    gene.finalize()
                    gene = parse_func(gene)
                    genes.append(gene)
                gene = GeneFeature(feature, seqrecord)
            elif feature.type in TranscriptFeature.types:
                gene.add_child(feature)
            elif feature.type in IntervalFeature.types:
                if len(gene.children) == 0:
                    gene.add_child(feature)
                else:
                    gene.children[-1].add_child(feature)

        # gene could be None if this record has no annotations
        if gene is not None and gene.has_children:
            gene.finalize()
            gene = parse_func(gene)
            genes.append(gene)

        if source is not None:
            source_qualifiers = source.qualifiers
        else:
            source_qualifiers = None

        annotation = AnnotationCollectionModel.Schema().load(
            dict(genes=genes, sequence_name=seqrecord.id, start=0, end=len(seqrecord), qualifiers=source_qualifiers)
        )
        yield ParsedAnnotationRecord(annotation=annotation, seqrecord=seqrecord)


def group_gene_records_by_locus_tag(
    record_iter: Iterable[SeqRecord],
    parse_func: Callable[[GeneFeature], Dict[str, Any]],
) -> Iterable[ParsedAnnotationRecord]:
    """Model 2: ``locus_tag`` defined GenBank.

    Args:
        record_iter: Iterable of SeqRecord objects.
        parse_func: Optional parse function implementation.

    Yields:
        :class:`ParsedAnnotationRecord`.
    """
    for seqrecord in record_iter:
        genes = []

        filtered_features = []
        source = None
        for f in seqrecord.features:
            if GenBankFeatures.has_value(f.type) and "locus_tag" in f.qualifiers:
                filtered_features.append(f)
            elif f.type == MetadataFeatures.SOURCE.value:
                source = f

        sorted_filtered_features = sorted(filtered_features, key=lambda f: f.qualifiers["locus_tag"])

        for locus_tag, gene_features in itertools.groupby(
            sorted_filtered_features, key=lambda f: f.qualifiers["locus_tag"][0]
        ):
            gene_features = list(gene_features)
            gene = gene_features[0]

            if gene.type not in GeneFeature.types:
                raise GenBankValidationError("Grouping by locus tag produced a mis-ordered interpretation")
            gene = GeneFeature(gene, seqrecord)

            for feature in gene_features[1:]:
                if feature.type in GeneFeature.types:
                    raise GenBankValidationError("Grouping by locus tag found two genes")
                elif feature.type in TranscriptFeature.types:
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

        annotation = AnnotationCollectionModel.Schema().load(
            dict(
                genes=genes,
                feature_collections=None,
                name=seqrecord.id,
                sequence_name=seqrecord.id,
                start=0,
                end=len(seqrecord),
                qualifiers=source_qualifiers,
            )
        )
        yield ParsedAnnotationRecord(annotation=annotation, seqrecord=seqrecord)
