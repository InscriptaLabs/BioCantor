"""
Parse GFF3 file by wrapping the library :mod:`gffutils`.

Functions that call :mod:`gffutils` directly require that the filepaths be local paths that exist, because gffutils
cannot handle remote streams. Other functions accept any type of open file handle.

This module contains the default parser function :meth:`default_parse_func()`. This function can be over-written
to write custom parsers.

Additionally, the lower-level interface to :mod:`gffutils` can be tweaked by adjusting the :class:`GffutilsParseArgs`
dataclass to adjust the arguments passed to :mod:`gffutils`.
"""
import logging
from dataclasses import dataclass
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Iterable, List, Optional, Callable, TextIO, Dict, Set, Any
import re
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils.feature import Feature
from gffutils.interface import FeatureDB
from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds_frame import CDSPhase
from inscripta.biocantor.io.gff3.constants import (
    GFF3Headers,
    BioCantorGFF3ReservedQualifiers,
    GFF3GeneFeatureTypes,
    BioCantorQualifiers,
    BIOCANTOR_QUALIFIERS_REGEX,
)
from inscripta.biocantor.io.features import extract_feature_types, extract_feature_name_id, merge_qualifiers
from inscripta.biocantor.io.gff3.exc import (
    GFF3FastaException,
    EmptyGFF3Exception,
    GFF3ChildParentMismatchError,
    GFF3LocusTagError,
)
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.location.strand import Strand

logger = logging.getLogger(__name__)


@dataclass
class GffutilsParseArgs:
    """These arguments are passed to gffutils directly."""

    id_spec: Optional[dict] = None
    merge_strategy: Optional[str] = "create_unique"


def filter_and_sort_qualifiers(qualifiers: Dict[str, List[str]]) -> Optional[Dict[str, List[str]]]:
    """Filter out the qualifiers for any terms we have extracted as BioCantor identifiers as well as any
    GFF3 special terms"""
    qualifiers = {
        key: sorted(vals) for key, vals in qualifiers.items() if not re.match(BIOCANTOR_QUALIFIERS_REGEX, key)
    }
    return qualifiers if qualifiers else None


def _parse_genes(chrom: str, db: FeatureDB) -> List[Dict]:
    """
    Parse canonical genes from this database.

    Args:
        chrom: A chromosome to parse.
        db: Database from :mod:`gffutils`.

    Returns:
        A list of nested dictionaries representing all genes on this chromosome.
    """
    parsed_genes = []
    for gene in db.region(
        seqid=chrom, featuretype=[GFF3GeneFeatureTypes.GENE.value, GFF3GeneFeatureTypes.PSEUDOGENE.value]
    ):
        gene_id = gene.attributes.get("gene_id", [None])[0]
        locus_tag = gene.attributes.get("locus_tag", [None])[0]
        gene_symbol = gene.attributes.get("gene_name", [gene.attributes.get("gene_symbol", None)])[0]
        gene_biotype = gene.attributes.get("gene_biotype", [gene.attributes.get("gene_type", None)])[0]
        gene_qualifiers = {x: y for x, y in gene.attributes.items() if not BioCantorGFF3ReservedQualifiers.has_value(x)}

        if Biotype.has_name(gene_biotype):
            gene_biotype = Biotype[gene_biotype]
        elif gene_biotype:
            gene_qualifiers["provided_biotype"] = [gene_biotype]
            gene_biotype = None

        transcripts = []
        for i, transcript in enumerate(db.children(gene, level=1)):

            transcript_id = transcript.attributes.get("transcript_id", [None])[0]
            transcript_symbol = transcript.attributes.get(
                "transcript_name", [gene.attributes.get("transcript_name", None)]
            )[0]
            transcript_qualifiers = {
                x: y for x, y in transcript.attributes.items() if not BioCantorGFF3ReservedQualifiers.has_value(x)
            }
            provided_transcript_biotype = gene.attributes.get(
                "transcript_biotype", [gene.attributes.get("transcript_type", None)]
            )[0]

            if Biotype.has_name(provided_transcript_biotype):
                transcript_biotype = Biotype[provided_transcript_biotype]
            else:
                # keep track of what they gave us, that did not match the enum
                if provided_transcript_biotype:
                    transcript_qualifiers["provided_transcript_biotype"] = provided_transcript_biotype
                # use the gene biotype
                transcript_biotype = gene_biotype

            if locus_tag is not None:
                if transcript_id is None:
                    transcript_id = locus_tag
                if transcript_symbol is None:
                    transcript_symbol = locus_tag

            exons = []
            cds = []
            for feature in db.children(transcript, level=1):
                if feature.featuretype == GFF3GeneFeatureTypes.EXON.value:
                    exons.append(feature)
                elif feature.featuretype == GFF3GeneFeatureTypes.CDS.value:
                    cds.append(feature)
                else:
                    logger.warning(f"Found non CDS/exon child of transcript in feature: {feature}")

            # This gene has only a CDS/exon feature as its direct child
            # therefore, we really have one interval here
            if len(exons) == 0:
                if transcript.featuretype not in [
                    GFF3GeneFeatureTypes.CDS.value,
                    GFF3GeneFeatureTypes.EXON.value,
                ]:
                    logger.warning(f"Gene child feature has type {transcript.featuretype}; skipping")
                    continue
                logger.info(f"gene {gene_id} had no transcript feature")
                if transcript.featuretype == GFF3GeneFeatureTypes.CDS.value:
                    exons = cds = [transcript]
                else:
                    exons = [transcript]

            exons = sorted(exons, key=lambda e: e.start)
            exon_starts = [x.start - 1 for x in exons]
            exon_ends = [x.end for x in exons]
            start = exon_starts[0]
            end = exon_ends[-1]
            assert start <= end
            strand = Strand.from_symbol(transcript.strand)

            if len(cds) == 0:
                cds_starts = cds_ends = cds_frames = None
                protein_id = product = None
            else:
                # sort by start and end in case two blocks start at the same position
                cds = sorted(cds, key=lambda c: (c.start, c.end))
                cds_starts = [x.start - 1 for x in cds]
                cds_ends = [x.end for x in cds]
                cds_frames = [CDSPhase.from_int(int(f.frame)).to_frame().name for f in cds]
                # NCBI encodes protein IDs and products on the CDS feature
                protein_id = cds[0].attributes.get("protein_id", [None])[0]
                product = cds[0].attributes.get("product", [None])[0]

            tx = dict(
                exon_starts=exon_starts,
                exon_ends=exon_ends,
                strand=strand.name,
                cds_starts=cds_starts,
                cds_ends=cds_ends,
                cds_frames=cds_frames,
                qualifiers=filter_and_sort_qualifiers(transcript_qualifiers),
                is_primary_tx=False,
                transcript_id=transcript_id,
                transcript_type=transcript_biotype.name if transcript_biotype else transcript_biotype,
                transcript_symbol=transcript_symbol,
                sequence_name=chrom,
                protein_id=protein_id,
                product=product,
            )
            transcripts.append(tx)

        if len(transcripts) == 0:
            # infer a transcript for a gene
            logger.info(f"Inferring a transcript for gene {gene_symbol}")
            tx = dict(
                exon_starts=[gene.start],
                exon_ends=[gene.end],
                strand=Strand.from_symbol(gene.strand).name,
                qualifiers=gene_qualifiers,
                transcript_type=gene_biotype.name if gene_biotype else gene_biotype,
                transcript_id=gene_id,
                sequence_name=gene.seqid,
            )
            transcripts.append(tx)

        gene = dict(
            transcripts=transcripts,
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            locus_tag=locus_tag,
            gene_type=gene_biotype.name if gene_biotype else gene_biotype,
            qualifiers=filter_and_sort_qualifiers(gene_qualifiers),
            sequence_name=chrom,
        )

        parsed_genes.append(gene)
    return parsed_genes


def _find_all_top_level_non_gene_features(chrom: str, db: FeatureDB, feature_types: List[str]) -> Iterable[Feature]:
    """
    Find all top-level non gene features. GFFutils lacks a way to do this directly, so we just iterate over everything.

    Args:
        chrom: A chromosome to parse.
        db: Database from :mod:`gffutils`.
        feature_types: A set of feature types that are in the database that are not genic.

    Yields:
        Iterable of ``Feature`` objects that are top-level.
    """
    for feature in db.region(seqid=chrom, featuretype=feature_types):
        try:
            _ = next(db.parents(feature.id))
        except StopIteration:
            yield feature


def _parse_child_features_to_feature_interval(
    features: List[Feature], locus_tag: Optional[str] = None
) -> Dict[str, Any]:
    """
    Extract values from a list of child features and produce a dictionary to build a
    :class:`~biocantor.io.models.FeatureIntervalModel` from.

    Can also be provided a ``locus_tag`` value from a parent, if applicable.

    This function combines all child features of a top-level non-gene feature
    """
    # extract feature types, including the base type
    feature_types = {feature.featuretype for feature in features}

    # loop over features to get names, identifiers, positions and qualifiers
    qualifiers = None
    feature_names = Counter()
    feature_ids = Counter()
    strand = None
    chrom = None
    interval_starts = []
    interval_ends = []
    for feature in features:
        if not chrom:
            chrom = feature.chrom
        elif chrom != feature.chrom:
            raise GFF3ChildParentMismatchError("Cannot have multiple child features on different sequences.")

        this_strand = Strand.from_symbol(feature.strand).name
        if not strand:
            strand = this_strand
        elif strand != this_strand:
            raise GFF3ChildParentMismatchError("Cannot have multiple child features on different strands.")

        interval_starts.append(feature.start - 1)
        interval_ends.append(feature.end)

        extract_feature_types(feature_types, feature.attributes)
        feature_name, feature_id = extract_feature_name_id(feature.attributes)
        feature_names[feature_name] += 1
        feature_ids[feature_id] += 1

        if BioCantorQualifiers.LOCUS_TAG.value in feature.attributes:
            this_locus_tag = feature.attributes[BioCantorQualifiers.LOCUS_TAG.value][0]
            if this_locus_tag != locus_tag:
                raise GFF3LocusTagError("Cannot have multiple child features with different locus tags.")

        if not qualifiers:
            qualifiers = feature.attributes
        else:
            qualifiers = merge_qualifiers(qualifiers, feature.attributes)

    if len(feature_names) > 0:
        feature_name = feature_names.most_common(1)[0][0]
    else:
        feature_name = locus_tag

    if len(feature_ids) > 0:
        feature_id = feature_ids.most_common(1)[0][0]
    else:
        feature_id = locus_tag

    return dict(
        interval_starts=sorted(interval_starts),
        interval_ends=sorted(interval_ends),
        strand=strand,
        qualifiers=filter_and_sort_qualifiers(qualifiers),
        feature_id=feature_id,
        feature_name=feature_name,
        feature_types=sorted(feature_types),
        sequence_name=chrom,
        is_primary_feature=False,
    )


def _parse_features(chrom: str, db: FeatureDB, feature_types: List[str]) -> List[Dict]:
    """
    Parse generic features from this database. These are anything that cannot be interpreted as a gene.

    If a feature is a top-level feature with no children, then infer a collection wrapper for it.

    Args:
        chrom: A chromosome to parse.
        db: Database from :mod:`gffutils`.
        feature_types: A set of feature types that are in the database that are not genic.

    Returns:
        A list of nested dictionaries representing non-gene features.
    """
    feature_collections = []
    for top_level_feature in _find_all_top_level_non_gene_features(chrom, db, feature_types):
        children = list(db.children(top_level_feature, level=1))

        # extract parent locus tag to compare to children
        locus_tag = None
        if BioCantorQualifiers.LOCUS_TAG.value in top_level_feature.attributes:
            locus_tag = top_level_feature.attributes[BioCantorQualifiers.LOCUS_TAG.value][0]

        if not children:
            # treat this isolated feature as both FeatureIntervalCollection and FeatureInterval
            feature = _parse_child_features_to_feature_interval([top_level_feature])
            # infer a FeatureCollection from the information on the FeatureInterval
            feature_collection = dict(
                feature_intervals=[feature],
                feature_collection_name=feature["feature_name"],
                feature_collection_id=feature["feature_id"],
                feature_collection_type=top_level_feature.featuretype,
                locus_tag=locus_tag,
                sequence_name=chrom,
                qualifiers=feature["qualifiers"],
            )
            # remove qualifiers from feature
            del feature["qualifiers"]
        else:
            # combine all children into a FeatureInterval
            feature = _parse_child_features_to_feature_interval(children, locus_tag=locus_tag)
            feature_collection_name, feature_collection_id = extract_feature_name_id(top_level_feature.attributes)

            feature_collection = dict(
                feature_intervals=[feature],
                feature_collection_name=feature_collection_name,
                feature_collection_id=feature_collection_id,
                feature_collection_type=top_level_feature.featuretype,
                locus_tag=locus_tag,
                sequence_name=chrom,
                qualifiers=filter_and_sort_qualifiers(top_level_feature.attributes),
            )

        feature_collections.append(feature_collection)
    return feature_collections


def _find_non_gene_feature_types(db: FeatureDB, feature_types_to_ignore: Optional[Set[str]] = None) -> List[str]:
    """Non-gene feature types are those that are not either a member of :class:`~biocantor.gene.biotype.Biotype`
    or :class:`~biocantor.io.gff3.constants.GFF3GeneFeatureTypes`. This combination of filters prevents genes being
    inadvertently pulled in from either of the two main styles of representing them.

    NCBI Style: {gene,pseudogene} -> {mRNA, tRNA, etc} -> {exon, cds}
    Ensembl/GENCODE Style: gene -> transcript -> {exon, cds}

    Args:
        db: Database from :mod:`gffutils`.
        feature_types_to_ignore: Feature types to ignore, if chosen. This is often used to ignore pointless features
            like chromosome representations.

    Returns:
        A list of strings representing the non-gene feature types found in the database.
    """
    non_gene_feature_types = []
    for feat_type in db.featuretypes():
        if GFF3GeneFeatureTypes.has_value(feat_type):
            continue
        elif feature_types_to_ignore and feat_type in feature_types_to_ignore:
            continue
        else:
            try:
                _ = Biotype[feat_type]
            except KeyError:
                non_gene_feature_types.append(feat_type)
    return non_gene_feature_types


def default_parse_func(db: FeatureDB, chroms: List[str]) -> Iterable[AnnotationCollectionModel]:
    """
    This is the default parser function. Mappings include:

    gene_id -> gene_id
    gene_name or if missing gene_symbol -> gene_symbol
    gene_biotype or if missing gene_type -> gene_biotype
    transcript_id -> transcript_id
    transcript_name or if missing transcript_name -> transcript_symbol
    transcript_biotype or if missing transcript_type -> transcript_biotype
    if no transcript_biotype or transcript_type, then gene result is used

    A list of chromosomes is required in order to allow there to be a specified order of data, otherwise they come back
    unordered from the database.

    Args:
        db: Database from :mod:`gffutils`.
        chroms: List of sequence names to iterate over.

    Yields:
        :class:`~biocantor.io.models.AnnotationCollectionModel`
    """
    non_gene_feature_types = _find_non_gene_feature_types(db)

    for chrom in chroms:
        parsed_genes = _parse_genes(chrom, db)
        if non_gene_feature_types:
            parsed_features = _parse_features(chrom, db, non_gene_feature_types)
        else:
            parsed_features = None

        annot = AnnotationCollectionModel.Schema().load(
            dict(genes=parsed_genes, feature_collections=parsed_features, sequence_name=chrom)
        )
        yield annot


def extract_seqrecords_from_gff3_fasta(gff3_with_fasta_handle: TextIO) -> List[SeqRecord]:
    """This function is **NOT** a function to apply FASTA information to a GFF3. This function is purely intended
    to extract the FASTA from a combined file and produce SeqRecords.

    Parses a GFF3 with a FASTA suffix. Will raise an exception if such a suffix is not found.

    Args:
        gff3_with_fasta_handle: Open file handle in text mode to a GFF3 file with a FASTA suffix.

    Raises:
        GFF3FastaException: if the GFF3 lacks a FASTA suffix.

    Returns:
        List of ``SeqRecord`` objects.
    """
    # first, find the header
    fasta_found = False
    for row in gff3_with_fasta_handle:
        if row.rstrip() == GFF3Headers.FASTA_HEADER.value:
            fasta_found = True
            break

    if fasta_found is False:
        raise GFF3FastaException("Did not find FASTA header in the GFF3 file.")

    data = StringIO(gff3_with_fasta_handle.read())
    recs = list(SeqIO.parse(data, format="fasta"))
    return recs


def parse_standard_gff3(
    gff: Path,
    gffutil_parse_args: Optional[GffutilsParseArgs] = GffutilsParseArgs(),
    parse_func: Optional[Callable[[FeatureDB, List[str]], Iterable[AnnotationCollectionModel]]] = default_parse_func,
    gffutil_transform_func: Optional[Callable[[Feature], Feature]] = None,
    db_fn: Optional[str] = ":memory:",
) -> Iterable[ParsedAnnotationRecord]:
    """Parses a GFF3 file using gffutils.

    The parameters parse_func, gffutil_parse_args are implemented separately for each data source. A default
    implementation exists in this module.

    Args:
        gff: Path to a GFF. Must be local or HTTPS.
        parse_func: Function that actually converts gffutils to BioCantor representation.
        gffutil_transform_func: Function that transforms feature keys. Can be necessary in cases where IDs are not
            unique.
        gffutil_parse_args: Parsing arguments to pass to gffutils.
        db_fn: Location to write a gffutils database. Defaults to `:memory:`, which means the database will be built
            transiently. This value can be set to a file location if memory is a concern, or if you want to retain
            the gffutils database. It will not be cleaned up.

    Yields:
        Iterable of ``ParsedAnnotationRecord`` objects.
    """
    db = gffutils.create_db(str(gff), db_fn, transform=gffutil_transform_func, **gffutil_parse_args.__dict__)
    if sum(db.count_features_of_type(i) for i in db.featuretypes()) == 0:
        raise EmptyGFF3Exception("Parsing this GFF3 led to zero features. Is it empty or corrupted?")
    logger.info(f"Parsed {gff}")
    for i in db.featuretypes():
        logger.info(f"Found feature type {i} with {db.count_features_of_type(i)} features")
    # get the sequences
    chrom_query = db.execute("SELECT DISTINCT seqid FROM features")
    chroms = [x["seqid"] for x in chrom_query]
    logger.info(f"Found {len(chroms)} sequences")
    for annot in parse_func(db, chroms):
        yield ParsedAnnotationRecord(annot)


def _produce_empty_records(
    seqrecords_dict: Dict[str, SeqRecord], seen_seqs: Set[str]
) -> Iterable[ParsedAnnotationRecord]:
    """
    Convenience function shared by :meth:`parse_gff3_embedded_fasta()` and :meth:`parse_gff3_fasta()` that appends
    empty ``ParsedAnnotationRecord`` objects to the end. This ensures that every sequence in the FASTA is still
    represented in the final object set, even if it has zero annotations.

    Args:
        seqrecords_dict: Dictionary mapping sequence names to SeqRecord objects.
        seen_seqs: Set of sequences that were found when parsing the GFF3.

    Yields:
        Iterable of ``ParsedAnnotationRecord`` objects with empty annotations.
    """
    for sequence_name in seqrecords_dict.keys() - seen_seqs:
        seqrecord = seqrecords_dict[sequence_name]
        annot = AnnotationCollectionModel.Schema().load(dict(sequence_name=seqrecord.id, start=0, end=len(seqrecord)))
        yield ParsedAnnotationRecord(annotation=annot, seqrecord=seqrecord)


def parse_gff3_embedded_fasta(
    gff3_with_fasta: Path,
    gffutil_parse_args: Optional[GffutilsParseArgs] = GffutilsParseArgs(),
    parse_func: Optional[Callable[[FeatureDB, List[str]], Iterable[AnnotationCollectionModel]]] = default_parse_func,
    gffutil_transform_func: Optional[Callable[[Feature], Feature]] = None,
    db_fn: Optional[str] = ":memory:",
) -> Iterable[ParsedAnnotationRecord]:
    """
    Parses a GFF3 with an embedded FASTA. Wraps :meth:`parse_gff()` to produce ``ParsedAnnotationRecord``.

    Args:
        gff3_with_fasta: Path to a GFF3 file with a FASTA suffix. Must be local or HTTPS.
        parse_func: Function that actually converts gffutils to BioCantor representation.
        gffutil_transform_func: Function that transforms feature keys. Can be necessary in
            cases where IDs are not unique.
        gffutil_parse_args: Parsing arguments to pass to gffutils.
        db_fn: Location to write a gffutils database. Defaults to `:memory:`, which means the database will be built
            transiently. This value can be set to a file location if memory is a concern, or if you want to retain
            the gffutils database. It will not be cleaned up.

    Raises:
        GFF3FastaException: if the GFF3 lacks a FASTA suffix.

    Yields:
        Iterable of ``ParsedAnnotationRecord`` objects.
    """
    with open(gff3_with_fasta, "r") as fh:
        seqrecords = extract_seqrecords_from_gff3_fasta(fh)
    seqrecords_dict = {x.id: x for x in seqrecords}

    # keep track of sequences we see in the GFF3 so we can append empty records (contigs with no annotations)
    seen_seqs = set()

    for annot_record in parse_standard_gff3(
        gff3_with_fasta, gffutil_parse_args, parse_func, gffutil_transform_func, db_fn
    ):
        annot_record.seqrecord = seqrecords_dict[annot_record.annotation.sequence_name]
        yield annot_record
        seen_seqs.add(annot_record.annotation.sequence_name)
    yield from _produce_empty_records(seqrecords_dict, seen_seqs)


def parse_gff3_fasta(
    gff3: Path,
    fasta: Path,
    gffutil_parse_args: Optional[GffutilsParseArgs] = GffutilsParseArgs(),
    parse_func: Optional[Callable[[FeatureDB, List[str]], Iterable[AnnotationCollectionModel]]] = default_parse_func,
    gffutil_transform_func: Optional[Callable[[Feature], Feature]] = None,
    db_fn: Optional[str] = ":memory:",
) -> Iterable[ParsedAnnotationRecord]:
    """
    Parses a GFF3 with a separate FASTA. Wraps :meth:`parse_gff()` to produce ``ParsedAnnotationRecord``.

    Args:
        gff3: Path to a GFF3 file. Must be local or HTTPS.
        fasta: Path to a FASTA file. Must be local or HTTPS.
        parse_func: Function that actually converts gffutils to BioCantor representation.
        gffutil_transform_func: Function that transforms feature keys. Can be necessary in
            cases where IDs are not unique.
        gffutil_parse_args: Parsing arguments to pass to gffutils.
        db_fn: Location to write a gffutils database. Defaults to ``:memory:``, which means the database will be built
            transiently. This value can be set to a file location if memory is a concern, or if you want to retain
            the gffutils database. It will not be cleaned up.

    Raises:
        GFF3FastaException: if the GFF3 lacks a FASTA suffix.

    Yields:
        Iterable of ``ParsedAnnotationRecord`` objects.
    """
    seqrecords_dict = {x.id: x for x in SeqIO.parse(fasta, format="fasta")}

    # keep track of sequences we see in the GFF3 so we can append empty records (contigs with no annotations)
    seen_seqs = set()

    for annot_record in parse_standard_gff3(gff3, gffutil_parse_args, parse_func, gffutil_transform_func, db_fn):
        try:
            seqrecord = seqrecords_dict[annot_record.annotation.sequence_name]
        except KeyError:
            logger.warning(
                f"Sequence symbol {annot_record.annotation.sequence_name} found in GFF3 but not in FASTA. "
                f"These annotation records will be ignored."
            )
            continue
        annot_record.seqrecord = seqrecord
        yield annot_record
        seen_seqs.add(annot_record.annotation.sequence_name)
    yield from _produce_empty_records(seqrecords_dict, seen_seqs)
