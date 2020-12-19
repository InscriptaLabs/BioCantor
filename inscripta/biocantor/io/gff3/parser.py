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
from io import StringIO
from typing import Iterable, List, Optional, Callable, TextIO

import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils.feature import Feature
from gffutils.interface import FeatureDB

from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.cds import CDSPhase
from inscripta.biocantor.gene.collections import AnnotationCollection
from inscripta.biocantor.io.exc import InvalidInputError
from inscripta.biocantor.io.fasta.fasta import fasta_to_parents, seq_to_parent
from inscripta.biocantor.io.gff3.constants import GFF3Headers, BioCantorGFF3ReservedQualifiers, GFF3FeatureTypes
from inscripta.biocantor.io.gff3.exc import GFF3FastaException, EmptyGff3Exception
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.io.parser import ParsedAnnotationRecord
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.sequence.alphabet import Alphabet

logger = logging.getLogger(__name__)


@dataclass
class GffutilsParseArgs:
    """These arguments are passed to gffutils directly."""

    id_spec: Optional[dict] = None
    merge_strategy: Optional[str] = "create_unique"


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
    for chrom in chroms:
        this_genes = []  # holds genes
        for gene in db.region(seqid=chrom, featuretype=["gene"]):
            gene_id = gene.attributes.get("gene_id", [None])[0]
            locus_tag = gene.attributes.get("locus_tag", [None])[0]
            gene_symbol = gene.attributes.get("gene_name", [gene.attributes.get("gene_symbol", None)])[0]
            gene_biotype = gene.attributes.get("gene_biotype", [gene.attributes.get("gene_type", None)])[0]
            gene_qualifiers = {
                x: y for x, y in gene.attributes.items() if not BioCantorGFF3ReservedQualifiers.has_value(x)
            }

            if gene_biotype is not None:
                try:
                    gene_biotype = Biotype[gene_biotype]
                except KeyError:
                    gene_biotype = Biotype["unknown"]
                    gene_qualifiers["provided_biotype"] = gene.attributes["gene_biotype"][0]
            else:
                gene_biotype = Biotype["unknown"]

            transcripts = []
            for i, transcript in enumerate(db.children(gene, level=1)):

                transcript_id = transcript.attributes.get("transcript_id", [None])[0]
                transcript_symbol = transcript.attributes.get(
                    "transcript_name", [gene.attributes.get("transcript_name", None)]
                )[0]
                transcript_qualifiers = {
                    x: y for x, y in transcript.attributes.items() if not BioCantorGFF3ReservedQualifiers.has_value(x)
                }
                transcript_biotype = gene.attributes.get(
                    "transcript_biotype", [gene.attributes.get("transcript_type", None)]
                )[0]

                if transcript_biotype is not None:
                    try:
                        transcript_biotype = Biotype[transcript_biotype]
                    except KeyError:
                        transcript_biotype = Biotype["unknown"]
                        transcript_qualifiers["provided_transcript_biotype"] = gene.attributes["transcript_biotype"][0]
                elif gene_biotype is not None:
                    transcript_biotype = gene_biotype

                if locus_tag is not None:
                    if transcript_id is None:
                        transcript_id = locus_tag
                    if transcript_symbol is None:
                        transcript_symbol = locus_tag

                exons = []
                cds = []
                for feature in db.children(transcript, level=1):
                    if feature.featuretype == GFF3FeatureTypes.EXON.value:
                        exons.append(feature)
                    elif feature.featuretype == GFF3FeatureTypes.CDS.value:
                        cds.append(feature)
                    else:
                        logger.warning(f"Found non CDS/exon child of transcript in feature: {feature}")

                # This gene has only a CDS/exon feature as its direct child
                # therefore, we really have one interval here
                if len(exons) == 0:
                    if transcript.featuretype not in [
                        GFF3FeatureTypes.CDS.value,
                        GFF3FeatureTypes.EXON.value,
                    ]:
                        logger.warning(f"Gene child feature has type {transcript.featuretype}; skipping")
                        continue
                    logger.info(f"gene {gene_id} had no transcript feature")
                    if transcript.featuretype == GFF3FeatureTypes.CDS.value:
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
                    protein_id = None
                else:
                    cds = sorted(cds, key=lambda c: c.start)
                    cds_starts = [x.start - 1 for x in cds]
                    cds_ends = [x.end for x in cds]
                    cds_frames = [CDSPhase.from_int(int(f.frame)).to_frame().name for f in cds]
                    # NCBI encodes protein IDs on the CDS feature
                    protein_id = cds[0].attributes.get("protein_id", [None])[0]

                tx = dict(
                    exon_starts=exon_starts,
                    exon_ends=exon_ends,
                    strand=strand.name,
                    cds_starts=cds_starts,
                    cds_ends=cds_ends,
                    cds_frames=cds_frames,
                    qualifiers=transcript_qualifiers,
                    is_primary_tx=False,
                    transcript_id=transcript_id,
                    transcript_type=transcript_biotype.name,
                    transcript_symbol=transcript_symbol,
                    sequence_name=chrom,
                    protein_id=protein_id,
                )
                transcripts.append(tx)

            gene = dict(
                transcripts=transcripts,
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                locus_tag=locus_tag,
                gene_type=gene_biotype.name,
                qualifiers=gene_qualifiers,
                sequence_name=chrom,
            )

            this_genes.append(gene)

        if len(this_genes) == 0:
            # empty sequence
            continue

        annot = AnnotationCollectionModel.Schema().load(dict(genes=this_genes, sequence_name=chrom))
        yield annot


def parse_standard_gff3(
    gff: str,
    gffutil_parse_args: Optional[GffutilsParseArgs] = GffutilsParseArgs(),
    parse_func: Optional[Callable[[FeatureDB, List[str]], Iterable[AnnotationCollectionModel]]] = default_parse_func,
    gffutil_transform_func: Optional[Callable[[Feature], Feature]] = None,
    db_fn: Optional[str] = ":memory:",
) -> Iterable[AnnotationCollectionModel]:
    """Parses a GFF3 file using gffutils.

    The parameters parse_func, gffutil_parse_args are implemented separately for each data source. A default
    implementation exists in this module.

    Args:
        gff: Path to a GFF. Must be local or HTTPS.
        parse_func: Function that actually converts gffutils to Inscripta representation.
        gffutil_transform_func: Function that transforms feature keys. Can be necessary in
            cases where IDs are not unique.
        gffutil_parse_args: Parsing arguments to pass to gffutils.
        db_fn: Location to write a gffutils database. Defaults to `:memory:`, which means the database will be built
            transiently. This value can be set to a file location if memory is a concern, or if you want to retain
            the gffutils database. It will not be cleaned up.

    Yields:
        :class:`~biocantor.io.models.AnnotationCollectionModel`
    """
    db = gffutils.create_db(gff, db_fn, transform=gffutil_transform_func, **gffutil_parse_args.__dict__)
    if sum(db.count_features_of_type(i) for i in db.featuretypes()) == 0:
        raise EmptyGff3Exception("Parsing this GFF3 led to zero features. Is it empty or corrupted?")
    logger.info(f"Parsed {gff}")
    for i in db.featuretypes():
        logger.info(f"Found feature type {i} with {db.count_features_of_type(i)} features")
    # get the sequences
    chrom_query = db.execute("SELECT DISTINCT seqid FROM features")
    chroms = [x["seqid"] for x in chrom_query]
    logger.info(f"Found {len(chroms)} sequences")
    yield from parse_func(db, chroms)


def extract_seqrecords_from_gff3_fasta(gff3_with_fasta_handle: TextIO) -> List[SeqRecord]:
    """This function is **NOT** a function to apply FASTA information to a GFF3. This function is purely intended
    to extract the FASTA from a combined file and produce SeqRecords. If you want to parse a GFF3 into memory with
    sequences from a combined file, you would do this:

    .. code:: python

        from inscripta.biocantor.io.gff3.parser import gff3_fasta_to_model, parse_gff
        q = list(gff3_fasta_to_model(parse_gff("tests/data/insO_frameshift.gff3"),
            gff3_with_fasta="tests/data/insO_frameshift.gff3"))

    Parses a GFF3 with a FASTA suffix. Will raise an exception if such a suffix is not found.

    This function is specifically for extracting the FASTA information from the GFF3; it does not parse the GFF3 itself.

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


def parse_gff3_embedded_fasta(
    gff3_with_fasta: str,
    gffutil_parse_args: Optional[GffutilsParseArgs] = GffutilsParseArgs(),
    parse_func: Optional[Callable[[FeatureDB, List[str]], Iterable[AnnotationCollectionModel]]] = default_parse_func,
    gffutil_transform_func: Optional[Callable[[Feature], Feature]] = None,
    db_fn: Optional[str] = ":memory:",
) -> Iterable[ParsedAnnotationRecord]:
    """
    Parses a GFF3 with an embedded FASTA. Wraps :meth:`parse_gff()` to produce ``ParsedAnnotationRecord``.

    Args:
        gff3_with_fasta: Path to a GFF3 file with a FASTA suffix.
        parse_func: Function that actually converts gffutils to Inscripta representation.
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

    for annot_collection in parse_standard_gff3(
        gff3_with_fasta, gffutil_parse_args, parse_func, gffutil_transform_func, db_fn
    ):
        yield ParsedAnnotationRecord(
            annotation=annot_collection, seqrecord=seqrecords_dict[annot_collection.sequence_name]
        )


def parse_gff3_fasta(
    gff3: str,
    fasta: str,
    gffutil_parse_args: Optional[GffutilsParseArgs] = GffutilsParseArgs(),
    parse_func: Optional[Callable[[FeatureDB, List[str]], Iterable[AnnotationCollectionModel]]] = default_parse_func,
    gffutil_transform_func: Optional[Callable[[Feature], Feature]] = None,
    db_fn: Optional[str] = ":memory:",
) -> Iterable[ParsedAnnotationRecord]:
    """
    Parses a GFF3 with FASTA. Wraps :meth:`parse_gff()` to produce ``ParsedAnnotationRecord``.

    Args:
        gff3: Path to a GFF3 file.
        fasta: Path to a FASTA file.
        parse_func: Function that actually converts gffutils to Inscripta representation.
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
    seqrecords_dict = {x.id: x for x in SeqIO.parse(fasta, format="fasta")}

    for annot_collection in parse_standard_gff3(gff3, gffutil_parse_args, parse_func, gffutil_transform_func, db_fn):
        try:
            seqrecord = seqrecords_dict[annot_collection.sequence_name]
        except KeyError:
            logger.warning(
                f"Sequence symbol {annot_collection.sequence_name} found in GFF3 but not in FASTA. "
                f"These annotation records will be ignored."
            )
            continue
        yield ParsedAnnotationRecord(annotation=annot_collection, seqrecord=seqrecord)


def gff3_fasta_to_model(
    annotations: Iterable[AnnotationCollectionModel],
    fasta_handle: Optional[TextIO] = None,
    gff3_with_fasta_handle: Optional[TextIO] = None,
    alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED,
) -> Iterable[AnnotationCollection]:
    """Convert a parsed GFF3 file to object model.

    Take a iterator of :class:`~biocantor.io.models.AnnotationCollectionModel` and yield an
    iterable of :class:`~biocantor.gene.collections.AnnotationCollection`.

    This incorporates sequence information on to each :class:`~biocantor.gene.transcript.TranscriptInterval`
    and :class:`~biocantor.gene.feature.FeatureInterval` object, if a ``fasta`` is provided.

    Args:
        annotations: Iterable that comes from a parser function.
        fasta_handle: Open handle to a FASTA file. Mutually exclusive with ``gff3_with_fasta_handle``.
        gff3_with_fasta_handle: Open handle to a GFF3 file with a FASTA suffix. Mutually exclusive with ``fasta``.
        alphabet: The alphabet to use.

    Yields:
        :class:`~biocantor.gene.collections.AnnotationCollection` with or without sequence information.
    """
    if fasta_handle and gff3_with_fasta_handle:
        raise InvalidInputError("Cannot pass both and gff3_with_fasta to the backend.")

    if fasta_handle:
        parents = fasta_to_parents(fasta_handle, alphabet)
    elif gff3_with_fasta_handle:
        parents = {}
        for rec in extract_seqrecords_from_gff3_fasta(gff3_with_fasta_handle):
            parents[rec.id] = seq_to_parent(str(rec.seq), alphabet=alphabet, seq_id=rec.id)
    else:
        parents = None

    # keep track of seqs we saw in order to handle empty sequences
    seen_seqs = set()
    for annotation in annotations:
        if parents:
            parent = parents[annotation.sequence_name]
            seen_seqs.add(annotation.sequence_name)
        else:
            parent = None
        annot_obj = annotation.to_annotation_collection(parent=parent)
        yield annot_obj

    # build empty objects to store sequences that have no annotations
    if parents:
        for sequence_name in parents.keys() - seen_seqs:
            parent = parents[sequence_name]
            annot_obj = AnnotationCollection(sequence_name=sequence_name, parent=parent)
            yield annot_obj
