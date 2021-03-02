"""
Writes GenBank formatted files out, given one or more :class:`~biocantor.gene.collections.AnnotationCollection`.

These models must have sequences to retrieve.

There are two flavors of GenBank supported here:

1. Prokaryotic -- each gene has a direct descendant, which is its intervals. In other words, all annotations come in
pairs where you have ``gene`` followed by ``[CDS, tRNA, rRNA, ...]``.
2. Eukaryotic -- a more standard gene model, where the top level is always called ``gene``, then there is a child
``[mRNA, tRNA, ...]`` and if the case where the child is ``mRNA``, then there are ``CDS`` features.
"""
import warnings
from typing import Iterable, List, Optional, TextIO, Dict, Hashable, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from inscripta.biocantor.gene.collections import AnnotationCollection, GeneInterval, FeatureIntervalCollection
from inscripta.biocantor.gene.feature import FeatureInterval
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.io.exc import StrandViolationWarning
from inscripta.biocantor.io.genbank.constants import (
    GeneFeatures,
    TranscriptFeatures,
    GeneIntervalFeatures,
    GenbankFlavor,
    FeatureCollectionFeatures,
    FeatureIntervalFeatures,
)
from inscripta.biocantor.io.genbank.exc import GenBankExportError
from inscripta.biocantor.location.strand import Strand


def collection_to_genbank(
    collections: Iterable[AnnotationCollection],
    genbank_file_handle_or_path: Union[TextIO, str],
    genbank_type: Optional[GenbankFlavor] = GenbankFlavor.PROKARYOTIC,
    force_strand: Optional[bool] = True,
    organism: Optional[str] = None,
    source: Optional[str] = None,
):
    """
    Take an instantiated :class:`~biocantor.gene.collections.AnnotationCollection` and produce a GenBank file.

    Args:
        collections: Iterable of AnnotationCollections. They must have sequences associated with them.∂
        genbank_file_handle_or_path: Open file handle or path to write GenBank file to.
        genbank_type: Are we writing an prokaryotic or eukaryotic style GenBank file?
        force_strand: Boolean flag; if ``True``, then strand on children is forced, if ``False``, then improper
            strands are instead skipped.
        organism: What string to put in the ORGANISM field? If not set, will be a period.
        source: What string to put in the SOURCE field? If not set, will be the basename of the GenBank path.
    """

    if organism is None:
        organism = "."

    seqrecords = []
    for collection in collections:

        if collection.sequence is None:
            raise GenBankExportError("Cannot export GenBank if collections do not have sequence information")

        seqrecord = SeqRecord(
            Seq(str(collection.sequence)),
            name=collection.sequence_name,
            id=collection.sequence_guid if collection.sequence_guid else collection.sequence_name,
            description="GenBank produced by BioCantor",
        )
        seqrecord.annotations["molecule_type"] = "DNA"
        seqrecord.annotations["source"] = source
        seqrecord.annotations["organism"] = organism

        for gene_or_feature in collection:
            seqrecord.features.extend(gene_to_feature(gene_or_feature, genbank_type, force_strand))

        seqrecords.append(seqrecord)

    SeqIO.write(seqrecords, genbank_file_handle_or_path, format="genbank")


def gene_to_feature(
    gene_or_feature: Union[GeneInterval, FeatureIntervalCollection], genbank_type: GenbankFlavor, force_strand: bool
) -> Iterable[SeqFeature]:
    """Converts either a :class:`~biocantor.gene.collections.GeneInterval` or a
    :class:`~biocantor.gene.collections.FeatureIntervalCollection` to a :class:`Bio.SeqFeature.SeqFeature`.

    :class:`Bio.SeqFeature.SeqFeature` are BioPython objects that will then be used to write to a GenBank file. There
    is one :class:`Bio.SeqFeature.SeqFeature` for every feature, or row group, in the output file. There will be one
    contiguous interval at the Gene level.

    While :class:`~biocantor.gene.collections.GeneInterval` always has its interval on the plus strand,
    GenBank files assume that a Gene has an explicit strand. Therefore, this function picks the most common strand
    and forces it on all of its children.

    Args:
        gene_or_feature: A :class:`~biocantor.gene.collections.GeneInterval` or
            :class:`~biocantor.gene.collections.FeatureIntervalCollection`.
        genbank_type: Are we writing an prokaryotic or eukaryotic style GenBank file?
        force_strand: Boolean flag; if ``True``, then strand on children is forced, if ``False``, then improper
            strands are instead skipped.

    Yields:
        ``SeqFeature``s, one for the gene, one for each child transcript, and one for each transcript's CDS if it
            exists.
    """
    location = gene_or_feature._location.to_biopython()
    # update the strand by picking the most common
    strands = [child.strand for child in gene_or_feature]
    strand = max(strands, key=strands.count)

    qualifiers = {key: list(vals) for key, vals in gene_or_feature.export_qualifiers().items()}

    # do our best to ensure there is a /gene tag
    symbol = None
    if isinstance(gene_or_feature, GeneInterval):
        if gene_or_feature.gene_symbol:
            symbol = gene_or_feature.gene_symbol
        elif gene_or_feature.gene_id:
            symbol = gene_or_feature.gene_id

        feature_type = GeneFeatures.GENE.value

    else:
        if gene_or_feature.feature_collection_name:
            symbol = gene_or_feature.feature_collection_name
        elif gene_or_feature.feature_collection_id:
            symbol = gene_or_feature.feature_collection_id

        feature_type = FeatureCollectionFeatures.FEATURE_COLLECTION.value

    qualifiers[feature_type] = [symbol]
    feature = SeqFeature(location, type=feature_type, strand=strand.value)
    feature.qualifiers = qualifiers

    yield feature

    if isinstance(gene_or_feature, GeneInterval):
        yield from transcripts_to_feature(
            gene_or_feature.transcripts, strand, genbank_type, force_strand, symbol, gene_or_feature.locus_tag
        )
    else:
        yield from feature_intervals_to_features(
            gene_or_feature.feature_intervals, strand, force_strand, symbol, gene_or_feature.locus_tag
        )


def transcripts_to_feature(
    transcripts: List[TranscriptInterval],
    strand: Strand,
    genbank_type: GenbankFlavor,
    force_strand: bool,
    gene_symbol: Optional[str] = None,
    locus_tag: Optional[str] = None,
) -> Iterable[SeqFeature]:
    """Converts a :class:`~biocantor.gene.transcripts.TranscriptInterval` to a :class:`Bio.SeqFeature.SeqFeature`.

    :class:`Bio.SeqFeature.SeqFeature` are BioPython objects that will then be used to write to a GenBank file. There
    is one :class:`Bio.SeqFeature.SeqFeature` for every feature, or row group, in the output file. There will be one
    joined interval at the transcript level representing the exonic structure.

    While transcript members of a gene can have different strands, for GenBank files that is not allowed. This function
    will explicitly force the strand and provide a warning that this is happening.

    In eukaryotic mode, this function will create mRNA features for coding genes, and biotype features for non-coding.
    Coding genes are then passed on to create CDS features.

    In prokaryotic mode, this function will only create biotype features for non-coding genes.

    Args:
        transcripts: A list of :class:`~biocantor.gene.transcript.TranscriptInterval`.
        strand: ``Strand`` that this gene lives on.
        genbank_type: Are we writing an prokaryotic or eukaryotic style GenBank file?
        force_strand: Boolean flag; if ``True``, then strand is forced, if ``False``, then improper strands are instead
            skipped.
        gene_symbol: An optional gene symbol.
        locus_tag: An optional locus tag.

    Yields:
        ``SeqFeature``s, one for each transcript and then one for each CDS of the transcript, if it exists.
    """
    for transcript in transcripts:
        location = transcript._location.to_biopython()

        transcript_qualifiers = {key: list(vals) for key, vals in transcript.export_qualifiers().items()}
        if gene_symbol is not None:
            transcript_qualifiers["gene"] = [gene_symbol]
        if locus_tag is not None:
            transcript_qualifiers["locus_tag"] = [locus_tag]

        if location.strand != strand.value:
            warn_str = f"Found strand mismatch between gene and transcript on transcript {transcript}. "
            if force_strand:
                warn_str += "Forcing this transcript to the gene orientation."
                warnings.warn(warn_str, StrandViolationWarning)
            else:
                warn_str += "Skipping this transcript."
                warnings.warn(warn_str, StrandViolationWarning)
                continue

        if transcript.transcript_type is not None and TranscriptFeatures.has_value(transcript.transcript_type.name):
            feat_type = TranscriptFeatures(transcript.transcript_type.name)
        # biotypes might be wrong, only trust the CDS interval
        elif transcript.is_coding:
            feat_type = TranscriptFeatures.CODING_TRANSCRIPT
        else:
            feat_type = TranscriptFeatures.MISC_RNA

        if feat_type == TranscriptFeatures.CODING_TRANSCRIPT and genbank_type == GenbankFlavor.PROKARYOTIC:
            # this is a coding gene in prokaryotic mode; skip straight to CDS
            yield add_cds_feature(transcript, transcript_qualifiers, strand)
        else:
            # build this feature; it could be a mRNA for eukaryotic, or non-coding for either prokaryotic or eukaryotic
            feature = SeqFeature(location, type=feat_type.value, strand=strand.value)
            feature.qualifiers = transcript_qualifiers.copy()

            # NCBI does not like protein_id on transcript level features
            if "protein_id" in feature.qualifiers:
                del feature.qualifiers["protein_id"]

            yield feature
            # only in eukaryotic mode for coding genes do we add a third layer
            if genbank_type == GenbankFlavor.EUKARYOTIC and feat_type == TranscriptFeatures.CODING_TRANSCRIPT:
                yield add_cds_feature(transcript, transcript_qualifiers, strand)


def add_cds_feature(
    transcript: TranscriptInterval,
    transcript_qualifiers: Dict[Hashable, List[Hashable]],
    strand: Strand,
) -> SeqFeature:
    """
    Converts a :class:`~biocantor.gene.transcript.TranscriptInterval` that has a CDS to a
    :class:`Bio.SeqFeature.SeqFeature`. that represents the spliced CDS interval.

    Args:
        transcript: A :class:`~biocantor.gene.transcript.TranscriptInterval`.
        strand: ``Strand`` that this transcript lives on.
        transcript_qualifiers: Qualifiers dictionary from the transcript level feature.

    Returns:
        ``SeqFeature`` for the CDS of this transcript.
    """
    location = transcript.cds._location.to_biopython()
    feature = SeqFeature(location, type=GeneIntervalFeatures.CDS.value, strand=strand.value)
    feature.qualifiers = transcript_qualifiers

    # if the sequence has N's, we cannot translate
    try:
        feature.qualifiers["translation"] = [str(transcript.get_protein_sequence())]
    except ValueError:
        pass

    return feature


def feature_intervals_to_features(
    features: List[FeatureInterval],
    strand: Strand,
    force_strand: bool,
    feature_name: Optional[str] = None,
    locus_tag: Optional[str] = None,
) -> Iterable[SeqFeature]:
    """Converts a :class:`~biocantor.gene.feature.FeatureInterval` to a :class:`Bio.SeqFeature.SeqFeature`.

    :class:`Bio.SeqFeature.SeqFeature` are BioPython objects that will then be used to write to a GenBank file. There
    is one :class:`Bio.SeqFeature.SeqFeature` for every feature, or row group, in the output file. There will be one
    joined interval at the transcript level representing the exonic structure.

    While transcript members of a gene can have different strands, for GenBank files that is not allowed. This function
    will explicitly force the strand and provide a warning that this is happening.

    Args:
        features: A list of :class:`~biocantor.gene.feature.TranscriptInterval`.
        strand: ``Strand`` that this gene lives on.
        force_strand: Boolean flag; if ``True``, then strand is forced, if ``False``, then improper strands are instead
            skipped.
        feature_name: An optional feature name.
        locus_tag: An optional locus tag.

    Yields:
        A ``SeqFeature``s for each feature.
    """
    for feature in features:
        location = feature._location.to_biopython()

        feature_qualifiers = {key: list(vals) for key, vals in feature.export_qualifiers().items()}
        if feature_name:
            feature_qualifiers["gene"] = [feature_name]
        if locus_tag:
            feature_qualifiers["locus_tag"] = [locus_tag]

        if location.strand != strand.value:
            warn_str = f"Found strand mismatch between gene and feature on feature {feature}. "
            if force_strand:
                warn_str += "Forcing this transcript to the gene orientation."
                warnings.warn(warn_str, StrandViolationWarning)
            else:
                warn_str += "Skipping this transcript."
                warnings.warn(warn_str, StrandViolationWarning)
                continue

        feature = SeqFeature(location, type=FeatureIntervalFeatures.FEATURE_INTERVAL.value, strand=strand.value)
        feature.qualifiers = feature_qualifiers.copy()

        yield feature
