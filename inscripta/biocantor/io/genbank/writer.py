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
from typing import Iterable, List, Optional, TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from inscripta.biocantor.gene.collections import AnnotationCollection, GeneInterval
from inscripta.biocantor.gene.transcript import TranscriptInterval
from inscripta.biocantor.io.genbank.constants import GeneFeatures, TranscriptFeatures, IntervalFeatures, GenbankFlavor
from inscripta.biocantor.io.genbank.exc import GenBankExportError
from inscripta.biocantor.location.strand import Strand


class StrandViolationWarning(UserWarning):
    pass


def collection_to_genbank(
    collections: Iterable[AnnotationCollection],
    genbank_file_handle: TextIO,
    genbank_type: Optional[GenbankFlavor] = GenbankFlavor.PROKARYOTIC,
    force_strand: Optional[bool] = True,
    organism: Optional[str] = None,
    source: Optional[str] = None,
):
    """
    Take an instantiated :class:`~biocantor.gene.collections.AnnotationCollection` and produce a GenBank file.

    Args:
        collections: Iterable of AnnotationCollections. They must have sequences associated with them.∂
        genbank_file_handle: Path to write GenBank file to.
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

        for gene in collection.genes:
            seqrecord.features.extend(gene_to_feature(gene, genbank_type, force_strand))

        seqrecords.append(seqrecord)

    SeqIO.write(seqrecords, genbank_file_handle, format="genbank")


def gene_to_feature(gene: GeneInterval, genbank_type: GenbankFlavor, force_strand: bool) -> Iterable[SeqFeature]:
    """Converts a :class:`~biocantor.gene.collections.GeneInterval` to a :class:`Bio.SeqFeature.SeqFeature`.

    :class:`Bio.SeqFeature.SeqFeature` are BioPython objects that will then be used to write to a GenBank file. There
    is one :class:`Bio.SeqFeature.SeqFeature` for every feature, or row group, in the output file. There will be one
    contiguous interval at the Gene level.

    While :class:`~biocantor.gene.collections.GeneInterval` always has its interval on the plus strand,
    GenBank files assume that a Gene has an explicit strand. Therefore, this function picks the most common strand
    and forces it on all of its children.

    Args:
        gene: A :class:`~biocantor.gene.collections.GeneInterval`.
        genbank_type: Are we writing an prokaryotic or eukaryotic style GenBank file?
        force_strand: Boolean flag; if ``True``, then strand on children is forced, if ``False``, then improper
            strands are instead skipped.

    Yields:
        ``SeqFeature``s, one for the gene, one for each child transcript, and one for each transcript's CDS if it
            exists.
    """
    location = gene.location.to_biopython()
    # update the strand by picking the most common
    strands = [tx.strand for tx in gene.transcripts]
    strand = max(strands, key=strands.count)
    feature = SeqFeature(location, type=GeneFeatures.GENE.value, strand=strand.value)

    feature.qualifiers.update({key: val for key, val in gene.qualifiers.items() if val is not None})
    feature.qualifiers.update({key: [val] for key, val in gene.identifiers_dict.items() if val is not None})

    if gene.gene_type:
        feature.qualifiers["gene_biotype"] = [gene.gene_type.name]

    gene_symbol = None
    if gene.gene_symbol is not None:
        gene_symbol = gene.gene_symbol
    elif gene.gene_id is not None:
        gene_symbol = gene.gene_id
    feature.qualifiers[GeneFeatures.GENE.value] = [gene_symbol]

    yield feature

    yield from transcripts_to_feature(gene.transcripts, strand, genbank_type, force_strand, gene_symbol, gene.locus_tag)


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
        transcripts: A list of :class:`~biocantor.gene.collections.TranscriptInterval`.
        strand: ``Strand`` that this gene lives on.
        genbank_type: Are we writing an prokaryotic or eukaryotic style GenBank file?
        force_strand: Boolean flag; if ``True``, then strand is forced, if ``False``, then improper strands are instead
            skipped.
        gene_symbol: An optional gene symbol. Required if ``iep_v1_compatible`` is set to ``True``.
        locus_tag: An optional locus tag.

    Yields:
        ``SeqFeature``s, one for each transcript and then one for each CDS of the transcript, if it exists.
    """
    for transcript in transcripts:
        location = transcript.location.to_biopython()

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
            yield add_cds_feature(transcript, strand)
        else:
            # build this feature; it could be a mRNA for eukaryotic, or non-coding for either prokaryotic or eukaryotic
            feature = SeqFeature(location, type=feat_type.value, strand=strand.value)
            feature.qualifiers.update({key: val for key, val in transcript.qualifiers.items() if val is not None})
            feature.qualifiers.update(
                {
                    key: [val]
                    for key, val in transcript.identifiers_dict.items()
                    if val is not None and key != "protein_id"
                }
            )

            if gene_symbol is not None:
                feature.qualifiers["gene"] = [gene_symbol]
            if locus_tag is not None:
                feature.qualifiers["locus_tag"] = [locus_tag]

            yield feature
            # only in eukaryotic mode for coding genes do we add a third layer
            if genbank_type == GenbankFlavor.EUKARYOTIC and feat_type == TranscriptFeatures.CODING_TRANSCRIPT:
                yield add_cds_feature(transcript, strand, gene_symbol, locus_tag)


def add_cds_feature(
    transcript: TranscriptInterval,
    strand: Strand,
    gene_symbol: Optional[str] = None,
    locus_tag: Optional[str] = None,
) -> SeqFeature:
    """
    Converts a :class:`~biocantor.gene.transcript.TranscriptInterval` that has a CDS to a
    :class:`Bio.SeqFeature.SeqFeature`. that represents the spliced CDS interval.

    Args:
        transcript: A :class:`~biocantor.gene.transcript.TranscriptInterval`.
        strand: ``Strand`` that this transcript lives on.
        gene_symbol: An optional gene symbol. Required if ``iep_v1_compatible`` is set to ``True``.
        locus_tag: An optional locus tag.

    Returns:
        ``SeqFeature`` for the CDS of this transcript.
    """
    quals = {key: val for key, val in transcript.qualifiers.items() if val is not None}
    quals.update({key: [val] for key, val in transcript.identifiers_dict.items() if val is not None})
    location = transcript.cds.location.to_biopython()
    feature = SeqFeature(location, type=IntervalFeatures.CDS.value, strand=strand.value)
    feature.qualifiers = quals

    # if the sequence has N's, we cannot translate
    try:
        feature.qualifiers["translation"] = [str(transcript.get_protein_sequence())]
    except ValueError:
        pass

    if gene_symbol is not None:
        feature.qualifiers["gene"] = [gene_symbol]
    if locus_tag is not None:
        feature.qualifiers["locus_tag"] = [locus_tag]

    return feature
