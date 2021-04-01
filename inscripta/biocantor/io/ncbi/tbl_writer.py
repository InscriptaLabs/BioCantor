"""
Write BioCantor data models to the NCBI .tbl format.

The .tbl format is used for NCBI genome submission, and can be validated with the tool ``tbl2asn``.
"""
import itertools
import random
import re
import warnings
from abc import ABC
from string import ascii_uppercase
from typing import Optional, TextIO, Iterable, Union, Dict, List, Set, Hashable

from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.codon import TranslationTable
from inscripta.biocantor.gene.collections import AnnotationCollection, GeneInterval, TranscriptInterval
from inscripta.biocantor.io.genbank.constants import (
    GeneFeatures,
    TranscriptFeatures,
    GeneIntervalFeatures,
    GenbankFlavor,
)
from inscripta.biocantor.io.ncbi.exc import TblExportException
from inscripta.biocantor.location import Location
from inscripta.biocantor.location.location_impl import CompoundInterval
from inscripta.biocantor.location.strand import Strand


def random_uppercase_str(size=10) -> str:
    """
    Generates a random uppercase string of size ``size``.
    Args:
        size: Size of string to produce

    Returns:
        A random string of size ``size``.
    """
    return "".join([random.choice(ascii_uppercase) for _ in range(size)])


class TblFeature(ABC):
    """
    Models one feature in a tbl file.

    tbl is a funky format with five tab delimited columns, separated with FASTA-like headers. It is probably best
    to show an example:

    .. code-block::
        >Feature gb|CM021127.1|
        <14406  14026   gene
                                gene    TDA8
                                locus_tag       GI527_G0000001
                                gene_synonym    YAL064C-A
                                gene_synonym    YAL065C-A
                                db_xref GeneID:851234
        <14406  14393   mRNA
        14390   14382
        14380   14026
                                product Tda8p
                                note    R64_transcript_id: NM_001180041.1
                                exception       low-quality sequence region
                                protein_id      gnl|WGS:JAAEAL|T0000001_1_prot
                                transcript_id   gnl|WGS:JAAEAL|T0000001
                                gene    TDA8
                                locus_tag       GI527_G0000001
                                gene_synonym    NM_001180041.1
                                gene_synonym    YAL065C-A
                                gene_synonym    YAL064C-A
                                db_xref GeneID:851234


    The rows that define intervals have the first 3 tab delimited columns populated, and the last 2 are empty.
    On the other hand, rows that define qualifiers have the first 3 tab delimited columns empty, and the last 2
    contain the key-value pairs.

    The caret on the mRNA region above says that the region is incomplete. This must be set if the start codon
    is invalid or the stop codon is invalid, depending on which direction we are translating in.

    The example gene here is on the negative strand, and this is signified with the intervals starting with a larger
    number than the end. The positions are 1 based inclusive. Another way to think of the intervals is that they are
    always 5' to 3'.
    """

    VALID_KEYS: Optional[Set] = None
    FEATURE_TYPE: Optional[Union[TranscriptFeatures, GeneFeatures, GeneIntervalFeatures]] = None
    children: Optional[List["TblFeature"]] = None
    chars_to_remove = re.compile(r"[\[\]\(\);]*")  # characters to remove from qualifier values

    def __init__(
        self,
        location: Location,
        start_is_incomplete: bool,
        end_is_complete: bool,
        is_pseudo: bool,
        qualifiers: Dict[str, List[str]],
        children: Optional[List["TblFeature"]] = None,
    ):
        self.location = location
        self.start_is_incomplete = start_is_incomplete
        self.end_is_complete = end_is_complete
        self.is_pseudo = is_pseudo
        self.qualifiers = qualifiers
        self.children = children if children else []

    def __str__(self):
        return f"{self._location_to_str()}\n{self._qualifiers_to_str()}"

    def __iter__(self):
        yield from self.iter_children()

    def iter_children(self):
        yield self
        for child in self.children:
            yield from child

    def _location_to_str(self) -> str:
        """Extract location blocks. Handle strand."""
        s = [[b.start + 1, b.end] for b in self.location.blocks]

        if self.location.strand == Strand.MINUS:
            # flip every block, then flip the set of blocks
            s = [b[::-1] for b in s][::-1]

        # handle the incomplete sequences now
        if self.start_is_incomplete:
            # 5p is incomplete
            s[0][0] = f"<{s[0][0]}"
        if self.end_is_complete:
            # 3p is incomplete
            s[-1][1] = f">{s[-1][1]}"

        # start building the string representation. Don't forget to include the feature type!
        r = []
        for i, (start, end) in enumerate(s):
            feature_type = self.FEATURE_TYPE.value if i == 0 else ""
            r.append(f"{start}\t{end}\t{feature_type}\t\t")

        return "\n".join(r)

    def _qualifiers_to_str(self) -> str:
        """
        Converts a qualifiers dictionary to TBL representation.

        Qualifiers are encoded as key-value pairs in the 4th and 5th columns of TBL rows, after the rows that
        represent the genomic interval of the feature. Keys with multiple values are represented by
        repeating the key on another line.

        NCBI does not like parenthesis, brackets or semicolons in the values of a qualifier, so these are removed.

        ``pseudo`` is a special qualifier key with no value, that marks the feature as being a pseudogene. This turns
        off many of the classifiers that ``tbl2asn`` applies to the gene, and so should be used carefully.
        """
        s = []
        for key, val in self.qualifiers.items():
            if not val or key not in self.VALID_KEYS:
                continue
            filtered_vals = [str(x) for x in val if x is not None]
            if not filtered_vals:
                continue
            for val in sorted(filtered_vals):
                val = re.sub(self.chars_to_remove, "", val)
                s.append(f"\t\t\t{key}\t{val}")
        if self.is_pseudo:
            s.append("\t\t\tpseudo\t")
        return "\n".join(s)

    @staticmethod
    def extract_dbxref_synonyms(
        parsed_qualifiers: Dict[Hashable, Set[str]],
        tbl_qualifiers: Dict[str, List[str]],
        gene_symbol: Optional[str] = None,
    ):
        """Update ``tbl_qualifiers`` with values from ``parsed_qualifiers`` if they are xrefs or synonyms."""
        for key, value in parsed_qualifiers.items():
            if "synonym" in key:
                if "gene_synonym" not in tbl_qualifiers:
                    tbl_qualifiers["gene_synonym"] = []
                for item in value:
                    if item != gene_symbol:
                        tbl_qualifiers["gene_synonym"].append(item)
            elif key == "db_xref":
                tbl_qualifiers["db_xref"] = list(value)


class GeneTblFeature(TblFeature):
    """
    A gene feature should have a single interval only. It also should have a fairly limited set of qualifiers.

    The BioCantor model has GeneIntervals always on the + strand to account for the possibility of mixed-strand
    children. This is not possible under the NCBI model, so we override the constructor of a generic table feature
    to set the strand.

    The signifiers for starts and ends being incomplete are not applied at the gene level.
    """

    FEATURE_TYPE = GeneFeatures.GENE
    VALID_KEYS = {"gene", "locus_tag", "gene_synonym", "db_xref", "note"}

    def __init__(
        self,
        gene: GeneInterval,
        locus_tag,
    ):
        self.locus_tag = locus_tag

        strands = [tx.strand for tx in gene.transcripts]
        strand = max(strands, key=strands.count)
        if len(set(strands)) != 1:
            warnings.warn(
                f"Mixed strands found for gene {gene.gene_symbol}. Assigning {strand.name} to all child features."
            )

        if gene.is_coding:
            # if this gene is coding, and there is an in-frame stop, flag it as a pseudogene.
            is_pseudo = any(tx.has_in_frame_stop for tx in gene.transcripts)
        else:
            is_pseudo = False

        location = gene.chromosome_location.reset_strand(strand)

        qualifiers = {"gene": [gene.gene_symbol], "locus_tag": [locus_tag], "note": []}

        if gene.locus_tag:
            qualifiers["note"].append(f"original locus_tag: {gene.locus_tag}")

        if qualifiers["gene"][0] is None and qualifiers["locus_tag"][0] is None:
            raise TblExportException("Must be able to provide locus_tag or gene.")
        # locus tag and gene can't be the same thing or tbl2asn reports LocusTagProblem
        elif qualifiers["gene"] == qualifiers["locus_tag"]:
            qualifiers["gene"][0] += "_gene"
        elif qualifiers["gene"][0] is None:
            qualifiers["gene"] = qualifiers["locus_tag"]

        # try to pull out the special qualifiers
        if gene.qualifiers:
            TblFeature.extract_dbxref_synonyms(gene.qualifiers, qualifiers, gene_symbol=qualifiers["gene"][0])

        # if we end up pulling out synonyms, but have no gene symbol, pick the first synonym as the new symbol
        if not gene.gene_symbol and "gene_synonym" in qualifiers:
            qualifiers["gene"] = [qualifiers["gene_synonym"].pop(0)]
            # if we only had one synonym, then remove the key from the qualifiers dictionary
            if not qualifiers["gene_synonym"]:
                del qualifiers["gene_synonym"]

        # incomplete markers are never put on gene features, so these are always false
        super().__init__(
            location, start_is_incomplete=False, end_is_complete=False, is_pseudo=is_pseudo, qualifiers=qualifiers
        )


class MRNATblFeature(TblFeature):
    """
    A mRNA feature.
    """

    FEATURE_TYPE = TranscriptFeatures.CODING_TRANSCRIPT
    VALID_KEYS = GeneTblFeature.VALID_KEYS | {"protein_id", "transcript_id", "product"}

    def __init__(
        self,
        transcript: TranscriptInterval,
        cds_feature: "CDSTblFeature",
    ):
        super().__init__(
            transcript.chromosome_location,
            start_is_incomplete=cds_feature.start_is_incomplete,
            end_is_complete=cds_feature.end_is_complete,
            is_pseudo=cds_feature.is_pseudo,
            qualifiers=cds_feature.qualifiers,
            children=[cds_feature],
        )


class CDSTblFeature(TblFeature):
    """
    A CDS feature.
    """

    FEATURE_TYPE = GeneIntervalFeatures.CDS
    VALID_KEYS = MRNATblFeature.VALID_KEYS | {"codon_start"}

    def __init__(
        self,
        transcript: TranscriptInterval,
        gene_feature: GeneTblFeature,
        submitter_lab_name: str,
        translation_table: TranslationTable,
    ):
        qualifiers = gene_feature.qualifiers.copy()

        if "product" in transcript.qualifiers:
            product = list(transcript.qualifiers["product"])[0]
            # NCBI does not allow underscores in product names
            product = product.replace("_", " ")
        else:
            # if we don't have a product, default to "hypothetical protein"
            product = "hypothetical protein"

        # NCBI also does not seem to like the term 'alpha' or 'alpha-1' as the product name
        # it considers it a hypothetical protein
        # NCBI also requires that the product contain string characters, cannot be only numbers
        if product in ["alpha", "alpha-1"] or not re.search("[A-Za-z-_]", product):
            qualifiers["note"].append(f"original product: {product}")
            product = "hypothetical protein"

        qualifiers["product"] = [product]

        # protein IDs as well as transcript IDs must be in the format gnl|dbname|string, where:
        # gnl: static string
        # dbname: submitter group's name
        # string: a unique string for this protein
        # transcript ID must not match protein ID
        # this code is not currently explicitly checking this, because the chances of a collision with 12 characters
        # is acceptably low.

        qualifiers["protein_id"] = [f"gnl|{submitter_lab_name}|{random_uppercase_str(size=12)}"]
        qualifiers["transcript_id"] = [f"gnl|{submitter_lab_name}|{random_uppercase_str(size=12)}"]

        # if we have protein or transcript IDs, retain them under the note field
        for key, val in [["protein_id", transcript.protein_id], ["transcript_id", transcript.transcript_id]]:
            if val:
                qualifiers["note"].append(f"original {key}: {val}")

        codon_start = next(transcript.cds.frame_iter()).value + 1
        qualifiers["codon_start"] = [codon_start]

        # try to pull out the special qualifiers
        if transcript.qualifiers:
            TblFeature.extract_dbxref_synonyms(
                transcript.qualifiers, qualifiers, gene_symbol=qualifiers.get("gene", [None])[0]
            )
            # apply these to the gene level also
            TblFeature.extract_dbxref_synonyms(
                transcript.qualifiers, gene_feature.qualifiers, gene_symbol=qualifiers.get("gene", [None])[0]
            )

        # hypothetical proteins cannot have /gene tags; items without /gene tags cannot have /gene_synonym tags
        if product == "hypothetical protein":
            for key, d in itertools.product(["gene", "gene_synonym"], [qualifiers, gene_feature.qualifiers]):
                if key in d:
                    del d[key]

        # start codon can look directly at the translation, because we have a codon_start value
        start_is_incomplete = not transcript.cds.has_start_codon_in_specific_translation_table(translation_table)

        # End completeness requires that there be a in-frame stop codon that is in the last mod3 position
        end_is_incomplete = len(transcript.cds) % 3 != (codon_start - 1) or not transcript.cds.has_valid_stop

        super().__init__(
            transcript.cds_location,
            start_is_incomplete=start_is_incomplete,
            end_is_complete=end_is_incomplete,
            is_pseudo=gene_feature.is_pseudo,
            qualifiers=qualifiers,
        )


class NcRNATblFeature(TblFeature):
    """
    A more specific ncRNA feature. Has a ncrna_class that must also be discerned.
    """

    FEATURE_TYPE = TranscriptFeatures.NONCODING_TRANSCRIPT
    VALID_KEYS = GeneTblFeature.VALID_KEYS | {"transcript_id", "ncRNA_class"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneTblFeature):
        qualifiers = gene_feature.qualifiers.copy()
        qualifiers["transcript_id"] = [transcript.transcript_id]
        qualifiers["ncRNA_class"] = [transcript.transcript_type.name]

        super().__init__(
            transcript.chromosome_location,
            start_is_incomplete=False,
            end_is_complete=False,
            is_pseudo=False,
            qualifiers=qualifiers,
        )


class MiscRNATblFeature(TblFeature):
    """
    A generic ncRNA feature. Also provides a shared constructor for all non-coding transcripts.
    """

    FEATURE_TYPE = TranscriptFeatures.MISC_RNA
    VALID_KEYS = GeneTblFeature.VALID_KEYS | {"transcript_id", "product"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneTblFeature):
        qualifiers = gene_feature.qualifiers.copy()
        qualifiers["transcript_id"] = [transcript.transcript_id]

        if "product" in transcript.qualifiers:
            qualifiers["product"] = transcript.qualifiers["product"]
        else:
            qualifiers["product"] = qualifiers["gene"]

        super().__init__(
            transcript.chromosome_location,
            start_is_incomplete=False,
            end_is_complete=False,
            is_pseudo=False,
            qualifiers=qualifiers,
        )


class TRNATblFeature(TblFeature):
    """
    A tRNA feature. Applies ``tRNA-Xxx`` to products if they are not correct.
    """

    FEATURE_TYPE = TranscriptFeatures.TRANSFER_RNA
    VALID_KEYS = GeneTblFeature.VALID_KEYS | {"transcript_id", "product"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneTblFeature):
        qualifiers = gene_feature.qualifiers.copy()
        qualifiers["transcript_id"] = [transcript.transcript_id]

        if "product" in transcript.qualifiers:
            product = str(list(transcript.qualifiers["product"])[0])
            if product.startswith("tRNA-"):
                qualifiers["product"] = [product]
        else:
            qualifiers["product"] = ["tRNA-Xxx"]

        super().__init__(
            transcript.chromosome_location,
            start_is_incomplete=False,
            end_is_complete=False,
            is_pseudo=False,
            qualifiers=qualifiers,
        )


class RRNATblFeature(TblFeature):
    """
    A rRNA feature. rRNA features must have a product, so these are
    """

    FEATURE_TYPE = TranscriptFeatures.RIBOSOMAL_RNA
    VALID_KEYS = GeneTblFeature.VALID_KEYS | {"transcript_id", "product"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneTblFeature):
        qualifiers = gene_feature.qualifiers.copy()
        qualifiers["transcript_id"] = [transcript.transcript_id]

        if "product" in transcript.qualifiers:
            # some tools encode rRNA with underscores, but tbl2asn does not like that
            # TODO: Try to find the complete list of rRNA names that they consider valid
            qualifiers["product"] = [transcript.qualifiers["product"].pop().replace("_", " ")]
        else:
            qualifiers["product"] = ["unknown ribosomal RNA"]

        super().__init__(
            transcript.chromosome_location,
            start_is_incomplete=False,
            end_is_complete=False,
            is_pseudo=False,
            qualifiers=qualifiers,
        )


class TblGene:
    """
    Container class that holds a gene and its descendant features.
    """

    def __init__(
        self,
        gene: GeneInterval,
        submitter_lab_name: str,
        locus_tag: Optional[str] = None,
        translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT,
    ):
        gene = GeneInterval.from_dict(gene.to_dict(), parent_or_seq_chunk_parent=gene.chunk_relative_location.parent)

        # the location of the transcripts and its CDS intervals must be merged because NCBI does not like
        # adjacent blocks. This must be performed before we instantiate GeneTblFeature because performing this process
        # may break translations. We want to capture broken translations in the /pseudo tag, even if the original
        # potentially overlapping representation was able to properly represent the ORF. This is an inherent limitation
        # of the TBL format (as well as the GenBank format).
        for tx in gene.transcripts:
            if isinstance(tx._location, CompoundInterval):
                tx._location = tx._location.optimize_and_combine_blocks()
            if gene.is_coding:
                if isinstance(tx.cds_location, CompoundInterval):
                    tx.cds = tx.cds.optimize_and_combine_blocks()
        self.gene = gene

        self.gene_tbl = GeneTblFeature(self.gene, locus_tag)

        # now build the transcript level feature(s). NCBI assumes that all transcripts have the same biotype
        # as the gene feature, and so we just apply that assumption here.
        # TODO: this will apply /pseudo to a gene with only some pseudogene isoforms. This will need to be fixed
        #    for mammalian genomes
        for tx in gene.transcripts:
            if gene.is_coding:
                # the CDS feature is built first so that the end-completeness status can be assessed
                cds_tbl = CDSTblFeature(tx, self.gene_tbl, submitter_lab_name, translation_table)
                tx_tbl = MRNATblFeature(tx, cds_tbl)
            # try to figure out if this biotype is one of the INSDC biotypes we support, otherwise fall back
            # to just a generic misc_RNA feature
            elif gene.gene_type == Biotype.rRNA:
                tx_tbl = RRNATblFeature(tx, self.gene_tbl)
            elif gene.gene_type == Biotype.tRNA:
                tx_tbl = TRNATblFeature(tx, self.gene_tbl)
            elif gene.gene_type == Biotype.misc_RNA:
                # NCBI has complained about misc_RNA being "not usual"
                tx_tbl = NcRNATblFeature(tx, self.gene_tbl)
            else:
                tx_tbl = NcRNATblFeature(tx, self.gene_tbl)
            self.gene_tbl.children.append(tx_tbl)

    def __iter__(self):
        yield from self.gene_tbl


def collection_to_tbl(
    collections: Iterable[AnnotationCollection],
    tbl_file_handle: TextIO,
    translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT,
    locus_tag_prefix: Optional[str] = None,
    genbank_flavor: Optional[GenbankFlavor] = GenbankFlavor.EUKARYOTIC,
    locus_tag_jump_size: Optional[int] = 5,
    submitter_lab_name: Optional[str] = None,
    random_seed: Optional[int] = None,
):
    """
    Take an iterable of :class:`~biocantor.gene.collections.AnnotationCollection` and produce a TBL file.

    Args:
        collections: Iterable of AnnotationCollections. They must have sequences associated with them.∂
        tbl_file_handle: Path to write TBL file to.
        translation_table: Translation table of this species. Used to determine if the ends of coding genes
            are considered complete or not.
        locus_tag_prefix: Locus tag prefix. If not set, the locus tag field present in the annotation set will be
            used. An exception will be raised if there is no locus tag on any GeneInterval.
        genbank_flavor: NCBI treats TBL files similarly to Genbank files. The primary distinction is that for
            prokaryotic genomes, the ``mRNA`` feature is not produced, and instead a ``CDS`` is a direct child
            of the ``gene`` feature.
        locus_tag_jump_size: NCBI likes to have locus tags jump so that they can be in sort order in subsequent
            genome iterations.
        submitter_lab_name: Name to put in the ``dbname`` section of unique identifiers. This is intended to be a
            string that uniquely labels the submitter lab. If not set, will be a random string.
        random_seed: A seed value for random string generation. Useful for reproducible runs of this function.
    """
    if random_seed:
        random.seed(random_seed)

    if not locus_tag_prefix:
        locus_tag_prefix = random_uppercase_str(size=8)
    if not submitter_lab_name:
        submitter_lab_name = random_uppercase_str(size=8)

    # keep track of the locus tag we last saw
    locus_tag_offset = 0

    for collection in collections:
        print(f">Features {collection.sequence_name}", file=tbl_file_handle)

        if collection.sequence_name is None:
            raise TblExportException("Must have a sequence name for tbl export.")

        for gene in collection.genes:

            locus_tag_offset += locus_tag_jump_size
            locus_tag = f"{locus_tag_prefix}_{locus_tag_offset}"

            tblgene = TblGene(gene, submitter_lab_name, locus_tag, translation_table)
            for obj in tblgene:
                if genbank_flavor == GenbankFlavor.PROKARYOTIC and type(obj) == MRNATblFeature:
                    continue
                print(str(obj), file=tbl_file_handle)
