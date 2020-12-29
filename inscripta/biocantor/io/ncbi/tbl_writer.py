"""
Write BioCantor data models to the NCBI .tbl format.

The .tbl format is used for NCBI genome submission, and can be validated with the tool ``tbl2asn``.
"""
import warnings
from abc import ABC
from typing import Optional, TextIO, Iterable, Union, Dict, List, Set

from inscripta.biocantor.gene.biotype import Biotype
from inscripta.biocantor.gene.codon import TranslationTable
from inscripta.biocantor.gene.collections import AnnotationCollection, GeneInterval, TranscriptInterval
from inscripta.biocantor.io.genbank.constants import GeneFeatures, TranscriptFeatures, IntervalFeatures, GenbankFlavor
from inscripta.biocantor.io.ncbi.exc import TblExportException, LocusTagException
from inscripta.biocantor.location import Location
from inscripta.biocantor.location.strand import Strand


class TblFeature(ABC):
    """
    Models one feature in a tbl file.

    tbl is a funky format with five tab delimited columns, separated with FASTA-like headers. It is probably best
    to show an example:

    ```
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
    ```

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
    FEATURE_TYPE: Optional[Union[TranscriptFeatures, GeneFeatures, IntervalFeatures]] = None
    children: Optional[List["TblFeature"]] = None

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
        s = []
        for key, val in self.qualifiers.items():
            if val is None or key not in self.VALID_KEYS:
                continue
            filtered_vals = [str(x) for x in val if x is not None]
            if len(filtered_vals) == 0:
                continue
            val = "; ".join(filtered_vals)
            s.append(f"\t\t\t{key}\t{val}")
        if self.is_pseudo:
            s.append("\t\t\tpseudo\t")
        return "\n".join(s)


class GeneFeature(TblFeature):
    """
    A gene feature should have a single interval only. It also should have a fairly limited set of qualifiers.

    The BioCantor model has GeneIntervals always on the + strand to account for the possibility of mixed-strand
    children. This is not possible under the NCBI model, so we override the constructor of a generic table feature
    to set the strand.

    The signifiers for starts and ends being incomplete are not applied at the gene level.
    """

    FEATURE_TYPE = GeneFeatures.GENE
    VALID_KEYS = {"gene", "locus_tag", "gene_synonym", "db_xref"}

    def __init__(
        self,
        gene: GeneInterval,
        locus_tag: Optional[str] = None,
    ):
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

        location = gene.location.reset_strand(strand)

        qualifiers = {"gene": [gene.gene_symbol], "locus_tag": [locus_tag] if locus_tag else [gene.locus_tag]}

        # try to pull out the special qualifiers
        if gene.qualifiers:
            for key, value in gene.qualifiers.items():
                if "synonym" in key:
                    if "gene_synonym" not in qualifiers:
                        qualifiers["gene_synonym"] = []
                    qualifiers["gene_synonym"].extend(value)
                elif key == "db_xref":
                    qualifiers["db_xref"] = value

        self.locus_tag = locus_tag

        # incomplete markers are never put on gene features, so these are always false
        super().__init__(
            location, start_is_incomplete=False, end_is_complete=False, is_pseudo=is_pseudo, qualifiers=qualifiers
        )


class MRnaFeature(TblFeature):
    """
    A mRNA feature.
    """

    FEATURE_TYPE = TranscriptFeatures.CODING_TRANSCRIPT
    VALID_KEYS = GeneFeature.VALID_KEYS | {"protein_id", "transcript_id"}

    def __init__(
        self,
        transcript: TranscriptInterval,
        cds_feature: "CDSFeature",
    ):
        super().__init__(
            transcript.location,
            start_is_incomplete=cds_feature.start_is_incomplete,
            end_is_complete=cds_feature.end_is_complete,
            is_pseudo=cds_feature.is_pseudo,
            qualifiers=cds_feature.qualifiers,
            children=[cds_feature],
        )


class CDSFeature(TblFeature):
    """
    A mRNA feature.
    """

    FEATURE_TYPE = IntervalFeatures.CDS
    VALID_KEYS = MRnaFeature.VALID_KEYS | {"product"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneFeature, translation_table: TranslationTable):
        qualifiers = gene_feature.qualifiers.copy()

        if "product" in transcript.qualifiers:
            qualifiers["product"] = transcript.qualifiers["product"]
        elif transcript.protein_id:
            qualifiers["product"] = qualifiers["protein_id"]
        else:
            qualifiers["product"] = qualifiers["gene"]

        if transcript.protein_id:
            qualifiers["protein_id"] = [transcript.protein_id]
        else:
            qualifiers["protein_id"] = qualifiers["gene"]

        if transcript.transcript_id:
            qualifiers["transcript_id"] = [transcript.transcript_id]

        codons = list(transcript.cds.scan_codons())
        start_is_incomplete = not codons[0].is_start_codon_in_specific_translation_table(translation_table)
        end_is_incomplete = not codons[-1].is_stop_codon
        super().__init__(
            transcript.cds.location,
            start_is_incomplete=start_is_incomplete,
            end_is_complete=end_is_incomplete,
            is_pseudo=gene_feature.is_pseudo,
            qualifiers=qualifiers,
        )


class NcRnaFeature(TblFeature):
    """
    A more specific ncRNA feature. Has a ncrna_class that must also be discerned.
    """

    FEATURE_TYPE = TranscriptFeatures.NONCODING_TRANSCRIPT
    VALID_KEYS = GeneFeature.VALID_KEYS | {"transcript_id", "ncRNA_class"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneFeature):
        qualifiers = gene_feature.qualifiers.copy()
        qualifiers["transcript_id"] = [transcript.transcript_id]
        qualifiers["ncRNA_class"] = [transcript.transcript_type.name]

        super().__init__(
            transcript.location,
            start_is_incomplete=False,
            end_is_complete=False,
            is_pseudo=False,
            qualifiers=qualifiers,
        )


class MiscRnaFeature(TblFeature):
    """
    A generic ncRNA feature. Also provides a shared constructor for all non-coding transcripts.
    """

    FEATURE_TYPE = TranscriptFeatures.MISC_RNA
    VALID_KEYS = GeneFeature.VALID_KEYS | {"transcript_id", "product"}

    def __init__(self, transcript: TranscriptInterval, gene_feature: GeneFeature):
        qualifiers = gene_feature.qualifiers.copy()
        qualifiers["transcript_id"] = [transcript.transcript_id]

        if "product" in transcript.qualifiers:
            qualifiers["product"] = transcript.qualifiers["product"]
        else:
            qualifiers["product"] = qualifiers["gene"]

        super().__init__(
            transcript.location,
            start_is_incomplete=False,
            end_is_complete=False,
            is_pseudo=False,
            qualifiers=qualifiers,
        )


class TRnaFeature(MiscRnaFeature):
    """
    A tRNA feature.
    """

    FEATURE_TYPE = TranscriptFeatures.TRANSFER_RNA


class RRnaFeature(MiscRnaFeature):
    """
    A rRNA feature.
    """

    FEATURE_TYPE = TranscriptFeatures.RIBOSOMAL_RNA


class TblGene:
    """
    Container class that holds a gene and its descendant features.
    """

    def __init__(
        self,
        gene: GeneInterval,
        locus_tag: Optional[str] = None,
        translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT,
    ):
        self.gene = gene
        self.gene_tbl = GeneFeature(self.gene, locus_tag)

        # now build the transcript level feature(s). NCBI assumes that all transcripts have the same biotype
        # as the gene feature, and so we just apply that assumption here.
        # TODO: this will apply /pseudo to a gene with only some pseudogene isoforms. This will need to be fixed
        #    for mammalian genomes
        for tx in gene.transcripts:
            if gene.is_coding:
                # the CDS feature is built first so that the end-completeness status can be assessed
                cds_tbl = CDSFeature(tx, self.gene_tbl, translation_table)
                tx_tbl = MRnaFeature(tx, cds_tbl)
            # try to figure out if this biotype is one of the INSDC biotypes we support, otherwise fall back
            # to just a generic misc_RNA feature
            elif gene.gene_type == Biotype.rRNA:
                tx_tbl = RRnaFeature(tx, self.gene_tbl)
            elif gene.gene_type == Biotype.tRNA:
                tx_tbl = TRnaFeature(tx, self.gene_tbl)
            else:
                tx_tbl = NcRnaFeature(tx, self.gene_tbl)
            self.gene_tbl.children.append(tx_tbl)

    def __iter__(self):
        yield from self.gene_tbl


def collection_to_tbl(
    collections: Iterable[AnnotationCollection],
    tbl_file_handle: TextIO,
    translation_table: Optional[TranslationTable] = TranslationTable.DEFAULT,
    locus_tag_prefix: Optional[str] = None,
    genbank_flavor: Optional[GenbankFlavor] = GenbankFlavor.EUKARYOTIC,
):
    """
    Take an instantiated :class:`~biocantor.gene.collections.AnnotationCollection` and produce a TBL file.

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
    """
    # keep track of the locus tag we last saw
    locus_tag_offset = 1

    for collection in collections:
        print(f">Features {collection.sequence_name}", file=tbl_file_handle)

        if collection.sequence_name is None:
            raise TblExportException("Must have a sequence name for tbl export.")

        for gene in collection.genes:

            if locus_tag_prefix:
                locus_tag = f"{locus_tag_prefix}_{locus_tag_offset:06}"
                locus_tag_offset += 1
            else:
                locus_tag = gene.locus_tag

            if not locus_tag:
                raise LocusTagException("Must explicitly assign locus tags to all genes if not providing a prefix.")

            tblgene = TblGene(gene, locus_tag, translation_table)
            for obj in tblgene:
                if genbank_flavor == GenbankFlavor.PROKARYOTIC and type(obj) == MRnaFeature:
                    continue
                print(str(obj), file=tbl_file_handle)
