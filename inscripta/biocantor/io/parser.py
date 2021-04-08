"""
Core parser functionality. Contains the dataclass :class:`ParsedAnnotationRecord` which wraps annotations produced
by any of the parser with optional sequence information.
"""
from dataclasses import dataclass
from typing import Optional, Iterable, TextIO, Union
from uuid import UUID

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from inscripta.biocantor.gene.collections import AnnotationCollection
from inscripta.biocantor.io.fasta.exc import FastaExportError
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.location.location_impl import SingleInterval
from inscripta.biocantor.location.strand import Strand
from inscripta.biocantor.parent import Parent, SequenceType
from inscripta.biocantor.sequence.alphabet import Alphabet
from inscripta.biocantor.sequence.sequence import Sequence


@dataclass
class ParsedAnnotationRecord:
    """Dataclass that wraps a :class:`~biocantor.io.models.AnnotationCollectionModel` along with an accompanying
    :class:`SeqRecord` to store sequence information.

    This is an intermediate that allows for sequence information to be applied to collection objects downstream. This
    can be done with :meth:`to_annotation_collection()`.
    """

    annotation: AnnotationCollectionModel
    seqrecord: Optional[SeqRecord] = None
    alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED

    def to_annotation_collection(self) -> AnnotationCollection:
        """Export to a final model. Will apply the sequence information, if it exists (there is a SeqRecord)."""
        if self.seqrecord:
            parent = seq_to_parent(str(self.seqrecord.seq), alphabet=self.alphabet, seq_id=str(self.seqrecord.id))
        else:
            parent = None
        return self.annotation.to_annotation_collection(parent)

    @staticmethod
    def parsed_annotation_records_to_model(
        annotations: Iterable["ParsedAnnotationRecord"],
    ) -> Iterable[AnnotationCollection]:
        """Convenience function for converting an iterable of ParsedAnnotationRecords file to object model.

        Take a iterator of class:`biocantor.io.parser.ParsedAnnotationRecord` and yield an
        iterable of :class:`~biocantor.gene.collections.AnnotationCollection`.

        This incorporates sequence information on to each :class:`~biocantor.gene.transcript.TranscriptInterval`
        and :class:`~biocantor.gene.feature.FeatureInterval` object.

        Args:
            annotations: Iterable that comes from a parser function.

        Yields:
            :class:`~biocantor.gene.collections.AnnotationCollection` with sequence information.
        """
        for annotation in annotations:
            yield annotation.to_annotation_collection()

    def to_fasta(self, fasta_file_handle: TextIO):
        """Convenience function that writes the associated SeqRecord in this record to FASTA.

        Args:
            fasta_file_handle: Open file handle to write to.

        Raises:
            ``FastaExportError`` if the associated ``SeqRecord`` is null.
        """
        if not self.seqrecord:
            raise FastaExportError("Cannot export FASTA without sequence information.")
        SeqIO.write(self.seqrecord, fasta_file_handle, format="fasta")


def seq_to_parent(
    seq: str,
    alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED,
    seq_id: Optional[str] = None,
    seq_type: Optional[str] = SequenceType.CHROMOSOME,
) -> Parent:
    """Convert a string into a Parent object. This is the intermediate that transfers a BioPython sequence object to
    a BioCantor sequence object.

    NOTE: This sequence is assumed to be the entire chromosome.

    Args:
        seq: String of sequence.
        alphabet: Alphabet this sequence is in.
        seq_id: ID to attach to the Parent.
        seq_type: Sequence type to attach to the Parent.

    Returns:
         A :class:`Parent` object.
    """
    return Parent(
        sequence=Sequence(seq, alphabet, type=seq_type, id=seq_id), location=SingleInterval(0, len(seq), Strand.PLUS)
    )


def seq_chunk_to_parent(
    seq: str,
    sequence_name: Union[UUID, str],
    start: int,
    end: int,
    strand: Optional[Strand] = Strand.PLUS,
    alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED,
) -> Parent:
    """Construct a sequence chunk parent from a sequence. This is used when an annotation collection is being
    instantiated with a subset of a genome sequence.

    NOTE: This sequence is assumed to be a subset of a chromosome. There is no way to validate that within this
    function.

    Args:
        seq: Sequence subset to use.
        sequence_name: The name of the sequence.
        start: The genomic start position of this sequence.
        end: The genomic end position of this sequence.
        strand: The strand this chunk is relative to the genome.
        alphabet: The alphabet the sequence is in.

    Returns:
        An instantiated Parent object ready to be passed to a constructor.
    """
    chunk_id = f"{sequence_name}:{start}-{end}"
    return Parent(
        id=chunk_id,
        sequence=Sequence(
            seq,
            alphabet,
            id=chunk_id,
            type=SequenceType.SEQUENCE_CHUNK,
            parent=Parent(
                location=SingleInterval(
                    start,
                    end,
                    strand,
                    parent=Parent(id=sequence_name, sequence_type=SequenceType.CHROMOSOME),
                )
            ),
        ),
    )
