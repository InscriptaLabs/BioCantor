"""
Core parser functionality. Contains the dataclass :class:`ParsedAnnotationRecord` which wraps annotations produced
by any of the parser with optional sequence information.
"""
from dataclasses import dataclass
from typing import Optional, Iterable

from Bio.SeqRecord import SeqRecord

from inscripta.biocantor.gene.collections import AnnotationCollection
from inscripta.biocantor.io.models import AnnotationCollectionModel
from inscripta.biocantor.parent import Parent
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


def seq_to_parent(
    seq: str,
    alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED,
    seq_id: Optional[str] = None,
    seq_type: Optional[str] = "chromosome",
) -> Parent:
    """Convert a string into a Parent object. This is the intermediate that transfers a BioPython sequence object to
    a BioCantor sequence object.

    Args:
        seq: String of sequence.
        alphabet: Alphabet this sequence is in.
        seq_id: ID to attach to the Parent.
        seq_type: Sequence type to attach to the Parent.

    Returns:
         A :class:`Parent` object.
    """
    return Parent(sequence=Sequence(seq, alphabet, type=seq_type, id=seq_id))
