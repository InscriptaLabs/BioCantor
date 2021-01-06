"""
Functions for handling FASTA files. Builds BioCantor sequence objects from various types of FASTA files.
"""
from typing import Optional, Dict, TextIO, Iterable

from Bio import SeqIO

from inscripta.biocantor.gene.collections import AnnotationCollection
from inscripta.biocantor.io.fasta.exc import FastaExportError
from inscripta.biocantor.io.parser import seq_to_parent
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence.alphabet import Alphabet


def fasta_to_parents(
    fasta_handle: TextIO, alphabet: Optional[Alphabet] = Alphabet.NT_EXTENDED_GAPPED
) -> Dict[str, Parent]:
    """Parser that converts a FASTA to a dictionary of Parents, with sequence symbols as keys.

    Args:
        fasta_handle: Open file handle in text mode.
        alphabet: Alphabet this file is in.

    Returns:
        Dictionary mapping the symbol of the sequence to a :class:`~biocantor.parent.parent.Parent` object.
    """
    parents = {}
    for rec in SeqIO.parse(fasta_handle, format="fasta"):
        parents[rec.id] = seq_to_parent(str(rec.seq), alphabet=alphabet, seq_id=rec.id)
    return parents


def collection_to_fasta(collections: Iterable[AnnotationCollection], fasta_file_handle: TextIO):
    """
    Take an instantiated :class:`~biocantor.gene.collections.AnnotationCollection` and produce a FASTA file.

    Args:
        collections: Iterable of instantiated :class:`~biocantor.gene.collections.AnnotationCollection` with
            attached sequence information.
        fasta_file_handle: Open file handle to write the FASTA to.

    Raises:
        FastaExportError if there is no sequence attached to any of the collections.
    """
    for collection in collections:
        if not collection.sequence:
            raise FastaExportError("Cannot export FASTA from collections without a sequence.")
        print(collection.sequence.to_fasta(), file=fasta_file_handle)
