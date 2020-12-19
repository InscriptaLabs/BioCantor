"""
Functions for writing GFF3. Capable of writing GFF3 + FASTA.
"""
from typing import Iterable, Optional, TextIO

from inscripta.biocantor.gene.collections import AnnotationCollection
from inscripta.biocantor.io.gff3.constants import GFF3Headers
from inscripta.biocantor.io.gff3.exc import GFF3ExportException


def collection_to_gff3(
    collections: Iterable[AnnotationCollection],
    gff3_handle: TextIO,
    add_sequences: Optional[bool] = False,
    ordered: Optional[bool] = True,
):
    """
    Take an instantiated :class:`~biocantor.gene.collections.AnnotationCollection` and produce a GFF3 file. Has the
    capability to write GFF3+FASTA format, but the necessary information must be present.

    Args:
        collections: Iterable of AnnotationCollections. If ``add_sequences`` is ``True``, then the collection must
            have instantiated sequences.
        gff3_handle: Open file handle to write GFF3 file to.
        add_sequences: If set to ``True``, the collections must have sequence information. The GFF3 written will
            have the FASTA sequence at the end of it.
        ordered: If set to ``True``, the output GFF3 will be sequence then position sorted.
    """
    if ordered is True:
        collections = sorted(collections, key=lambda c: c.sequence_name)

    print(GFF3Headers.HEADER.value, file=gff3_handle)

    if add_sequences:
        for collection in collections:
            if collection.sequence is None:
                raise GFF3ExportException("Cannot export FASTA in GFF3 if collection has no associated sequence")
            print(
                GFF3Headers.SEQUENCE_HEADER.value.format(
                    symbol=collection.sequence_name, length=len(collection.sequence)
                ),
                file=gff3_handle,
            )

    for collection in collections:
        for item in collection.to_gff(ordered=ordered):
            print(item, file=gff3_handle)

    if add_sequences:
        print(GFF3Headers.FASTA_HEADER.value, file=gff3_handle)
        for collection in collections:
            print(collection.sequence.to_fasta(60), file=gff3_handle)
