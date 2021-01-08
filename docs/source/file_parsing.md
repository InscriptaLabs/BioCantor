# Parsing annotation files with BioCantor


BioCantor supports two common annotation file formats, GenBank and GFF3. GenBank parsing is backed by the parser
present in BioPython, while GFF3 parsing relies on the library `gffutils`.

Because annotation file format specifications are often interpreted differently, BioCantor is designed to allow for
custom parsing functions to be written for all input formats. There exists default implementations in the library,
but these can be replaced with your own methods in order to best map the data in your annotation file on to the values
in the BioCantor annotation model.

## GFF3

BioCantor supports import and export of GFF3 with or without an embedded FASTA file. Embedded FASTA files are a method
of combining sequence and annotation information. The annotation section is standard GFF3 up top, delimited with the
header line `##gff-version 3`. The FASTA section comes at the end, and is delimited by the break line `##FASTA`.

This version of GFF3 could be constructed from two separate files with the command:

```
(cat ${GFF3}; echo -e "##FASTA\n"; cat ${FASTA}) > gff3_with_fasta.gff3
```

To enable GFF3 parsing, BioCantor leverages the library [gffutils](http://daler.github.io/gffutils/). 
This library builds a sqlite database of the input file, and has a lot of flexibility that allows for parsing of
the many ways that the GFF3 spec can be interpreted. BioCantor takes the resulting database and interprets it into
the data model. Users can pass commands down to `gffutils` in order to tweak how it interprets the files.

## GenBank

BioCantor relies on `BioPython` to perform the core parsing of GenBank files. The resulting data structures are then
ran through the default BioCantor parsing function to build the BioCantor annotation model.

Because GenBank files do not have an explicit hierarchical structure to annotations like GFF3 files do, the hierarchy
must be inferred. BioCantor offers two ways to do this. The first is to assume that the file is sorted, and the second
is to use the key qualifier `locus_tag` to group together objects that belong to the same gene or feature.

Using the `locus_tag` method works for NCBI GenBank files produced for bacterial genomes. However, if the annotation
in question has more than one isoform for a gene, this will not work, as it is ambiguous which transcript level feature
maps on to which CDS feature. In this case, the `sorted` method is the right choice, after ensuring that your GenBank
file is correctly sorted.
 
The GenBank parser understands two types of sorted hierarchical relationships, one for prokaryotic genomes and one
for eukaryotic genomes. These do not need to be specified on parsing, and will be inferred based on the data found.
The key difference between these two are the presence or absence of certain types of features.

In eukaryotic GenBank files, the parser assumes that the sorted order looks like:

```
Coding genes: gene -> mRNA -> CDS
Non-coding genes: gene -> tRNA/rRNA/etc
```

In prokaryotic GenBank files, the parser assumes that the sorted order looks like:

```
Coding genes: gene -> CDS
Non-coding genes: gene -> tRNA/rRNA/etc
```
