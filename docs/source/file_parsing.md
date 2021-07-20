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

Because GenBank files do not have an explicit hierarchical structure to annotations like 
GFF3 files do, the hierarchy must be inferred. BioCantor offers two ways to do this. 
The first is to assume that the file is sorted, and the second 
is to use the key qualifier `locus_tag` to group together objects that belong to the 
same gene or feature.

## Parser Implementations


### LocusTag (default parser)

The `LocusTag` parser implementation groups features in the GenBank file based on the value
of the `locus_tag` qualifier. 

Within each `locus_tag` group, genes are identified by looking for the presence of a `gene`
feature. If a `locus_tag` group exists with no `gene` features, an exception will be raised.

Any features without a `locus_tag` will be interpreted as a `FeatureIntervalCollection` with
one `FeatureInterval`.

While the GenBank parser defaults to the `LocusTag` parser, two other implementations 
also exist that can be chosen by specifying a different `GenBankParserType` to the
`gbk_type` keyword argument in the parser function.

### Sorted

This parser implementation relies entirely on the sort order of the GenBank file. Genes are identified by looking
at the ordering of features within the file, splitting them up into groups every time a `gene` feature is
found. This parser will also handle isolated `ncRNA` or `CDS` features that have no `gene` feature, and a `gene`
feature will be inferred.

### Hybrid

The hybrid parser combines both `LocusTag` and `Sorted` parsing together. The initial parsing pass uses `LocusTag`
parsing to group features together by `locus_tag`. After this step, the remaining features are passed to the
`Sorted` parser.
