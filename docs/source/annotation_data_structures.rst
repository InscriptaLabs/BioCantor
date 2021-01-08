BioCantor Annotation Data Structures
====================================

The above object paradigm is used to construct a hierarchical set of objects that encapsulate common genomic
concepts, and associated interval arithmetic operations.

The core data structure is :class:`biocantor.gene.feature.AbstractFeatureInterval`, which has two child classes:

1. :class:`biocantor.gene.transcript.TranscriptInterval`.
2. :class:`biocantor.gene.feature.FeatureInterval`.

`TranscriptInterval`s are used to model *transcribed features*. They can be coding or non-coding.
In contrast, `FeatureInterval`s are intended to be generic intervals, which may be transcribed or non-transcribed.
Both of these intervals implement a variety of methods that allow for operations such as:

1. Coordinate translations between genomic, CDS and transcript coordinate systems. These operations can be performed
    as either points or as intervals.
2. Sequence extraction in all possible coordinate systems.
3. Export to GFF3, BED and JSON.


Both transcripts and features have **container classes**, called :class:`biocantor.gene.collections.GeneInterval`
and :class:`biocantor.gene.collections.FeatureIntervalCollection` respectively. These classes contain one or more
of their constituent types, and provide methods for iterating over the members.

The final core data model is :class:`biocantor.gene.collections.AnnotationCollection`. This object contains one or more
:class:`~biocantor.gene.collections.GeneInterval` and :class:`~biocantor.gene.collections.FeatureIntervalCollection`.
This is intended to model an arbitrary genomic interval *on a single sequence/chromosome*. This object contains
methods that allow a user to subquery the interval by position or by identifiers of constituent objects. This object
is produced when parsing annotation file formats.
