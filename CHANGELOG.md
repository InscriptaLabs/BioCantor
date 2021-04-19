# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- Handle genbank files with broken intervals gracefully.
- Fix interval parsing for negative strand features.


## [0.4.3]
### Fixed
- `AnnotationCollection._subset_parent()` now uses `seq_chunk_to_parent` and pulls out the chromosome ID from the chromosome record.
- `CDSInterval.from_dict()` now passes along the parent provided.

### Added
- `strict_parent_compare` parameter for binary set theory operations.
- `AnnotationCollection.query_by_position()` has a new boolean flag `expand_location_to_children` that defaults to False, but if set to True 
will expand the interval to contain the transcripts. When False, it may be the case that transcripts will have their underlying location objects 
sliced down from their original coordinates. The original coordinates are still retained as integer members. 
If the query position is entirely intronic for an isoform, this isoform will have a `EmptyLocation` chunk relative location,
but will still retain a `chromosome_location`.

### Changed
- Added a parent-level sequence identifier to the output of `biocantor.io.parser.seq_chunk_to_parent()`.
- Added a `strand` argument to `biocantor.io.parser.seq_chunk_to_parent()` that allows for the sequence chunk to be strand-referenced.
- `Location.parent_to_relative_location` and `Location.location_relative_to` now has a `optimize_blocks` flag that defaults to True. 
If this flag is False, then these operations will not collapse adjacent or overlapping blocks. 


## [0.4.2]
### Added
- ``Location.union_preserve_overlaps()`` function added. This function produces the union of intervals, while preserving all overlaps.

## [0.4.1]
### Added
- Improved docstrings on interval objects.
- Location objects now have a `full_span` optional flag on all `intersction`, `overlaps` and `contains` functions. This flag has compound intervals be treated as their full span, i.e. from start to end, regardless of compound structure. This flag defaults to `False` in all cases. When two `CompoundInterval` are compared, they are both always compared in their full spans when this flag is `True`.
- `Interval` and `IntervalCollection` objects now are capable of being lifted to arbitrary coordinate systems, returning a new copy. These operations rely on first lifting to a shared chromosomal coordinate system.

### Changed
- New `SequenceType` enum stores whether interval sequences are `chromosome` or `chunk_relative`.
- All objects that accept `SequenceType` information accepts either the `SequenceType` enum OR raw strings.
- `AnnotationCollection` will look at the provided `parent_or_seq_chunk_parent` to see if the bounds of the object can be inferred from the parent object. This is only performed if no `start`/`end` are explicitly provided. If neither are provided, the bounds of the collection are the bounds of its children.
- Refactored `CDSInterval` to be based on `AbstractFeatureInterval`. Moved `CDSPhase` and `CDSFrame` to accomodate the circular import this introduced.
- All `Interval` objects are allowed to have *chromosome* parents without sequence information.    
- Removed versioneer in favor of hard coded versions.


### Fixed
- Some functions on Interval objects were not operating in chromosome coordinates
- `AnnotationCollection.query_by_position()` was not returning valid results if the parent was a sequence chunk.
- GFF3 parser was not inferring transcripts for a gene feature with no children.
- Fixed a bug with missing gene biotypes in GFF3 parsing.


## [0.4.0]
### Added
- All Interval objects now have the ability to be built from subsets of genome sequence (called `sequence_chunk`).
- Querying `AnnotationCollection` objects by coordinates produce new objects with sliced sequences with chunk-relative coordinates. 
- Interval objects built from sequence subsets can be exported in chunk-relative coordinates to GFF3/GenBank.
- Interval objects have new coordinate translation methods that operate in chunk-relative space. Coordinate methods that operate in genomic coordinate space were retained.
- Non-transcribed feature identifier parsing looks in the `note` special field for identifiers.

### Changed
- All Interval objects now must be built directly from coordinates, and do not accept Location objects.
- All Interval objects now hide their Location member. This is to avoid confusion about what coordinate system the Location may be on.
- All interval collections have `__iter__` functions that call `__iter_children()` functions.
- All Interval objects have their core `._location` object hidden, and offer two accessors -- `.chromosome_location` and `.chunk_relative_location`. Note that `.chromosome_location` will not have sequence information attached to it if a sequence chunk was used. Generally, it is advised to not access `.location` objects directly.


## [0.3.1]
### Fixed
- Feature interval identifier regex should exactly match qualifier keys

### Added
- Unified API for identifiers on all interval objects with new property methods `.id` and .`name`.


## [0.3.0]
### Fixed
- `Biotype` enum improperly mapped `protein_coding` and `protein-coding` to different values. Added `mRNA` as another synonym for this type.
- GFF3, BED and GenBank export from Interval objects now raise an exception when the sequence name field is null.

### Added
- Parse `FeatureInterval` and `FeatureIntervalCollection` from GFF3 or GenBank, and write back as well.

### Changed
- `FeatureInterval` now has multiple types, stored as sets. `FeatureIntervalCollection` stores the union of these types, in addition to optionally having its own type.


## [0.2.0] 2021-01-06
### Fixed
- `CompoundInterval.relative_interval_to_parent_location()` in the case of overlapping blocks. Had previously been double counting overlap region.
- `CompoundInterval.gap_list()` in the case of overlapping blocks. Had been raising an error in that case.

### Added
- `CDSInterval.scan_codon_locations()` method. Returns an iterator over codon locations.
- Implement `__hash__()` for `CompoundInterval` and `CDSInterval`
- Implemented data structures for `TranscriptInterval` and `FeatureInterval` that model transcribed and non-transcribed genomic features.
- Implemented data structures for `GeneInterval` and `FeatureIntervalCollection`, that model groups of intervals as genes or generic feature groups.
- Implemented a wrapper data structure `AnnotationCollection` that contains groups of genes and non-transcribed feature groups.
- Implemented the ability to build BioCantor gene models from GFF3 and GenBank files.
- Implemented the ability to export the above data structures as GFF3, GenBank, BED and NCBI TBL formats.
- Implemented Marshmallow dataclasses that allow for serialization and deserialization of the above data structures.
- Copied the bins implementation from gffutils to avoid needing the full dependency set in a minimal install.
- Added a Biotype enumeration that tracks known biotypes. 
- Added caching of sequence retrieval to `Interval` objects.

### Changed
- Migrated sphinx documentation from `automodapi` to `autoapi`. 
- Performance upgrades to interval arithmetic operations.

### Removed
- `CDSInterval.intersect()` method. Frame math was incorrect for complex CDSs and was deemed too difficult to implement correctly.


## [0.1.1] 2020-11-17
### Fixed
- Remove duplicated code that had been created by a merge issue
- Don't do modular arithmetic with CDSFrame.NONE (value -1)


## [0.1.0] 2020-10-29
### Added
- Initial release
