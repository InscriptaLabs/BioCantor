# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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
