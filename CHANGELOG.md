# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [0.19.0] 2022-10-21
### Added
- `AA_EXTENDED`, `AA_STRICT_GAPPED`, `AA_EXTENDED_GAPPED`, and `AA_STRICT_UNKNOWN` alphabets.
- ASV tests created for benchmarking.
- `extended_gencode` dictionary in constants, holding codon translations for N-containing codons.

### Changed
- `AnnotationCollection.query_by_position()` will now perform faster queries when `cgranges` is installed.
- `AnnotationCollection.query_by_interval_guids()`, `query_by_transcript_interval_guids()`, and `query_by_feature_interval_guids`() now operate on a cached map that will improve performance.
- `CDSInterval().translate()` can now perform translations on unknown codons, using `X` for the unknown amino acid code.
- `Codons` object now promoted from Enum to full object, and can handle ambiguous IUPAC sequences.
- `Codons` object can now translate some N-containing codons if the N is a wobble base.

### Fixed
- `AnnotationCollection.query_by_transcript_interval_guids` and `AnnotationCollection.query_by_feature_interval_guids` now properly returns only the requested TranscriptInterval or FeatureInterval
as child of the GeneInterval or FeatureCollectionInterval objects. Fix bug introduced in 0.18.0 that was returning all children.
- Any of the identifier queries on `AnnotationCollection` would raise an exception if the collection was chunk-relative and the bounds of the returned genes exceeded the bounds of the chunk.

## [0.18.1] 2022-08-29
### Fixed
- ParentModel.to_parent() is no longer case sensitive for type parameter, and will create sequence chunk parents correctly regardless of string casing.

## [0.18.0] 2022-06-29
### Added
- Implemented `VariantInterval` and `VariantIntervalCollection`
- Implemented VCF parser
- Added functionality to `AnnotationCollection` to automatically associated Variants with Intervals, generating alternative haplotypes to use
- Added functionality to all implementations of `AbstractInterval` to take a Variant and convert to a new haplotype
- Adding missing `from_chunk_relative_location` on `FeatureInterval`
- If the library `cgranges` is installed, intersection operations will make use if it to improve runtime
- Added ability to expand window when scanning codons on `CSInterval` rather than the default behavior (to only include complete codons)

### Changed
- Moved CDS coordinate conversion logic from `TranscriptInterval` to `CDSInterval`. The original functions on `TranscriptInterval` still exist, and call the new child functions.
- Default value for all children of `AnnotationCollectionModel` are now an empty list instead of None
- `UUID` is no longer a valid `sequence_name` of `ParentModel`

### Fixed
- `FeatureInterval.from_location` was not checking for chunk-relativity
- `CompondInterval._single_intervals` caching was not working as intended
- Serializing `AnnotationCollection` to disk via `AnnotationCollectionModel.Schema().dump()` was not serializing sequence information

### Dependency notes
- Python3.10 cannot be supported until `pysam` supports python3.10


## [0.17.0] 2022-05-24
### Fixed
- `CDSInterval._prepare_single_exon_window_for_scan_codon_locations()` would improperly raise exceptions when operating in genome-relative coordinates
-  Fixed jupyter notebooks in documentation

### Added
- More documentation

## [0.16.1] 2022-05-10
### Fixed
- GenBank writer module was not propagating locus tags to child features if the locus tag was derived from the gene symbol because the feature lacked a locus tag


## [0.16.0] 2022-04-13
### Fixed
- Duplicate `/locus_tag` detection in GenBank parser was not working in all cases
- Serialized and simplified the implementation of `CompoundLocation.relative_interval_to_parent_location`

### Changed
- Allow any version of python3 above 3.7

### Added
- GFF3 parsing now supports providing a pre-built gffutils database

## [0.15.0] 2022-03-3
### Changed
- Added flag `allow_duplicate_sequence_identifiers` to GenBank parser
- GenBank parser now handles multiple isoforms of coding and non-coding genes


## [0.14.0] 2022-02-18
### Added
- Alphabet module now supports `NT_STRICT_UNKNOWN` (`ATGCN`)

### Changed
- GenBank parser was overhauled
- GenBank parsing now errs on the side of warnings instead of exceptions when encountering bad data
- Hybrid GenBank parsing now handles duplicate locus tags by shunting them to sorted mode


## [0.13.1] 2022-02-04
### Fixed
- Both Sorted and LocusTag GenBank parser modes did not properly discard mixed-strand annotations if the strand-mixing was present on the `CDS` or `mRNA` feature


## [0.13.0] 2022-01-28
## Added
- `AnnotationCollections` can now export their `Parent` objects in a dictionary representation
- `AnnotationCollectionModel` have a new optional `Parent` member that represents a `Parent` object that can be serialized to disk
- `CDSInterval.scan_chromosome_codon_locations` and `CDSInterval.scan_chunk_relative_codon_locations` now have optional arguments to restrict the iteration to section of the CDS based on *chromosomal* coordinates

### Changed
- Improved documentation
- GenBank parser now handles `CDSInterval`s that exceed the bounds of their Exons by raising a warning and truncating the CDS interval
- GenBank parser now handles duplicate `TranscriptInterval` or `FeatureInterval` objects by removing the duplicates
- GenBank export will no longer propagate `/translation` tags to `mRNA` features when generating eukaryotic style GenBank files

### Fixed
- `AnnotationCollection.query_by_guids` did not accept a single GUID as an argument
- `FeatureInterval`, `TranscriptInterval`, and `CDSInterval` did not take strand into account when generating a GUID
- GenBank parser was inadvertently merging alternative isoforms in all parser modes instead of keeping them separate


## [0.12.0] 2021-10-22
### Fixed

- Fixed `CDSInterval._scan_codon_locations_multi_exon` to properly handle chunk-relative CDS.
- Fixed `CDSInterval.chunk_relative_frames` to properly handle chunk-relative CDS.
- GenBank parsing now handles the intervals of overlapping CDS, although it cannot still infer the correct frame downstream of the overlap.
- Fixed `CompoundInterval.minus` to handle subtraction of overlapping intervals that should lead to a `EmptyLocation`.
- Fixed `CompoundInterval.optimize_blocks` and `CompoundInterval.optimize_and_combine_blocks` to handle the possibility of the resulting `Location` being a `EmptyLocation`.

### Added
- `CDSInterval` now has methods `sequence_pos_to_cds` and `sequence_pos_to_amino_acid` to convert sequence positions to amino acid or CDS positions.
- Support for python3.8.


## [0.11.1] 2021-09-17
### Changed
- Use the optimized version of `scan_codon_locations` in the backwards compatible function.

## [0.11.0] 2021-09-17
### Changed
- Re-added `CDSInterval.scan_codon_locations` as a backwards compatible function. Raises a deprecation warning.


## [0.10.0] 2021-09-17
### Changed
- `CDSInterval` object now has methods to access the number of codons and codon locations in both chunk-relative and chromosome coordinates. Chromosome accessors will always return the full original CDS.
- If duplicate sequence identifiers are found when parsing GenBank/FASTA files, an exception is raised.
- The `scan_codon_locations` methods on `CDSInterval` now operate on two algorithms, one simpler algorithm for canonical transcripts (no programmed frameshifts and no offset frames) and the original more robust algorithm otherwise.
- `_parse_genes` now handles features that have a CDS with no gene parent at all, and ensures that they are marked as protein coding.  

## [0.9.0] 2021-09-03
### Fixed
- `LocusTag` parser is now able to handle CDS-only features, and have more informative error reporting for locus tag collisions. 
- GenBank parser no longer raises an exception for multi-stranded feature; instead it warns and moves on


## [0.8.0]
### Fixed
- Do not trust feature type annotations to define coding vs. non-coding when parsing GenBank files; only rely on the presence/absence of CDS intervals associated with the transcript.
- Setup requirements, build tests, and sphinx config updated to allow building of documentation without installing the package.
- GenBank parser was not handling `exon` features as direct descendants of `gene` correctly.

### Changed
- `Sorted` parser now sorts features by position, then gene/mRNA/CDS/other. This helps deal with genbank files that are oddly ordered.
- Introduced new `Hybrid` GenBank parser mode that does both `LocusTag` and `Sorted` parsing at the same time.
- All `.chromosome_location` accessor of `Interval` objects always return the full length `Location`, even if the `Interval` itself is chunk relative such that the underlying `Location` object cannot represent the full length. As a result of this, the `.chromosome_location` of a chunk-relative location cannot have associated sequence information.

### Added
- `TranscriptInterval` and `FeatureInterval` now have accessor methods to get `Location` objects for their introns/gaps and full span.

## [0.7.0]
### Changed
- GenBank position-sorted parser can now handle CDS records that are not directly following a gene record.
- Refactor `Location`, `Parent` and `Sequence` to have base classes `AbstractLocation`, `AbstractParent` and `AbstractSequence` that are in the base of the `inscripta.biocantor.location` module. This greatly helps with resolving circular imports.
- Optimized checking `sequence` and `location` members to explicitly check for `None`. This avoids a call to `__len__`.
- `CompoundInterval._single_intervals` is now lazily evaluated, because it is expensive to generate many `SingleInterval` objects.
- `CompoundInterval` now stores the positions as two sorted integer lists.
- `CompoundInterval` constructor accepts tuples in addition to lists of integer values to avoid list construction overhead.
- `CompoundInterval.is_overlapping` and `CompoundInterval.is_contiguous` are lazily evaluated.
- `CompoundInterval._combine_blocks` now always removes empty blocks. The new implementation also avoids producing a new interval if the result is identical to the start.
- `unique_value_or_none` was pulled out of `Parent` into its own separate function with an associated cache. This function was optimized to use sets.
- Added `__slots__` to all child classes of `AbstractLocation`, `AbstractSequence` and `AbstractParent`.
- Removed unnecessary call to `strip_location_info()` in `Sequence` constructor.
- Removed all unnecessary instances of constructing lists, replacing them with iterators and tuples.
- GenBank export now defaults to not updating `/translation` tag in order to save execution time. The original behavior can be restored by setting `update_translations=True`.

### Fixed
- GenBank parser was not properly handling 0bp intervals, which can be sometimes seen as a way to represent insertions.
- GenBank parser was not capturing CDS qualifiers when parsing eukaryotic style GenBank files that have mRNA level features


## [0.6.0]
### Changed
- Added `raise_on_reserved_attributes` flag to GFF3 export that controls whether reserved attributes lead to warnings or exceptions.
- Added more top-level imports to simplify imports
- Try more common identifiers when parsing gene symbols from GFF3 files
- Attempt to infer frame from GFF3 files with null Phase columns on CDS records
- Update Tox tests to have a separate formatting case


## [0.5.0]
### Changed
- Added ability to parse non-transcribed features from GenBank records without a parent /gene record in the position-sorted parser.
- Added ability to export `SeqRecord` annotations when writing to GenBank.
- Added methods to `FeatureInterval` that mirror `TranscriptInterval`.
- Added support for translating with non-standard codon tables.

## [0.4.5]
### Changed
- Remove contributor license agreement, which is superseded by the MIT license.

## [0.4.4]
### Fixed
- Handle genbank files with broken intervals gracefully.
- Fix interval parsing for negative strand features.

### Changed
- The tag `Name` can now be used to identify a feature interval in a GFF3/GenBank file.


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
