# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- `CompoundInterval.relative_interval_to_parent_location()` in the case of overlapping blocks
- `CompoundInterval.gap_list()` in the case of overlapping blocks
### Added
- `CDSInterval.scan_codon_locations()` method
### Removed
- `CDSInterval.intersect()` method

## [0.1.1] 2020-11-17
### Fixed
- Remove duplicated code that had been created by a merge issue
- Don't do modular arithmetic with CDSFrame.NONE (value -1)

## [0.1.0] 2020-10-29
### Added
- Initial release
