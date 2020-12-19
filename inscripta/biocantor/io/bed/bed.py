"""
BED comes in a handful of flavors: BED3, BED6, BED12, and BED12+.

The 12 defined columns are:

1. ``chrom``: Sequence.
2. ``start``: 0-based start.
3. ``end``: 0-based exclusive end.
4. ``name``: A name.
5. ``score``: A score. Should be an integer between 0 and 1000.
6. ``strand``: A string. Any of ``[+, -, .]``.
7. ``thickStart``: Translation start site.
8. ``thickEnd``: Translation end site.
9. ``itemRgb``: String encoded tuple for an RGB value to display.
10. ``blockCount``: Number of spliced blocks.
11. ``blockSizes``: Length of each block.
12. ``BlockStarts``: Start of each block relative to ``start``.

BED3 is the first 3 columns, BED6 is the first 6, and BED12 is the full thing.

BED12+ basically is just people tacking on extra self-defined columns.

These representations **cannot be present in the same file**. For this reason, we have separate subclasses
for each type that only return their own type.
"""

from dataclasses import dataclass, astuple
from typing import List

from inscripta.biocantor.location import Strand


@dataclass(frozen=True)
class RGB:
    r: int = 0
    g: int = 0
    b: int = 0

    def __str__(self) -> str:
        return ",".join(str(color) for color in astuple(self))


@dataclass
class BED3:
    """BED3 includes basic interval information; the simplest type of interval"""

    chrom: str
    start: int
    end: int

    def __str__(self) -> str:
        return ",".join(str(col) for col in astuple(self))


@dataclass
class BED6(BED3):
    """BED6 includes name, score and strand information. Cannot be spliced."""

    name: str
    score: int
    strand: Strand

    def __str__(self) -> str:
        return "\t".join(map(str, [self.chrom, self.start, self.end, self.name, self.score, self.strand.to_symbol()]))


@dataclass
class BED12(BED6):
    """BED12 contains splicing and CDS information. It does not contain frame/phase information."""

    thick_start: int
    thick_end: int
    item_rgb: RGB
    block_count: int
    block_sizes: List[int]
    block_starts: List[int]

    def __str__(self) -> str:
        return "\t".join(
            (
                str(x)
                for x in [
                    self.chrom,
                    self.start,
                    self.end,
                    self.name,
                    self.score,
                    self.strand.to_symbol(),
                    self.thick_start,
                    self.thick_end,
                    self.item_rgb,
                    self.block_count,
                    ",".join(map(str, self.block_sizes)),
                    ",".join(map(str, self.block_starts)),
                ]
            )
        )
