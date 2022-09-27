# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from inscripta.biocantor.io.genbank.parser import parse_genbank, ParsedAnnotationRecord
from pathlib import Path


class ParseGenBank:
    params = ["tests/data/INSC1003.gbk", "tests/data/INSC1006_chrI.gbff", "tests/data/MG1655_subset.gbff", "tests/data/R64_subset.gbff"]

    def parse_genbank(self, gb):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(gb)))
