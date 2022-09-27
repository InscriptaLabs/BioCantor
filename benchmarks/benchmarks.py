# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
from pathlib import Path

from inscripta.biocantor.io.genbank.parser import parse_genbank, ParsedAnnotationRecord

DATA_DIR = Path(__file__).parent.parent / "tests/data"


class ParseGenBank:
    params = [
        "INSC1003.gbk",
        "INSC1006_chrI.gbff",
        "MG1655_subset.gbff",
        "R64_subset.gbff",
    ]

    def time_parse_genbank(self, gb):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))

    def time_parse_genbank_translate(self, gb):
        recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))
        for rec in recs:
            for gene in rec.genes:
                _ = gene.get_primary_protein()
