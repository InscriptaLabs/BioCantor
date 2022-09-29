from pathlib import Path

from inscripta.biocantor.io.genbank.parser import parse_genbank, ParsedAnnotationRecord
from inscripta.biocantor.io.gff3.parser import parse_standard_gff3

DATA_DIR = Path(__file__).parent.parent / "tests/data"


GENBANKS = [
    "INSC1003.gbk",
    "INSC1006_chrI.gbff",
    "MG1655_subset.gbff",
    "R64_subset.gbff",
]


GFF3 = [
    "INSC1003.gff3",
    "INSC1006_chrI.gff3",
]


class ParseGenBank:
    params = GENBANKS
    repeat = (1, 5, 30.0)

    def time_parse_genbank(self, gb):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))

    def mem_parse_genbank(self, gb):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))


class ParseGFF3:
    params = GFF3
    repeat = (1, 5, 30.0)

    def time_parse_gff3(self, gff3):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / gff3)))

    def mem_parse_gff3(self, gff3):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / gff3)))


class GenBankSequenceExtraction:
    params = GENBANKS
    repeat = (1, 10, 10.0)

    def setup(self, gb):
        self.recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))

    def time_get_primary_protein(self, gb):
        for rec in self.recs:
            for gene in rec.genes:
                _ = gene.get_primary_protein()

    def time_get_primary_transcript_sequence(self, gb):
        for rec in self.recs:
            for gene in rec.genes:
                _ = gene.get_primary_transcript_sequence()


class AnnotationCollectionIntervalQuery:
    params = [True, False]
    param_names = "completely_within"
    repeat = (1, 10, 10.0)

    def setup_cache(self):
        parsed_gb = list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "MG1655_subset.gbff"))
        )[0]
        return parsed_gb

    def time_interval_query_50000_100000(self, parsed_gb, completely_within):
        start, end = 50000, 100000
        _ = parsed_gb.query_by_position(start, end, completely_within)
