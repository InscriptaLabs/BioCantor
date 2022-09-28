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

    def time_parse_genbank(self, gb):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))

    def mem_parse_genbank(self, gb):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / gb)))


class ParseGFF3:
    params = GFF3

    def time_parse_genbank(self, gff3):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / gff3)))

    def mem_parse_genbank(self, gff3):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / gff3)))


class GenBankSequenceExtraction:
    params = GENBANKS

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

    def time_scan_codons(self, gb):
        for rec in self.recs:
            for gene in rec.genes:
                if not gene.is_coding:
                    continue
                _ = list(gene.get_primary_transcript().cds.scan_codon_locations)


class AnnotationCollectionIntervalQuery:
    params = [
        [0, 100000, True],
        [50000, 100000, True],
        [0, 50000, True],
        [20000, 50000, True],
        [0, 100000, False],
        [50000, 100000, False],
        [0, 50000, False],
        [20000, 50000, False],
    ]

    def setup_cache(self):
        return list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "MG1655_subset.gbff"))
        )[0]

    def time_interval_query(self, range, rec):
        start, end, completely_within = range
        _ = rec.query_by_position(start, end, completely_within)
