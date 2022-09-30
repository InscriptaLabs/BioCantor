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
    repeat = (1, 5, 20.0)

    def time_parse_genbank_INSC1003(self):
        list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "INSC1003.gbk")))

    def mem_parse_genbank_INSC1003(self):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "INSC1003.gbk")))

    def time_parse_genbank_INSC1006_chrI(self):
        list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "INSC1006_chrI.gbff")))

    def mem_parse_genbank_INSC1006_chrI(self):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "INSC1006_chrI.gbff")))

    def time_parse_genbank_MG1655_subset(self):
        list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "MG1655_subset.gbff")))

    def mem_parse_genbank_MG1655_subset(self):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "MG1655_subset.gbff")))

    def time_parse_genbank_R64_subset(self):
        list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "R64_subset.gbff")))

    def mem_parse_genbank_R64_subset(self):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "R64_subset.gbff")))


class ParseGFF3:
    repeat = (1, 5, 20.0)

    def time_parse_gff3_INSC1003(self):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / "INSC1003.gff3")))

    def mem_parse_gff3_INSC1003(self):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / "INSC1003.gff3")))

    def time_parse_gff3_INSC1006_chrI(self):
        _ = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / "INSC1006_chrI.gff3")))

    def mem_parse_gff3_INSC1006_chrI(self):
        return list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_standard_gff3(DATA_DIR / "INSC1006_chrI.gff3")))


class GenBankSequenceExtraction:
    repeat = (1, 10, 10.0)

    def setup(self, *args, **kwargs):
        self.recs = list(ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "R64_subset.gbff")))

    def time_get_primary_protein(self):
        for rec in self.recs:
            for gene in rec.genes:
                gene.get_primary_protein()

    def time_get_primary_transcript_sequence(self):
        for rec in self.recs:
            for gene in rec.genes:
                gene.get_primary_transcript_sequence()

    def time_get_transcript_sequence(self):
        self.recs[0].genes[20].transcripts[0].get_transcript_sequence()

    def time_scan_codon_locations(self):
        [x for x in self.recs[0].genes[20].transcripts[0].cds.scan_codon_locations()]

    def time_chromosome_codon_locations(self):
        self.recs[0].genes[20].transcripts[0].cds.chromosome_codon_locations

    def time_chunk_relative_codon_locations(self):
        self.recs[0].genes[20].transcripts[0].cds.chunk_relative_codon_locations


class AnnotationCollectionIntervalQuery:
    repeat = (1, 10, 10.0)

    def setup(self, *args, **kwargs):
        gb_record = list(
            ParsedAnnotationRecord.parsed_annotation_records_to_model(parse_genbank(DATA_DIR / "MG1655_subset.gbff"))
        )[0]
        self.gb_record = gb_record
        self.transcript_guids = [next(iter(c.guid_map)) for c in gb_record.genes]
        self.feature_guids = [next(iter(c.guid_map)) for c in gb_record.feature_collections]

    def time_query_by_position(self):
        start, end = 5000, 100000
        self.gb_record.query_by_position(start, end)

    def time_query_by_position_cached(self):
        start, end = 5000, 100000
        self.gb_record.query_by_position(start, end)
        # query twice for caching
        self.gb_record.query_by_position(start, end)

    def time_query_by_interval_guids(self):
        all_intervals = self.transcript_guids + self.feature_guids
        self.gb_record.query_by_interval_guids(all_intervals)

    def time_query_by_interval_guids_cached(self):
        all_intervals = self.transcript_guids + self.feature_guids
        self.gb_record.query_by_interval_guids(all_intervals)
        # query twice for caching
        self.gb_record.query_by_interval_guids(all_intervals)

    def time_query_by_transcript_interval_guids(self):
        self.gb_record.query_by_transcript_interval_guids(self.transcript_guids)

    def time_query_by_transcript_interval_guids_cached(self):
        self.gb_record.query_by_transcript_interval_guids(self.transcript_guids)
        # query twice for caching
        self.gb_record.query_by_transcript_interval_guids(self.transcript_guids)

    def time_query_by_feature_interval_guids(self):
        self.gb_record.query_by_feature_interval_guids(self.feature_guids)

    def time_query_by_feature_interval_guids_cached(self):
        self.gb_record.query_by_feature_interval_guids(self.feature_guids)
        # query twice for caching
        self.gb_record.query_by_feature_interval_guids(self.feature_guids)
