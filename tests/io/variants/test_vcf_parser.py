import json

import pytest

from inscripta.biocantor.io.genbank.parser import ParsedAnnotationRecord, parse_genbank
from inscripta.biocantor.io.models import VariantIntervalCollectionModel
from inscripta.biocantor.io.vcf.parser import parse_vcf_file


@pytest.mark.parametrize(
    "vcf_file,expected_variants_json",
    [
        ["INSC1006_chrI.simulated_variants.vcf", "INSC1006_chrI.simulated_variants.json"],
    ],
)
def test_parser_as_path(test_data_dir, vcf_file, expected_variants_json):
    recs = parse_vcf_file(test_data_dir / vcf_file)

    with open(test_data_dir / expected_variants_json) as fh:
        j = json.load(fh)
        expected_recs = {key: VariantIntervalCollectionModel.Schema().load(j[key], many=True) for key in j.keys()}

    assert expected_recs == recs


@pytest.mark.parametrize(
    "vcf_file,expected_variants_json",
    [
        ["INSC1006_chrI.simulated_variants.vcf", "INSC1006_chrI.simulated_variants.json"],
    ],
)
def test_parser_as_handle(test_data_dir, vcf_file, expected_variants_json):
    with open(test_data_dir / vcf_file) as fh:
        recs = parse_vcf_file(fh)

    with open(test_data_dir / expected_variants_json) as fh:
        j = json.load(fh)
        expected_recs = {key: VariantIntervalCollectionModel.Schema().load(j[key], many=True) for key in j.keys()}

    assert expected_recs == recs


@pytest.mark.parametrize(
    "vcf_file,genbank_file",
    [
        ["INSC1006_chrI.simulated_variants.vcf", "INSC1006_chrI.gbff"],
    ],
)
def test_parser_with_genbank(test_data_dir, vcf_file, genbank_file):
    variant_recs = parse_vcf_file(test_data_dir / vcf_file)
    recs = list(
        ParsedAnnotationRecord.parsed_annotation_records_to_model(
            parse_genbank(test_data_dir / genbank_file, parsed_variants=variant_recs)
        )
    )
    assert recs[0].variant_collections
    recs = list(
        ParsedAnnotationRecord.parsed_annotation_records_to_model(
            parse_genbank(test_data_dir / genbank_file, variant_handle_or_path=test_data_dir / vcf_file)
        )
    )
    assert recs[0].variant_collections
