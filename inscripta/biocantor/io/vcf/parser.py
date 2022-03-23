"""
Generic VCF parser. This is the default implementation.

This implementation uses random UUIDs, so that every `Variant` is unique and maps 1-1 to the `VariantInterval`
produced.
"""
import itertools
import warnings
from typing import List

import vcf.model

from inscripta.biocantor.io.models import VariantIntervalModel


def parse_vcf_records(recs: List[vcf.model._Record]) -> List[VariantIntervalModel]:
    """
    Default parser for VCF files. Converts VCF records into `VariantIntervalModel`.
    """
    variants = []
    for seq_id, seq_variants in itertools.groupby(recs, key=lambda v: v.CHROM):
        for seq_variant in seq_variants:
            if len(seq_variant.samples) > 1:
                warnings.warn(f"Variant at {seq_variant.POS} has more than one call. Extra calls will be ignored.")
            start = seq_variant.affected_start
            end = seq_variant.affected_end
            sample = seq_variant.samples[0]

            for alt in seq_variant.ALT:
                variant_interval = dict(
                    start=start,
                    end=end,
                    sequence=alt.sequence,
                    type=alt.type,
                )
                if hasattr(sample.data, "PS"):
                    variant_interval["phase_block"] = sample.data.PS
                variant_interval = VariantIntervalModel.Schema().load(variant_interval)
                variants.append(variant_interval)
    return variants
