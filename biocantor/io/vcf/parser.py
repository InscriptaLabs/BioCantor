"""
Generic VCF parser. This is the default implementation.

This implementation uses random UUIDs, so that every `Variant` is unique and maps 1-1 to the `VariantInterval`
produced.
"""
import itertools
import warnings
import vcf
from typing import List, Dict, Optional, TextIO, Union
import pathlib

import vcf.model

from biocantor.io.models import VariantIntervalCollectionModel


def convert_vcf_records_to_model(recs: List[vcf.model._Record]) -> Dict[str, List[VariantIntervalCollectionModel]]:
    """
    Default parser for VCF files. Converts VCF records into `VariantIntervalModel`.
    """
    variants = {}
    for seq_id, seq_variants in itertools.groupby(recs, key=lambda v: v.CHROM):
        these_variants = []
        for seq_variant in seq_variants:
            if len(seq_variant.samples) > 1:
                warnings.warn(f"Variant at {seq_variant.POS} has more than one call. Extra calls will be ignored.")
            start = seq_variant.affected_start
            end = seq_variant.affected_end
            # variants always exist on an interval in BioCantor language
            if start == end:
                end += 1
            sample = seq_variant.samples[0]

            # construct one VariantInterval for each alternative haplotype
            for alt in seq_variant.ALT:
                variant_interval = dict(
                    start=start,
                    end=end,
                    sequence=alt.sequence,
                    variant_type=alt.type,
                )
                if hasattr(sample.data, "PS"):
                    variant_interval["phase_block"] = sample.data.PS
                these_variants.append(variant_interval)

        # group VariantInterval's into VariantIntervalCollection based on their phase block
        # by VCF convention, phase blocks are positive integers so we sort unphased variants using a -1 sentinel value
        sorted_variants = sorted(these_variants, key=lambda x: x.get("phase_block", -1))
        grouped_variants = []
        for phase_block, phased_variants in itertools.groupby(sorted_variants, key=lambda x: x.get("phase_block")):
            phased_variants = list(phased_variants)
            if phase_block is not None:
                grouped_variants.append(
                    VariantIntervalCollectionModel.Schema().load(
                        dict(
                            variant_intervals=phased_variants,
                            variant_collection_id=str(phase_block),
                            sequence_name=seq_id,
                        )
                    )
                )
            else:
                for phased_variant in phased_variants:
                    grouped_variants.append(
                        VariantIntervalCollectionModel.Schema().load(
                            dict(
                                variant_intervals=[phased_variant],
                                sequence_name=seq_id,
                            )
                        )
                    )
        variants[seq_id] = grouped_variants

    return variants


def parse_vcf_file(
    variant_handle_or_path: Optional[Union[TextIO, str, pathlib.Path]] = None
) -> Dict[str, List[VariantIntervalCollectionModel]]:
    """
    Wrapper for :meth:`convert_vcf_records_to_model()` that handles file opening. Primary VCF parsing function.
    """
    if isinstance(variant_handle_or_path, str):
        handle = vcf.Reader(filename=variant_handle_or_path)
    elif isinstance(variant_handle_or_path, pathlib.Path):
        handle = vcf.Reader(filename=str(variant_handle_or_path))
    else:
        handle = vcf.Reader(variant_handle_or_path)
    recs = list(handle)
    return convert_vcf_records_to_model(recs)
