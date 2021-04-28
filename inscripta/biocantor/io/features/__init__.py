"""
This module contains code shared between the Genbank and GFF3 parser for extracting feature information.

These regexes, enums and functions are used to identify primary name and ID values, as well as pull
out possible feature types.
"""
import itertools
import re
import string
from collections import defaultdict
from enum import IntEnum
from typing import Set, Dict, Tuple, Hashable, List


class FeatureIntervalNameQualifiers(IntEnum):
    """Keys that will be looked for to find human readable identifier(s) for a feature.
    Will be checked in a case-insensitive fashion.

    The ordering of this enumeration matters, because the key-value pair found with the smallest value is the most
    important.

    In order to be future proofed, the values of this enum are separated by 10, so that up to 9 new items can be
    inserted at every level without changing the values.

    NOTE: If names are added to or changed in this enum, you must also change FEATURE_INTERVAL_NAME_QUALIFIERS.
    """

    FEATURE_NAME = 0
    STANDARD_NAME = 10
    NAME = 15
    GENE = 20
    GENE_NAME = 30
    LABEL = 40
    OPERON = 50


class FeatureIntervalIDQualifiers(IntEnum):
    """Keys that will be looked for to find identifiers for a feature.

    The ordering of this enumeration matters, because the key-value pair found with the smallest value is the most
    important.

    In order to be future proofed, the values of this enum are separated by 10, so that up to 9 new items can be
    inserted at every level without changing the values.

    NOTE: If names are added to or changd in this enum, you must also change FEATURE_INTERVAL_ID_QUALIFIERS.
    """

    FEATURE_ID = 0
    ID = 255


FEATURE_INTERVAL_NAME_QUALIFIERS = {"feature_name", "name", "standard_name", "gene", "gene_name", "label", "operon"}
FEATURE_INTERVAL_NAME_QUALIFIERS_REGEX = re.compile(
    r"({})".format("|".join(f"^{k}$" for k in FEATURE_INTERVAL_NAME_QUALIFIERS)), re.IGNORECASE
)

FEATURE_INTERVAL_ID_QUALIFIERS = {"feature_id", "id"}
FEATURE_INTERVAL_ID_QUALIFIERS_REGEX = re.compile(
    r"({})".format("|".join(f"^{k}$" for k in FEATURE_INTERVAL_ID_QUALIFIERS)), re.IGNORECASE
)


# FEATURE_TYPE_IDENTIFIERS are case-insensitive substrings to match for identifying feature types to include.
FEATURE_TYPE_IDENTIFIERS = {"_class", "gbkey", "_type"}
FEATURE_TYPE_IDENTIFIERS_REGEX = re.compile(
    r"({})".format("|".join(k for k in FEATURE_TYPE_IDENTIFIERS)), re.IGNORECASE
)


def extract_feature_types(feature_types: Set[str], feature_qualifiers: Dict[str, List[str]]):
    """
    Extracts feature types from a qualifiers dictionary.

    Args:
        feature_types: Set of feature types. Starts with the primary feature type.
        feature_qualifiers: Qualifiers associated with a record from a genbank parsing event.
    """
    for key, vals in feature_qualifiers.items():
        if re.search(FEATURE_TYPE_IDENTIFIERS_REGEX, key):
            feature_types.update(vals)


def extract_feature_name_id(feature_qualifiers: Dict[str, List[str]]) -> Tuple[str, str]:
    """
    Extract primary feature name and ID from a qualifiers dictionary based on the hierarchy in
    FeatureIntervalNameKeys.

    If there is more than one value for a key, the first is chosen.

    Args:
        feature_qualifiers: Qualifiers associated with a record from a genbank parsing event.

    Returns:
        The primary identifier name and ID.
    """
    feature_name = None
    feature_key = None
    feature_id = None
    feature_id_key = None
    for qualifier, vals in feature_qualifiers.items():
        # exact case insensitive match
        if re.match(FEATURE_INTERVAL_NAME_QUALIFIERS_REGEX, qualifier):
            this_feature_key = FeatureIntervalNameQualifiers[qualifier.upper()]
            if not feature_key or this_feature_key < feature_key:
                feature_name = vals[0]
                feature_key = this_feature_key
        elif re.match(FEATURE_INTERVAL_ID_QUALIFIERS_REGEX, qualifier):
            this_feature_id_key = FeatureIntervalIDQualifiers[qualifier.upper()]
            if not feature_id_key or this_feature_id_key < feature_id_key:
                feature_id = vals[0]
                feature_id_key = this_feature_id_key

    # if no standard identifiers are present, try using /note
    if not feature_name and not feature_id and "note" in feature_qualifiers:
        try:
            feature_name = feature_id = feature_qualifiers["note"][0].split()[0].strip(string.punctuation)
        except IndexError:
            pass

    return feature_name, feature_id


def merge_qualifiers(
    qualifiers: Dict[Hashable, List[str]], other_qualifiers: Dict[Hashable, List[str]]
) -> Dict[Hashable, List[str]]:
    """Merges two dicts of lists using sets.

    Could be made more efficient, probably.
    """
    merged = defaultdict(set)
    for key, vals in itertools.chain(qualifiers.items(), other_qualifiers.items()):
        merged[key].update(vals)
    return {key: sorted(vals) for key, vals in merged.items()}
