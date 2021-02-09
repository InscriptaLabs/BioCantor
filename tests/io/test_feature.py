"""
The inscripta.biocantor.io.feature module contains some functions and enums in the __init__.py that are shared between
GFF3 and GenBank parsing. These functions try to extract FeatureInterval information.
"""
from inscripta.biocantor.io.features import extract_feature_types, extract_feature_name_id, merge_qualifiers
import pytest


@pytest.mark.parametrize(
    "q1,q2,expected",
    [
        ({"key1": ["a", "1"], "key2": ["c"]}, {"key1": ["b"]}, {"key1": ["1", "a", "b"], "key2": ["c"]}),
        ({"key2": ["c"]}, {"key1": ["b"]}, {"key1": ["b"], "key2": ["c"]}),
    ],
)
def test_merge_qualifiers(q1, q2, expected):
    assert merge_qualifiers(q1, q2) == expected


@pytest.mark.parametrize(
    "feature_types,feature_qualifiers,expected",
    [
        (set(), {"gbkey": ["abc"], "regulatory_class": ["cool"]}, {"abc", "cool"}),
        (set(), {"gbkey": ["abc"], "Regulatory_Class": ["cool"]}, {"abc", "cool"}),
        (set(), {"gbkey": ["abc"], "Regulator_Class": ["cool"]}, {"abc", "cool"}),
        (set(), {"gbkey": ["abc"], "Regulator_lass": ["cool"]}, {"abc"}),
        (set(), {"type": ["abc"]}, set()),
        (set(), {"feature_type": ["abc"]}, {"abc"}),
        ({"feature"}, {"feature_type": ["abc"]}, {"abc", "feature"}),
    ],
)
def test_extract_feature_types(feature_types, feature_qualifiers, expected):
    # feature_types is normally a pre-populated set that gets updated
    extract_feature_types(feature_types, feature_qualifiers)
    assert feature_types == expected


@pytest.mark.parametrize(
    "feature_qualifiers,feature_name,feature_id",
    [
        ({"feature_name": ["abc"], "feature_id": ["123"]}, "abc", "123"),
        ({"feature_name": ["abc"]}, "abc", None),
        ({"feature_id": ["123"]}, None, "123"),
        # pick first in list
        ({"feature_id": ["123", "abc"]}, None, "123"),
        # resolve multiple based on hierarchy
        ({"feature_id": ["123"], "ID": ["234"]}, None, "234"),
        ({"feature_id": ["123"], "ID": ["234"], "standard_name": ["cool"], "gene": ["notcool"]}, "cool", "234"),
        # case insensitive
        ({"feature_ID": ["123"], "ID": ["234"], "Standard_name": ["cool"], "Gene": ["notcool"]}, "cool", "234"),
        # must be an exact match
        ({"gene_synonym": ["123"]}, None, None),
        ({"my_gene": ["123"]}, None, None),
        ({"feature_id_id": ["123"]}, None, None),
        ({"my_feature_id": ["123"]}, None, None),
    ],
)
def test_extract_feature_name_id(feature_qualifiers, feature_name, feature_id):
    assert extract_feature_name_id(feature_qualifiers) == (feature_name, feature_id)
