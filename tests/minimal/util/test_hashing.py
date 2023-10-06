"""
Prove consistent hashing across instances even with unordered datatypes.
"""
import pytest
from biocantor.util.hashing import digest_object, _encode_object_for_digest
from uuid import UUID


@pytest.mark.parametrize(
    "val,str_rep,uuid",
    [
        ({"a", "b", "c"}, ["['a', 'b', 'c']"], UUID("eea45728-5a61-f212-e4bb-aaf890263ab4")),
        ({"a", "c", "b"}, ["['a', 'b', 'c']"], UUID("eea45728-5a61-f212-e4bb-aaf890263ab4")),
        ({"a", "b", 1}, ["['1', 'a', 'b']"], UUID("b4fb5de6-eb7f-5fe3-9a86-c9cb8d7e89cf")),
        ({"a", 1, "b"}, ["['1', 'a', 'b']"], UUID("b4fb5de6-eb7f-5fe3-9a86-c9cb8d7e89cf")),
    ],
)
def test_sets(val, str_rep, uuid):
    assert str_rep == list(_encode_object_for_digest(val))
    assert digest_object(val) == uuid


@pytest.mark.parametrize(
    "val,str_rep,uuid",
    [
        ({"key1": {"a", "b", "c"}}, ["key1", "['a', 'b', 'c']"], UUID("0586f397-e249-5a02-d9e2-f2c4a27d8e44")),
        ({"key1": {"c", "b", "a"}}, ["key1", "['a', 'b', 'c']"], UUID("0586f397-e249-5a02-d9e2-f2c4a27d8e44")),
        ({"key1": {"a", "b", 1}}, ["key1", "['1', 'a', 'b']"], UUID("099ab1fd-be7f-ba53-1a9d-d35a76995d0e")),
    ],
)
def test_dicts_of_sets(val, str_rep, uuid):
    assert list(_encode_object_for_digest(val)) == str_rep
    assert digest_object(val) == uuid


@pytest.mark.parametrize(
    "val,str_rep,uuid",
    [
        (
            {"key1": {"inner": {"a", "b", "c"}}},
            ["key1", "inner", "['a', 'b', 'c']"],
            UUID("f8338949-0b25-5211-3f2d-8aefdabfe8a7"),
        ),
        (
            {"key1": {"inner": {"a", "c", "b"}}},
            ["key1", "inner", "['a', 'b', 'c']"],
            UUID("f8338949-0b25-5211-3f2d-8aefdabfe8a7"),
        ),
        (
            {"key1": {"inner": {"a", "c", "b"}, "key2": "abc"}},
            ["key1", "inner", "['a', 'b', 'c']", "key2", "abc"],
            UUID("e1ecf098-5eff-ced3-7b42-156913463792"),
        ),
    ],
)
def test_nested(val, str_rep, uuid):
    assert list(_encode_object_for_digest(val)) == str_rep
    assert digest_object(val) == uuid


@pytest.mark.parametrize(
    "kwargs,str_rep,uuid",
    [
        (
            {"key1": {"inner": {"a", "b", "c"}}},
            ["key1", "inner", "['a', 'b', 'c']"],
            UUID("f8338949-0b25-5211-3f2d-8aefdabfe8a7"),
        ),
        (
            {"key1": {"inner": {"a", "c", "b"}, "key2": "abc"}},
            ["key1", "inner", "['a', 'b', 'c']", "key2", "abc"],
            UUID("e1ecf098-5eff-ced3-7b42-156913463792"),
        ),
    ],
)
def test_nested_kwargs(kwargs, str_rep, uuid):
    assert list(_encode_object_for_digest(**kwargs)) == str_rep
    assert digest_object(**kwargs) == uuid


@pytest.mark.parametrize(
    "args,kwargs,str_rep,uuid",
    [
        (
            [],
            {"key1": {"inner": {"a", "b", "c"}}},
            ["key1", "inner", "['a', 'b', 'c']"],
            UUID("f8338949-0b25-5211-3f2d-8aefdabfe8a7"),
        ),
        (
            [],
            {"key1": {"inner": {"a", "c", "b"}, "key2": "abc"}},
            ["key1", "inner", "['a', 'b', 'c']", "key2", "abc"],
            UUID("e1ecf098-5eff-ced3-7b42-156913463792"),
        ),
        # order of args does matter
        (
            [12, "ab"],
            {"key1": {"inner": {"a", "b", "c"}}},
            ["12", "ab", "key1", "inner", "['a', 'b', 'c']"],
            UUID("901cb8ed-0030-7e7c-b167-e2f6f990ca47"),
        ),
        (
            ["ab", 12],
            {"key1": {"inner": {"a", "c", "b"}, "key2": "abc"}},
            ["ab", "12", "key1", "inner", "['a', 'b', 'c']", "key2", "abc"],
            UUID("21e41837-5aab-64f0-3fc9-b8b3bce1aab1"),
        ),
    ],
)
def test_nested_args_kwargs(args, kwargs, str_rep, uuid):
    assert list(_encode_object_for_digest(*args, **kwargs)) == str_rep
    assert digest_object(*args, **kwargs) == uuid


@pytest.mark.parametrize(
    "args,kwargs,str_rep,uuid",
    [
        (
            [],
            {"key1": {"inner": {"a", "b", "c"}}, "key2": {"inner": {"a", "b", "c"}}},
            ["key1", "inner", "['a', 'b', 'c']", "key2", "inner", "['a', 'b', 'c']"],
            UUID("fd76d879-6f7c-2aaa-b2ad-e0f5febbb4b0"),
        ),
        (
            [],
            {"key1": {"inner": {"a", "b", "c"}}, "key3": {"inner": {"a", "b", "c"}}},
            ["key1", "inner", "['a', 'b', 'c']", "key3", "inner", "['a', 'b', 'c']"],
            UUID("d694fdb5-f4ba-998a-8f76-0b04b0c1c1f6"),
        ),
        (
            [{"key1": {"inner": {"a", "b", "c"}}, "key2": {"inner": {"a", "b", "c"}}}],
            {},
            ["key1", "inner", "['a', 'b', 'c']", "key2", "inner", "['a', 'b', 'c']"],
            UUID("fd76d879-6f7c-2aaa-b2ad-e0f5febbb4b0"),
        ),
        (
            [{"key1": {"inner": {"a", "b", "c"}}, "key3": {"inner": {"a", "b", "c"}}}],
            {},
            ["key1", "inner", "['a', 'b', 'c']", "key3", "inner", "['a', 'b', 'c']"],
            UUID("d694fdb5-f4ba-998a-8f76-0b04b0c1c1f6"),
        ),
    ],
)
def test_same_values_different_keys(args, kwargs, str_rep, uuid):
    assert list(_encode_object_for_digest(*args, **kwargs)) == str_rep
    assert digest_object(*args, **kwargs) == uuid
