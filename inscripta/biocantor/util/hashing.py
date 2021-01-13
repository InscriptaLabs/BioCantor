"""
Perform MD5 digests of arbitrary objects in memory in python. All objects are unpacked and their string representations
are digested.
"""
import hashlib
from typing import Dict, Hashable, Any, List, Set, Iterable
from uuid import UUID


def _order_set(s: Set[Hashable]) -> List[str]:
    """
    Lexicographically orders a set and converts all values to strings to enable order comparison.

    Args:
        s: A set of hashable items

    Returns:
        An lexicographically ordered list of strings
    """
    return sorted(str(x) for x in s)


def _order_dict_of_possible_sets(d: Dict[Hashable, Any]) -> Iterable[str]:
    """
    Recursively search a dictionary for sets and order them.

    Args:
        d: Any dictionary

    Returns:
        Ordered list of tuples
    """
    for key in sorted(d):
        val = d[key]
        yield str(key)
        if isinstance(val, dict):
            yield from _order_dict_of_possible_sets(val)
        elif isinstance(val, set):
            yield str(_order_set(val))
        else:
            yield str(d[key])


def _encode_object_for_digest(*args, **kwargs) -> Iterable[str]:
    """
    Inner function for :meth:`digest_object()` that produces the string representations. This helps with debugging.
    """
    if args:
        for member in args:
            if isinstance(member, dict):
                yield from _order_dict_of_possible_sets(member)
            elif isinstance(member, set):
                yield str(_order_set(member))
            else:
                yield str(member)
    yield from _order_dict_of_possible_sets(kwargs)


def digest_object(*args, **kwargs) -> UUID:
    """MD5 digest of any arbitrary set of python objects. Must be utf-8 encodeable.

    Because the string representation of sets are not guaranteed to be ordered across different
    python interpreters, we recursively hunt for sets in args and kwargs and then sort them. This sorting requires
    conversion to string representation for set members. All items in args and kwargs are also utf-8 encoded
    before hashing.

    Args can be any set of objects with stable string representations, as well as sets.
    Kwargs can be any set of objects with stable string representations, including nested dicts of sets. The argument
    name in the kwargs dictionary are part of the hash produced.
    """
    hasher = hashlib.md5()
    for val in _encode_object_for_digest(*args, **kwargs):
        hasher.update(val.encode("utf-8"))
    return UUID(hasher.hexdigest())
