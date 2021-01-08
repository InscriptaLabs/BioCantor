"""
Perform MD5 digests of arbitrary objects in memory in python. All objects are unpacked and their string representations
are digested.
"""
import hashlib
from uuid import UUID


def digest_object(*args, **kwargs) -> UUID:
    """MD5 digest of any arbitrary set of python objects. Must be utf-8 encodeable."""
    hasher = hashlib.md5()
    for member in args:
        hasher.update(str(member).encode("utf-8"))
    for val in kwargs.values():
        hasher.update(str(val).encode("utf-8"))
    return UUID(hasher.hexdigest())
