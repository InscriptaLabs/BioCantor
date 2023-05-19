"""
Test GFF3 attribute export.
"""
import pytest
from biocantor.io.gff3.rows import GFFAttributes, GFF3ExportException, ReservedKeyWarning


class TestAttributes:
    @pytest.mark.parametrize(
        "id,qualifiers,name,parent,expected",
        [
            (
                # qualifiers have a single item and a paired item
                "abc",
                {"key1": {"val1"}, "key2": {"a", 12}},
                None,
                None,
                "ID=abc;key1=val1;key2=12,a",
            ),
            (
                "abc",
                {"key1": {"val1"}, "key2": {"a", 12}},
                "myname",
                "parent1",
                "ID=abc;Parent=parent1;Name=myname;key1=val1;key2=12,a",
            ),
            (
                # escape special characters in values
                # commas are not escaped if they join values or if they exist in input data
                "abc",
                {"key1": {"semi;colon", "a=1", "b,2", "space space", "tab\ttab", "newline\nnewline", "a>a"}},
                "myname",
                "parent1",
                "ID=abc;Parent=parent1;Name=myname;"
                "key1=a%3D1,a%3Ea,b,2,newline%0Anewline,semi%3Bcolon,space%20space,tab%09tab",
            ),
            (
                # commas must be escaped in special key-value pairs like Parent, who can only have one value
                "abc,123",
                {"key1": {"val1"}},
                "myname,myname2",
                "parent1,parent2",
                "ID=abc%2C123;Parent=parent1%2Cparent2;Name=myname%2Cmyname2;key1=val1",
            ),
        ],
    )
    def test_construction(self, id, qualifiers, name, parent, expected):
        attrs = GFFAttributes(id=id, qualifiers=qualifiers, name=name, parent=parent)
        assert str(attrs) == expected

    @pytest.mark.parametrize(
        "id,qualifiers,name,parent,expected_exception",
        [
            ("abc", {"key1": ["123"]}, None, None, GFF3ExportException),
            ("abc", {"Name": {"123"}}, None, None, GFF3ExportException),
        ],
    )
    def test_exceptions(self, id, qualifiers, name, parent, expected_exception):
        with pytest.raises(expected_exception):
            _ = str(GFFAttributes(id=id, qualifiers=qualifiers, name=name, parent=parent))

    @pytest.mark.parametrize(
        "id,qualifiers,name,parent,expected_warning",
        [
            ("abc", {"Target": {"123"}}, None, None, ReservedKeyWarning),
            ("abc", {"Alias": {"123"}}, None, None, ReservedKeyWarning),
        ],
    )
    def test_warnings(self, id, qualifiers, name, parent, expected_warning):
        with pytest.warns(expected_warning):
            _ = str(GFFAttributes(id=id, qualifiers=qualifiers, name=name, parent=parent))
