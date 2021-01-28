from inscripta.biocantor.gene.biotype import Biotype
import pytest


@pytest.mark.parametrize(
    "a,b",
    (
        [Biotype["protein_coding"], Biotype["protein-coding"]],
        [Biotype.misc_RNA == Biotype.miscRNA],
        [Biotype.pseudogene == Biotype.pseudo],
        [Biotype.lncRNA == Biotype.lnc_RNA],
    ),
)
def test_redundant_biotypes(a, b):
    assert a == b
