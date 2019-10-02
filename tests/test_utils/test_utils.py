from primer3plus.utils import anneal
from primer3plus.utils import reverse_complement


def test_reverse_complement():
    s = "ACGTGATGTCGTGAGTAGTA"
    print(reverse_complement(s))


def test_anneal():
    s = "ACGTGTATGTGATGATGTGCGTGTGTCGTGTAGCTTATTATATGCGGAGTCGTTGATGCTGTGAGT"
    fwd, rev = anneal(s, ["AAAAAGTGCGTGTGTCGTGTAG"])
    match = list(fwd)[0]
    assert match["start"] == 16
    assert match["end"] == 33
    assert match["anneal"] == "GTGCGTGTGTCGTGTAG"
    assert match["overhang"] == "AAAAA"
    print(match)


def test_anneal_end():
    s = "ACGTGTATGTGATGATGTGCGTGTGTCGTGTAGCTTATTATATGCGGAGTCGTTGATGCTGTGAGT"
    fwd, rev = anneal(s, [s[-16:]])
    match = list(fwd)[0]
    assert match["start"] == len(s) - 16
    assert match["end"] == len(s)
