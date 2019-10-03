from primer3plus.utils import anneal
from primer3plus.utils import reverse_complement


def test_reverse_complement():
    s = "ACGTGATGTCGTGAGTAGTA"
    s2 = reverse_complement(s)
    assert s2 == "TGCACTACAGCACTCATCAT"[::-1]


def test_anneal_fwd():
    s = "ACGTGTATGTGATGATGTGCGTGTGTCGTGTAGCTTATTATATGCGGAGTCGTTGATGCTGTGAGT"
    fwd, rev = list(anneal(s, ["NNNNN" + s[20:40]]))
    assert fwd[0]["start"] == 20
    assert fwd[0]["end"] == 40
    assert fwd[0]["anneal"] == s[20:40]
    assert fwd[0]["overhang"] == "NNNNN"
    assert fwd[0]["strand"] == 1


def test_anneal_rev():
    s = "ACGTGTATGTGATGATGTGCGTGTGTCGTGTAGCTTATTATATGCGGAGTCGTTGATGCTGTGAGT"
    fwd, rev = list(anneal(s, ["ANNNN" + reverse_complement(s[20:40])]))
    print(rev[0])
    assert rev[0]["start"] == 40
    assert rev[0]["end"] == 20
    assert rev[0]["anneal"] == reverse_complement(s[20:40])
    assert rev[0]["overhang"] == "ANNNN"
    assert rev[0]["strand"] == -1


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
