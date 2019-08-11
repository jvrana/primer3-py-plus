from primer3plus.utils import reverse_complement, initial_matches


def test_reverse_complement():
    s = "ACGTGATGTCGTGAGTAGTA"
    print(reverse_complement(s))


def test_init_matches():
    s = "ACGTGTATGTGATGATGTGCGTGTGTCGTGTAGCTTATTATATGCGGAGTCGTTGATGCTGTGAGT"
    match = list(initial_matches(s, ["AAAAAGTGCGTGTGTCGTGTAG"]))[0]
    assert match["start"] == 16
    assert match["end"] == 33
    assert match["anneal"] == "GTGCGTGTGTCGTGTAG"
    assert match["overhang"] == "AAAAA"
