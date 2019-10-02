import re
import typing
import warnings
from itertools import product

from Bio.Seq import Seq

warnings.simplefilter("ignore", PendingDeprecationWarning)


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def _extend_match(seq, primer, length, end):
    anneal = primer[-length:]
    if anneal not in seq:
        raise ValueError
    try_anneal = anneal
    expected_seq = seq[end - length : end]
    while try_anneal == expected_seq and length <= len(primer):
        anneal = try_anneal
        length += 1
        try_anneal = primer[-length:]
        expected_seq = seq[end - length : end]
    return anneal


def _iter_anneal(seq: str, primer_list: typing.List[typing.Tuple[str, str]], n_bases=9):
    li = iter(
        (i, i + n_bases, seq[i : i + n_bases]) for i in range(len(seq) - n_bases + 1)
    )

    for (start, end, anneal), p in product(li, primer_list):
        if isinstance(p, str):
            name = None
        else:
            p, name = p
        length = end - start
        if anneal == p[-length:]:
            anneal = _extend_match(seq, p, length, end)
            yield {
                "name": name,
                "anneal": anneal,
                "overhang": p[: -len(anneal)],
                "primer": p,
                "end": end,
                "start": end - len(anneal),
            }


def anneal(seq: str, primer_list: typing.List[typing.Tuple[str, str]], n_bases=9):
    fwd = list(_iter_anneal(seq, primer_list, n_bases))
    for f in fwd:
        f["strand"] = 1

    rev = list(_iter_anneal(reverse_complement(seq), primer_list, n_bases))
    for r in rev:
        r["strand"] = -1
        r["end"] = len(seq) - r["start"]
        r["start"] = len(seq) - r["end"]

    return fwd, rev
