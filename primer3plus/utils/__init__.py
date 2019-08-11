import warnings
import typing
import re
from itertools import product

warnings.simplefilter("ignore", PendingDeprecationWarning)
from Bio.Seq import Seq


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


def initial_matches(
    seq: str, primer_list: typing.List[typing.Tuple[str, str]], n_bases=9
):
    li = iter(
        (i, i + n_bases, seq[i : i + n_bases]) for i in range(len(seq) - n_bases + 1)
    )

    for (start, end, anneal), p in product(li, primer_list):
        if isinstance(p, str):
            name = None
        else:
            p, name = p
        x = [end - len(p), end + 1]
        length = end - start
        if anneal == p[-length:]:
            anneal = _extend_match(seq, p, length, end)
            yield {
                "name": name,
                "anneal": anneal,
                "overhang": p[: -len(anneal)],
                "end": end,
                "start": end - len(anneal),
            }
