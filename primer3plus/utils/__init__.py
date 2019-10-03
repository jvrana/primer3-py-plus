import warnings
from itertools import product
from typing import Dict
from typing import Iterator
from typing import List
from typing import Tuple
from typing import Union

from Bio.Seq import Seq

warnings.simplefilter("ignore", PendingDeprecationWarning)


def reverse_complement(seq: str):
    return str(Seq(seq).reverse_complement())


def _extend_match(seq: str, primer: str, length: int, end: int):
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


def _iter_anneal(
    seq: str, primer_list: List[Union[str, Tuple[str, str]]], n_bases=9
) -> Iterator[Dict[str, Union[str, int]]]:
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
                "length": len(anneal),
            }


def anneal_iter(
    seq: str, primer_list: List[Union[str, Tuple[str, str]]], n_bases=10
) -> Tuple[Iterator[Dict[str, Union[str, int]]], Iterator[Dict[str, Union[str, int]]]]:
    """
    Anneal a list of primers to the sequence. Note the position conventions
    between fwd and reverse are different in order to match the conventions used
    by primer3. It is often more intuitive to use the 'slice_index' key.

    ::
        0123456789
        |--------------------------------|
          [<-----)      start=2  end=9  length=7  strand = 1
         (----->]       start=8  end=1  length=7  strand = -1

        However, the slice_index in both cases would be [2,9)

    :param seq: the template sequence
    :param primer_list: a list of either Tuples[name, bases] or just bases
    :param n_bases: number of bases for seed matching (default: 10)
    :return: Two iterators of dictionary results.
    """
    fwd = list(_iter_anneal(seq, primer_list, n_bases))
    for f in fwd:
        f["strand"] = 1

    rev = list(_iter_anneal(reverse_complement(seq), primer_list, n_bases))
    for r in rev:
        r["strand"] = -1
        s, e = r["start"], r["end"]
        r["end"] = len(seq) - e
        r["start"] = len(seq) - s

    return fwd, rev


def anneal(
    seq: str, primer_list: List[Union[str, Tuple[str, str]]], n_bases=10
) -> Tuple[List[Dict[str, Union[str, int]]], List[Dict[str, Union[str, int]]]]:
    """
    Anneal a list of primers to the sequence. Note the position conventions
    between fwd and reverse are equivalent, only the `strand` key changes.

    ::
        0123456789
        |--------------------------------|
          [<-----)      start=2  end=9  strand = 1
          [----->)      start=2  end=9  strand = -1

    :param seq: the template sequence
    :param primer_list: a list of either Tuples[name, bases] or just bases
    :param n_bases: number of bases for seed matching (default: 10)
    :return: Tuple of list of dictionary results.
    """
    fwd, rev = anneal_iter(seq, primer_list, n_bases)
    return list(fwd), list(rev)
