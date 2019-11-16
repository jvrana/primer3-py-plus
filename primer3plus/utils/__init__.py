import warnings
from itertools import product
from typing import Dict
from typing import Iterator
from typing import List
from typing import Tuple
from typing import Union

warnings.simplefilter("ignore", PendingDeprecationWarning)
rcdict = dict(zip("agtcn ", "tcagn "))
rcdict.update({k.upper(): v.upper() for k, v in rcdict.items()})


def reverse_complement(seq: str):
    return "".join(rcdict[x] for x in seq[::-1])


def _extend_match(
    seq: str, primer: str, length: int, end: int, ignore_case: bool = True
):
    anneal = primer[-length:]
    if ignore_case:
        if anneal.lower() not in seq.lower():
            raise ValueError
    elif anneal not in seq:
        raise ValueError
    try_anneal = anneal
    expected_seq = seq[end - length : end]

    if ignore_case:

        def is_eq(s1, s2):
            return s1.lower() == s2.lower()

    else:

        def is_eq(s1, s2):
            return s1 == s2

    while is_eq(try_anneal, expected_seq) and length <= len(primer):
        anneal = try_anneal
        length += 1
        try_anneal = primer[-length:]
        expected_seq = seq[end - length : end]
    return anneal


def _iter_anneal(
    seq: str,
    primer_list: List[Union[str, Tuple[str, str]]],
    n_bases=9,
    ignore_case: bool = True,
) -> Iterator[Dict[str, Union[str, int]]]:
    li = iter(
        (i, i + n_bases, seq[i : i + n_bases]) for i in range(len(seq) - n_bases + 1)
    )

    if ignore_case:

        def is_eq(s1, s2):
            return s1.lower() == s2.lower()

    else:

        def is_eq(s1, s2):
            return s1 == s2

    for (start, end, anneal), p in product(li, primer_list):
        if isinstance(p, str):
            name = None
        else:
            p, name = p
        length = end - start
        if is_eq(anneal, p[-length:]):
            anneal = _extend_match(seq, p, length, end, ignore_case)
            yield {
                "name": name,
                "anneal": anneal,
                "overhang": p[: -len(anneal)],
                "primer": p,
                "start": end - len(anneal),
                "length": len(anneal),
                "top_strand_slice": (end - len(anneal), end),
            }


def anneal_iter(
    seq: str,
    primer_list: List[Union[str, Tuple[str, str]]],
    n_bases=10,
    ignore_case=True,
) -> Tuple[Iterator[Dict[str, Union[str, int]]], Iterator[Dict[str, Union[str, int]]]]:
    """Anneal a list of primers to the sequence. Returns two iterables with
    elements of the forms:

    .. code-block:: JSON

        {
            "name": "<str> the name of the primer annealed",
            "anneal": "<str> the 5'->3' sequence of the annealing portion of the
                        primer"
            "overhang": "<str> the 5'->3' sequence of the remaining portion of the
                        primer"
            "start": "<int> the inclusive starting position of the primer"
            "length": "<int> the length of the annealing portion of the primer",
            "top_strand_slice": "tuple[int, int] the location on the top strand of the
                                 annealing portion of the primer."
        }

    Note the position conventions
    between fwd and reverse are different in order to match the conventions used
    by primer3. It is often more intuitive to use the 'top_strand_slice' key.

    ::
        0123456789
        |--------------------------------|
          [<-----)      start=2  length=7  strand = 1
         (----->]       start=8  length=7  strand = -1

        However, the top_strand_slice in both cases would be [2,9)

    :param seq: the template sequence
    :param primer_list: a list of either Tuples[name, bases] or just bases
    :param n_bases: number of bases for seed matching (default: 10)
    :return: Two iterators of dictionary results.
    """
    if isinstance(primer_list, str):
        raise TypeError("Expected a list of primer sequences, not a str")
    fwd = list(_iter_anneal(seq, primer_list, n_bases, ignore_case))
    for f in fwd:
        f["strand"] = 1

    rev = list(_iter_anneal(reverse_complement(seq), primer_list, n_bases, ignore_case))
    for r in rev:
        r["strand"] = -1
        s = r["start"]
        e = s + r["length"]
        r["top_strand_slice"] = (len(seq) - e, len(seq) - s)
        r["start"] = len(seq) - s - 1

    return fwd, rev


def anneal(
    seq: str,
    primer_list: List[Union[str, Tuple[str, str]]],
    n_bases=10,
    ignore_case: bool = True,
) -> Tuple[List[Dict[str, Union[str, int]]], List[Dict[str, Union[str, int]]]]:
    """
    Anneal a list of primers to the sequence. Note the position conventions
    between fwd and reverse are equivalent, only the `strand` key changes.

    ::
        0123456789
        |--------------------------------|
          [<-----)      start=2  end=9  strand = 1
          [----->)      start=9  end=2  strand = -1

    :param seq: the template sequence
    :param primer_list: a list of either Tuples[name, bases] or just bases
    :param n_bases: number of bases for seed matching (default: 10)
    :return: Tuple of list of dictionary results.
    """
    fwd, rev = anneal_iter(seq, primer_list, n_bases, ignore_case)
    return list(fwd), list(rev)


_warned_list = []


def depreciated_warning(msg, only_once=True):
    if only_once is False or msg not in _warned_list:
        warnings.warn(msg, DeprecationWarning)
        _warned_list.append(msg)
