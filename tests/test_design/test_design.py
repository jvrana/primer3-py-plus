from primer3plus.design import Design
import random
from primer3plus.utils import reverse_complement, iter_anneal
from itertools import chain


def test_init():
    design = Design()


def iter_random_primer(n, seq, l):

    for i in range(n):
        start = random.randint(0, len(seq) - l)
        end = start + l
        yield seq[start:end]


def test_set(gfp):
    design = Design()
    design.set.template(gfp)
    design.set.left_sequence(gfp[60:85])
    design.set.as_generic_task()
    results = design.run()
    print(results)
    assert results


def check_primers(gfp, primerlist):

    fwd, rev = iter_anneal(gfp, primerlist, n_bases=16)

    span_region_ok = []
    for primer in fwd:
        s = primer["start"]
        l = len(primer["anneal"])
        span_region_ok.append([s, l, -1, -1])
    for primer in rev:
        s = len(gfp) - primer["end"]
        l = len(primer["anneal"])
        span_region_ok.append([-1, -1, s, l])
    return span_region_ok


def test_gfp(gfp):
    fwd_primers = list(iter_random_primer(10, gfp, 16))
    rev_primers = list(iter_random_primer(10, reverse_complement(gfp), 16))
    primers = fwd_primers + rev_primers

    region_ok = check_primers(gfp, primers)
    # print(region_ok[:1])
    design = Design()
    design.set.template(gfp)
    print(region_ok)
    design.set.pair_region_list(region_ok[:10])
    design.set.product_size([50, 2000])
    pairs, explain = design.run()
    print(explain)
    print(pairs)
