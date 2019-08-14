from primer3plus.design import Design
import random
from primer3plus.utils import reverse_complement, anneal


def test_init():
    design = Design()


def test_set(gfp):
    design = Design()
    design.set.template(gfp)
    design.set.left_sequence(gfp[60:85])
    design.set.as_generic_task()
    results = design.run()
    print(results)
    assert results


def check_primers(gfp, primerlist):

    fwd, rev = anneal(gfp, primerlist, n_bases=16)

    span_region_ok = []
    for primer in fwd:
        s = primer["start"]
        l = len(primer["anneal"])
        span_region_ok.append([s, l, -1, -1])
    for primer in rev:
        s = primer["start"]
        l = len(primer["anneal"])
        span_region_ok.append([-1, -1, s, l])
    return span_region_ok


def test_gfp(gfp, iter_random_primer):
    fwd_primers = list(iter_random_primer(25, gfp[:250], 16))
    rev_primers = list(iter_random_primer(25, reverse_complement(gfp)[:250], 16))

    primers = fwd_primers + rev_primers

    # from itertools import product
    #
    # for f, r in product(fwd_primers,  rev_primers):
    #     design = Design()
    #     design.set.left_sequence(f)
    #     design.set.right_sequence(r)
    #     design.set.task('check_primers')
    #     design.set.template(gfp)
    #     design.run()

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
