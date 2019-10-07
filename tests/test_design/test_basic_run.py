import json
import random

import pytest

from primer3plus.design import Design
from primer3plus.exceptions import Primer3PlusRunTimeError
from primer3plus.utils import anneal
from primer3plus.utils import reverse_complement as rc


def pprint(pairs):
    print(json.dumps(pairs, indent=1))


def test_init():
    design = Design()


def test_set(gfp):
    design = Design()
    design.settings.template(gfp)
    design.settings.left_sequence(gfp[60:85])
    design.settings.as_generic_task()
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
    rev_primers = list(iter_random_primer(25, rc(gfp)[:250], 16))

    primers = fwd_primers + rev_primers

    # from itertools import product
    #
    # for f, r in product(fwd_primers,  rev_primers):
    #     design = Design()
    #     design.settings.left_sequence(f)
    #     design.settings.right_sequence(r)
    #     design.settings.task('check_primers')
    #     design.settings.template(gfp)
    #     design.run()

    region_ok = check_primers(gfp, primers)
    # print(region_ok[:1])
    design = Design()
    design.settings.template(gfp)
    print(region_ok)
    design.settings.pair_region_list(region_ok[:10])
    design.settings.product_size([50, 2000])
    pairs, explain = design.run()
    print(explain)
    print(pairs)


def test_product_size_list(gfp):

    design = Design()
    design.settings.template(gfp)
    design.settings.product_size([(50, 100), (200, 300)])
    design.run()


class TestExclude:
    def test_as_tuple(self, gfp):
        design = Design()
        design.settings.template(gfp).excluded((100, 300))
        pairs, explain = design.run()
        assert pairs

    def test_as_list(self, gfp):
        design = Design()
        design.settings.template(gfp).excluded([100, 300])
        print(design.SEQUENCE_EXCLUDED_REGION)
        pairs, explain = design.run()
        assert pairs

    def test_as_list_of_tuples(self, gfp):
        design = Design()
        design.settings.template(gfp).excluded([(100, 300)])
        pairs, explain = design.run()
        assert pairs

    def test_as_str(self, gfp):
        design = Design()
        design.settings.template(gfp).excluded(gfp[100:300])
        pairs, explain = design.run()
        assert pairs


class TestIncluded:
    def test_as_tuple(self, gfp):
        design = Design()
        design.settings.template(gfp).included((100, 100))
        print(design.SEQUENCE_INCLUDED_REGION)
        pairs, explain = design.run()
        assert pairs

    def test_as_list(self, gfp):
        design = Design()
        design.settings.template(gfp).included([100, 100])
        print(design.SEQUENCE_INCLUDED_REGION)
        pairs, explain = design.run()
        assert pairs

    def test_as_str(self, gfp):
        design = Design()
        design.settings.template(gfp).included(gfp[100:300])
        print(design.SEQUENCE_INCLUDED_REGION)
        pairs, explain = design.run()
        assert pairs

    def test_raises_value_error(self, gfp):
        design = Design()
        design.settings.template(gfp).included([None, 100])
        print(design.SEQUENCE_INCLUDED_REGION)
        with pytest.raises(TypeError):
            design.run()

    def test_raises_value_error(self, gfp):
        design = Design()
        with pytest.raises(TypeError):
            design.settings.template(gfp).included([100, 100, 100])


class TestOverhangs:
    """Tests for resolving primer overhangs."""

    # TODO: add overhangs to results
    def test_left_overhang(self, gfp):
        design = Design()
        design.settings.template(gfp)
        design.settings.left_sequence("AGGCGGCTGA" + gfp[0:20])
        design.settings.use_overhangs()
        pairs, explain = design.run_and_optimize(5)
        print(explain)
        assert pairs

    def test_left_overhang_default_behavior(self, gfp):
        """Without setting overhangs, should raise run time error"""
        design = Design()
        design.settings.template(gfp)
        design.settings.left_sequence("AGGCGGCTGA" + gfp[20:40])
        with pytest.raises(Primer3PlusRunTimeError):
            design.run()

    def test_right_overhang(self, gfp):
        design = Design()
        design.settings.template(gfp)

        rseq = "AGGCGGCTGA" + rc(gfp[200:225])
        design.settings.right_sequence(rseq)
        design.settings.use_overhangs()
        pairs, explain = design.run_and_optimize(5)
        assert pairs

    def test_right_overhang_default_behavior(self, gfp):
        """Without setting overhangs, should raise run time error"""
        design = Design()
        design.settings.template(gfp)

        rseq = "AAAAAA" + rc(gfp[200:220])
        design.settings.right_sequence(rseq)
        with pytest.raises(Primer3PlusRunTimeError):
            design.run()

    def test_both_overhangs(self, gfp):
        design = Design()
        design.settings.template(gfp)
        lseq = "TTTTTT" + gfp[10:30]
        rseq = rc(gfp[200:220] + "AAAAAAA")

        design.settings.left_sequence(lseq)
        design.settings.right_sequence(rseq)
        design.settings.use_overhangs()
        pairs, explain = design.run_and_optimize(5)
        assert pairs


class TestLongPrimers:
    def test_long_left_primer(self, gfp):

        design = Design()
        design.settings.template(gfp)
        design.settings.left_sequence(gfp[0:50])
        design.PRIMER_MAX_TM.value = 75.0
        design.PRIMER_MAX_SIZE = 35
        design.PRIMER_PICK_ANYWAY = 1
        design.settings.long_ok()
        pairs, explain = design.run_and_optimize(10)
        print(explain)
        assert pairs

    def test_long_left_primer_default_behavior(self, gfp):

        design = Design()
        design.settings.template(gfp)
        design.settings.left_sequence(gfp[0:50])
        design.PRIMER_MAX_TM.value = 75.0
        design.PRIMER_MAX_SIZE = 35
        design.PRIMER_PICK_ANYWAY = 1
        with pytest.raises(Primer3PlusRunTimeError):
            design.run_and_optimize(10)

    def test_long_right_primer(self, gfp):

        design = Design()
        design.settings.template(gfp)
        design.settings.right_sequence(rc(gfp[-50:]))
        design.PRIMER_MAX_TM.value = 75.0
        design.PRIMER_MAX_SIZE = 35
        design.PRIMER_PICK_ANYWAY = 1
        design.settings.long_ok()
        pairs, explain = design.run_and_optimize(10)
        print(explain)
        assert pairs

    def test_long_right_primer_default_behavior(self, gfp):

        design = Design()
        design.settings.template(gfp)
        design.settings.right_sequence(rc(gfp[-50:]))
        design.PRIMER_MAX_TM.value = 75.0
        design.PRIMER_MAX_SIZE = 35
        design.PRIMER_PICK_ANYWAY = 1
        with pytest.raises(Primer3PlusRunTimeError):
            design.run_and_optimize(10)


class TestOverOrigin:
    def test_included_over_origin(self, gfp):

        design = Design()
        design.settings.template(gfp)
        design.settings.included((len(gfp) - 100, 110))
        design.run()


def test_long_right_primer_with_overhangs(gfp):

    design = Design()
    design.settings.template(gfp)
    design.settings.right_sequence(rc(gfp[-50:]))
    design.settings.product_size((len(gfp), len(gfp)))

    loverhang = "TTTTTTT"
    roverhang = "GGGAGAG"

    design.settings.left_overhang(loverhang)
    design.settings.right_overhang(roverhang)
    design.settings.long_ok()
    design.settings.use_overhangs()
    design.settings.pick_anyway()
    pairs, explain = design.run()

    pprint(pairs)
    assert pairs
    for pair in pairs.values():
        assert pair["PAIR"]["PRODUCT_SIZE"] == len(gfp)
        assert pair["RIGHT"]["OVERHANG"] == roverhang
        assert pair["LEFT"]["OVERHANG"] == loverhang
        assert pair["RIGHT"]["SEQUENCE"] == rc(gfp[-50:])

        lloc = pair["LEFT"]["location"]
        rloc = pair["RIGHT"]["location"]

        assert lloc[0] == 0

        assert rloc[0] == len(gfp) - 1
        assert rloc[1] == 50
    assert pairs


def test_long_left_primer_with_overhangs(gfp):

    design = Design()
    design.settings.template(gfp)
    design.settings.left_sequence(gfp[:50])
    design.settings.product_size((len(gfp), len(gfp)))

    loverhang = "TTTTTTT"
    roverhang = "GGGAGAG"

    design.settings.left_overhang(loverhang)
    design.settings.right_overhang(roverhang)
    design.settings.long_ok()
    design.settings.use_overhangs()
    design.settings.pick_anyway()
    pairs, explain = design.run()

    pprint(pairs)
    assert pairs
    for pair in pairs.values():
        assert pair["PAIR"]["PRODUCT_SIZE"] == len(gfp)
        assert pair["RIGHT"]["OVERHANG"] == roverhang
        assert pair["LEFT"]["OVERHANG"] == loverhang
        assert pair["LEFT"]["SEQUENCE"] == gfp[:50]

        lloc = pair["LEFT"]["location"]
        rloc = pair["RIGHT"]["location"]

        assert lloc[0] == 0
        assert lloc[1] == 50

        assert rloc[0] == len(gfp) - 1
    assert pairs


def test_set_long_overhang(gfp):
    design = Design()
    design.settings.template(gfp)
    design.settings.left_sequence(gfp[0:50])
    design.PRIMER_MAX_TM.value = 75.0
    design.PRIMER_MAX_SIZE = 35
    design.PRIMER_PICK_ANYWAY = 1
    design.settings.left_overhang("AAAAAAAAA")
    design.settings.long_ok()
    design.settings.use_overhangs()
    design.settings.pick_anyway()
    pairs, explain = design.run()

    print(explain)
    import json

    print(json.dumps(pairs[0]["LEFT"], indent=1))
    assert pairs
