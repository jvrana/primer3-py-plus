from primer3plus.design import Design
import random
from primer3plus.utils import reverse_complement, anneal


def test_run_and_optimize_basic(gfp, iter_random_primer):
    fwd = gfp[:25]
    rev = reverse_complement(gfp[-25:])
    design = Design()
    design.set.template(gfp)
    design.set.left_sequence(fwd)
    design.set.right_sequence(rev)
    design.set.task("check_primers")
    pairs, explain = design.run()
    assert not pairs

    pairs, explain = design.run_and_optimize(15)
    print(explain)
    assert pairs
