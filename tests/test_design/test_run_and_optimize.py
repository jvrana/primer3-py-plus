import random

from primer3plus.design import Design
from primer3plus.utils import anneal
from primer3plus.utils import reverse_complement


def test_run_and_optimize_basic(gfp, iter_random_primer):
    fwd = gfp[:25]
    rev = reverse_complement(gfp[-25:])
    design = Design()
    design.settings.template(gfp)
    design.settings.left_sequence(fwd)
    design.settings.right_sequence(rev)
    design.settings.task("check_primers")
    pairs, explain = design.run()
    assert not pairs

    pairs, explain = design.run_and_optimize(15)
    print(explain)
    assert pairs
