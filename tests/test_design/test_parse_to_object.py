from primer3plus import Design
from primer3plus.design.results import to_pair_result


def test_as_str(gfp):
    design = Design()
    design.presets.template(gfp)
    design.presets.included(gfp[100:300])
    pairs, explain = design.run()
    assert pairs
    for v in pairs.values():
        pair = to_pair_result(v)
        print(pair)
        print(pair.thermo_tm())
