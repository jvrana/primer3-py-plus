import pytest
import json
from primer3plus.p3p import Primer3Params, Primer3Design


def test_primer3_params_init():
    params = Primer3Params()
    # print(params.GLOBAL)
    print(params.GLOBAL.PRIMER_EXPLAIN_FLAG)


def test_set():
    params = Primer3Params()
    params["SEQUENCE_TEMPLATE"] = "agtacaga"
    assert params["SEQUENCE_TEMPLATE"] == "agtacaga"


def test_update_values():
    params = Primer3Params()
    params.update(
        {"SEQUENCE_TEMPLATE": "gggaggagag", "SEQUENCE_INCLUDED_REGION": [500, 700]}
    )

    assert params["SEQUENCE_TEMPLATE"] == "gggaggagag"
    assert params["SEQUENCE_INCLUDED_REGION"] == [500, 700]


def test_df():
    params = Primer3Params()
    params.df()


@pytest.mark.parametrize("parse", [True, False])
def test_basic_design(parse, gfp):
    designer = Primer3Design()
    print(designer.params["PRIMER_OPT_SIZE"])
    results = designer.design(
        {"SEQUENCE_TEMPLATE": gfp, "SEQUENCE_INCLUDED_REGION": [100, 300]},
        {"PRIMER_PICK_ANYWAY": False},
        parse=parse,
    )
    print(results)


def test_check_primers():
    designer = Primer3Design()
    p1 = "atgatgaagttgccgccct"
    p2 = "tcccaattcttgttgaattagatggtgat"
    pairs, other = designer.check_primers(p1, p2)
    print(pairs)


def test_pick_cloning_primers(gfp):
    designer = Primer3Design()
    pairs, other = designer.pick_cloning_primers(gfp, max_iterations=15)
    print(pairs)
    print(other)
    assert pairs


def test_pick_pcr_primers(gfp):
    designer = Primer3Design()
    pairs, other = designer.pick_pcr_primers(
        gfp, (300, 400), (200, 1000), max_iterations=15
    )
    for p in pairs:
        print(pairs[p])
    assert pairs


def test_pick_right_pcr_primer(gfp):
    designer = Primer3Design()
    pairs, other = designer.pick_sequencing_primers(gfp, "actggagttgtcccaattc")
    for p in pairs:
        print(pairs[p])
    assert pairs


def test_pick_left_pcr_primer(gfp):
    designer = Primer3Design()
    right_primer = "ctatttgtatagttcat"
    pairs, other = designer.pick_left_pcr_primer(
        gfp, right_primer, [], (100, 1000), max_iterations=15
    )
    for p in pairs:
        print(pairs[p])


def test_pick_sequencing_primer(gfp):
    designer = Primer3Design()
    pairs, other = designer.pick_sequencing_primers(
        gfp, "gagttgtcccaattcttgttgaattagat"
    )
    print(pairs)


def test_pick_left_sequencing_primer(gfp):
    designer = Primer3Design()
    pairs, other = designer.pick_sequencing_primers(
        gfp, "gatgaagttgccgccctcg", left_only=True
    )
    print(other)
    print(pairs)


def test_pick_right_sequencing_primer(gfp):
    designer = Primer3Design()
    pairs, other = designer.pick_sequencing_primers(
        gfp, "gagttgtcccaattcttgttgaattagat", right_only=True
    )
    print(pairs)
    print(other)


def test_pick_right_pcr_primers(gfp):
    designer = Primer3Design()

    primers = [gfp[20:40], gfp[60:85], gfp[120:145]]

    pairs, other = designer.pick_right_pcr_primer(gfp, primers, size_range=(200, 300))
    print(pairs)
    print(other)


def test_check_many_primers(gfp):
    designer = Primer3Design()

    lprimers = [gfp[20:40], gfp[60:85], gfp[120:145]]

    rprimers = ["gtagtgacaagtgttggccatgga"]
    results = designer.check_primers(lprimers, rprimers)
    print(results)


def test_check_pcr_primers(gfp):
    designer = Primer3Design()

    lprimers = [gfp[20:40], gfp[60:85], gfp[120:145]]

    rprimers = ["tcaacaagaattgggac", "gtagtgacaagtgttggccatgga"]
    pairs, other = designer.check_pcr_primers(gfp, lprimers, rprimers, size_range=(40, 400))
    print(json.dumps(pairs, indent=2))
    print(json.dumps(other, indent=2))
