import pytest

from primer3plus.p3p import Primer3Design


@pytest.mark.parametrize("parse", [True, False])
def test_benchmark_basic_design(benchmark, parse, gfp):
    designer = Primer3Design()
    print(designer.params["PRIMER_OPT_SIZE"])

    benchmark(
        designer.design,
        {"SEQUENCE_TEMPLATE": gfp, "SEQUENCE_INCLUDED_REGION": [100, 300]},
        {"PRIMER_PICK_ANYWAY": False},
        parse=parse,
    )
