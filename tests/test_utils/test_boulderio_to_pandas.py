import sys

import pytest

from primer3plus import Design

if sys.version_info[0] == 3 and sys.version_info[1] <= 6:
    pytest.skip("test requires version 3.6", allow_module_level=True)


def test_to_pandas():
    import pandas as pd

    boulderio = Design.DEFAULT_PARAMS
    rows = []
    for k, v in boulderio._params.items():
        rows.append(
            {
                "name": v.name,
                "type": str(v.ptype.type),
                "value": v.value,
                "default": v.ptype.default,
            }
        )
    return pd.DataFrame(rows)
