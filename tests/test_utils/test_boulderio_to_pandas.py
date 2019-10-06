import pandas as pd

from primer3plus import Design


def test_to_pandas():
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
