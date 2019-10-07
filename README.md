# Primer3Plus

Primer3Plus is a Python DNA primer design tool based off of Primer3 and the
Python primer3 wrapper (https://pypi.org/project/primer3-py/).

```python
import json

design = Design()
template = """
TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCAC
ATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAG
GTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTAT
TTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGC
ATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGC
TCGAAGGCTTTAATTTGC
"""
template = template.replace('\n', '').replace(' ', '')
design.settings.template(template)
design.settings.as_cloning_task()
design.settings.primer_num_return(1)
results, explain = design.run()

print(json.dumps(results, indent=1))
print(json.dumps(explain, indent=1))
```

```json
{
 "0": {
  "PAIR": {
   "PENALTY": 11.204301707622733,
   "COMPL_ANY_TH": 0.0,
   "COMPL_END_TH": 0.0,
   "PRODUCT_SIZE": 248
  },
  "LEFT": {
   "PENALTY": 9.027129166714644,
   "SEQUENCE": "TCATGTAATTAGTTATGTCACGCTTAC",
   "location": [
    0,
    27
   ],
   "TM": 57.972870833285356,
   "GC_PERCENT": 33.333333333333336,
   "SELF_ANY_TH": 0.0,
   "SELF_END_TH": 0.0,
   "HAIRPIN_TH": 0.0,
   "END_STABILITY": 2.34
  },
  "RIGHT": {
   "PENALTY": 2.1771725409080886,
   "SEQUENCE": "GCAAATTAAAGCCTTCGAGCG",
   "location": [
    247,
    21
   ],
   "TM": 58.82282745909191,
   "GC_PERCENT": 47.61904761904762,
   "SELF_ANY_TH": 0.0,
   "SELF_END_TH": 0.0,
   "HAIRPIN_TH": 38.006257959698985,
   "END_STABILITY": 5.03
  }
 }
}
{
 "PRIMER_LEFT_EXPLAIN": "considered 10, low tm 9, ok 1",
 "PRIMER_RIGHT_EXPLAIN": "considered 10, low tm 3, high tm 4, ok 3",
 "PRIMER_PAIR_EXPLAIN": "considered 1, ok 1",
 "PRIMER_LEFT_NUM_RETURNED": 1,
 "PRIMER_RIGHT_NUM_RETURNED": 1,
 "PRIMER_INTERNAL_NUM_RETURNED": 0,
 "PRIMER_PAIR_NUM_RETURNED": 1
}
```
## Installation

```
pip install primer3plus -U
```

## Requirements

python >= 3.5