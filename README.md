# Primer3Plus

Primer3Plus is a Python DNA primer design tool based off of Primer3 [1]_ [2]_ and the
Python primer3 wrapper (https://pypi.org/project/primer3-py/).

Documentation: https://jvrana.github.io/primer3-py-plus/

```python
design = Design()
design.presets.template("AGCGTCGTGTATGGTAGTGTATTTGCGTTGACGTTGCTGACGTCGTTGAGTCGT")
design.presets.as_cloning_task()
design.presets.left_sequence("AGCGTCGTGTATGGTAGTG")
design.run_and_optimize(5)
```