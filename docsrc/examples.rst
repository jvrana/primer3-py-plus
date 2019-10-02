Examples
========

setting parameters
------------------

The preferred way to set parameters is to use the primer3plus.Design.presets() property.

using :class:`DesignPresets <primer3plus.design.DesignPresets>`
****************************************************************

.. testsetup::

    from primer3plus import Design

.. testcode::

    design = Design()
    design.presets.template("ACGGGGAGTTGTCTGTAGGTTGATTATGTGTGTCGTGTGTGTATATGGGTCTGA")
    print(design.get('SEQUENCE_TEMPLATE').value)

.. testoutput::

    ACGGGGAGTTGTCTGTAGGTTGATTATGTGTGTCGTGTGTGTATATGGGTCTGA

from a single key-value pair
********************************

.. code-block::

    # preferred
    design.set('SEQUENCE_TEMPLATE', 'foo')

    # alternative ways of setting parameters using a known parameter name
    design.SEQUENCE_TEMPLATE.value = 'foo'
    design.get('SEQUENCE_TEMPLATE').value = 'foo'
    design.params['SEQUENCE_TEMPLATE'] = 'foo'


from a dictionary
************************

.. code-block::

    design.update({
        'SEQUENCE_TEMPLATE': "AGGGGTAGTAGTATGTGAAGGGGTAGTAGTATGTGAAGGGGTAGTAGTATGTGAAGGGGTAGTAGTATGTGA",
        'LEFT_SEQUENCE': 'TAGTAGTATGTGAAGG'
    })

Design cloning primers
------------------------

.. testcode::

    import json

    design = Design()
    design.presets.template('TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGC')
    design.presets.as_cloning_task()
    design.presets.primer_num_return(1)
    results, explain = design.run()

    print(json.dumps(results, indent=1))
    print(json.dumps(explain, indent=1))

.. testoutput::

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

Design primers that target the region

Designing the right primer only
------------------------------

.. code-block::

    design = Design()
    design.presets.template("TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGC")
    design.presets.left_sequence('GTTATGTCACGCTTACATTCACG')
    design.presets.as_cloning_task()
    design.run()

Design primers targeting interval
---------------------------------

.. code-block::

    design = Design()
    design.presets.template("TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGC")
    design.presets.target((50, 150))
    design.run()

Relaxing parameters
-------------------

In this example, the parameter conditions are too strict to find a primer pair
the first time around:

.. testcode::

    design = Design()
    design.presets.template("TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGC")
    design.presets.target((25, 150))
    res, explain = design.run()
    print("Results: ", json.dumps(res, indent=1))
    print("Explain: ", json.dumps(explain, indent=1))

.. testoutput::

    Results:  {}
    Explain:  {
     "PRIMER_LEFT_EXPLAIN": "considered 36, low tm 36, ok 0",
     "PRIMER_RIGHT_EXPLAIN": "considered 515, low tm 238, high tm 104, ok 173",
     "PRIMER_PAIR_EXPLAIN": "considered 0, ok 0",
     "PRIMER_LEFT_NUM_RETURNED": 0,
     "PRIMER_RIGHT_NUM_RETURNED": 0,
     "PRIMER_INTERNAL_NUM_RETURNED": 0,
     "PRIMER_PAIR_NUM_RETURNED": 0
    }

We can run the relaxation procedure using :meth:`run_and_optimize <primer3plus.Design.run_and_optimize>`:

.. testcode::

    design = Design()
    design.presets.template("TCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTATTTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTTTTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGTTTTGGGACGCTCGAAGGCTTTAATTTGC")
    design.presets.target((25, 150))
    design.presets.primer_num_return(1)
    res, explain = design.run_and_optimize(5)
    print("Gradient used: ", design.DEFAULT_GRADIENT)
    print("Results: ", json.dumps(res, indent=1))
    print("Explain: ", json.dumps(explain, indent=1))

.. testoutput::

    Gradient used:  {'PRIMER_MAX_SIZE': (1, 27, 36), 'PRIMER_MIN_SIZE': (-1, 16, 27), 'PRIMER_MAX_TM': (1, 27, 80), 'PRIMER_MIN_TM': (-1, 48, 57.0), 'PRIMER_MAX_HAIRPIN_TH': (1, 47.0, 60)}
    Results:  {
     "0": {
      "PAIR": {
       "PENALTY": 7.892226720976964,
       "COMPL_ANY_TH": 0.0,
       "COMPL_END_TH": 0.0,
       "PRODUCT_SIZE": 235
      },
      "LEFT": {
       "PENALTY": 7.713186819997588,
       "SEQUENCE": "TCATGTAATTAGTTATGTCACGCT",
       "location": [
        0,
        24
       ],
       "TM": 56.28681318000241,
       "GC_PERCENT": 33.333333333333336,
       "SELF_ANY_TH": 0.0,
       "SELF_END_TH": 0.0,
       "HAIRPIN_TH": 0.0,
       "END_STABILITY": 5.07
      },
      "RIGHT": {
       "PENALTY": 0.17903990097937594,
       "SEQUENCE": "TTCGAGCGTCCCAAAACCTT",
       "location": [
        234,
        20
       ],
       "TM": 60.179039900979376,
       "GC_PERCENT": 50.0,
       "SELF_ANY_TH": 0.0,
       "SELF_END_TH": 0.0,
       "HAIRPIN_TH": 0.0,
       "END_STABILITY": 3.5
      }
     }
    }
    Explain:  {
     "PRIMER_LEFT_EXPLAIN": "considered 45, low tm 43, ok 2",
     "PRIMER_RIGHT_EXPLAIN": "considered 618, GC content failed 1, low tm 250, high tm 99, ok 268",
     "PRIMER_PAIR_EXPLAIN": "considered 1, ok 1",
     "PRIMER_LEFT_NUM_RETURNED": 1,
     "PRIMER_RIGHT_NUM_RETURNED": 1,
     "PRIMER_INTERNAL_NUM_RETURNED": 0,
     "PRIMER_PAIR_NUM_RETURNED": 1
    }