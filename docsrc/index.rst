.. DASi documentation master file, created by
   sphinx-quickstart on Sun Nov 19 22:18:51 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Primer3Plus
===========

Release v\ |version|. (:doc:`Changelog <developer/changelog>`)

Github |gitpage|

Primer3Plus is a Python DNA primer design tool based off of primer3 [1]_ [2]_ and the
Python primer3 wrapper https://pypi.org/project/primer3-py/

.. code-block:: python

   design = Design()
   design.set.template("AGCGTCGTGTATGGTAGTGTATTTGCGTTGACGTTGCTGACGTCGTTGAGTCGT")
   design.set.as_cloning_task()
   design.run()

.. code-block:: python

   design = Design()
   design.set.template("AGCGTCGTGTATGGTAGTGTATTTGCGTTGACGTTGCTGACGTCGTTGAGTCGT")
   design.set.as_cloning_task()
   design.set.left_sequence("AGCGTCGTGTATGGTAGTG")
   design.run_and_optimize(5)

Installation
------------

.. code-block::

   pip install primer3plus -U

References
----------

.. [1] Untergasser A, Cutcutache I, Koressaar T, et al. Primer3--new capabilities and
interfaces. Nucleic Acids Res. 2012;40(15):e115. doi:10.1093/nar/gks596

.. [2] http://primer3.ut.ee/primer3web_help.htm

Table of Contents
=================

.. toctree::
   :maxdepth: 3

   examples
   developer/api_reference
   developer/changelog
