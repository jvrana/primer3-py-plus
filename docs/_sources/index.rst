:github_url: |homepage|


Primer3Plus
===========
v\ |version|. (:doc:`Changelog <developer/changelog>`)

Github: |homepage|

Primer3Plus is a Python DNA primer design tool based off of primer3 [1]_ [2]_ and the
Python primer3 wrapper https://pypi.org/project/primer3-py/.

.. code-block:: python

   design = Design()
   design.presets.template("AGCGTCGTGTATGGTAGTGTATTTGCGTTGACGTTGCTGACGTCGTTGAGTCGT")
   design.presets.as_cloning_task()
   design.run()

.. code-block:: python

   design = Design()
   design.presets.template("AGCGTCGTGTATGGTAGTGTATTTGCGTTGACGTTGCTGACGTCGTTGAGTCGT")
   design.presets.as_cloning_task()
   design.presets.left_sequence("AGCGTCGTGTATGGTAGTG")
   design.run_and_optimize(5)

For more info, take a look at the :ref:`API Docs <api>`

For a list of design parameters available, take a look at the
:ref:`BoulderIO Parameters <api_default_parameters>`

For more examples, see the :doc:`Usage Docs <usage>` and some of the following:

- :ref:`Cloning primers <cloning_primers>`
- :ref:`Setting primer sequences <setting_primers>`
- :ref:`Handling primers with overhangs <handle_overhangs>`
- :ref:`Handling long sequences <handle_long_sequences>`
- :ref:`Auto-relaxation of design params <autorelax>`

Getting started
---------------

Installation
************

.. code-block::

   pip install primer3plus -U

Making a primer design task
***************************

The preferred way to start a design is to use the :meth:`primer3plus.Design.presets`
property, which provides many helper methods for setting parameters.
If you're using an interactive notebook or terminal, hitting `TAB` after typing
`design.presets` will provide all of the standard design settings

.. code-block::

   design = Design()
   design.presets.left_sequence('AGCGTCGTGTATGGTAGTG')
   design.presets.template('AGGGGCGGAGGTGTAGTCGTCGTTAGCGTTAGTCTA')

You can also interact with the parameters more directly. To get the parameters, use
the :meth:`get <primer3plus.Design.get>` method. From there you can set the value,
print parameter information, set defaults

.. code-block::

   param = design.get('LEFT_SEQUENCE')

   # print parameter information
   print(param)

   # print help url
   print(param.help())

   # set a new value
   param.value = "AGGTAGTAT"

   # set back to the default value
   param.set_default()

You can set the value of the parameter using `value` or using the
:meth:`set <primer3plus.Design.set` method:

.. code-block::

   design.set('LEFT_SEQUENCE', 'AGGTATTAGTATATGATAT')

You can also access descriptors of the parameters. Since there
are many parameters, this is useful in interactive environments:

.. code-block:: python

   design.SEQUENCE_ID.help()
   design.SEQUENCE_ID.value = "Seq11123123"

To view all of the parameters as a dictionary, use:

.. code-block::

   print(dict(design.params))

Running design tasks
********************

To run, simply use the :meth:`run <primer3plus.Design.run>` method:

.. code-block::

   design.run()

Provided is also a primer3 relaxation procedure, which will gradually relax
parameters such as melting temperature and primer size until primer pairs are found:

.. code-block::

   design.run_and_optimize(5)  # run for max of 5 iterations

Table of Contents
=================

.. toctree::
   :maxdepth: 3

   usage
   developer/api_reference
   developer/changelog


References
==========

.. [1] Untergasser A, Cutcutache I, Koressaar T, et al. Primer3--new capabilities and
interfaces. Nucleic Acids Res. 2012;40(15):e115. doi:10.1093/nar/gks596

.. [2] http://primer3.ut.ee/primer3web_help.htm