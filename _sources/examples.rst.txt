=========
Examples
=========

The `alignparse <https://jbloomlab.github.io/alignparse/>`_ package is designed to aid in the analysis of long-read sequencing of user-defined amplicons. Major functionaly includes:
    - The :class:`alignparse.targets.Targets` class for defining reference sequence `Targets`, aligning reads, and parsing features from aligned reads.
    - The :mod:`alignparse.consensus` for analyzing mutations and sequencing accuracy, especially for sequences grouped by shared barcodes.

Here are some examples that demonstrate how to utilize this functionality and other aspects of the `alignparse <https://jbloomlab.github.io/alignparse/>`_ package:

.. toctree::
   :maxdepth: 1

   recA_DMS
   lasv_pilot
   flu_virus_seq_example

The actual notebooks and input files used in these examples are `available here <https://github.com/jbloomlab/alignparse/tree/master/notebooks>`_.

The above examples can also be run as interactive Jupyter notebooks on `mybinder <https://mybinder.readthedocs.io>`_ by going to the `following link <https://mybinder.org/v2/gh/jbloomlab/alignparse/master?filepath=notebooks>`_ (it may take a minute to load) and then opening the notebook you want to run.
