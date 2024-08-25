``alignparse`` documentation
===================================================

`alignparse <https://jbloomlab.github.io/alignparse/>`_ is a Python package written by `the Bloom lab <https://jbloomlab.org>`_.
It is designed to align long sequencing reads (such as those from PacBio circular consensus sequencing) to `targets`, filter these alignments based on user-provided specifications, and parse out user-defined sequence `features`. 

For each read that passes the user-defined filters, the user can define what information about each feature (e.g. sequence, mutations, and sequencing accuracy) should be retained for further analyses. The `alignparse.consensus` module provides tools for such further analyses, including analyzing mutations identified in sequence features and defining consensus sequences from barcoded sequencing reads.

The `alignparse source code <https://github.com/jbloomlab/alignparse>`_ is freely available on GitHub.

If you use `alignparse <https://jbloomlab.github.io/alignparse/>`_ for your work, please cite the references in `Acknowledgements <https://jbloomlab.github.io/alignparse/acknowledgements.html>`_.

See below for information and examples of how to use this package.

Contents
----------
.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 2

   installation
   examples
   alignparse
   package_index
   acknowledgements
