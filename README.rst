===============================
alignparse
===============================

.. image:: https://img.shields.io/pypi/v/alignparse.svg
        :target: https://pypi.python.org/pypi/alignparse

.. image:: https://img.shields.io/travis/jbloomlab/alignparse.svg
        :target: https://travis-ci.com/jbloomlab/alignparse

.. image:: https://mybinder.org/badge_logo.svg
        :target: https://mybinder.org/v2/gh/jbloomlab/alignparse/master?filepath=notebooks

.. image:: https://zenodo.org/badge/194140958.svg
   :target: https://zenodo.org/badge/latestdoi/194140958

.. image:: https://joss.theoj.org/papers/10.21105/joss.01915/status.svg
   :target: https://doi.org/10.21105/joss.01915

``alignparse`` is a Python package written by `the Bloom lab <https://research.fhcrc.org/bloom/en.html>`_. 
It is designed to align long sequencing reads (such as those from PacBio circular consensus sequencing) to targets, filter these alignments based on user-provided specifications, and parse out user-defined sequence features.
For each read that passes the filters, information about the features (e.g. accuracy, sequence, mutations) is retained for further analyses. 

See the `alignparse documentation <https://jbloomlab.github.io/alignparse>`_ for details on how to install and use ``alignparse``.

Please `cite alignparse <https://jbloomlab.github.io/alignparse/acknowledgements.html>`_ if you use it in your work.

The source code is `on GitHub <https://github.com/jbloomlab/alignparse>`_.

To contribute to this package, read the instructions in `CONTRIBUTING.rst <CONTRIBUTING.rst>`_.

To report issues or seek support, please use the `GitHub issue tracker <https://github.com/jbloomlab/alignparse/issues>`_.
