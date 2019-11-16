===============================
alignparse
===============================

.. image:: https://img.shields.io/pypi/v/alignparse.svg
        :target: https://pypi.python.org/pypi/alignparse

.. image:: https://img.shields.io/travis/jbloomlab/alignparse.svg
        :target: https://travis-ci.org/jbloomlab/alignparse

.. image:: https://mybinder.org/badge_logo.svg
        :target: https://mybinder.org/v2/gh/jbloomlab/alignparse/master?filepath=notebooks

Align sequences and then parse features.

``alignparse`` is a Python package written by `the Bloom lab <https://research.fhcrc.org/bloom/en.html>`_. 

This package is designed to align long sequencing reads (such as those from PacBio circular consensus sequencing) to `targets`, filter these alignments based on user-provided specifications, and parse out user-defined sequence `features`. For each read that passes the user-defined filters, information about the features (e.g. accuracy, mutations) is retained for further analyses.

The source code is `on GitHub <https://github.com/jbloomlab/alignparse>`_.

See the `alignparse documentation <https://jbloomlab.github.io/alignparse>`_ for details on how to install and use ``alignparse``.

To contribute to this package, read the instructions in `CONTRIBUTING.rst <CONTRIBUTING.rst>`_.
