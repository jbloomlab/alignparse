Installation
--------------

``alignparse`` requires Python 3.6 or 3.7 and a Mac OS X or Linux operating system.

The easiest way to install ``alignparse`` is from `PyPI <https://pypi.org/>`_ using `pip <https://pip.pypa.io>`_ with::

    pip install alignparse

Note that installing the dependency ``pysam`` may require additional packages on the host system.
See the `pysam documentation <https://pysam.readthedocs.io/en/latest/installation.html#pypi-installation>`_ for more details.

The source code for ``alignparse`` is available on GitHub at https://github.com/jbloomlab/alignparse.

In order to use the :mod:`alignparse.minimap2` module, you need to install the `minimap2 <https://github.com/lh3/minimap2>`_ executable.
You can do that via the `bioconda recipe <https://bioconda.github.io/recipes/minimap2/README.html>`_ or from the program's release page `as described here <https://github.com/lh3/minimap2#install>`_.
