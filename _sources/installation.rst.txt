Installation
--------------

`alignparse <https://jbloomlab.github.io/alignparse/>`_ requires Python 3.8 or higher and is currently tested on Python 3.11.
A Mac OS X or Linux operating system is also required.

The easiest way to install `alignparse <https://jbloomlab.github.io/alignparse/>`_ is from `PyPI <https://pypi.org/>`_ using `pip <https://pip.pypa.io>`_ with::

    pip install alignparse

Note that installing the dependency `pysam <https://pysam.readthedocs.io/en/latest/api.html>`_ may require additional packages on the host system.
See the `pysam documentation <https://pysam.readthedocs.io/en/latest/installation.html#pypi-installation>`_ for more details.

In order to use the :mod:`alignparse.minimap2` module, you need to install the `minimap2 <https://github.com/lh3/minimap2>`_ executable (version 2.17 or higher).
You can do that via the `bioconda recipe <https://bioconda.github.io/recipes/minimap2/README.html>`_ or from the program's release page `as described here <https://github.com/lh3/minimap2#install>`_.

An ``environment.yml`` file is provided with the source code. This can be used to create a `conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ that installs `alignparse <https://jbloomlab.github.io/alignparse/>`_ and its dependencies using::

    conda env create -f environment.yml
