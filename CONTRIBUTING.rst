=====================================
How to contribute to this package
=====================================

This document describes how to edit the package, run the tests, build the docs, put tagged versions on PyPI_, etc.

Editing the project
---------------------

Package structure
++++++++++++++++++
- Package code is in `alignparse <alignparse>`_
- Docs are in docs_
- `Jupyter notebooks`_ are in notebooks_
- Unit tests are in tests_

Modify the code via pull requests
+++++++++++++++++++++++++++++++++++
To make changes to the code, you should make a branch or fork, make your changes, and then submit a pull request.
If you aren't sure about pull requests:

 - A general description of pull requests: https://help.github.com/en/articles/about-pull-requests

 - How to create a pull request: https://help.github.com/en/articles/creating-a-pull-request

 - How to structure your pull requests (one conceptually distinct change per pull request): https://medium.com/@fagnerbrack/one-pull-request-one-concern-e84a27dfe9f1

Tests and documentation
+++++++++++++++++++++++
You should document your code clearly with `numpy style documentation`_.
You may also want to write sphinx_ documentation / examples in docs_ or the notebooks_ to demonstrate large-scale functionality.

You should add tests.
For simple things, these can be `doctests <https://docs.python.org/3/library/doctest.html>`_ in the code.
For more elaborate functionality, put unit tests in tests_.
Note also that the `Jupyter notebooks`_ in notebooks_ are tested via nbval_.
If you are getting errors on these notebook tests due to testing cells that output objects that can't be properly tested (such as widgets), see the *nbval-ignore-output* tag option discussed in the nbval_ docs.

Versions and CHANGELOG
++++++++++++++++++++++
The version is `single sourced <https://packaging.python.org/guides/single-sourcing-package-version/>`_ in `__init__.py`_.
When modifying a tagged version (e.g., ``0.1.0``), indicate you are working on a development version by adding a ``dev`` (e.g., ``0.1.dev1``).
See `here <https://www.python.org/dev/peps/pep-0440/>`_ for more information on version numbers.

Conceptual descriptions of changes should also be tracked in the CHANGELOG_.

Adding dependencies
+++++++++++++++++++++
When you add code that uses a new package that is not in the standard python library, you should add it to the dependencies specified under the `install_requires` option in `setup.py <setup.py>`_.
`See here <https://packaging.python.org/discussions/install-requires-vs-requirements/>`_ for information on how to do this, and how to specify minimal required versions.
As described in the above link, you should **not** pin exact versions in `install_requires` in `setup.py <setup.py>`_ unless absolutely necessary.

Notebooks on mybinder
-----------------------
The `Jupyter notebooks`_ in notebooks_ can be run interactively on mybinder_ by going to the following link:
https://mybinder.org/v2/gh/jbloomlab/alignparse/master?filepath=notebooks

In order for this to work, you need to keep the `environment.yml <environment.yml>`_ configuration file up to date with the dependencies for running these notebooks as `described here <https://mybinder.readthedocs.io/en/latest/config_files.html>`_.
Note that unlike for the `install_requires` in `setup.py <setup.py>`_, you may want to pin exact versions here to get reproducible installations.
Look into the `pip freeze <https://pip.pypa.io/en/stable/reference/pip_freeze/>`_ and `conda env export <https://packaging.python.org/discussions/install-requires-vs-requirements>`_ commands on how to automatically create such a configuration file.

Testing
---------

Adding tests
++++++++++++++
As you add new codes, you should create tests to make sure it is working correctly.
These can include:

  - doctests in the code

  - unit tests in the `./tests/ <tests>`_ subdirectory

Running the tests locally
++++++++++++++++++++++++++
After you make changes, you should run two sets of tests.
To run the tests, go to the top-level package directory.
Then make sure that you have installed the packages listed in `test_requirements.txt <test_requirements.txt>`_.
If these are not installed, install them with::

    pip install -r test_requirements.txt

Then use flake8_ to `lint the code <https://en.wikipedia.org/wiki/Lint_%28software%29>`_ by running::

    flake8

If you need to change the flake8_ configuration, edit the `.flake8 <.flake8>`_ file.

Then run the tests with pytest_ by running::

    pytest

If you need to change the pytest_ configuration, edit the `pytest.ini <pytest.ini>`_ file.

Automated testing on Travis
+++++++++++++++++++++++++++
The aforementioned flake8_ and pytest_ tests will be run automatically by the Travis_ continuous integration system as specified in the `.travis.yml <.travis.yml>`_ file.
Note that running the Travis_ tests requires you to register the project with Travis_.

If the tests are passing, you will see this on the Travis_ badge on GitHub repo main page.

Slack notifications of test results
+++++++++++++++++++++++++++++++++++++
You can configure Travis_ to provide automatic Slack notifications of the test results.
To do that, follow the `instructions here <https://docs.travis-ci.com/user/notifications/#configuring-slack-notifications>`_.


Building documentation
------------------------
See `docs/README.rst <docs/README.rst>`_ for information on how to build the documentation.

Tagging versions and putting on PyPI
-------------------------------------
When you have a new stable release, you will want to tag it and put it on PyPI_ where it can be installed with pip_.
First, make sure the version number is up-to-date in `__init__.py`_ and the CHANGELOG_.
Then commit the code to GitHub if you haven't already done so.
Next tag the version, as in::

    git tag -a 0.1.0 -m 'version 0.1.0'

and then push the tag to GitHub with::

    git push --tags

Finally, upload to PyPI_ with twine_ as `described here <https://github.com/pypa/twine>`_.
Note that this requires you to have registered the package on PyPI_ if this is the first version of the package there.

.. _pytest: https://docs.pytest.org
.. _flake8: http://flake8.pycqa.org
.. _Travis: https://docs.travis-ci.com
.. _PyPI: https://pypi.org/
.. _pip: https://pip.pypa.io
.. _sphinx: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _tests: tests
.. _docs: docs
.. _notebooks: notebooks
.. _`Jupyter notebooks`: https://jupyter.org/
.. _`__init__.py`: alignparse/__init__.py
.. _CHANGELOG: CHANGELOG.rst
.. _twine: https://github.com/pypa/twine
.. _`numpy style documentation`: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
.. _nbval: https://nbval.readthedocs.io
.. _mybinder: https://mybinder.readthedocs.io
