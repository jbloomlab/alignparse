===========================
Documentation
===========================

This subdirectory contains the reStructuredText_ documentation that can be built with sphinx_.

Build the documentation
-----------------------------

To build the documentation, you need ensure that you have installed the packages in `doc_requirements.txt <doc_requirements.txt>`_.
If these packages are not installed, you can install them with::

    pip install -r doc_requirements.txt

You also may need to install ``pandoc``:

    mamba install pandoc

Then simply type::

    make html

and the HTML documentation will be built in ``_build/html/``.
You can open the HTML files in that directory to see how the docs look.


Jupyter notebooks as examples
------------------------------
The documentation is set up to include `Jupyter notebooks`_ via nbsphinx_.
Assuming you have them in the `../notebooks/ <../notebooks/>`_ subdirectory, you need to use nbsphinx-link_ to include them.
For instance, if you have a notebook at the path ``../notebooks/example_notebook.ipynb``, then you would first create a file in this docs directory called ``example_notebook.nblink`` with the following contents::

    {
        "path": "../notebooks/example_notebook.ipynb"
    }

You would then simply include the notebook by adding ``example_notebook`` to one of your ``*.rst`` files (e.g., `examples.rst <examples.rst>`_) the same way that you would link to another ``*.rst`` file.

Pushing docs to GitHub pages
------------------------------
The docs are hosted on `GitHub pages`_ at https://jbloomlab.github.io/alignparse

After you have built the documentation as described above, here is how you push it to `GitHub pages`_:

Building the documentation for the first time
+++++++++++++++++++++++++++++++++++++++++++++++
If you are building the docs for the **very first time**, then there will **not** be any branch called ``gh-pages`` on the project's GitHub repo.
In that case, you need to create that branch in the ``_build/html/`` subdirectory.

First, clone the repo into ``_build/html/``::

    mkdir _build
    git clone https://github.com/jbloomlab/alignparse _build/html

Then create a new branch called ``gh-pages`` in this subdirectory::

    cd _build/html
    git checkout -b gh-pages

Now go back to the top directory of the docs, remove all the existing files in ``_build/html``, and build the sphinx_ docs into this directory::

    cd ../../
    make html

Now got into ``_build/html`` and commit this initial version of your docs::

    cd _build/html
    git add .
    git commit -m 'initial build of `gh-pages` docs'

Finally, push the ``gh-pages`` branch to GitHub::

    git push -u origin gh-pages

After a few minutes, your docs should now show up at https://jbloomlab.github.io/alignparse

Building a new version of the docs
++++++++++++++++++++++++++++++++++
If you have already built and earlier version of the docs that ``gh-pages`` already exists, then you simply need to commit a new version.

First, remove any existing docs::

    make clean

Now clone the existing ``gh-pages`` branch in ``_build/html``::

    git clone -b gh-pages https://github.com/jbloomlab/alignparse _build/html

Next, build the docs::

    make html

Finally, commit and push this new version of the docs::

    cd _build/html
    git add --all
    git commit -m 'new version of docs'
    git push

.. _reStructuredText: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _`GitHub pages`: https://help.github.com/en/articles/what-is-github-pages>
.. _sphinx: http://www.sphinx-doc.org
.. _nbsphinx: https://nbsphinx.readthedocs.io
.. _nbsphinx-link: https://github.com/vidartf/nbsphinx-link
.. _`Jupyter notebooks`: https://jupyter.org/
