=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com>`_.

0.7.0
------
Added
+++++
* Added ``allow_arbitrary_numbers`` to ``MutationRenumber`` to allow arbitrary strings as site numbers.

0.6.3
-----

Fixed
+++++
* Fix bug in handling ``minimap2`` errors ([see this issue](https://github.com/jbloomlab/alignparse/issues/99))
* Pass formatting with new ``black`` version
* Pass tests with new ``pandas`` version.
* Fixed ``simple_mut_consensus`` for newer versions of ``pandas`` when goruping by just one variable.

Changed
+++++++
* Change code linting to ``ruff`` rather than ``flake8``.
* Test with GitHub Actions rather than Travis CI.
* Remove ``mybinder`` examples.
* Test on Python 3.11 rather than 3.9.
* Don't allow ``pysam`` version 0.22.1 as it was causing some type of OPENSSL import error.
* Test with ``minimap2`` version 2.22

0.6.2
-----

Fixed
+++++
* Require python 3.8 or greater.


0.6.1
-----

Fixed
+++++
* Pass tests on newer dependency versions, test on python 3.9.

0.6.0
-----

Added
+++++
* Can parse reports from ``pbccs`` version 6.0.

0.5.0
-----

Added
+++++
* ``MutationRenumber`` accepts parameter ``allow_letter_suffixed_numbers`` to allow sites with lowercase letter suffixes (e.g, ``214a``)

0.4.1
-----

Fixed
+++++
* Remove deprecated ``pandas.DataFrame.append``

0.4.0
------

Added
++++++
* ``MutationRenumber`` accepts gap (``-``) and stop codo (``*``) characters.

* Started requring ``black`` formatting

0.3.0
------

Added
+++++
* ``InFrameDeletionsToSubs`` class.

* ``simple_mutconsensus`` takes ``None`` as parameter for ``max_sub_diffs`` and ``max_indel_diffs``

* ``simple_mutconsensus`` takes arguments ``min_support`` and ``max_minor_greater_or_equal``

0.2.6
-----

Fixed
+++++
* Fixed bug in ``add_mut_info_cols`` when index is not unique.

0.2.5
-----

Added
+++++
* ``merge_dels`` function that merges consecutive deletion strings.

0.2.4
-----

Fixed
+++++
* ``consensus.simple_mutconsensus`` works with just one grouping column.

0.2.3
-----

Fixed
+++++
* ``consensus`` module handles negative site numbers.

0.2.2
-----

Fixed
++++++
* ``sort_mutations`` and ``MutationRenumber`` handle negative site numbers.

0.2.1
-----

Added
+++++
* ``sort_mutations`` function to sort mutations by site, and optionally concatenate them.

0.2.0
------

Added
+++++
* ``MutationRenumber`` class to enable re-numbering of mutations.

0.1.6
------

Fixed
++++++
* Removed ``--for-only`` from  ``OPTIONS_VIRUS_W_DEL``

0.1.5
-----

Added
+++++
* ``Summaries`` now parses ``np`` (number of passes) tags from ``ccs`` version 5.0 FASTQs.

0.1.4
-----

Added
+++++
* ``Summaries`` now handles summaries from ``ccs`` version 5.0

0.1.3
------

Added
+++++
* Added ``hspace`` as option to ``Targets.plot``.

0.1.2
-----

Added
+++++
* Added ``select_target_names`` option to ``Targets``.

0.1.1
-----

Fixed
+++++
* Fixed DataFrame querying bug in ``./alignparse/ccs.py``

0.1.0
-----
Initial release

