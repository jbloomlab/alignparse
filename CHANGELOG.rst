=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com>`_.

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

