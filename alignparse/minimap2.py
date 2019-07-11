"""
===============
minimap2
===============

Runs `minimap2 <https://lh3.github.io/minimap2/>`_ aligner.

"""


import gzip
import os
import subprocess
import tempfile

import packaging.version


OPTIONS_CODON_DMS = ('-A2',
                     '-B4',
                     '-O12',
                     '-E2',
                     '--end-bonus=13',
                     '--score-N=4',
                     '--secondary=no',
                     '--cs',
                     )
"""tuple: ``minimap2`` options for codon-mutant libraries.

Note
----
These options work well for sequences expected to have codon-level mutations
(3 nucleotides in a row) and not expected to have many indels.

Options have the following meaning / rationale:
  - ``-A2`` : score for a match
  - ``-B4`` : penalty for a mismatch
  - ``-O12`` : gap opening penalty (high, because few indels expected)
  - ``-E2`` : gap extension penalty
  - ``--end-bonus=13`` : favor codon mutation over gap at termini.
  - ``--score-N=4`` : mismatch with ambiguous base as bad as other mismatches.
  - ``--secondary=no`` : do not output secondary alignments.
  - ``--cs`` : add short ``cs`` tag ((see https://github.com/lh3/minimap2#cs)

"""

OPTIONS_VIRUS_W_DEL = ('-xsplice',
                       '-un',
                       '-C0',
                       '--splice-flank=no',
                       '-M=1',
                       '--for-only',
                       '--end-seed-pen=2',
                       '--end-bonus=1',
                       '--secondary=no',
                       '--cs',
                       )
"""tuple: ``minimap2`` options for viral genes.

Note
----
These options work well for sequences expected to have nucleotide mutations
and long internal deletions. It handles long deletions by having ``minimap2``
treat them like introns. This means that in the output, introns should be
parsed like deletions.

Options have the following meaning / rationale:
  - ``-xsplice`` : pre-set options to long-read spliced alignment.
  - ``-un`` : do not try to find canonical splice sites.
  - ``-C0`` : no cost for non-canonical splice sites.
  - ``--splice-flank=no`` : no assumptions about base next to splice donor.
  - ``-M=1`` : mark secondary chain that overlaps with another by this much.
  - ``--end-seed-pen=2`` : described as helping avoid tiny terminal exons.
  - ``--end-bonus=1`` : bonus for extending to end of alignment.
  - ``--secondary=no`` : do not output secondary alignments.
  - ``--cs`` : add short ``cs`` tag ((see https://github.com/lh3/minimap2#cs)

"""


class Mapper:
    """Run `minimap2 <https://lh3.github.io/minimap2/>`_.

    Parameters
    ----------
    options : tuple or list
        Command line options to ``minimap2``. For instance, see
        :data:`OPTIONS_CODON_DMS` or :data:`OPTIONS_VIRUS_W_DEL`.
    prog : str
        Path to ``minimap2`` executable.
    min_version : str
        Minimum version of ``minimap2``. :class:`Mapper` has only
        been tested with versions >= 2.12.
    check_cs : bool
        Check that `options` includes ``--cs`` to output the short ``cs`` tag.

    Attributes
    ----------
    options : list
        Options to ``minimap2``.
    prog : str
        Path to ``minimap2`` executable.
    version : str
        Version of ``minimap2``.

    """

    def __init__(self, options, *, prog='minimap2', min_version='2.12',
                 check_cs=True):
        """See main :class:`Mapper` doc string."""
        try:
            version = subprocess.check_output([prog, '--version'])
        except Exception:
            raise ValueError(f"Can't execute `prog` {prog}. Is it installed?")
        self.version = version.strip().decode('utf-8')
        min_version = packaging.version.parse(min_version)
        if packaging.version.parse(self.version) < min_version:
            raise ValueError(f"You have `minimap2` version {0}, need >= {1}")
        self.prog = prog
        self.options = list(options)
        if check_cs and not any('--cs' == opt for opt in self.options):
            raise ValueError(f"`options` do not include `--cs`:\n{options}")

    def map_to_sam(self, targetfile, queryfile, samfile):
        """Map query sequences to targets and output SAM file.

        Parameters
        ----------
        targetfile : str
            Alignment targets (reference), typically FASTA (can be gzipped).
        queryfile : str
            Queries to align, typically FASTA or FASTQ (can be gzipped).
        samfile : str
            Name of created SAM file. If it has the extension ``.gz`` then
            the created file is gzipped.

        """
        for fname, f in [('target', targetfile), ('query', queryfile)]:
            if not os.path.isfile(f):
                raise IOError(f"cannot find `{fname}file` {fname}")

        if not any('-a' == opt for opt in self.options):
            options = self.options + '-a'

        cmds = [self.prog] + options + [targetfile, queryfile]

        if os.path.splitext(samfile)[0] == '.gz':
            sam = gzip.open(samfile, 'w')
        else:
            sam = open(samfile, 'w')
        err = tempfile.TemporaryFile(mode='w')
        try:
            _ = subprocess.check_call(cmds, stdout=sam, stderr=err)
        except Exception as exc:
            exc_type = type(exc).__name__
            sam.close()
            os.remove(samfile)
            err.seek(0)
            raise RuntimeError(f"Error running these commands:\n{cmds}\n\n"
                               f"Exception message:\n{exc_type}: {exc}\n\n"
                               f"Standard error message:\n{err.read()}\n\n"
                               f"Removing the output SAM file {samfile}")
        finally:
            err.close()
        sam.close()
        assert os.path.isfile(samfile)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
