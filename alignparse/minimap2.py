"""
===============
minimap2
===============

Runs `minimap2 <https://lh3.github.io/minimap2/>`_ aligner.
``minimap2`` must be version 2.17 or higher.

"""

import contextlib
import os
import re
import subprocess
import tempfile
import textwrap  # noqa: F401

import packaging.version

import pysam


OPTIONS_CODON_DMS = ('-A2',
                     '-B4',
                     '-O12',
                     '-E2',
                     '--end-bonus=13',
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
  - ``--secondary=no`` : do not output secondary alignments.
  - ``--cs`` : add short ``cs`` tag ((see https://github.com/lh3/minimap2#cs)

"""

OPTIONS_VIRUS_W_DEL = ('-xsplice:hq',
                       '-un',
                       '-C0',
                       '--splice-flank=no',
                       '-M=1',
                       '--end-seed-pen=2',
                       '--end-bonus=2',
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
  - ``-xsplice:hq`` : preset for long-read high-quality spliced alignment.
  - ``-un`` : do not try to find canonical splice sites.
  - ``-C0`` : no cost for non-canonical splice sites.
  - ``--splice-flank=no`` : no assumptions about base next to splice donor.
  - ``-M=1`` : mark secondary chain that overlaps with another by this much.
  - ``--end-seed-pen=2`` : described as helping avoid tiny terminal exons.
  - ``--end-bonus=2`` : bonus for extending to end of alignment.
  - ``--secondary=no`` : do not output secondary alignments.
  - ``--cs`` : add short ``cs`` tag ((see https://github.com/lh3/minimap2#cs)

"""


class Mapper:
    r"""Run `minimap2 <https://lh3.github.io/minimap2/>`_.

    Parameters
    ----------
    options : tuple or list
        Command line options to ``minimap2``. For instance, see
        :data:`OPTIONS_CODON_DMS` or :data:`OPTIONS_VIRUS_W_DEL`.
    prog : str
        Path to ``minimap2`` executable.
    min_version : str
        Minimum version of ``minimap2``. :class:`Mapper` has only
        been tested with versions >= 2.17.
    check_cs : bool
        Check that `options` includes ``--cs`` to output the short ``cs`` tag.
    retain_tags : None or list
        Retain these tags in query sequences in output alignments.

    Attributes
    ----------
    options : list
        Options to ``minimap2``.
    prog : str
        Path to ``minimap2`` executable.
    version : str
        Version of ``minimap2``.
    retain_tags : None or list
        List of tags to retain.

    Note
    -----
    If `retain_tags` is set, then :meth:`Mapper.map_to_sam` searches for tags
    in query FASTQ headers like the `np`` tags here::

        @m54228_181120_212724/4194377/ccs   np:i:45

    and retains them in the `samfile` alignment output. Such tags
    are present in FASTQ files created by PacBio ``ccs`` version 4.0. For
    the above header, you'd use `retain_tags=['np']`.

    Examples
    --------
    First, create toy example target and query file:

    >>> tempdir = tempfile.TemporaryDirectory()
    >>> targetfile = os.path.join(tempdir.name, 'targets.fasta')
    >>> queryfile = os.path.join(tempdir.name, 'queries.fastq')
    >>> with open(targetfile, 'w') as f:
    ...     _ = f.write(textwrap.dedent('''
    ...         >refseq
    ...         ATGCAANNNGATGCATAGTATTAGCATAAATAGGATAGCCATAAGGTTACTGCATAAGGGTAT
    ...         ''').strip())
    >>> with open(queryfile, 'w') as f:
    ...     _ = f.write(textwrap.dedent('''
    ...         @m54228_181120_212724/4194376/ccs   np:i:127
    ...         ATGCAAAATGATGCATAGTATTAGCATAAATAGGATAGCCATAAGGTTACTGCATAAGAGTAT
    ...         +
    ...         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ...         ''').strip())

    Now map **not** retaining tags from the FASTQ (this is default behavior):

    >>> samfile = os.path.join(tempdir.name, 'alignments.sam')
    >>> mapper = Mapper(OPTIONS_CODON_DMS)
    >>> mapper.map_to_sam(targetfile, queryfile, samfile)
    >>> for a in pysam.AlignmentFile(samfile):
    ...     tag_names = [tup[0] for tup in a.get_tags()]
    ...     print(a.tostring())  # doctest: +NORMALIZE_WHITESPACE
    m54228_181120_212724/4194376/ccs 0 refseq 1 1 63M * 0 0
    ATGCAAAATGATGCATAGTATTAGCATAAATAGGATAGCCATAAGGTTACTGCATAAGAGTAT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NM:i:4 ms:i:111 AS:i:111 nn:i:3 tp:A:P cm:i:7 s1:i:41 s2:i:0 de:f:0.0167
    cs:Z::6*na*na*nt:49*ga:4 rl:i:0
    >>> print(tag_names)
    ['NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 's2', 'de', 'cs', 'rl']

    Now map retaining the `np` tag in the FASTQ:

    >>> samfile_tags = os.path.join(tempdir.name, 'alignments_tags.sam')
    >>> mapper_tags = Mapper(OPTIONS_CODON_DMS, retain_tags=['np'])
    >>> mapper_tags.map_to_sam(targetfile, queryfile, samfile_tags)
    >>> for a in pysam.AlignmentFile(samfile_tags):
    ...     tag_names = [tup[0] for tup in a.get_tags()]
    ...     print(a.tostring())  # doctest: +NORMALIZE_WHITESPACE
    m54228_181120_212724/4194376/ccs 0 refseq 1 1 63M * 0 0
    ATGCAAAATGATGCATAGTATTAGCATAAATAGGATAGCCATAAGGTTACTGCATAAGAGTAT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NM:i:4 ms:i:111 AS:i:111 nn:i:3 tp:A:P cm:i:7 s1:i:41 s2:i:0 de:f:0.0167
    cs:Z::6*na*na*nt:49*ga:4 rl:i:0 np:i:127
    >>> print(tag_names)  # doctest: +NORMALIZE_WHITESPACE
    ['NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 's2', 'de', 'cs', 'rl', 'np']

    Remove the temporary directory for the example:

    >>> tempdir.cleanup()

    """

    def __init__(self, options, *, prog='minimap2', min_version='2.17',
                 check_cs=True, retain_tags=None):
        """See main :class:`Mapper` doc string."""
        try:
            version = subprocess.check_output([prog, '--version'])
        except Exception:
            raise ValueError(f"Can't execute `prog` {prog}. Is it installed?")
        self.version = version.strip().decode('utf-8')
        min_version = packaging.version.parse(min_version)
        if packaging.version.parse(self.version) < min_version:
            raise ValueError(f"You have `minimap2` version {self.version}, "
                             f"but need >= {min_version}")
        self.prog = prog
        self.options = list(options)
        if check_cs and not any('--cs' == opt for opt in self.options):
            raise ValueError(f"`options` do not include `--cs`:\n{options}")
        if retain_tags:
            if isinstance(retain_tags, str):
                self.retain_tags = [retain_tags]
            else:
                self.retain_tags = list(retain_tags)
        else:
            self.retain_tags = None

    def map_to_sam(self, targetfile, queryfile, samfile):
        """Map query sequences to targets and output SAM file.

        Parameters
        ----------
        targetfile : str
            Alignment targets (reference), typically FASTA (can be gzipped).
        queryfile : str
            Queries to align, typically FASTA or FASTQ (can be gzipped).
        samfile : str
            Name of created SAM file (should have suffix ``.sam``).

        """
        for fname, f in [('target', targetfile), ('query', queryfile)]:
            if not os.path.isfile(f):
                raise IOError(f"cannot find `{fname}file` {f}")

        if os.path.splitext(samfile)[1] != '.sam':
            raise ValueError(f"`samfile` lacks extension '.sam': {samfile}")

        if not any('-a' == opt for opt in self.options):
            options = self.options + ['-a']

        cmds = [self.prog] + options + [targetfile, queryfile]

        with tempfile.TemporaryDirectory() as tempdir:
            # if retaining tags, we initially align to temporary file
            # then add tags to final alignment file
            if self.retain_tags:
                untagged_samfile = os.path.join(tempdir,
                                                os.path.basename(samfile))
            else:
                untagged_samfile = samfile

            # do the alignment
            with contextlib.ExitStack() as stack:
                sam = stack.enter_context(open(untagged_samfile, 'w'))
                err = stack.enter_context(tempfile.TemporaryFile(mode='w'))
                try:
                    _ = subprocess.check_call(cmds, stdout=sam, stderr=err)
                except Exception as exc:
                    if os.path.isfile(untagged_samfile):
                        os.remove(untagged_samfile)
                    exc_type = type(exc).__name__
                    err.seek(0)
                    raise RuntimeError(f"Error running commands:\n{cmds}\n\n"
                                       f"Exception:\n{exc_type}: {exc}\n\n"
                                       f"Error message:\n{err.read()}\n\n")

            # if needed, add retained tags from queries FASTQ file
            if self.retain_tags:
                matchers = {tag: re.compile(
                                    r'(?:^|\s)'  # start or space
                                    rf"{tag}:(?P<valtype>\w):(?P<val>\S+)"
                                    r'(?:\s|$)'  # end or space
                                    )
                            for tag in self.retain_tags}
                with contextlib.ExitStack() as tagstack:
                    untagged_sam = tagstack.enter_context(
                                    pysam.AlignmentFile(untagged_samfile))
                    tagged_sam = tagstack.enter_context(
                                    pysam.AlignmentFile(samfile, mode='w',
                                                        template=untagged_sam))
                    queries = tagstack.enter_context(
                                    pysam.FastxFile(queryfile))

                    query = next(queries)
                    for a in untagged_sam:
                        name = a.query_name
                        try:
                            while query.name != name:
                                query = next(queries)
                            if self.retain_tags and not query.comment:
                                raise ValueError('specified `retain_tags` but '
                                                 'no tags in `queryfile` '
                                                 f"{queryfile}")
                            for tag in self.retain_tags:
                                m = matchers[tag].search(query.comment)
                                if not m:
                                    raise ValueError(f"no tag {tag}:\n{query}")
                                valtype = m.group('valtype')
                                if valtype == 'i':
                                    val = int(m.group('val'))
                                elif valtype == 'f':
                                    val = float(m.group('val'))
                                elif valtype in {'A', 'Z'}:
                                    val = m.group('val')
                                else:
                                    raise ValueError(f"bad tag type {valtype} "
                                                     f"for tag {tag}")
                                a.set_tag(tag, val, valtype)
                        except StopIteration:
                            raise ValueError(f"No entry in {queryfile} for "
                                             f"query {name}\nMaybe alignment "
                                             f"is not sorted in query order?")
                        tagged_sam.write(a)

        assert os.path.isfile(samfile)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
