"""
=======
cs_tag
=======

Parse SAM entries with ``cs`` short tag from
`minimap2 <https://lh3.github.io/minimap2/>`_ output.

See `here <https://lh3.github.io/minimap2/minimap2.html>`_ for details on
short ``cs`` tag format.

"""


import functools

import numpy

import regex

import alignparse.utils

_CS_OPS = {
    'identity': ':[0-9]+',
    'substitution': r'\*[acgtn][acgtn]',
    'insertion': r'\+[acgtn]+',
    'deletion': r'\-[acgtn]+'
    }
"""dict: Short ``cs`` tag operation regular expression matches."""

_INTRON_OP = r'\~[acgtn]{2}\d+[acgtn]{2}'
"""str: Short ``cs`` tag operation regular expression for intron."""

_INTRON_OP_REGEX = regex.compile(_INTRON_OP)
"""regex.Regex: matches short ``cs`` tag operation for intron."""

_CS_OPS_W_INTRON = {key: val for key, val in
                    list(_CS_OPS.items()) + [('intron', _INTRON_OP)]}
"""dict: ``cs`` tag operation regular expression matches including introns."""

_CS_STR_REGEX = regex.compile('(' + '|'.join(list(_CS_OPS.values())) + ')*')
"""regex.Regex: matches full-length short ``cs`` tags."""

_CS_STR_REGEX_W_INTRON = regex.compile(
            '(' + '|'.join(list(_CS_OPS_W_INTRON.values())) + ')*')
"""regex.Regex: matches full-length short ``cs`` tags including introns."""

_CS_OP_REGEX = regex.compile('|'.join(f"(?P<{op_name}>{op_str})" for
                                      op_name, op_str in _CS_OPS.items()))
"""regex.Regex: matches single ``cs`` operation, group name is operation."""

_CS_OP_REGEX_W_INTRON = regex.compile(
            '|'.join(f"(?P<{op_name}>{op_str})" for
                     op_name, op_str in _CS_OPS_W_INTRON.items()))
"""regex.Regex: matches single ``cs`` operation including introns,
group name is operation."""


@functools.lru_cache(maxsize=16384)
def split_cs(cs_string, *, invalid='raise', allow_intron=False):
    """Split a short ``cs`` tag into its constituent operations.

    Parameters
    ----------
    cs_string : str
        The short ``cs`` tag.
    invalid : {'raise', 'ignore'}
        If `cs_string` is not a valid string, raise an error or ignore it
        and return `None`.
    allow_intron : bool
        Are introns allowed as ``cs`` operations?

    Return
    ------
    tuple or None
        Tuple of the individual ``cs`` operations, or `None` if invalid
        `cs_string` and `invalid` is 'ignore'.

    Example
    -------
    >>> split_cs(':32*nt*na:10-gga:5+aaa:10')
    (':32', '*nt', '*na', ':10', '-gga', ':5', '+aaa', ':10')

    >>> split_cs('bad:32*nt*na:10-gga:5', invalid='ignore') is None
    True

    >>> split_cs('bad:32*nt*na:10-gga:5')
    Traceback (most recent call last):
    ...
    ValueError: invalid `cs_string` of bad:32*nt*na:10-gga:5

    """
    if allow_intron:
        m = _CS_STR_REGEX_W_INTRON.fullmatch(cs_string)
    else:
        m = _CS_STR_REGEX.fullmatch(cs_string)
    if m is None:
        if invalid == 'ignore':
            return None
        elif invalid == 'raise':
            raise ValueError(f"invalid `cs_string` of {cs_string}")
        else:
            raise ValueError(f"invalid `invalid` of {invalid}")
    else:
        return tuple(m.captures(1))


@functools.lru_cache(maxsize=16384)
def cs_op_type(cs_op, *, invalid='raise'):
    """Get type of ``cs`` operation.

    Parameters
    ----------
    cs_op : str
        A **single** operation in a short ``cs`` tag.
    invalid : {'raise', 'ignore'}
        If `cs_string` is not a valid string, raise an error or ignore it
        and return `None`.

    Returns
    -------
    {'substitution', 'insertion', 'deletion', 'identity', None}
        Type of ``cs`` operation, or `None` if `cs_string` invalid and
        `invalid` is 'ignore'.

    Example
    -------
    >>> cs_op_type('*nt')
    'substitution'
    >>> cs_op_type(':45')
    'identity'
    >>> cs_op_type('*nt:45')
    Traceback (most recent call last):
    ...
    ValueError: invalid `cs_op` of *nt:45
    >>> cs_op_type('*nt:45', invalid='ignore') is None
    True

    """
    m = _CS_OP_REGEX.fullmatch(cs_op)
    if m is None:
        if invalid == 'ignore':
            return None
        elif invalid == 'raise':
            raise ValueError(f"invalid `cs_op` of {cs_op}")
        else:
            raise ValueError(f"invalid `invalid` of {invalid}")
    else:
        return m.lastgroup


@functools.lru_cache(maxsize=16384)
def cs_op_len_target(cs_op, *, invalid='raise'):
    """Get length of valid ``cs`` operation.

    Parameters
    ----------
    cs_op : str
        A **single** operation in a short ``cs`` tag.
    invalid : {'raise', 'ignore'}
        If `cs_string` is not a valid string, raise an error or ignore it
        and return `None`.

    Returns
    -------
    int or None
        Length of given cs_op. This length is based on the target sequence,
        so insertions in the query have length 0 and deletions are the length
        of the target sequence deleted from the query. 'None' if `cs_op`
        is invalid and `invalid` is `ignore`.

    Example
    -------
    >>> cs_op_len_target('*nt')
    1
    >>> cs_op_len_target(':45')
    45
    >>> cs_op_len_target('-at')
    2
    >>> cs_op_len_target('+gc')
    0
    >>> cs_op_len_target('*nt:45')
    Traceback (most recent call last):
    ...
    ValueError: invalid `cs_op` of *nt:45
    >>> cs_op_len_target('*nt:45', invalid='ignore') is None
    True

    """
    op_type = cs_op_type(cs_op, invalid=invalid)
    if op_type == 'identity':
        return int(cs_op[1:])
    elif op_type == 'substitution':
        return 1
    elif op_type == 'deletion':
        return len(cs_op) - 1
    elif op_type == 'insertion':
        return 0
    elif op_type is None:
        return None
    else:
        raise ValueError(f"invalid `op_type` of {op_type}")


@functools.lru_cache(maxsize=16384)
def cs_introns_to_deletions(cs, targetseq):
    """Convert introns to deletions in ``cs`` tag.

    Parameters
    ----------
    cs : str
        Short-format ``cs`` tag.
    targetseq : str
        Region of target sequenced covered by `cs` alignment.

    Returns
    -------
    str
        Version of `cs` where all introns have been converted to deletions.

    Examples
    --------
    >>> cs_introns_to_deletions(':3-ggaac:2', 'ATGGGAACAT')
    ':3-ggaac:2'
    >>> cs_introns_to_deletions(':3~gg5ac:2', 'ATGGGAACAT')
    ':3-ggaac:2'
    >>> cs_introns_to_deletions(':2*ga~gg6ta:2-gc:2~at4cg:3',
    ...                  'ATGGGCCTATTGCTAATCGAAA')
    ':2*ga-ggccta:2-gc:2-atcg:3'

    """
    if not _INTRON_OP_REGEX.search(cs):
        return cs
    itarget = 0
    new_cs = []
    for op in split_cs(cs, allow_intron=True):
        if _INTRON_OP_REGEX.fullmatch(op):
            op_len = int(op[3: -2])
            target_subseq = _ambiguous_to_n(
                        targetseq[itarget: itarget + op_len]).lower()
            itarget += op_len
            new_cs.append(f"-{target_subseq}")
            assert target_subseq[: 2] == op[1: 3], "{target_subseq}\n{op}"
            assert target_subseq[-2:] == op[-2:], "{target_subseq}\n{op}"
        else:
            new_cs.append(op)
            itarget += cs_op_len_target(op)
    return ''.join(new_cs)


def _ambiguous_to_n(seq):
    """Convert all ambiguous nucleotides to 'N'.

    Parameters
    ----------
    seq : str
        Sequence.

    Returns
    -------
    str
        Version of `seq` where all non-N IUPAC ambiguous nucleotides have
        been converted to 'N'.

    Example
    -------
    >>> _ambiguous_to_n('ATGYCAkac')
    'ATGNCANac'

    """
    return regex.sub('[MmRrWwSsYyKkVvHhDdBb]', 'N', seq)


class Alignment:
    """Process a SAM alignment with a ``cs`` tag to extract features.

    Parameters
    ----------
    sam_alignment : pysam.AlignedSegment
        Aligned segment from `pysam <https://pysam.readthedocs.io>`_,
        must have a short format ``cs`` tag (see
        https://lh3.github.io/minimap2/minimap2.html). This must be
        a mapped read, you will get an error if unmapped.
    introns_to_deletions : bool
        Convert all introns in the ``cs`` tag to deletions.
    target_seqs : dict
        Required if `introns_to_deletions` is `True`. Is keyed by target
        names with values target sequences as str.

    Attributes
    ----------
    query_name : str
        Name of query in alignment.
    target_name : str
        Name of alignment target.
    cs : str
        The ``cs`` tag.
    query_clip5 : int
        Length at 5' end of query not in alignment.
    query_clip3 : int
        Length at 3' end of query not in alignment.
    target_clip5 : int
        Length at 5' end of target not in alignment.
    target_lastpos : int
        Last position of alignment in target (exclusive).
    orientation :  {'+', '-'}
        Does query align to the target (+) or reverse complement (-).

    """

    def __init__(self, sam_alignment, *,
                 introns_to_deletions=False, target_seqs=None):
        """See main class docstring."""
        if sam_alignment.is_unmapped:
            raise ValueError(f"`sam_alignment` {sam_alignment.query_name} "
                             'is unmapped')

        self.query_name = sam_alignment.query_name
        self.target_name = sam_alignment.reference_name
        self.query_clip5 = sam_alignment.query_alignment_start
        self.query_clip3 = (sam_alignment.query_length -
                            sam_alignment.query_alignment_end)
        self.target_clip5 = sam_alignment.reference_start
        self.target_lastpos = sam_alignment.reference_end
        if sam_alignment.is_reverse:
            self.orientation = '-'
        else:
            self.orientation = '+'
        self.cs = str(sam_alignment.get_tag('cs'))

        if introns_to_deletions:
            if target_seqs is None:
                raise ValueError('must set `target_seqs`')
            targetseq = target_seqs[self.target_name][self.target_clip5:
                                                      self.target_lastpos]
            self.cs = cs_introns_to_deletions(self.cs, targetseq)

        self._cs_ops = split_cs(self.cs)
        self._cs_ops_lengths_target = numpy.array([cs_op_len_target(op)
                                                   for op in self._cs_ops])

        # sites are 0-indexed and exclusive
        self._cs_ops_ends = (self.target_clip5 +
                             numpy.cumsum(self._cs_ops_lengths_target))
        self._cs_ops_starts = numpy.append(numpy.array(self.target_clip5),
                                           self._cs_ops_ends[:-1])

        self._nops = len(self._cs_ops)
        assert self._nops == len(self._cs_ops_lengths_target)
        assert self._nops == len(self._cs_ops_ends)
        assert self._nops == len(self._cs_ops_starts)
        assert (self._cs_ops_ends - self._cs_ops_starts ==
                self._cs_ops_lengths_target).all()

        self._sam_alignment = sam_alignment

    def get_accuracy(self, targetstart, targetend):
        """Get accuracy of part of query aligned to a region of the target.

        Parameters
        ----------
        targetstart : int
            Start of region in target in 0, 1, ... numbering.
        targetend : int
            End of region in target (not inclusive of this site).

        Returns
        -------
        float or NaN
            The accuracy of the region, computed as the average of the
            Q-scores for all aligned sites in the query, or `NaN` if
            there are no query sites aligned to this region.

        """
        if targetstart >= targetend:
            raise ValueError('`targetstart` must be < `targetend`')

        # if function hasn't yet been called, set up necessary variables
        if not hasattr(self, '_qs'):
            # array of all Q values in query
            self._qs = numpy.asarray(self._sam_alignment.query_qualities,
                                     dtype='int')
            # arrays of query and target sites that are aligned
            self._aligned_query, self._aligned_target = map(
                    numpy.array,
                    zip(*self._sam_alignment
                        .get_aligned_pairs(matches_only=True))
                    )
            assert len(self._aligned_query) == len(self._aligned_target)

        # get index of first aligned site >= targetstart, and last
        # aligned site <= targetend
        istart = numpy.searchsorted(self._aligned_target, targetstart)
        iend = numpy.searchsorted(self._aligned_target, targetend)
        n = len(self._aligned_query) - 1
        istart = min(istart, n)
        iend = min(iend, n)
        assert istart >= 0
        assert istart <= iend

        # compute accuracy
        return alignparse.utils.qvals_to_accuracy(
                self._qs[self._aligned_query[istart]:
                         self._aligned_query[iend]])

    def extract_cs(self, start, end):
        """Extract ``cs`` tag corresponding to feature in target.

        Parameters
        ----------
        start : int
            Start of feature in target in 0, 1, ... numbering.
        end : int
            End of feature in target (not inclusive of this site).

        Returns
        -------
        tuple or None
            If feature is not present in alignment, return `None`. Otherwise,
            return `(cs, clip5, clip3)` where `cs` is the ``cs`` string for
            the aligned portion of the feature, and `clip5` and `clip3` are
            the lengths that must be clipped off the end of the feature to
            get to that alignment.

        Note
        ----
        If an insertion is at the boundary of two features, it is assigned
        as being at the end of the first feature.

        """
        if start < 0:
            raise ValueError(f"invalid `start` of {start}")
        if end <= start:
            raise ValueError(f"`end` {end} not > `start` {start}")
        if start >= numpy.amax(self._cs_ops_ends):
            # feature starts at or after end of last cs op
            return None
        if end <= numpy.amin(self._cs_ops_starts):
            # feature ends at or before beginning of first cs op
            return None

        clip5 = clip3 = 0

        # Get `start_idx` as index of cs op that contains feature start
        # add to `feature_cs` overlapping part of first cs op
        start_idx = numpy.searchsorted(self._cs_ops_ends, start, side='right')
        start_op_start = self._cs_ops_starts[start_idx]
        start_op = self._cs_ops[start_idx]
        assert start_idx < self._nops
        feature_cs = []
        if start_op_start > start:
            # feature starts before first cs op
            assert start_idx == 0, "5' clip not at 5' end."
            clip5 = start_op_start - start
            feature_cs.append(start_op)
        else:
            # feature starts at or within specific cs op
            start_op_end = self._cs_ops_ends[start_idx]
            assert start < start_op_end
            start_op_type = cs_op_type(start_op)
            if start_op_start == start and end >= start_op_end:
                feature_cs.append(start_op)
            elif start_op_type == 'identity':
                feature_cs.append(f":{start_op_end - start}")
            elif start_op_type == 'deletion':
                feature_cs.append(f"-{start_op[start - start_op_end:]}")
            elif start_op_type == 'insertion':
                raise RuntimeError('insertion should not be feature start')
            else:
                raise RuntimeError(f"unrecognized op type of {start_op_type}")

        # Get `end_idx` as index of cs op that contains feature end, and
        # make `feat_cs_end` the overlapping part of this last cs op
        end_idx = max(0, numpy.searchsorted(self._cs_ops_ends, end,
                                            side='right') - 1)
        while ((end > self._cs_ops_ends[end_idx]) and
               (end_idx + 1 < self._nops)):
            end_idx += 1
        end_op_start = self._cs_ops_starts[end_idx]
        end_op_end = self._cs_ops_ends[end_idx]
        end_op = self._cs_ops[end_idx]
        assert start_idx <= end_idx <= self._nops
        assert end <= end_op_end or end_idx == self._nops - 1
        assert end >= end_op_start
        if end > end_op_end:
            assert end_idx == self._nops - 1, 'clip3 not at end'
            clip3 = end - end_op_end
            feat_cs_end = end_op
        elif end == end_op_end:
            # feature ends at a specific cs op
            feat_cs_end = end_op
        else:
            # feature ends within specific cs op
            end_overlap = end - end_op_start
            end_op_type = cs_op_type(end_op)
            if end_op_type == 'identity':
                feat_cs_end = f":{end_overlap}"
            elif end_op_type == 'deletion':
                feat_cs_end = end_op[: end_overlap + 1]
            elif end_op_type == 'insertion':
                raise RuntimeError('should not get here as end == end_op_end')
            else:
                raise RuntimeError(f"unrecognized op type of {end_op_type}")

        if start_idx == end_idx:
            # avoid double-counting feature, clip properly
            assert start_op == end_op
            op_type = cs_op_type(start_op)
            if op_type == 'identity':
                feature_cs = f":{end - start - clip5 - clip3}"
            elif op_type == 'substitution':
                feature_cs = end_op
            elif op_type == 'deletion':
                del_start = max(0, start - start_op_start)
                del_end = end_op_end - end_op_start - max(0, end_op_end - end)
                feature_cs = '-' + end_op[del_start + 1: del_end + 1]
            elif op_type == 'insertion':
                raise RuntimeError('start_idx != end_idx for insertion')
            else:
                raise RuntimeError(f"unrecognized op type of {op_type}")
        else:
            feature_cs.extend(self._cs_ops[start_idx + 1: end_idx])
            feature_cs.append(feat_cs_end)
            feature_cs = ''.join(feature_cs)

        # this next assert might be costly, so maybe remove eventually
        assert (sum(cs_op_len_target(op) for op in split_cs(feature_cs)) +
                clip5 + clip3 == end - start), f"{feature_cs},{clip5},{clip3}"

        return (feature_cs, clip5, clip3)


@functools.lru_cache(maxsize=16384)
def cs_to_sequence(cs, seq):
    """Convert ``cs`` tag to a sequence.

    Parameters
    ----------
    cs : str
        `cs` string
    seq : str
        Sequence of target for region corresponding to `cs` string.

    Returns
    -------
    str
        Nucleotide sequence generated by applying `cs` to `seq`.

    Example
    -------
    >>> cs_to_sequence(':4*nt-tc:2+g:2', 'CGGANTCCAAT')
    'CGGATCAGAT'

    """
    seq_loc = 0
    seq_list = []
    for cs_op in split_cs(cs):
        op_type = cs_op_type(cs_op)
        if op_type == 'identity':
            op_len = cs_op_len_target(cs_op)
            seq_list.append(seq[seq_loc: seq_loc + op_len])
            seq_loc += op_len
        elif op_type == 'substitution':
            seq_list.append(cs_op[2])
            seq_loc += 1
        elif op_type == 'insertion':
            seq_list.append(cs_op[1:])
        elif op_type == 'deletion':
            seq_loc += len(cs_op) - 1
        else:
            raise ValueError(f"Invalid cs `op_type` of {op_type}")

    return ''.join(seq_list).upper()


@functools.lru_cache(maxsize=16384)
def cs_to_mutation_str(cs, offset=0):
    """Convert ``cs`` tag to a descriptive string of mutations.

    Parameters
    ----------
    cs : str
        A ``cs`` tag.
    offset : int
        Mutations are numbered in 1, ... numbering **plus** this offset.

    Returns
    -------
    str
        Space-delimited string of form 'A5T G86A ins7ACG del19to24'
        for all mutations specified in `cs`.

    Example
    -------
    >>> cs_to_mutation_str(':4*nt-tc:2+ga:6')
    'del6to7 ins10GA'
    >>> cs_to_mutation_str(':4*at-tc:2+ga:6')
    'A5T del6to7 ins10GA'
    >>> cs_to_mutation_str(':4*at-tc:2+ga:6', offset=2)
    'A7T del8to9 ins12GA'
    >>> cs_to_mutation_str(':45')
    ''

    Note
    ----
    Mutation strings use "human readable" indexing, so the first nucleotide of
    the sequence is 1 and deletions are inclusive of the last number.

    Changes from ambiguous nucleotides to any other identity are **not**
    considered mutations in the returned strings.

    """
    seq_loc = 1 + offset
    mut_strs_list = []
    for cs_op in split_cs(cs):
        op_type = cs_op_type(cs_op)
        if op_type == 'identity':
            seq_loc += cs_op_len_target(cs_op)
        elif op_type == 'substitution':
            if cs_op[1] != 'n':
                sub = ''.join([cs_op[1], str(seq_loc), cs_op[2]]).upper()
                mut_strs_list.append(sub)
            seq_loc += 1
        elif op_type == 'insertion':
            ins = ''.join(['ins', str(seq_loc), cs_op[1:].upper()])
            mut_strs_list.append(ins)
        elif op_type == 'deletion':
            deletion = ''.join(['del', str(seq_loc), 'to',
                                str(seq_loc+len(cs_op)-2)])
            mut_strs_list.append(deletion)
            seq_loc += len(cs_op) - 1
        else:
            raise ValueError(f"Invalid cs `op_type` of {op_type}")

    return ' '.join(mut_strs_list)


@functools.lru_cache(maxsize=16384)
def cs_to_nt_mutation_count(cs):
    """Count the number of nucleotide mutations in ``cs`` tag.

    Parameters
    ----------
    cs : str
        `cs` string

    Returns
    -------
    int
        Number of nucleotides that are mutated. Insertions / deletions
        are counted as the number of nucleotides in the indel.
        Changes from an ambiguous nucleotide to are **not** considered
        mutations.

    Example
    -------
    >>> cs_to_nt_mutation_count(':4*nt-tc:2+g')
    3
    >>> cs_to_nt_mutation_count(':4*gt-tc:2+g')
    4

    """
    nt_mut_count = 0
    for cs_op in split_cs(cs):
        op_type = cs_op_type(cs_op)
        if op_type == 'substitution':
            if cs_op[1] != 'n':
                nt_mut_count += 1
        elif op_type == 'insertion' or op_type == 'deletion':
            nt_mut_count += len(cs_op) - 1
        elif op_type != 'identity':
            raise ValueError(f'Invalid cs `op_type` of {op_type}.')

    return nt_mut_count


@functools.lru_cache(maxsize=16384)
def cs_to_op_mutation_count(cs):
    """Count the number of mutation operations in ``cs`` tag.

    Parameters
    ----------
    cs : str
        The ``cs`` tag.

    Returns
    -------
    int
        Number of mutation operations in the query sequence. Each indel or
        substitution is counted as a single mutation operation regardless
        of how many mutations it contains. Changes from ambiguous nucleotides
        to another nucleotide are **not** counted.

    Example
    -------
    >>> cs_to_op_mutation_count(':4*nt-tc:2+g')
    2
    >>> cs_to_op_mutation_count(':4*gt-tc:2+g')
    3

    """
    op_mut_count = 0
    for cs_op in split_cs(cs):
        op_type = cs_op_type(cs_op)
        if op_type == 'substitution':
            if cs_op[1] != 'n':
                op_mut_count += 1
        elif op_type == 'insertion' or op_type == 'deletion':
            op_mut_count += 1
        elif op_type != 'identity':
            raise ValueError(f'Invalid cs `op_type` of {op_type}.')

    return op_mut_count


if __name__ == '__main__':
    import doctest
    doctest.testmod()
