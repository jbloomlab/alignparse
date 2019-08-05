"""
=======
cs_tag
=======

Parse SAM entries with ``cs`` short tag from ``minimap2`` output.

See here for details on short ``cs`` tag format:
https://lh3.github.io/minimap2/minimap2.html

"""

import numpy

import regex

_CS_OPS = {
    'identity': ':[0-9]+',
    'substitution': r'\*[acgtn][acgtn]',
    'insertion': r'\+[acgtn]+',
    'deletion': r'\-[acgtn]+',
    }
"""dict: Short ``cs`` tag operation regular expression matches."""

_CS_STR_REGEX = regex.compile('(' + '|'.join(list(_CS_OPS.values())) + ')+')
"""regex.Regex: matches full-length short ``cs`` tags."""

_CS_OP_REGEX = regex.compile('|'.join(f"(?P<{op_name}>{op_str})" for
                                      op_name, op_str in _CS_OPS.items()))
"""regex.Regex: matches single ``cs`` operation, group name is operation."""


def split_cs(cs_string, *, invalid='raise'):
    """Split a short ``cs`` tag into its constituent operations.

    Parameters
    ----------
    cs_string : str
        The short ``cs`` tag.
    invalid : {'raise', 'ignore'}
        If `cs_string` is not a valid string, raise an error or ignore it
        and return `None`.

    Return
    ------
    list or None
        List of the individual ``cs`` operations, or `None` if invalid
        `cs_string` and `invalid` is 'ignore'.

    Example
    -------
    >>> split_cs(':32*nt*na:10-gga:5+aaa:10')
    [':32', '*nt', '*na', ':10', '-gga', ':5', '+aaa', ':10']

    >>> split_cs('bad:32*nt*na:10-gga:5', invalid='ignore') is None
    True

    >>> split_cs('bad:32*nt*na:10-gga:5')
    Traceback (most recent call last):
    ...
    ValueError: invalid `cs_string` of bad:32*nt*na:10-gga:5

    """
    m = _CS_STR_REGEX.fullmatch(cs_string)
    if m is None:
        if invalid == 'ignore':
            return None
        elif invalid == 'raise':
            raise ValueError(f"invalid `cs_string` of {cs_string}")
        else:
            raise ValueError(f"invalid `invalid` of {invalid}")
    else:
        return m.captures(1)


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


class Alignment:
    """Process a SAM alignment with a ``cs`` tag to extract features.

    Parameters
    ----------
    sam_alignment : pysam.AlignedSegment
        Aligned segment from `pysam <https://pysam.readthedocs.io>`_,
        must have a short format ``cs`` tag (see
        https://lh3.github.io/minimap2/minimap2.html). This must be
        a mapped read, you will get an error if unmapped.

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

    def __init__(self, sam_alignment):
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

        self._cs_ops = split_cs(self.cs)
        self._cs_ops_lengths_target = numpy.array([cs_op_len_target(op)
                                                   for op in self._cs_ops])

        # sites are 0-indexed and exclusive
        self._cs_ops_ends = (self.target_clip5 +
                             numpy.cumsum(self._cs_ops_lengths_target))
        self._cs_ops_starts = numpy.append(numpy.array(self.target_clip5),
                                           self._cs_ops_ends[:-1])

        nops = len(self._cs_ops)
        assert nops == len(self._cs_ops_lengths_target)
        assert nops == len(self._cs_ops_ends)
        assert nops == len(self._cs_ops_starts)
        assert (self._cs_ops_ends - self._cs_ops_starts ==
                self._cs_ops_lengths_target).all()

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

        """
        if start >= numpy.amax(self._cs_ops_ends):
            # feature starts at or after end of last cs op
            return None
        
        if end <= numpy.amin(self._cs_ops_starts):
            # feature ends at or before beginning of first cs op
            return None

        clip5 = 0
        clip3 = 0

        # Get `start_idx` as index of cs op that contains feature start, and
        # add to `feature_cs` any partial first cs op.
        start_idx = numpy.searchsorted(self._cs_ops_ends, start)
        feature_cs = []
        start_overlap = self._cs_ops_ends[start_idx] - start
        if self._cs_ops_starts[start_idx] > start:
            assert start_idx == 0, "5' clip not at 5' end."
            clip5 = self._cs_ops_starts[start_idx] - start
            feature_cs.append(self._cs_ops[start_idx])
        elif self._cs_ops_starts[start_idx] == start:
            feature_cs.append(self._cs_ops[start_idx])
        elif cs_op_type(self._cs_ops[start_idx]) == 'identity':
            feature_cs.append(f":{start_overlap}")
        elif cs_op_type(self._cs_ops[start_idx]) == 'deletion':
            feature_cs.append(f"-{self._cs_ops[start_idx][-start_overlap:]}")
        elif cs_op_type(self._cs_ops[start_idx]) == 'insertion':
            raise RuntimeError('Splitting insertions not implemented.')
        else:
            raise RuntimeError('Getting `cs` start failed.')

        # Get `end_idx` as index of cs op that contains feature end, and
        # add set as `feat_cs_end`
        end_idx = numpy.searchsorted(self._cs_ops_ends, end)
        feat_cs_end = []
        end_overlap = end - self._cs_ops_ends[end_idx-1]
        if end > numpy.amax(self._cs_ops_ends):
            assert end_idx == len(self._cs_ops_ends), "3' clip, not at 3' end"
            clip3 = end - numpy.amax(self._cs_ops_ends)
        elif self._cs_ops_ends[end_idx] == end:
            feat_cs_end.append(self._cs_ops[end_idx])
        elif cs_op_type(self._cs_ops[end_idx]) == 'identity':
            feat_cs_end.append(f":{end_overlap}")
        elif cs_op_type(self._cs_ops[end_idx]) == 'deletion':
            feat_cs_end.append(f"-{self._cs_ops[end_idx][:end_overlap]}")
        elif cs_op_type(self._cs_ops[end_idx]) == 'insertion':
            raise RuntimeError('Splitting insertions not implemented.')
        else:
            raise RuntimeError('Getting `cs` end failed.')

        if start_idx == end_idx:
            if cs_op_type(self._cs_ops[start_idx]) == 'identity':
                feature_cs = [f":{end-start}"]
            else:  # entire feature deleted in query
                return None
        elif end_idx - start_idx == 1:
            feature_cs.extend(feat_cs_end)
        else:
            feature_cs.extend(self._cs_ops[start_idx+1: end_idx])
            feature_cs.extend(feat_cs_end)

        return (''.join(feature_cs), clip5, clip3)

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
