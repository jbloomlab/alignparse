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
    'clip': '<clip[0-9]+>'
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

    >>> split_cs('<clip8>:32*nt*na:10-gga:5')
    ['<clip8>', ':32', '*nt', '*na', ':10', '-gga', ':5']

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
    >>> cs_op_type('<clip8>')
    'clip'
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
    >>> cs_op_len_target('<clip8>')
    8
    >>> cs_op_len_target('<clip46>')
    46
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
    elif op_type == 'clip':
        return int(regex.findall(r'\d+', cs_op)[0])
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

        self._nops = len(self._cs_ops)
        assert self._nops == len(self._cs_ops_lengths_target)
        assert self._nops == len(self._cs_ops_ends)
        assert self._nops == len(self._cs_ops_starts)
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
        start_op_end = self._cs_ops_ends[start_idx]
        start_op = self._cs_ops[start_idx]
        assert start_idx < self._nops
        assert start < start_op_end
        feature_cs = []
        if start_op_start > start:
            # feature starts before first cs op
            assert start_idx == 0, "5' clip not at 5' end."
            clip5 = start_op_start - start
            feature_cs.append(start_op)
        else:
            # feature starts at or within specific cs op
            start_overlap = start_op_end - start
            assert start_overlap >= 0
            start_op_type = cs_op_type(start_op)
            if start_op_start == start and end >= start_op_end:
                feature_cs.append(start_op)
            elif start_op_type == 'identity':
                feature_cs.append(f":{start_overlap}")
            elif start_op_type == 'deletion':
                feature_cs.append(f"-{start_op[-start_overlap:]}")
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
            assert end_overlap > 0
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

def cs_to_sequence(cs, feat_name, feat_seq, *, custom_cs=False):
    """Convert `cs` string for a feature in `feat_seq` to nt sequence.

    Paramters
    ---------
    cs : str
        `cs` string for feature
    feat_name : str
        Name of feature for which the cs string is being converted to seq.
    custom_cs : bool
        If `True`, will process custom `cs` strings that include clip amounts
        on the 5' and/or 3' ends. Default is False to require users to
        to consider if they want to process `cs` strings with clipping.

    Returns
    -------
    sequence : str
        Nucleotide sequence for specified feature in the query.
    """
    raise RuntimeError('not yet implemented')

def cs_to_mutation_str(cs, feat_name, *, custom_cs=False):
    """Convert `cs` string for a feature in `feat_seq` to mutation string.
    
    Paramters
    ---------
    cs : str
        `cs` string for feature
    feat_name : str
        Name of feature for which the cs string is being converted to a
        mutation string.
    custom_cs : bool
        If `True`, will process custom `cs` strings that include clip amounts
        on the 5' and/or 3' ends. Default is False to require users to
        to consider if they want to process `cs` strings with clipping.

    Returns
    -------
    mut_str : str
        Mutation string of form 'A56T G86A' for all mutations in the feature
        in the query compared to the target sequence.
    """
    raise RuntimeError('not yet implemented')

def cs_to_mutation_count(cs):
    """Count the number of nucleotide mutations in `cs` string.
    
    Paramters
    ---------
    cs : str
        `cs` string for feature

    There is no `custom_cs` parameter becuase this function does ignores
    clipping.

    Returns
    -------
    mut_count : int
        Number of nucleotides that are mutated in the query sequence. All
        substituted, inserted, or deleted nucleotides are counted. Clipped
        nucelotides are not included.

    Example
    -------
    >>> cs_to_mutation_count(':4*nt-tc:2+g')
    4

    >>> cs_to_mutation_count('<clip4>:4*nt-tc:2+g')
    4
    """
    mut_count = 0
    cs_list = split_cs(cs)
    for cs_op in cs_list:
        op_type = cs_op_type(cs_op)
        if op_type == 'substitution':
            mut_count += 1
        elif op_type == 'insertion' or op_type == 'deletion':
            mut_count += len(cs_op) - 1
        else:
            mut_count += 0

    return mut_count


def cs_to_clip_count(cs, *, custom_cs=True):
    """Count the number of clipped nucleotides in `cs` string.
    
    Paramters
    ---------
    cs : str
        `cs` string for feature
    custom_cs : bool
        If `True`, will process custom `cs` strings that include clip amounts
        on the 5' and/or 3' ends. Default is `True` as this function uses the
        information specified in the custom `cs` strings that specify 5' and/
        or 3' clipping.

    Returns
    -------
    clip_count : int
        Number of nucleotides that are clipped from either or both ends of the
        feature in the query sequence.
    """
    raise RuntimeError('not yet implemented')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
