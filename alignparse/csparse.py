"""
=======
csparse
=======

Parses ``cs`` tag from ``minimap2`` output.
"""

import re

CS_OPERATORS = (':', '*', '+', '-', '~', '=')

"""tuple: characters used as operators in the ``cs`` tag.

From https://lh3.github.io/minimap2/minimap2.html#10

"""


def get_ambiguous(ambiguousDNA):
    """
    Given a ``cs`` tag of ambiguous reference nucleotides aligned to
    the query sequence, return the query sequence.

    Parameters
    ----------
    ambiguousDNA : str
        ``cs`` tag containing a single amibiguous nucleotide

    Returns
    -------
    queryseq : str
        query sequence for ambiguous nucleotide

    >>> get_ambiguous('*ng')
    'g'

    >>> get_ambiguous('ctat')
    Traceback (most recent call last):
    ...
    ValueError: No ambiguous nucleotide.

    """
    if '*n' == ambiguousDNA[:2]:
        return ambiguousDNA[2]
    else:
        raise ValueError('No ambiguous nucleotide.')


def cs_list(cs, *, cs_splits=CS_OPERATORS):
    """
    Turn a ``cs`` tag string into a list of component tags.

    Parameters
    ----------
    cs : str
        full ``cs`` tag string
    cs_splits : list
        list of operators by which to split the ``cs`` tag into
        component tags

    Returns
    -------
    cs_list : list
        ``cs`` tag in list format with each item corresponding to
        one tag fragment (or chunk of DNA)

    """
    cs_list = []
    split = ''
    i = 0
    while i < len(cs):
        if cs[i] in cs_splits:
            split = ''
            split += cs[i]
            i += 1
        while i < len(cs) and cs[i] not in cs_splits:
            split += cs[i]
            i += 1
        cs_list.append(split)
    return cs_list


def convert_len(tag):
    """
    Convert length of portion of ``cs`` tag into corresponding ref length.

    Parameters
    ----------
    tag : str
        single ``cs`` tag fragment

    Returns
    -------
    int
        corresponding length in reference sequence

    >>> convert_len(':8')
    8

    >>> convert_len('*tg')
    1

    >>> convert_len('-gtc')
    3

    >>> convert_len('+agct')
    0

    >>> convert_len(':8*ng')
    Traceback (most recent call last):
    ...
    ValueError: Tag :8*ng not supported

    """
    if sum(tag.count(x) for x in CS_OPERATORS) > 1:
        raise ValueError(f'Tag {tag} not supported')

    if re.match(r'^:', tag):
        return int(tag[1:])
    elif re.match(r'^\*', tag):
        return 1
    elif re.match(r'^\+', tag):
        return 0
    elif re.match(r'^\-', tag):
        if re.match(r'^\-\d+', tag):
            return int(re.findall(r'\d+', tag)[0])
        else:
            return (len(tag)-1)
    elif re.match(r'^\=', tag):
        return(len(tag)-1)
    elif re.match(r'^\~', tag):
        raise RuntimeError('Handling splice sites not yet implemented.')
    else:
        raise ValueError(f'Tag {tag} not supported.')


def split(cs, alignstart, alignend, featstart, featend):
    """
    Give ``cs`` tag for region of ref specified by ref_indices.

    Parameters
    ----------
    cs : str
        ``cs`` tag for read alignment
    featstart : int
        start index for region to extract
    featend: int
        end index for region to extract

    Returns
    -------
    feat_cs : str
        ``cs`` tag for region specified by featstart and featend

    """
    # return None if feature coordinates not in alignment
    if featend <= alignstart:
        return None
    elif featstart >= alignend:
        return None

    else:
        feat_cs = ''
        ref_loc = alignstart
        cslist = cs_list(cs)
        csindex = 0

        # if feature starts at beginning of alignment
        if featstart == ref_loc:
            ref_loc += convert_len(cslist[csindex])
            while ref_loc < featend:
                feat_cs += cslist[csindex]
                csindex += 1
                ref_loc += convert_len(cslist[csindex])
            if ref_loc == featend:
                feat_cs += cslist[csindex]
            elif cslist[csindex][0] == ':':
                feat_cs = feat_cs + ':' + str(featend - (ref_loc -
                                                convert_len(cslist[csindex])))
            elif cslist[csindex][0] == '+' or cslist[csindex][0] == '-':
                raise RuntimeError('Breaks at indels not yet implemented.')

        # if feature starts before alignment, but overlaps with alignment
        elif featstart < ref_loc and featend >= ref_loc:
            feat_cs = '-{0}'.format(ref_loc - featstart)
            ref_loc += convert_len(cslist[csindex])
            while ref_loc < featend:
                feat_cs += cslist[csindex]
                csindex += 1
                ref_loc += convert_len(cslist[csindex])
            if ref_loc == featend:
                feat_cs += cslist[csindex]
            elif cslist[csindex][0] == ':':
                feat_cs = feat_cs + ':' + str(featend - (ref_loc -
                                                convert_len(cslist[csindex])))

        # if feature starts after alignment starts
        elif featstart > ref_loc:
            ref_loc += convert_len(cslist[csindex])
            while ref_loc < featstart:
                csindex += 1
                ref_loc += convert_len(cslist[csindex])
            # if feature starts between cs tags
            if featstart == ref_loc:
                csindex += 1
                ref_loc += convert_len(cslist[csindex])
                while ref_loc < featend:
                    feat_cs += cslist[csindex]
                    csindex += 1
                    ref_loc += convert_len(cslist[csindex])
                if ref_loc == featend:
                    feat_cs += cslist[csindex]
                if cslist[csindex][0] == ':':
                    feat_cs = feat_cs + ':' + str(featend - (ref_loc -
                                                convert_len(cslist[csindex])))
            # if feature starts in middle of match, get length of match first
            # THIS IS WRONG
            elif ref_loc > featstart:
                csindex += 1
                if cslist[csindex-1][0] == ':':
                    feat_cs += feat_cs + ':' + str(ref_loc - featstart)
                while ref_loc <= featend:
                    feat_cs += cslist[csindex]
                    csindex += 1
                    ref_loc += convert_len(cslist[csindex])
                if cslist[csindex][0] == ':':
                    feat_cs = feat_cs + ':' + str(featend - (ref_loc -
                                                convert_len(cslist[csindex])))

            # Need to implement if feat start and end are in last tag.

    feat_list = cs_list(feat_cs)
    feat_len = 0
    for feat in feat_list:
        feat_len += convert_len(feat)

    if feat_len != (featend - featstart):
        raise ValueError('`feat_cs` incorrect length.')

    return feat_cs


if __name__ == '__main__':
    import doctest
    doctest.testmod()
