"""
=======
csparse
=======

Parses ``cs`` tag from ``minimap2`` output. 
"""

import re
import alignparse.constants


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
    AssertionError: No ambiguous nucleotide
    """
    assert '*n' == ambiguousDNA[:2], 'No ambiguous nucleotide'
    return ambiguousDNA[2]

def cs_list(cs, cs_splits):
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
    AssertionError: Tag :8*ng not supported
    """

    assert sum(tag.count(x) for x in alignparse.constants.CS_OPERATORS) <= 1, \
            f'Tag {tag} not supported'

    try:
        if re.match(r'^:', tag):
            return int(tag[1:])
        elif re.match(r'^\*', tag):
            return 1
        elif re.match(r'^\+', tag):
            return 0
        elif re.match(r'^\-', tag):
            return (len(tag)-1)
        elif re.match(r'^\=', tag):
            return(len(tag)-1)
        elif re.match(r'^\~', tag):
            raise RuntimeError('Handling splice sites not yet implemented.')
    except:
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
    split_tag : str
        ``cs`` tag for region specified by ref_indices 
    """

    #return None if feature coordinates not in alignment
    if featstart <= alignstart and featend <= alignstart:
        return None
    elif featstart >= alignend and featend >= alignend:
        return None
    #else:



if __name__ == '__main__':
    import doctest
    doctest.testmod()
