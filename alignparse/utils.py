"""
=====
utils
=====

"""


import math
import numbers
import re

import numpy

import pandas as pd  # noqa: F401

import alignparse.consensus


def qvals_to_accuracy(qvals, encoding='numbers'):
    r"""Convert set of quality scores into average accuracy.

    Parameters
    ----------
    qvals : numpy.array, number, or str
        Q-values, for how they are encoded see `encoding`.
    encoding : {'numbers', 'sanger'}
        If 'numbers' then `qvals` should be a numpy.array of Q-values
        or a number giving a single Q-value. If 'sanger', then `qvals`
        is a string, with the Q-value being the ASCII value minus 33.

    Returns
    -------
    float or nan
        The average accuracy if the Q-values.
        `nan` if `qvals` is empty.

    Note
    ----
    The probability :math:`p` of an error at a given site is related to
    the Q-value :math:`Q` by :math:`Q = -10 \log_{10} p`. The accuracy
    is one minus the average error rate.

    Example
    -------

    >>> qvals = numpy.array([13, 77, 93])
    >>> round(qvals_to_accuracy(qvals), 3)
    0.983
    >>> round(qvals_to_accuracy(qvals[1 : ]), 3)
    1.0
    >>> qvals_to_accuracy(numpy.array([]))
    nan

    >>> qvals_str = '.n~'
    >>> round(qvals_to_accuracy(qvals_str, encoding='sanger'), 3)
    0.983

    >>> round(qvals_to_accuracy(15), 3)
    0.968

    """
    if encoding == 'numbers':
        if isinstance(qvals, numbers.Number):
            qvals = numpy.array([qvals])
        elif isinstance(qvals, list):
            qvals = numpy.array(qvals)

    if len(qvals) == 0:
        return math.nan

    if encoding == 'numbers':
        pass
    if encoding == 'sanger':
        qvals = numpy.array([ord(q) - 33 for q in qvals])
    elif encoding != 'numbers':
        raise ValueError(f"invalid `encoding`: {encoding}")

    return (1 - 10**(qvals / -10)).sum() / len(qvals)


def sort_mutations(mut_strs):
    """Sort mutation string by site, and combine multiple mutation strings.

    Parameters
    ----------
    mut_strs : str or list
        A single mutation string or a list of such strings.

    Returns
    -------
    str
        A single mutation string with all mutations sorted by site.

    Example
    -------
    Sort a single mutation string:

    >>> sort_mutations('ins7GC A5C del2to3')
    'del2to3 A5C ins7GC'

    Sort a list of two mutation strings, including a negative site:

    >>> sort_mutations(['ins7GC', 'A-5C del2to3'])
    'A-5C del2to3 ins7GC'

    """
    if isinstance(mut_strs, str):
        mut_strs = [mut_strs]
    decorated_list = []
    for mut_str in mut_strs:
        for mut in mut_str.split():
            m = re.fullmatch(r'ins(\-?\d+)[A-Z]+|'
                             r'[A-Z](\-?\d+)[A-Z]|'
                             r'del(\-?\d+)to\-?\d+',
                             mut)
            if not m:
                raise ValueError(f"failed to match {mut} in:\n{mut_str}")
            site = [site for site in m.groups() if site is not None]
            assert len(site) == 1
            site = int(site[0])
            decorated_list.append((site, mut))
    return ' '.join(mut for _, mut in sorted(decorated_list))


def merge_dels(s):
    """Merge consecutive deletions

    Parameters
    ----------
    s : str
        A single string of mutations.

    Returns
    -------
    str
        A mutation strings where consecutive deletions have been merged,
        and all mutations are sorted by site.

    Example
    -------
    Merge consecutive deletions:

    >>> merge_dels('del12to15 del21to30 del210to300 del16to20 '
    ...            'del1702to1909 del1910to1930 G885T G85T')
    'del12to30 G85T del210to300 G885T del1702to1930'

    """
    # parse deletions
    subs, deletions, insertions = alignparse.consensus.process_mut_str(s)

    # if no deletions are available return current mutation
    if not deletions:
        return sort_mutations(s)

    else:
        # extract position and from:to deletion list
        mut_list = []
        for deletion in deletions:
            mut_sites = re.fullmatch(r'del(?P<start>\-?\d+)to(?P<end>\-?\d+)',
                                     deletion)
            if not mut_sites:
                raise ValueError(f"cannot match deletion {deletion}")
            mut_list.append([int(mut_sites.group('start')),
                             int(mut_sites.group('end'))])

        # merge consecutive ranges
        mut_list.sort()
        new_ranges = []
        left, right = mut_list[0]
        for del_range in mut_list[1:]:
            next_left, next_right = del_range
            if right + 1 < next_left:
                new_ranges.append([left, right])
                left, right = next_left, next_right
            else:
                right = max(right, next_right)
        new_ranges.append([left, right])

        # create range-corrected (RC) deletion list
        deletion_RC = [f"del{left}to{right}" for left, right in new_ranges]

        # overwrite mutation tuple with new range
        return sort_mutations([*subs, *deletion_RC, *insertions])


class MutationRenumber:
    """Re-number mutations.

    Arguments
    ----------
    number_mapping : pandas.DataFrame
        Data frame giving mapping from old to new numbering scheme.
    old_num_col : str
        Column in `number_mapping` giving old site number.
    new_num_col : str
        Column in `number_mapping` giving new site number.
    wt_nt_col : str or None
        Column in `number_mapping` giving wildtype nucleotide at each site,
        or `None` to not check identity.
    err_suffix : str
        Append this message to any errors raised about invalid sites or
        mutation strings. Can be useful for debugging.

    Attributes
    ----------
    old_to_new_site : dict
        Maps old site number to new one.
    old_to_wt : dict or None
        Maps old site number to wildtype nucleotide if using `wt_nt_col`.

    Example
    --------
    >>> number_mapping = pd.DataFrame({'old': [1, 2, 3],
    ...                                'new': [5, 6, 7],
    ...                                'wt_nt': ['A', 'C', 'G']})
    >>> renumberer = MutationRenumber(number_mapping=number_mapping,
    ...                               old_num_col='old',
    ...                               new_num_col='new',
    ...                               wt_nt_col='wt_nt')
    >>> renumberer.old_to_new_site
    {1: 5, 2: 6, 3: 7}
    >>> renumberer.old_to_wt
    {1: 'A', 2: 'C', 3: 'G'}
    >>> renumberer.renumber_muts('A1C del2to3 ins3GC')
    'A5C del6to7 ins7GC'

    """

    def __init__(self, number_mapping, old_num_col, new_num_col, wt_nt_col,
                 *, err_suffix=''):
        """See main class docstring."""
        self._err_suffix = err_suffix
        for col in [old_num_col, new_num_col]:
            if col not in number_mapping.columns:
                raise ValueError(f"`number_mapping` lacks column {col}" +
                                 self._err_suffix)
            if number_mapping[col].dtype != int:
                raise ValueError(f"`number_mapping` column {col} not integer" +
                                 self._err_suffix)
        self.old_to_new_site = (number_mapping
                                .set_index(old_num_col)
                                [new_num_col]
                                .to_dict()
                                )
        if not (len(self.old_to_new_site) ==
                len(set(self.old_to_new_site.values())) ==
                len(number_mapping)):
            raise ValueError('site numbers not unique' + self._err_suffix)

        self._old_to_new_site_str = {str(old): str(new) for old, new
                                     in self.old_to_new_site.items()}

        if wt_nt_col:
            if wt_nt_col not in number_mapping.columns:
                raise ValueError(f"`number_mapping` lacks column {col}" +
                                 self._err_suffix)
            if not all(isinstance(nt, str) and len(nt) == 1
                       for nt in number_mapping[wt_nt_col]):
                raise ValueError(f"`number_mapping` column {col} not letters" +
                                 self._err_suffix)
            self.old_to_wt = (number_mapping
                              .set_index(old_num_col)
                              [wt_nt_col]
                              .to_dict()
                              )
            self._old_to_wt_str = {str(old): wt for old, wt
                                   in self.old_to_wt.items()}
        else:
            self.old_to_wt = None

    def renumber_muts(self, mut_str):
        """Get re-numbered mutation string.

        Parameters
        ----------
        mut_str : str
            Mutations in format 'A1C del2to3 ins3GG'.

        Returns
        -------
        str
            A version of `mut_str` where sites have been renumbered.

        """
        new_muts = []
        for mut in mut_str.split():
            try:
                # try to match substitutions
                m = re.fullmatch(
                        r'(?P<wt>[A-Z])(?P<site>\-?\d+)(?P<mut>[A-Z])',
                        mut)
                if m:
                    site = m.group('site')
                    if self.old_to_wt is not None:
                        if self._old_to_wt_str[site] != m.group('wt'):
                            expected_wt = self._old_to_wt_str[m.group('site')]
                            raise ValueError(f"Mutation {mut} invalid wt, "
                                             f"expected {expected_wt}" +
                                             self._err_suffix)
                    new_muts.append(
                            m.group('wt') +
                            self._old_to_new_site_str[site] +
                            m.group('mut')
                            )
                    continue
                # try to match insertion
                m = re.fullmatch(r'ins(?P<site>\-?\d+)(?P<insertion>[A-Z]+)',
                                 mut)
                if m:
                    new_muts.append(
                            'ins' +
                            self._old_to_new_site_str[m.group('site')] +
                            m.group('insertion')
                            )
                    continue
                # try to match deletion
                m = re.fullmatch(r'del(?P<site1>\-?\d+)to(?P<site2>\-?\d+)',
                                 mut)
                if m:
                    new_muts.append(
                            'del' +
                            self._old_to_new_site_str[m.group('site1')] +
                            'to' +
                            self._old_to_new_site_str[m.group('site2')]
                            )
                    continue
                # problem if we made it here, couldn't match anything
                raise ValueError(f"Cannot match {mut} in {mut_str}" +
                                 self._err_suffix)
            except KeyError:
                raise ValueError(f"Mutation {mut} site out of numbering range"
                                 + self._err_suffix)
        return ' '.join(new_muts)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
