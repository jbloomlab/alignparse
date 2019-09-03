"""
=====
utils
=====

"""


import math
import numbers

import numpy


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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
