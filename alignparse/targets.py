"""
=======
targets
=======

Defines :class:`Targets`, which holds a several :class:`Target` objects that
define the alignment targets. Each :class:`Target` in a collection of
:class:`Targets` has some :class:`Feature` regions.

"""


class Feature:
    """A sequence feature within a :class:`Target` sequence.

    Parameters
    ----------
    name : str
    seq : str
    start: int
    end : int

    Attributes
    ----------
    name : str
        Name of feature.
    seq : str
        Sequence of feature in strand sense of :class:`Target`.
    start : int
        Feature start in :class:`Target`, using Python-like 0, ... numbering.
    end : int
        Feature end in :class:`Target` using Python-like 0, ... numbering.
    length: int
        Length of feature.

    """

    def __init__(self, *, name, seq, start, end):
        """See main class docstring."""
        self.name = name
        self.seq = seq
        if end - start != len(seq):
            raise ValueError('length of `seq` not equal to `end` - `start`')
        self.end = end
        self.start = start
        self.length = end - start


class Target:
    """A single target sequence.

    Parameters
    ----------
    seqrecord : Bio.SeqRecord.SeqRecord
        BioPython sequence record of target. Must have `seq`, `name`,
        and `features` attributes.
    req_features : set or other iterable
        Required features in `seqrecord`.
    opt_features: set of other iterable
        Optional features in `seqrecord`.
    allow_extra_features : bool
        Can `seqrecord` have features not in `req_features` or `opt_features`?

    Attributes
    ----------
    seq : str
        Full sequence of target.
    features : list
        List of all features as :class:`Feature` objects.
    """

    def __init__(self, *, seqrecord, req_features=set(), opt_features=set(),
                 allow_extra_features=False):
        """See main class docstring."""
        raise RuntimeError('not yet implemented')

    def get_feature(self, name):
        """Get :class:`Feature` by name.

        Parameters
        ----------
        name : str
            Name of :class:`Feature`.

        Returns
        -------
        :class:`Feature`
            Returns the feature, or raises `ValueError` if no such feature.

        """
        raise RuntimeError('not yet implemented')

    def draw(self):
        """Draws the feature."""
        raise RuntimeError('not yet implemented')


class Targets:
    """Collection of :class:`Target` sequences.

    Parameters
    ----------
    seqsfile : str or list
        Name of file specifying the targets, or list of such files. So
        if multiple targets they can all be in one file or in separate files.
    req_features : set or other iterable
        Required features for each target in `seqsfile`.
    opt_features: set of other iterable
        Optional features for each target in `seqsfile`.
    allow_extra_features : bool
        Can targets have features not in `req_features` or `opt_features`?
    seqsfileformat : {'genbank'}
        Format of `seqsfile`.

    Attributes
    ----------
    targets : list
        List of all :class:`Target` objects.

    """

    def __init__(self, *, seqsfile, req_features=set(), opt_features=set(),
                 allow_extra_features=False, seqsfileformat='genbank'):
        """See main class docstring."""
        raise RuntimeError('not yet implemented')

    def get_target(self, name):
        """Get :class:`Target` by name.

        Parameters
        ----------
        name : str
            Name of :class:`Target`.

        Returns
        -------
        :class:`Target`
            Returns the target, or raises `ValueError` if no such target.

        """
        raise RuntimeError('not yet implemented')

    def write_fasta(self, fastafile):
        """Write all targets to a FASTA file.

        Parameters
        ----------
        filename : str or file-like object.
            Name of created FASTA file, or file-like object to write to.
        
        """
        raise RuntimeError('not yet implemented')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
