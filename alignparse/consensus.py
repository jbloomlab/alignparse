"""
=========
consensus
=========

Get consensus of sequences / mutations. Useful for grouping sequences
within a barcode.

"""


import collections
import io  # noqa: F401
import itertools
import math
import re
import textwrap  # noqa: F401

import numpy

import pandas as pd

import scipy.optimize
import scipy.special


Mutations = collections.namedtuple(
                        'Mutations',
                        ['substitutions', 'deletions', 'insertions'])


_MUT_REGEX = {
    'substitution': re.compile(r'[ACGTN](?P<start>\-?\d+)[ACGTN\-]'),
    'deletion': re.compile(r'del(?P<start>\-?\d+)to(?P<end>\-?\d+)'),
    'insertion': re.compile(r'ins(?P<start>\-?\d+)(?:len\d+|[ACGTN]+)'),
    }
"""dict: Mutation regular expression matches."""


def process_mut_str(s):
    """Process a string of mutations.

    Parameters
    ----------
    s : str
        Space-delimited mutations. Substitutions in form 'A6T'. Deletions
        in form 'del5to7'. Insertions in form 'ins6len2' or 'ins6TA'.
        Negative site numbers are allowed.

    Returns
    -------
    namedtuple
        The tuple is named `Mutations` and has the following elements in order:
          - `substitutions`: list of substitutions in order of sequence.
          - `deletions`: list of deletions in order of sequence.
          - `insertions`: list of insertions in order of sequence.

    Example
    -------
    >>> s = 'A1T del5to7 G9C ins6len2 del12to12'
    >>> process_mut_str(s)  # doctest: +NORMALIZE_WHITESPACE
    Mutations(substitutions=['A1T', 'G9C'],
              deletions=['del5to7', 'del12to12'],
              insertions=['ins6len2'])
    >>> s = 'A11T del5to7 G-9C ins-6GA del12to15 ins13AAT'
    >>> process_mut_str(s)  # doctest: +NORMALIZE_WHITESPACE
    Mutations(substitutions=['G-9C', 'A11T'],
              deletions=['del5to7', 'del12to15'],
              insertions=['ins-6GA', 'ins13AAT'])

    """
    mut_lists = collections.defaultdict(list)
    for mut_str in s.split():
        for mut_type, regex in _MUT_REGEX.items():
            match = regex.fullmatch(mut_str)
            if match:
                start = int(match.group('start'))
                mut_lists[mut_type].append((start, mut_str))
                break
        else:
            raise ValueError(f"cannot match mutation {mut_str}\n(in {s})")

    if any(len(mlist) != len(set(mlist)) for mlist in mut_lists.values()):
        raise ValueError(f"duplicate mutation in: {s}")

    return Mutations(
            substitutions=[m for _, m in sorted(mut_lists['substitution'])],
            deletions=[m for _, m in sorted(mut_lists['deletion'])],
            insertions=[m for _, m in sorted(mut_lists['insertion'])],
            )


def add_mut_info_cols(df,
                      *,
                      mutation_col='mutations',
                      sub_str_col=None,
                      del_str_col=None,
                      ins_str_col=None,
                      indel_str_col=None,
                      n_sub_col=None,
                      n_del_col=None,
                      n_ins_col=None,
                      n_indel_col=None,
                      overwrite_cols=False,
                      ):
    """Expand information about mutations in a data frame.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame with each row giving information about a variant.
    mutation_col : str
        Name of column in `df` with str containing space-delimited
        list of mutations in format parseable by :func:`process_mut_str`.
    sub_str_col : str or None
        Name of added column with string giving substitutions.
    del_str_col : str or None
        Name of added column with string giving deletions.
    ins_str_col : str or None
        Name of added column with string giving insertions.
    indel_str_col : str or None
        Name of added column with string giving deletions and insertions.
    n_sub_col : str or None
        Name of added column with number of substitutions.
    n_del_col : str or None
        Name of added column with number of deletions.
    n_ins_col : str or None
        Name of added column with number of insertions.
    n_indel_col : str or None
        Name of added column with number of deletions and insertions.
    overwrite_cols : bool
        If column to be created is already in `df`, overwrite it?

    Returns
    -------
    pandas.DataFrame
        A **copy** of `df` with the indicated columns added.

    Example
    -------
    Data frame for which we get mutation info:

    >>> df = pd.DataFrame({
    ...       'name': ['seq1', 'seq2', 'seq3'],
    ...       'mutations': ['A6G del7to8 T3A', '', 'T5A ins5AAT del8to9'],
    ...       })

    Get info on substitutions:

    >>> add_mut_info_cols(df, sub_str_col='subs', n_sub_col='nsubs')
       name            mutations     subs  nsubs
    0  seq1      A6G del7to8 T3A  T3A A6G      2
    1  seq2                                    0
    2  seq3  T5A ins5AAT del8to9      T5A      1

    Get info on substitutions and indels:

    >>> pd.set_option('max_columns', 10)
    >>> add_mut_info_cols(df, sub_str_col='subs', n_sub_col='nsubs',
    ...                   indel_str_col='indels', n_indel_col='nindels')
       name            mutations     subs           indels  nsubs  nindels
    0  seq1      A6G del7to8 T3A  T3A A6G          del7to8      2        1
    1  seq2                                                     0        0
    2  seq3  T5A ins5AAT del8to9      T5A  del8to9 ins5AAT      1        2

    Just get counts for all types of mutations:

    >>> add_mut_info_cols(df, n_sub_col='nsubs', n_del_col='ndels',
    ...                   n_ins_col='nins', n_indel_col='nindels')
       name            mutations  nsubs  ndels  nins  nindels
    0  seq1      A6G del7to8 T3A      2      1     0        1
    1  seq2                           0      0     0        0
    2  seq3  T5A ins5AAT del8to9      1      1     1        2

    """
    df = df.reset_index(drop=True)

    if mutation_col not in df.columns:
        raise ValueError(f"`df` lacks `mutation_col` {mutation_col}")

    new_cols = collections.OrderedDict(
                [(key, val)
                 for key, val in [('substitutions_str', sub_str_col),
                                  ('deletions_str', del_str_col),
                                  ('insertions_str', ins_str_col),
                                  ('indel_str', indel_str_col),
                                  ('substitutions_n', n_sub_col),
                                  ('deletions_n', n_del_col),
                                  ('insertions_n', n_ins_col),
                                  ('indel_n', n_indel_col)]
                 if val is not None])
    if mutation_col in set(new_cols.values()):
        raise ValueError('`mutation_col` also name of col to add')
    already_has = set(new_cols.values()).intersection(set(df.columns))
    if already_has and not overwrite_cols:
        raise ValueError(f"`df` already has these columns: {already_has}")

    def _mut_info(s):
        m = process_mut_str(s)
        returnlist = []
        for col in new_cols.keys():
            mut_type, col_type = col.split('_')
            if mut_type == 'indel':
                mlist = m.deletions + m.insertions
            else:
                mlist = getattr(m, mut_type)
            if col_type == 'str':
                returnlist.append(' '.join(mlist))
            else:
                assert col_type == 'n'
                returnlist.append(len(mlist))
        return returnlist

    mut_info_df = (pd.DataFrame(
                        [_mut_info(m) for m in df[mutation_col].values],
                        columns=new_cols.values())
                   .reindex(index=df.index)
                   )

    return (df
            .drop(columns=new_cols.values(), errors='ignore')
            .join(mut_info_df)
            )


class _LnL_error_rate:
    """Log likelihood of sequences per group as function of error rate.

    Note
    ----
    Used as utility function for computing the empirical accuracy by
    :func:`empirical_accuracy`. See the docs for that function to understand
    the notation / nomenclature used here.

    Parameters
    ----------
    df : pandas.DataFrame
        A different row for each group of sequences.
    n_col : str
        Column in `df` with number of sequences per group.
    u_col : str
        Column in `df` with number of unique sequences per group.
    count_col : str
        Column in `df` with number of groups with this value of :math:`n`
        and :math:`u`.

    Example
    -------
    >>> df = pd.DataFrame(
    ...         [[4, 3, 1],
    ...          [4, 2, 2],
    ...          [4, 1, 1],
    ...          ],
    ...         columns=['n', 'u', 'count'],
    ...         )
    >>> lnl = _LnL_error_rate(df, n_col='n', u_col='u', count_col='count')
    >>> round(lnl.maxlik_eps(), 3)
    0.25

    """

    def __init__(self, df, *, n_col, u_col, count_col):
        """See main class docstring."""
        self._df = (
            df
            .assign(
                n=lambda x: x[n_col],
                u=lambda x: x[u_col],
                count=lambda x: x[count_col],
                binom=lambda x: scipy.special.binom(x['n'], x['u'] - 1),
                delta_un=lambda x: (x['n'] == x['u']).astype(int),
                )
            )

    def lnlik(self, eps):
        """Log likelihood for error rate `eps`."""
        return sum(self._df['count'] *
                   numpy.log(self._df['binom'] *
                             (1 - eps)**(self._df['n'] - self._df['u'] + 1) *
                             eps**(self._df['u'] - 1) +
                             self._df['delta_un'] * eps**self._df['n']
                             )
                   )

    def neg_lnlik(self, eps):
        """Negative log likelihood for error rate `epsilon`."""
        return -self.lnlik(eps)

    def maxlik_eps(self):
        """Maximum likelihood value of error rate `epsilon`."""
        res = scipy.optimize.minimize_scalar(self.neg_lnlik,
                                             bounds=(1e-8, 1 - 1e-8),
                                             method='bounded')
        if not res.success:
            raise RuntimeError(f"optimization failed:\n{res}")
        return res.x


def empirical_accuracy(df,
                       *,
                       group_cols='barcode',
                       upstream_group_cols='library',
                       mutation_col='mutations',
                       accuracy_col='accuracy',
                       sort_mutations=True,
                       ):
    r"""Accuracy from number of identical sequences in a group (i.e., barcode).

    Note
    ----
    The accuracy :math:`1 - \epsilon` is calculated as follows:

    Let :math:`\epsilon` be the probability a sequence has an error.
    Consider a group of :math:`n` sequences that should be identical in
    the absence of errors (i.e., they all have the same barcode), and
    assume errors only make sequences **different** (the chance of two
    sequences having the same error is very low). We would like to
    calculate the probability of observing :math:`u` unique sequences in
    the group given :math:`\epsilon` and :math:`n`.

    First, consider the case where there are two sequences in the group
    (:math:`n = 2`). The probability of having :math:`u = 1` unique
    sequences in the group is the probability that neither have an error:

    .. math::

       \Pr\left(u=1 | n=2, \epsilon\right) = \left(1 - \epsilon\right)^2.

    The probability of having :math:`u = 2` unique sequences is the
    the probability that either one or both have errors:

    .. math::

       \Pr\left(u=2 | n=2, \epsilon\right) = 2\left(1 - \epsilon\right)\epsilon
       + \epsilon^2.

    Generalizing to arbitrary :math:`n`, we have:

    .. math::

       \Pr\left(u | n, \epsilon\right) = \binom{n}{u-1}
       \left(1 - \epsilon\right)^{n - u + 1} \epsilon^{u - 1} +
       \delta_{un} \epsilon^n

    where :math:`\binom{n}{u-1}` is the binomial coefficient and
    :math:`\delta_{un}` is the Kronecker delta.

    Let there be :math:`G` groups of sequences, with the size and number
    of unique sequences in group :math:`g` be :math:`n_g` and :math:`u_g`
    (where :math:`g = 1, \ldots, G`).
    The overall likelihood :math:`L` of the observed numbers of unique
    sequences given the group sizes and error rate is:

    .. math::

       L &=& \Pr\left(\left\{u_g\right\}|\left\{n_g\right\},\epsilon\right) \\
       &=& \prod_g \Pr\left(u_g | n_g, \epsilon\right).

    To find the maximum likelihood of error rate, we simply use numerical
    optimization to find the value of :math:`\epsilon` that maximizes
    :math:`L`. In practice, we actually maximize :math:`\ln\left(L\right)`.

    Parameters
    ----------
    df : pandas.DataFrame
        Each row gives data for a different sequence.
    group_cols : str or list
        Column(s) in `df` indicating how we group sequences before computing
        fraction within group that is identical.
    upstream_group_cols : str or None or list
        Column(s) in `df` that we use to group **before** doing the grouping
        for computing fraction identical. So a different average fraction
        will be returned for each of these upstream groups.
    mutation_col : str
        Column in `df` that we compare when determining if sequences are
        identical. Typically list of mutations or sequences themselves.
    accuracy_col : str
        Name of column in returned data frame that gives empirical accuracy.
    sort_mutations : bool
        Useful if you have strings of space-delimited mutations not guaranteed
        to be consistently ordered. If `True`, sort such space-delimited
        mutations.

    Returns
    -------
    pandas.DataFrame
        Has all columns in `upstream_group_cols` plus the new column with
        name given by `accuracy_col`.

    Example
    -------
    >>> df = pd.DataFrame({
    ...         'barcode': ['AT', 'AT', 'TG', 'TA', 'TA', 'TA'],
    ...         'mutations': ['A1G', 'A1G', '', 'T2A C6G', 'T5C', 'C6G T2A']})
    >>> with pd.option_context('precision', 4):
    ...     empirical_accuracy(df, upstream_group_cols=None)
       accuracy
    0    0.7692

    If we do the same without sorting the mutations, we get lower values
    as 'T2A C6G' and 'C6G T2A' are then considered as different:

    >>> with pd.option_context('precision', 4):
    ...     empirical_accuracy(df, upstream_group_cols=None,
    ...                        sort_mutations=False)
       accuracy
    0    0.4766

    Now another example with two libraries and non-standard column
    names. Note that we only get results for the two libraries
    with barcodes found multiple times:

    >>> df2 = pd.DataFrame({
    ...     'bc'  :['AT', 'AT', 'TG', 'TA', 'TA', 'TA', 'TA'],
    ...     'var' :['v1', 'v1', 'v2', 'v3', 'v4', 'v3', 'v3'],
    ...     'lib' :['s1', 's1', 's2', 's3', 's3', 's3', 's4']})
    >>> with pd.option_context('precision', 4):
    ...     empirical_accuracy(df2, upstream_group_cols='lib', group_cols='bc',
    ...                        mutation_col='var', sort_mutations=False)
      lib  accuracy
    0  s1    1.0000
    1  s3    0.6667

    """
    reserved_cols = ['_n', '_u', '_ngroups', '_dummy']
    for col in reserved_cols:
        if col in df.columns:
            raise ValueError(f"`df` cannot have column named {col}")

    if accuracy_col in df.columns:
        raise ValueError(f"`df` has column `accuracy_col` {accuracy_col}")

    if isinstance(group_cols, str):
        group_cols = [group_cols]
    if (upstream_group_cols is None) or upstream_group_cols == []:
        drop_upstream_col = True
        upstream_group_cols = '_dummy'
        df = df.assign(**{upstream_group_cols: 'dummy'})
    else:
        drop_upstream_col = False
    if isinstance(upstream_group_cols, str):
        upstream_group_cols = [upstream_group_cols]

    cols = group_cols + upstream_group_cols
    if len(set(cols)) != len(cols):
        raise ValueError('duplicate in `group_cols` / `upstream_group_cols`')
    if mutation_col in cols:
        raise ValueError(f"`mutation_col` {mutation_col} also grouping column")
    for col in cols + [mutation_col]:
        if col not in df.columns:
            raise ValueError(f"no column {col} in `df`: {list(df.columns)}")

    result = (
        df
        # number of sequences in each group
        .assign(_n=lambda x: x.groupby(cols)[mutation_col].transform('count'))
        # only retain groups with >1 sequence
        .query('_n > 1')
        # sort mutations
        .assign(**{mutation_col: (lambda x:
                                  x[mutation_col]
                                  .map(lambda s: ' '.join(sorted(s.split())))
                                  if sort_mutations else x[mutation_col]
                                  )}
                )
        # number of unique sequences in each group
        .assign(_u=lambda x: (x.groupby(cols)[mutation_col]
                              .transform('nunique'))
                )
        # number of groups with each combination of n and u
        .groupby(upstream_group_cols + ['_n', '_u'])
        .size()
        .rename('_ngroups')
        .reset_index()
        # get error rate
        .groupby(upstream_group_cols)
        .apply(lambda x: 1 - _LnL_error_rate(x,
                                             n_col='_n',
                                             u_col='_u',
                                             count_col='_ngroups'
                                             ).maxlik_eps()
               )
        .rename(accuracy_col)
        .reset_index()
        )

    if drop_upstream_col:
        assert len(upstream_group_cols) == 1
        result = result.drop(upstream_group_cols, axis='columns')

    return result


def simple_mutconsensus(df,
                        *,
                        group_cols=('library', 'barcode'),
                        mutation_col='mutations',
                        max_sub_diffs=1,
                        max_indel_diffs=2,
                        max_minor_sub_frac=0.1,
                        max_minor_indel_frac=0.25,
                        max_minor_greater_or_equal=False,
                        min_support=1,
                        support_col='variant_call_support',
                        ):
    """Get simple consensus of mutations with group (i.e., barcode).

    Parameters
    ----------
    df : pandas.DataFrame
        Each row gives data on mutations for a different sequence.
    group_cols : list, tuple, or str
        Group all sequences that share values in these column(s) of `df`.
    mutation_col : str
        Column in `df` with mutations in form that can be processed by
        :func:`process_mut_str`.
    max_sub_diffs : int or None
        Drop groups where any variant differs from all other variants
        by more than this many substitution (point mutation) differences.
        Set to ``None`` if no limit.
    max_indel_diffs : int or None
        Drop groups where any variant differs from all other variants
        by more than this many indel differences. Set to ``None`` if no limit.
    max_minor_sub_frac : float
        Drop groups with a minor (non-consensus) substitution in > the
        **ceiling** of this fraction times the number of sequences in group.
    max_minor_indel_frac : float
        Drop groups with a minor (non-consensus) indel in > the
        **ceiling** of this fraction times the number of sequences in group.
    max_minor_greater_or_equal : bool
        For ``max_minor_sub_frac`` and ``max_minor_indel_frac``, use >= rather
        than >. This is may help if you want to be conservative in
        avoiding bad consensuses. But don't reject if zero differences.
    min_support : int
        Require at least this many supporting sequences to call consensus.
    support_col : str
          Name of column in returned `consensus` data frame with number
          of sequences supporting the consensus call.

    Note
    ----
    The rationale behind this consensus calling scheme is that we want
    to build a consensus except in the following two cases, each of which
    indicate a likely problem **beyond** simple rare sequencing errors (such
    as strand exchange or barcode collisions):

      1. A sequence within the group differs by too much from others in
         group (given by `max_sub_diffs` and `max_indel_diffs`).

      2. There are "minor" (non-consensus) mutations present at too high
         frequency (given by `max_minor_indel_frac` and `max_minor_sub_frac`).

    Returns
    -------
    tuple
        The returned tuple `(consensus, dropped)` consists of two
        pandas data frames:

          - `consensus`: Gives consensus for all groups for which
            this can be constructed. Columns are `group_cols`, column
            with name given by `mutation_col` that gives consensus mutations,
            and a column with name given by `support_col` that gives the
            number of variants used to generate the consensus.

          - `dropped`: Gives each dropped group, the reason it was dropped,
            and the number of variants (sequences) in that group.

    Example
    -------
    Create a dummy CSV file with our variants and mutations:

    >>> f = io.StringIO()
    >>> _ = f.write(textwrap.dedent('''
    ...        library,barcode,mutations
    ...        lib_1,AG,A2C del5to7
    ...        lib_1,AG,A2C
    ...        lib_1,TA,G3A ins4len3
    ...        lib_1,TA,G3A ins4len3
    ...        lib_1,TA,G3A ins4len3
    ...        lib_2,TA,C5A T-6C
    ...        lib_2,TA,ins5len1 T-6C
    ...        lib_2,TA,T-6C
    ...        lib_2,TA,A4G T-6C
    ...        lib_2,TG,T6A
    ...        lib_2,TG,A2G
    ...        lib_2,GG,del1to4
    ...        lib_2,GG,A1C
    ...        lib_2,AA,
    ...        lib_2,AA,
    ...        lib_2,AA,T6C
    ...        lib_2,AA,T6C
    ...        lib_3,AA,ins1len1 del1to2 T6G
    ...        lib_3,AA,del5to7 T6G
    ...        '''))
    >>> _ = f.seek(0)

    Read into data frame. Note use of `na_filter=False` so empty strings
    are read as such rather than as `na` values:

    >>> df = pd.read_csv(f, na_filter=False)

    Get the consensus mutations, and the dropped libraries / barcodes:

    >>> consensus, dropped = simple_mutconsensus(df)
    >>> consensus
      library barcode     mutations  variant_call_support
    0   lib_1      AG           A2C                     2
    1   lib_1      TA  G3A ins4len3                     3
    2   lib_2      GG                                   2
    3   lib_2      TA          T-6C                     4
    >>> dropped
      library barcode              drop_reason  nseqs
    0   lib_2      AA  minor subs too frequent      4
    1   lib_2      TG      subs diff too large      2
    2   lib_3      AA    indels diff too large      2

    Get the consensus just for library 1:

    >>> lib1_consensus, _ = simple_mutconsensus(df.query('library == "lib_1"'),
    ...                                         group_cols='barcode')
    >>> lib1_consensus
      barcode     mutations  variant_call_support
    0      AG           A2C                     2
    1      TA  G3A ins4len3                     3

    Set ``max_sub_diffs`` to None:

    >>> consensus, dropped = simple_mutconsensus(df, max_sub_diffs=None)
    >>> consensus
      library barcode     mutations  variant_call_support
    0   lib_1      AG           A2C                     2
    1   lib_1      TA  G3A ins4len3                     3
    2   lib_2      GG                                   2
    3   lib_2      TA          T-6C                     4
    4   lib_2      TG                                   2
    >>> dropped
      library barcode              drop_reason  nseqs
    0   lib_2      AA  minor subs too frequent      4
    1   lib_3      AA    indels diff too large      2

    Set ``max_minor_greater_or_equal`` to True:

    >>> consensus, dropped = simple_mutconsensus(
    ...         df, max_minor_greater_or_equal=True)
    >>> consensus
      library barcode     mutations  variant_call_support
    0   lib_1      TA  G3A ins4len3                     3
    >>> dropped
      library barcode                drop_reason  nseqs
    0   lib_1      AG  minor indels too frequent      2
    1   lib_2      AA    minor subs too frequent      4
    2   lib_2      GG    minor subs too frequent      2
    3   lib_2      TA    minor subs too frequent      4
    4   lib_2      TG        subs diff too large      2
    5   lib_3      AA      indels diff too large      2

    Use ``min_support``:

    >>> consensus, dropped = simple_mutconsensus(df, min_support=3)
    >>> consensus
      library barcode     mutations  variant_call_support
    0   lib_1      TA  G3A ins4len3                     3
    1   lib_2      TA          T-6C                     4
    >>> dropped
      library barcode              drop_reason  nseqs
    0   lib_1      AG        too few sequences      2
    1   lib_2      AA  minor subs too frequent      4
    2   lib_2      GG        too few sequences      2
    3   lib_2      TG        too few sequences      2
    4   lib_3      AA        too few sequences      2

    """
    if isinstance(group_cols, str):
        group_cols = [group_cols]
    else:
        group_cols = list(group_cols)
    if set(group_cols) > set(df.columns):
        raise ValueError(f"`df` lacks `group_cols`: {group_cols}")
    if mutation_col in group_cols:
        raise ValueError(f"`mutation_col` {mutation_col} in `group_cols`")
    if mutation_col not in df.columns:
        raise ValueError(f"`df` lacks `mutation_col`: {mutation_col}")

    df = df[group_cols + [mutation_col]]
    if df.isnull().values.any():
        raise ValueError('`df` contains `na` (null) entries. It is possible '
                         'that you read a CSV without using `na_filter=False` '
                         'so that empty mutation strings were read as `na` '
                         'rather than empty strings.')

    dropped = []
    consensus = []
    for g, g_df in df.groupby(group_cols, observed=True)[mutation_col]:

        if len(group_cols) == 1:
            g = [g]

        nseqs = len(g_df)
        half_nseqs = 0.5 * nseqs
        assert nseqs > 0

        if nseqs < min_support:
            dropped.append((*g, 'too few sequences', nseqs))
            continue

        mutations = [process_mut_str(s) for s in g_df.values]

        drop_reason = None
        g_consensus = []

        mutlists = {'subs': [m.substitutions for m in mutations],
                    'indels': [m.deletions + m.insertions for m in mutations],
                    }

        for mtype, maxd, max_frac in [
                        ('subs', max_sub_diffs, max_minor_sub_frac),
                        ('indels', max_indel_diffs, max_minor_indel_frac),
                        ]:
            if maxd is not None:
                # is max_sub_diffs or max_indel_diffs satisfied?
                ndiffs_by_seq = collections.defaultdict(list)
                for (i1, m1set), (i2, m2set) in itertools.combinations(
                         enumerate(mutlists[mtype]), 2):
                    ndiffs = len(numpy.setxor1d(m1set, m2set,
                                                assume_unique=True))
                    ndiffs_by_seq[i1].append(ndiffs)
                    ndiffs_by_seq[i2].append(ndiffs)
                if any(min(ndifflist) > maxd for ndifflist
                       in ndiffs_by_seq.values()):
                    drop_reason = f"{mtype} diff too large"
                    break

            # see if max_minor_mut_frac is satisfied
            max_muts = math.ceil(max_frac * nseqs)
            nseqs_minus_max_muts = nseqs - max_muts
            counts = collections.Counter(itertools.chain.from_iterable(
                            mutlists[mtype]))
            if max_minor_greater_or_equal:
                if any((max_muts <= count <= nseqs_minus_max_muts)
                       and (count != nseqs)
                       for count in counts.values()):
                    drop_reason = f"minor {mtype} too frequent"
                    break
            elif any(max_muts < count < nseqs_minus_max_muts for count
                     in counts.values()):
                drop_reason = f"minor {mtype} too frequent"
                break

            # consensus for mutation type
            g_consensus += [m for m, c in counts.items() if c > half_nseqs]

        if drop_reason is not None:
            dropped.append((*g, drop_reason, nseqs))
        else:
            consensus.append((*g, ' '.join(g_consensus), nseqs))

    consensus = pd.DataFrame(consensus,
                             columns=group_cols + [mutation_col, support_col])
    dropped = pd.DataFrame(dropped,
                           columns=group_cols + ['drop_reason', 'nseqs'])

    return consensus, dropped


if __name__ == '__main__':
    import doctest
    doctest.testmod()
