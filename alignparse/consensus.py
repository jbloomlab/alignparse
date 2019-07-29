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


Mutations = collections.namedtuple(
                        'Mutations',
                        ['substitutions', 'deletions', 'insertions'])


_MUT_REGEX = {
    'substitution': re.compile(r'[ACGTN](?P<start>\d+)[ACGTN]'),
    'deletion': re.compile(r'del(?P<start>\d+)to\d+'),
    'insertion': re.compile(r'ins(?P<start>\d+)(?:len\d+|[ACGTN]+)'),
    }
"""dict: Mutation regular expression matches."""


def process_mut_str(s):
    """Process a string of mutations.

    Parameters
    ----------
    s : str
        Space-delimited mutations. Substitutions in form 'A6T'. Deletions
        in form 'del5to7'. Insertions in form 'ins6len2' or 'ins6TA'.

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
    >>> s = 'A11T del5to7 G9C ins6GA del12to15 ins13AAT'
    >>> process_mut_str(s)  # doctest: +NORMALIZE_WHITESPACE
    Mutations(substitutions=['G9C', 'A11T'],
              deletions=['del5to7', 'del12to15'],
              insertions=['ins6GA', 'ins13AAT'])

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


def fracident(df,
              *,
              group_cols='barcode',
              upstream_group_cols='library',
              mutation_col='mutations',
              sort_mutations=True,
              ):
    """Fraction of sequences with identical mutations in group (i.e., barcode).

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
        identical.
    sort_mutations : bool
        Useful if you have strings of space-delimited mutations not guaranteed
        to be consistently ordered. If `True`, sort such space-delimited
        mutations.

    Returns
    -------
    pandas.DataFrame
        Gives the fraction identical and implied accuracy of individual
        sequences. The columns are all all columns in `upstream_group_cols`,
        plus:

         - 'fraction_identical': average fraction of pairs of sequences in
           same group that are identical. This is the Simpson diversity index
           without replacement: https://en.wikipedia.org/wiki/Diversity_index

         - 'accuracy': the square root of the fraction identical, which is
           the accuracy of each sequence under the assumption that errors
           always make sequences different.

    Example
    -------
    >>> df = pd.DataFrame({
    ...         'barcode': ['AT', 'AT', 'TG', 'TA', 'TA', 'TA'],
    ...         'mutations': ['A1G', 'A1G', '', 'T2A C6G', 'T5C', 'C6G T2A']})
    >>> fracident(df, upstream_group_cols=None)
       fraction_identical  accuracy
    0                 0.5  0.707107

    If we do the same without sorting the mutations, we get lower values
    as 'T2A C6G' and 'C6G T2A' are then considered as different:

    >>> fracident(df, upstream_group_cols=None, sort_mutations=False)
       fraction_identical  accuracy
    0                0.25       0.5

    Now another example with two libraries and non-standard column
    names. Note that we only get results for the two libraries
    with barcodes found multiple times:

    >>> df2 = pd.DataFrame({
    ...     'bc'  :['AT', 'AT', 'TG', 'TA', 'TA', 'TA', 'TA'],
    ...     'var' :['v1', 'v1', 'v2', 'v3', 'v4', 'v3', 'v3'],
    ...     'lib' :['s1', 's1', 's2', 's3', 's3', 's3', 's4']})
    >>> fracident(df2, upstream_group_cols='lib', group_cols='bc',
    ...           mutation_col='var', sort_mutations=False)
      lib  fraction_identical  accuracy
    0  s1            1.000000   1.00000
    1  s3            0.333333   0.57735

    """
    # column names used in calculations below
    reserved_cols = ['_groupcounts', '_sequencecounts', '_npair',
                     '_simpson_diversity', '_weight_diversity', '_dummy',
                     'accuracy', 'fraction_identical']
    for col in reserved_cols:
        if col in df.columns:
            raise ValueError(f"`df` cannot have column named {col}")

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
        # get just sequences that have a barcode found multiple times
        .assign(_groupcounts=1)
        .assign(_groupcounts=lambda x: x.groupby(cols).transform('count'))
        .query('_groupcounts > 1')
        # sort mutations
        .assign(**{mutation_col: (lambda x:
                                  x[mutation_col]
                                  .map(lambda s: ' '.join(sorted(s.split())))
                                  if sort_mutations else x[mutation_col]
                                  )}
                )
        # within each barcode, count number of sequences of each mutation combo
        .groupby(cols + ['_groupcounts', mutation_col])
        .size()
        .rename('_sequencecounts')
        .reset_index()
        # compute Simpson diversity without replacement for each barcode
        .groupby(cols + ['_groupcounts'])
        .apply(lambda x: ((x['_sequencecounts'] * (x['_sequencecounts'] - 1)) /
                          (x['_groupcounts'] * (x['_groupcounts'] - 1))).sum())
        .reset_index(name='_simpson_diversity')
        # compute weighted average of fraction identical across all pairs
        .assign(
            _npair=lambda x: x['_groupcounts'] * (x['_groupcounts'] - 1) / 2,
            _weight_diversity=lambda x: x['_npair'] * x['_simpson_diversity'])
        .groupby(upstream_group_cols)
        .apply(lambda x: x['_weight_diversity'].sum() / x['_npair'].sum())
        .reset_index(name='fraction_identical')
        # estimate accuracy as square root of fraction identical
        .assign(accuracy=lambda x: numpy.sqrt(x['fraction_identical']))
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
        :func:`processe_mut_str`.
    max_sub_diffs : int
        Drop any group where any variant differs from all other variants
        by more than this many substitution (point mutation) differences.
    max_indel_diffs : int
        Drop any group where any variant differs from all other variants
        by more than this many indel differences.
    max_minor_sub_frac : float
        Drop groups with a minor (non-consensus) substitution in more than the
        **ceiling** of this fraction times the number of sequences in group.
    max_minor_indel_frac : float
        Drop groups with a minor (non-consensus) indel in more than the
        **ceiling** of this fraction times the number of sequences in group.
    support_col : str
          Name of column in returned `consensus` data frame with number
          of sequences supporting the consensus call.

    Note
    ----
    The rationale behind this consensus calling scheme is that we want
    to build a consensus except in the following two cases, each of which
    indicate a likely problem **beyond** simple rare sequencing errors:

      1. Any two sequences within the group differ by too much
         (given by `max_sub_diffs` and `max_indel_diffs`).

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
    ...        lib_2,TA,C5A T6C
    ...        lib_2,TA,ins5len1 T6C
    ...        lib_2,TA,T6C
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
    1   lib_1      TA  G3A ins4len3                     1
    2   lib_2      GG                                   2
    3   lib_2      TA           T6C                     3
    >>> dropped
      library barcode              drop_reason  nseqs
    0   lib_2      AA  minor subs too frequent      4
    1   lib_2      TG      subs diff too large      2
    2   lib_3      AA    indels diff too large      2

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

        nseqs = len(g_df)
        half_nseqs = 0.5 * nseqs
        assert nseqs > 0

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

            # are max_sub_diffs and max_indel_diffs satisfied?
            for m1set, m2set in itertools.combinations(mutlists[mtype], 2):
                ndiffs = len(numpy.setxor1d(m1set, m2set, assume_unique=True))
                if ndiffs > maxd:
                    drop_reason = f"{mtype} diff too large"
                    break
            if drop_reason is not None:
                break

            # see if max_minor_mut_frac is satisfied
            max_muts = math.ceil(max_frac * nseqs)
            nseqs_minus_max_muts = nseqs - max_muts
            counts = collections.Counter(itertools.chain.from_iterable(
                            mutlists[mtype]))
            if any(max_muts < count < nseqs_minus_max_muts for count
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
