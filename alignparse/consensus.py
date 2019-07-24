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
