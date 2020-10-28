"""
=======
targets
=======

Defines :class:`Targets`, which holds :class:`Target` objects that define
alignment targets. Each :class:`Target` has some :class:`Feature` regions.

"""


import contextlib
import copy
import itertools
import os
import re
import tempfile

import Bio.SeqIO

import dna_features_viewer

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt

import pandas as pd

import pathos

import pysam

import yaml

from alignparse.constants import CBPALETTE
from alignparse.cs_tag import (Alignment,
                               cs_to_mutation_str,
                               cs_to_nt_mutation_count,
                               cs_to_op_mutation_count,
                               cs_to_sequence,
                               )


class Feature:
    """A sequence feature within a :class:`Target` sequence.

    Parameters
    ----------
    name : str
        Name of feature.
    seq : str
        Sequence of feature.
    start: int
        Feature start in :class:`Target`, using Python-like 0, ... numbering.
    end : int
        Feature end in :class:`Target` using Python-like 0, ... numbering.

    Attributes
    ----------
    name : str
        Name of feature.
    seq : str
        Sequence of feature.
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
        if ',' in self.name:
            raise ValueError(f"comma not allowed in feature name: {self.name}")
        self.seq = seq
        if end - start != len(seq):
            raise ValueError('length of `seq` not equal to `end` - `start`')
        self.end = end
        self.start = start
        self.length = end - start

    def __repr__(self):
        """Get string representation."""
        return (f"{self.__class__.__name__}(name={self.name}, seq={self.seq}, "
                f"start={self.start}, end={self.end})")


class Target:
    """A single target sequence.

    Parameters
    ----------
    seqrecord : Bio.SeqRecord.SeqRecord
        BioPython sequence record of target. Must have `seq`, `name`,
        and `features` attributes. Currently only handles + strand features.
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
    name : str
        Name of target.
    length : str
        Length of sequence.
    features : list
        List of all features as :class:`Feature` objects.
    feature_names : list
        List of names of all features.

    """

    def __repr__(self):
        """Get string representation."""
        return (f"{self.__class__.__name__}(name={self.name}, seq={self.seq}, "
                f"features={self.features})")

    def __init__(self, *, seqrecord, req_features=frozenset(),
                 opt_features=frozenset(), allow_extra_features=False):
        """See main class docstring."""
        self.name = self.get_name(seqrecord)
        if ',' in self.name:
            raise ValueError(f"comma not allowed in target name: {self.name}")
        if not hasattr(seqrecord, 'seq'):
            raise ValueError('`seqrecord` does not define a seq')
        self.seq = str(seqrecord.seq)

        self.length = len(self.seq)

        allow_features = set(req_features) | set(opt_features)

        self._features_dict = {}
        self.features = []
        self.feature_names = []
        for bio_feature in seqrecord.features:
            feature_name = bio_feature.type
            if feature_name in self._features_dict:
                raise ValueError(f"duplicate feature {feature_name} when "
                                 f"creating Target {self.name}")
            if not (allow_extra_features or (feature_name in allow_features)):
                raise ValueError(f"feature {feature_name} not allowed feature")
            if bio_feature.strand != 1:
                raise ValueError(f"feature {feature_name} of {self.name} is - "
                                 'strand, but only + strand features handled')
            feature_seq = str(bio_feature.location.extract(seqrecord).seq)
            feature = Feature(name=feature_name,
                              seq=feature_seq,
                              start=bio_feature.location.start,
                              end=bio_feature.location.end,
                              )
            self.features.append(feature)
            self.feature_names.append(feature_name)
            self._features_dict[feature_name] = feature

        missing_features = set(req_features) - set(self._features_dict)
        if missing_features:
            raise ValueError(f"{self.name} lacks features: {missing_features}")

    @classmethod
    def get_name(cls, seqrecord):
        """Get name of target from sequence record.

        Parameters
        ----------
        seqrecord : Bio.SeqRecord.SeqRecord
            Sequence record as passed to :class:`Target`.

        Returns
        -------
        str
            Name parsed from `seqrecord`.

        """
        if not hasattr(seqrecord, 'name'):
            raise ValueError('`seqrecord` does not define a name')
        else:
            return seqrecord.name

    def has_feature(self, name):
        """Check if a feature is defined for this target.

        Parameters
        ----------
        name : str
            Name of :class:`Feature`.

        Returns
        -------
        bool
            `True` if target has feature of this name, `False` otherwise.

        """
        return (name in self._features_dict)

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
        if self.has_feature(name):
            return self._features_dict[name]
        else:
            raise ValueError(f"Target {self.name} has no feature {name}")

    def image(self, *, color_map=None, feature_labels=None,
              plots_indexing='genbank'):
        """Get image of the target.

        Parameters
        ----------
        color_map : None or dict
            To specify colors for each feature, provide a dict mapping
            feature names to colors. Otherwise automatically chosen.
        feature_labels : None or dict
            Map feature names to text labels shown on plot. Otherwise
            features just labeled by name.
        plots_indexing : {'biopython', 'genbank'}
            Does image use 0-based ('biopython') or 1-based ('genbank')
            indexing of nucleotide sites?

        Returns
        -------
        dna_features_viewer.GraphicRecord.GraphicRecord
            Image of target, which has `.plot` and `.plot_with_bokeh` methods:
            https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer

        """
        if color_map is None:
            if len(self.features) < len(CBPALETTE):
                color_map = {feature.name: CBPALETTE[i + 1] for
                             i, feature in enumerate(self.features)}
            else:
                cmap = matplotlib.cm.jet
                color_map = {feature.name: matplotlib.colors.to_hex(
                                           cmap(i / len(self.features)))
                             for i, feature in enumerate(self.features)}
        else:
            missing_colors = [feature.name for feature in self.features
                              if feature.name not in color_map]
            if missing_colors:
                raise ValueError(f"no `color_map` entry for {missing_colors}")

        if feature_labels is None:
            feature_labels = {}
        for feature in self.features:
            if feature.name not in feature_labels:
                feature_labels[feature.name] = feature.name

        graph_features = []
        for feature in self.features:
            graph_features.append(
                dna_features_viewer.GraphicFeature(
                    start=feature.start,
                    end=feature.end,
                    label=feature_labels[feature.name],
                    color=color_map[feature.name],
                    strand=1,
                    )
                )

        graph_record = dna_features_viewer.GraphicRecord(
                sequence_length=self.length,
                features=graph_features,
                sequence=self.seq,
                plots_indexing=plots_indexing,
                )

        return graph_record


class Targets:
    """Collection of :class:`Target` sequences.

    Parameters
    ----------
    seqsfile : str or list
        Name of file specifying the targets, or list of such files. So
        if multiple targets they can all be in one file or in separate files.
    feature_parse_specs : dict or str
        How :meth:`Targets.parse_alignment` parses alignments. Specify dict
        or name of YAML file. Keyed by names of targets, values target-level
        dicts keyed by feature names. The feature-level dicts have two keys:

          - 'filter': dict keyed by 'clip5', 'clip3', 'mutation_nt_count',
            and 'mutation_op_count' giving max clipping at each end, number
            of nucleotide mutations, and number of ``cs`` tag mutation
            operations allowed for feature. If 'filter' itself or any of
            the keys are missing, the value is set to zero. If the value
            is `None` ('null' in YAML notation), then no filter is applied.

          - 'return': str or list of strings indicating what to return for this
            feature. If 'returns' is absent or the value is `None` ('null' in
            YAML notation), nothing is returned for this feature. Otherwise
            list one or more of 'sequence', 'mutations', 'accuracy', 'cs',
            'clip5', and 'clip3' to get the sequence, mutation string, ``cs``
            tag, or number of clipped nucleotides from each end.

        In addition, target-level dicts should have keys 'query_clip5' and
        'query_clip3' which give the max amount that can be clipped from
        each end of the query prior to the alignment. Use a value of
        `None` ('null' in YAML notation) to have no filter on this clipping.
        Filters will be applied in the order the features appear in the
        `feature_parse_specs`.

    allow_extra_features : bool
        Can targets have features not in `feature_parse_specs`?
    seqsfileformat : {'genbank'}
        Format of `seqsfile`. Currently, 'genbank' is the only supported
        option. The GenBank Flat File format is described `here
        <https://www.ncbi.nlm.nih.gov/genbank/samplerecord/>`_, but not all
        fields are required. The documentation includes `examples
        <https://jbloomlab.github.io/alignparse/examples.html>`_ that show
        what fields should typically be included. GenBank files can be readily
        generated using several sequence editing programs, such as `ApE
        <https://jorgensen.biology.utah.edu/wayned/ape/>`_ or `Benchling
        <https://www.benchling.com/>`_.
    allow_clipped_muts_seqs : bool
        Returning sequence or mutations for features where non-zero
        clipping is allowed is dangerous, since as described in
        :meth:`Targets.parse_alignment` these will only be for unclipped
        region and so are easy to mis-interpret. So you must explicitly
        set this option to `True` in order to allow return of mutations /
        sequences for features with clipping allowed; otherwise you'll get
        an error if you try to recover such sequences / mutations.
    ignore_feature_parse_specs_keys : None or list
        Ignore these target-level keys in `feature_parse_specs`. Useful for
        YAML with default keys that don't represent actual targets.
    select_target_names : None or list
        If `None`, the created object is for all sequences in `seqsfile`.
        Otherwise pass a list with names of just the sequences of interest.

    Attributes
    ----------
    targets : list
        List of all :class:`Target` objects.
    target_names : list
        List of names of all targets.
    target_seqs : dict
        Keyed by target name, value is sequence as str.

    """

    def __repr__(self):
        """Get string representation."""
        return f"{self.__class__.__name__}(targets={self.targets})"

    def __init__(self, *, seqsfile, feature_parse_specs,
                 allow_extra_features=False, seqsfileformat='genbank',
                 allow_clipped_muts_seqs=False,
                 ignore_feature_parse_specs_keys=None,
                 select_target_names=None,
                 ):
        """See main class docstring."""
        # read feature_parse_specs
        if isinstance(feature_parse_specs, str):
            with open(feature_parse_specs) as f:
                self._feature_parse_specs = yaml.safe_load(f)
        else:
            self._feature_parse_specs = copy.deepcopy(feature_parse_specs)

        if ignore_feature_parse_specs_keys:
            for key in ignore_feature_parse_specs_keys:
                if key in self._feature_parse_specs:
                    del self._feature_parse_specs[key]
                else:
                    raise KeyError(f"`feature_parse_specs` lacks key {key} "
                                   'in `ignore_feature_parse_specs_keys`')

        # names of parse alignment columns with clipping
        self._clip_cols = ['query_clip5', 'query_clip3']

        # reserved columns for parsing, cannot be name of a feature
        self._reserved_cols = ['query_name'] + self._clip_cols

        # suffixes in feature columns returned parse_alignment
        self._return_suffixes = ['_mutations', '_sequence', '_accuracy', '_cs',
                                 '_clip5', '_clip3']

        # valid filtering keys
        self._filterkeys = ['clip5', 'clip3', 'mutation_nt_count',
                            'mutation_op_count']

        # get targets from seqsfile
        if select_target_names is not None:
            if not (isinstance(select_target_names, list) and
                    len(select_target_names) >= 1):
                raise ValueError('`select_target_names` must be none or '
                                 'non-empty list')
        if isinstance(seqsfile, str):
            seqrecords = list(Bio.SeqIO.parse(seqsfile, format=seqsfileformat))
        else:
            seqrecords = []
            for f in seqsfile:
                seqrecords += list(Bio.SeqIO.parse(f, format=seqsfileformat))
        self.targets = []
        self._target_dict = {}
        for seqrecord in seqrecords:
            tname = Target.get_name(seqrecord)
            if select_target_names and (tname not in select_target_names):
                continue
            target = Target(seqrecord=seqrecord,
                            req_features=self.features_to_parse(tname, 'name'),
                            allow_extra_features=allow_extra_features,
                            )
            if target.name in self._target_dict:
                raise ValueError(f"duplicate target name of {target.name}")
            self.targets.append(target)
            self._target_dict[target.name] = target
            # ensure feature names are not reserved or have reserved suffix
            for feature in target.features:
                if feature.name in self._reserved_cols:
                    raise ValueError(f"feature cannot be named {feature.name}")
                if re.search('|'.join(s + '$' for s in self._return_suffixes),
                             feature.name):
                    raise ValueError('feature name cannot end in ' +
                                     str(self._return_suffixes))
        self.target_names = [target.name for target in self.targets]
        self.target_seqs = {target.name: target.seq for target in self.targets}
        if not self.targets:
            raise ValueError('no targets found')

        # check needed for `to_csv` option of `parse_alignment`.
        if len(self.target_names) != len({tname.replace(' ', '_') for
                                         tname in self.target_names}):
            raise ValueError('target names must be unique even after '
                             'replacing spaces with underscores.')

        # make sure we have all targets to parse
        extra_targets = set(self._feature_parse_specs) - set(self.target_names)
        if extra_targets:
            raise ValueError('`feature_parse_specs` includes non-existent '
                             f"targets {extra_targets}")

        self._set_feature_parse_specs_defaults()

        self._set_parse_filters()

        # check we are not set to return sequence / mutations for clipped
        # features unless flag to do this explicitly set
        if not allow_clipped_muts_seqs:
            for t in self.target_names:
                for f in self.features_to_parse(t, 'name'):
                    for return_name in ['sequence', 'mutations']:
                        if return_name in self._parse_returnvals(t, f):
                            filt = self._parse_filters[t][f]
                            if any(map(lambda c: c not in filt or filt[c] > 0,
                                       ['clip5', 'clip3'])):
                                raise ValueError(
                                        f"You asked to return {return_name} "
                                        f"for {t}, {f}, but clipping is not "
                                        '0 for this feature. To do this, set'
                                        '`allow_clipped_muts_seqs` to `True`')

    def _set_parse_filters(self):
        """Set `_parse_filters` attribute.

        This is dict keyed by targetname, then keyed by feature name,
        then keyed by each parameter to filter on with value being
        max allowed. If the value is `None` (no filter), it is not in dict.

        This is a simpler-to-access version of information in
        `feature_parse_specs`.

        """
        self._parse_filters = {}
        for tname in self.target_names:
            self._parse_filters[tname] = {}
            for fname in self.features_to_parse(tname, 'name'):
                self._parse_filters[tname][fname] = {}
                filterspecs = self._feature_parse_specs[tname][fname]['filter']
                for k, v in filterspecs.items():
                    if v is not None:
                        if not isinstance(v, int):
                            raise ValueError('`feature_parse_spec` filter for'
                                             f" {tname}, {fname}, {k} is not "
                                             f"`None` or an integer: {v}")
                        self._parse_filters[tname][fname][k] = v

    def _set_feature_parse_specs_defaults(self):
        """Set missing values in `feature_parse_specs` to defaults.

        These defaults are described in the main :class:`Targets` docs.

        """
        for tname, tspecs in self._feature_parse_specs.items():
            if set(self._clip_cols) - set(tspecs):
                raise ValueError(f"`feature_parse_specs` for {tname} "
                                 f"lacks {self._clip_cols}")
            for fname, fdict in tspecs.items():
                if fname in self._clip_cols:
                    continue
                if set(fdict.keys()) - {'return', 'filter'}:
                    raise ValueError(f"`feature_parse_specs` for {tname} "
                                     f"{fname} has extra keys: only 'return' "
                                     "and 'filter' are allowed.")
                if 'return' not in fdict:
                    fdict['return'] = []
                else:
                    if isinstance(fdict['return'], str):
                        fdict['return'] = [fdict['return']]
                    for returnval in fdict['return']:
                        if returnval not in [suffix[1:] for suffix in
                                             self._return_suffixes]:
                            raise ValueError(
                                    f"`feature_parse_specs` of {tname} {fname}"
                                    f" has invalid return type {returnval}")
                if 'filter' not in fdict:
                    fdict['filter'] = {}
                if set(fdict['filter'].keys()) - set(self._filterkeys):
                    raise ValueError(f"`feature_parse_specs` for {tname} "
                                     f"{fname} has invalid filter type. Only "
                                     f"{self._filterkeys} are allowed.")
                for filterkey in self._filterkeys:
                    if filterkey not in fdict['filter']:
                        fdict['filter'][filterkey] = 0

    def features_to_parse(self, targetname, feature_or_name='feature'):
        """Features to parse for a target.

        Parameters
        ----------
        targetname : str
            Name of target.
        feature_or_name : {'feature', 'name'}
            Get the :class:`Feature` objects themselves or their names.

        Returns
        -------
        list
            Features to parse for this target, as specified in
            :meth:`Targets.feature_parse_specs`.

        """
        if feature_or_name == 'name':
            if not hasattr(self, '_fnames_to_parse'):
                self._fnames_to_parse = {}
                for tname, tdict in self._feature_parse_specs.items():
                    self._fnames_to_parse[tname] = [f for f in tdict if
                                                    f not in self._clip_cols]
            try:
                return self._fnames_to_parse[targetname]
            except KeyError:
                raise ValueError(f"target {targetname} not in "
                                 '`feature_parse_specs`')
        elif feature_or_name == 'feature':
            if not hasattr(self, '_features_to_parse'):
                self._features_to_parse = {}
                for tname, tdict in self._feature_parse_specs.items():
                    target = self.get_target(tname)
                    self._features_to_parse[tname] = [target.get_feature(f)
                                                      for f in tdict if
                                                      f not in self._clip_cols]
            try:
                return self._features_to_parse[targetname]
            except KeyError:
                raise ValueError(f"target {targetname} not in "
                                 '`feature_parse_specs`')
        else:
            raise ValueError(f"invalid `feature_or_name` {feature_or_name}")

    def feature_parse_specs(self, returntype):
        """Get the feature parsing specs.

        Note
        ----
        Filters will be applied in the order they are listed in the
        `feature_parse_specs` `yaml` file or `dict`. Once a read fails a
        filter, other filters will not be applied. As such, it is recommended
        to have features with filters for 5' and 3' clipping listed first.

        Parameters
        ----------
        returntype : {'dict', 'yaml'}
            Return a Python `dict` or a YAML string representation.

        Returns
        -------
        dict or str
            The feature parsing specs set by the `feature_parse_specs` at
            :class:`Targets` initialization, but with any missing default
            values explicitly filled in.

        """
        if returntype == 'dict':
            return copy.deepcopy(self._feature_parse_specs)
        elif returntype == 'yaml':
            return yaml.dump(self._feature_parse_specs,
                             default_flow_style=False,
                             sort_keys=False)
        else:
            raise ValueError(f"invalid `returntype` of {returntype}")

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
        if name in self._target_dict:
            return self._target_dict[name]
        else:
            raise ValueError(f"no target named {name}")

    def write_fasta(self, fastafile):
        """Write all targets to a FASTA file.

        Parameters
        ----------
        filename : str or writable file-like object
            Write targets to this file.

        """
        try:
            for target in self.targets:
                fastafile.write(f">{target.name}\n{target.seq}\n")
            fastafile.flush()
        except AttributeError:
            with open(fastafile, 'w') as f:
                for target in self.targets:
                    f.write(f">{target.name}\n{target.seq}\n")

    def plot(self,
             *,
             sharex=True,
             ax_width=5,
             ax_height=3,
             hspace=0.4,
             **kwargs):
        """Plot all the targets.

        Note
        ----
        For more customizable plots, call :meth:`Target.image` for individual
        targets.

        Parameters
        ----------
        sharex : bool
            Share x-axis among plots for each target?
        ax_width : float
            Width of each axis in inches.
        ax_height : float
            Height of each axis in inches.
        hspace : float
            Vertical space between axes as fraction of `ax_height`.
        ``**kwargs``
            Keyword arguments passed to :meth:`Target.image`.

        Returns
        -------
        matplotlib.pyplot.figure
            Figure showing all targets.

        """
        fig, axes = plt.subplots(nrows=len(self.targets),
                                 ncols=1,
                                 sharex=False,
                                 squeeze=False,
                                 gridspec_kw={'hspace': hspace},
                                 figsize=(ax_width,
                                          len(self.targets) * ax_height),
                                 )
        for ax, target in zip(axes.ravel(), self.targets):
            image = target.image(**kwargs)
            image.plot(ax=ax)
            ax.set_title(target.name)

        return fig

    def align(self, queryfile, alignmentfile, mapper):
        """Align query sequences to targets.

        Parameters
        ----------
        queryfile : str
            The query sequences to align (FASTQ or FASTA, can be gzipped).
        alignmentfile : str
            SAM file created by `mapper` with alignments of queries to the
            target sequences within this :class:`Targets` object.
        mapper : :class:`alignparse.minimap2.Mapper`
            Mapper that runs ``minimap2``. Alignment options set when creating
            this mapper.

        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa') as targetfile:
            self.write_fasta(targetfile)
            mapper.map_to_sam(targetfile.name, queryfile, alignmentfile)

    def align_and_parse(self,
                        df,
                        mapper,
                        outdir,
                        *,
                        name_col='name',
                        queryfile_col='queryfile',
                        group_cols=None,
                        to_csv=False,
                        overwrite=False,
                        multi_align='primary',
                        filtered_cs=False,
                        skip_sups=True,
                        ncpus=-1
                        ):
        """Align query sequences and then parse alignments.

        Note
        ----
        This is a convenience method to run :meth:`Targets.align` and
        :meth:`Targets.parse_alignment` on multiple queries and collate
        the results.

        It also allows multiple queries to be handled simultaneously using
        multiprocessing.

        Parameters
        ----------
        df : pandas.DataFrame
            Data frame with information on queries to align.
        mapper : :class:`alignparse.minimap2.Mapper`
            Mapper that runs ``minimap2``. Alignment options set when creating
            this mapper.
        outdir : str
            Name of directory with created alignments and parsing files.
            Created if it does not exist.
        name_col : str
            Column in `df` with the name of each set of queries.
        queryfile_col :str
            Column in `df` with FASTQ file with queries.
        group_cols : `None` or str or list
            Columns in `df` used to "group" results. These columns are
            in all created data frames. For instance, might specify
            different libraries or samples.
        to_csv : bool
            Write CSV files rather than return data frames. Useful to
            avoid reading large data frames into memory.
        overwrite : bool
            If some of the created output files already exist, do we
            overwrite them or raise and error?
        multi_align : {'primary'}
            How to handle multiple alignments. Currently only option is
            'primary', which ignores all secondary alignments.
        filtered_cs : bool
            Add `cs` tag that failed the filter to filtered dataframe along
            with filter reason. Allows for more easily investigating why
            reads are failing the filters.
        skip_sups : bool
            Whether or not to skip supplementary alignments when parsing.
            Supplementary alignments are additional possible alignments for
            a read due to the read potentially being a chimeric. The default
            is to skip these alignments and *not* parse them.
        ncpus : int
            Number of CPUs to use; -1 means all available.

        Returns
        -------
        (readstats, aligned, filtered) : tuple
            Same meaning as for :meth:`Targets.parse_alignment` except
            the data frames / CSV files all have additional columns indicating
            name of each query set (`name_cols`) as well as any `group_cols`.

        """
        # check columns in `df`
        if not group_cols:
            addtl_cols = [name_col]
        else:
            if isinstance(group_cols, str):
                group_cols = [group_cols]
            addtl_cols = group_cols + [name_col]
            if len(set(addtl_cols)) != len(addtl_cols):
                raise ValueError('`name_col` and `group_cols` have redundant '
                                 f"entries: {addtl_cols}")
        expected_cols = addtl_cols + [queryfile_col]
        if not set(df.columns).issuperset(set(expected_cols)):
            raise ValueError('`df` does not contain all expected columns: ' +
                             str(expected_cols))

        os.makedirs(outdir, exist_ok=True)

        # For each query we create a "full name" that includes name +
        # grouping cols, a subdirectory that holds all files
        # for this query, and the samfile for the alignments.
        reserved_cols = {'fullname', 'subdir', 'samfile'}
        if set(expected_cols).intersection(reserved_cols):
            raise ValueError(f"`df` cannot have columns: {reserved_cols}")
        df = (
            df
            [expected_cols]
            .astype(str)
            .assign(
                fullname=lambda x: x.apply(lambda r: '_'.join(r[c] for c
                                                              in addtl_cols),
                                           axis=1),
                subdir=lambda x: x['fullname'].map(lambda n: os.path.join(
                                                   outdir, n)),
                samfile=lambda x: x['subdir'].map(lambda d: os.path.join(
                                                  d, 'alignments.sam')),
                )
            .reset_index(drop=True)
            )
        if len(df) != df['fullname'].nunique():
            raise ValueError('Names the queries are not unique even after '
                             'adding grouping columns.')
        if any(df[col].str.contains(',').any() for col in df.columns):
            raise ValueError('`name_col` and `group_cols` entry contains ","')

        # set up multiprocessing pool
        if ncpus == -1:
            ncpus = pathos.multiprocessing.cpu_count()
        else:
            ncpus = min(pathos.multiprocessing.cpu_count(), ncpus)
        if ncpus < 1:
            raise ValueError('`ncpus` must be >= 1')
        if ncpus > 1:
            pool = pathos.pools.ProcessPool(ncpus)
            map_func = pool.map
        else:
            def map_func(f, *args):
                return [f(*argtup) for argtup in zip(*args)]

        # now make the alignments
        for tup in df.itertuples():
            os.makedirs(tup.subdir, exist_ok=True)
            if os.path.isfile(tup.samfile):
                if overwrite:
                    os.remove(tup.samfile)
                else:
                    raise IOError(f"file {tup.samfile} already exists")
        _ = map_func(self.align,
                     df[queryfile_col],
                     df['samfile'],
                     itertools.repeat(mapper))
        assert all(os.path.isfile(f) for f in df['samfile'])

        # now parse the alignments
        parse_results = map_func(self.parse_alignment,
                                 df['samfile'],
                                 itertools.repeat(multi_align),
                                 itertools.repeat(True),
                                 df['subdir'],
                                 itertools.repeat(overwrite),
                                 itertools.repeat(filtered_cs),
                                 itertools.repeat(skip_sups)
                                 )

        # close, clear pool: https://github.com/uqfoundation/pathos/issues/111
        if ncpus > 1:
            pool.close()
            pool.join()
            pool.clear()

        # Set up to gather overall readstats, aligned, and filtered by
        # getting and checking column names:
        readstatcols = ['category', 'count']
        alignedcols = {t: self._parse_returnvals(t) for t in self.target_names}
        filteredcols = ['query_name', 'filter_reason']
        disallowedcols = readstatcols + filteredcols
        for tc in alignedcols.values():
            disallowedcols += tc
        if set(addtl_cols).intersection(set(disallowedcols)):
            raise ValueError('`name_col`, `group_cols` cannot have any of ' +
                             str(disallowedcols))
        alignedcols = {t: addtl_cols + tc for t, tc in alignedcols.items()}
        filteredcols = addtl_cols + filteredcols

        # set up data frames or names of files
        readstats = pd.DataFrame([], columns=addtl_cols + readstatcols)
        filtered = {t: os.path.join(outdir, t.replace(' ', '_') +
                                    '_filtered.csv')
                    for t in self.target_names}
        aligned = {t: os.path.join(outdir, t.replace(' ', '_') +
                                   '_aligned.csv')
                   for t in self.target_names}
        for f in list(filtered.values()) + list(aligned.values()):
            if os.path.isfile(f):
                if not overwrite:
                    raise IOError(f"file {f} already exists.")
                else:
                    os.remove(f)

        # collect read stats for all runs
        for i, (ireadstats, _, _) in enumerate(parse_results):
            readstats = readstats.append(
                    (ireadstats
                     .assign(**{c: df.at[i, c] for c in addtl_cols})
                     [readstats.columns]
                     ),
                    ignore_index=True,
                    sort=False,
                    )
        readstats = readstats.assign(count=lambda x: x['count'].astype(int))

        # collect aligned and filtered for all runs
        for t in self.target_names:
            with contextlib.ExitStack() as stack:
                # Define callback to delete CSV files on error. See here:
                # https://docs.python.org/3/library/contextlib.html#replacing-any-use-of-try-finally-and-flag-variables
                @stack.callback
                def delete_files_on_err():
                    for fname in [aligned[t], filtered[t]]:
                        if os.path.isfile(fname):
                            os.remove(fname)

                alignedfile = stack.enter_context(open(aligned[t], 'w'))
                filteredfile = stack.enter_context(open(filtered[t], 'w'))
                alignedfile.write(','.join(alignedcols[t]) + '\n')
                filteredfile.write(','.join(filteredcols) + '\n')

                for i, (_, ialigned, ifiltered) in enumerate(parse_results):
                    addtl_text = ','.join(df.at[i, c] for c in addtl_cols)
                    for fin_name, fout, cols in [
                            (ialigned[t], alignedfile, alignedcols[t]),
                            (ifiltered[t], filteredfile, filteredcols),
                            ]:
                        with open(fin_name) as fin:
                            firstline = fin.readline()
                            assert (firstline[: -1].split(',') ==
                                    cols[len(addtl_cols):]), fin_name
                            for line in fin:
                                fout.write(addtl_text)
                                fout.write(',')
                                fout.write(line)
                            fout.flush()

                stack.pop_all()  # no callback to delete files if reached here

        if not to_csv:
            aligned = {t: pd.read_csv(f).fillna('')
                       for t, f in aligned.items()}
            filtered = {t: pd.read_csv(f) for t, f in filtered.items()}

        return readstats, aligned, filtered

    def parse_alignment(self, samfile, multi_align='primary',
                        to_csv=False, csv_dir=None, overwrite_csv=False,
                        filtered_cs=False, skip_sups=True):
        """Parse alignment features as specified in `feature_parse_specs`.

        Parameters
        ----------
        samfile : str
            SAM file with ``minimap2`` alignments with ``cs`` tag, typically
            created by :meth:`Targets.align`.
        multi_align : {'primary'}
            How to handle multiple alignments. Currently only option is
            'primary', which ignores all secondary alignments.
        to_csv : bool
            Return CSV file names rather than return data frames. Useful to
            avoid reading large data frames into memory.
        csv_dir : None or str
            If `to_csv` is `True`, name of directory to which we
            write CSV files (created if needed). If `None`, write
            to current directory.
        overwrite_csv : bool
            If using `to_csv`, do we overwrite existing CSV files or
            raise an error if they already exist?
        filtered_cs : bool
            Add `cs` tag that failed the filter to filtered dataframe along
            with filter reason. Allows for more easily investigating why
            reads are failing the filters.
        skip_sups : bool
            Whether or not to skip supplementary alignments when parsing.
            Supplementary alignments are additional possible alignments for
            a read due to the read potentially being a chimeric. The default
            is to skip these alignments and *not* parse them.

        Returns
        -------
        (readstats, aligned, filtered) : tuple

            - `readstats` is `pandas.DataFrame` with numbers of unmapped reads,
              and for each target the number of mapped reads that are validly
              aligned and that fail filters in `feature_parse_specs`.

            - `aligned` is a dict keyed by name of each target. Entries are
              `pandas.DataFrame` with rows for each validly aligned read. Rows
              give query name, query clipping at each end of alignment, and any
              feature-level info specified for return in `feature_parse_specs`
              in columns with names equal to feature suffixed by '_sequence',
              '_mutations', '_accuracy', '_cs', '_clip5', and '_clip3'.
              and '_clip3'.

            - `filtered` is a dict keyed by name of each target. Entries are
              `pandas.DataFrame` with a row for each filtered aligned read
              giving the query name and the reason it was filtered.
              If `filtered_cs` is `True` then, add a column to the "filtered"
              `pandas.DataFrame`s with the `cs` tag that failed the filter.

            If `to_csv` is `True`, then `aligned` and `filtered` give
            names of CSV files holding data frames.

        Note
        ----
        The ``cs`` tags are in the short format returned by ``minimap2``;
        see here for details: https://lh3.github.io/minimap2/minimap2.html

        When parsing features, if an insertion occurs between two features,
        it is assigned to the end of the first feature.

        Returned sequences, mutation strings, and ``cs`` tags are only for
        for the portion of the feature that aligns, and do **not** indicate
        clipping, which you instead get in the '_clip*' columns. The
        sequences are simply what the ``cs`` tag implies, indels / mutations
        are not indicated in this column. Mutation strings are space-delimited
        with these operations in **1-based** (1, 2, ...) numbering from start
        of the feature:

          - 'A2G' : substitution at site 2 from A to G

          - 'ins5TAA' : insertion of 'TAA' starting at site 5

          - 'del5to6' : deletion of sites 5 to 6, inclusive

        The returned accuracy is the average accuracy of the aligned
        query sites as calculated from the Q-values, and is `nan` if
        there are no aligned query sites.

        """
        if multi_align == 'primary':
            primary_only = True
        else:
            raise ValueError(f"invalid `multi_align` {multi_align}")

        if csv_dir is None:
            csv_dir = ''
        elif csv_dir:
            os.makedirs(csv_dir, exist_ok=True)

        unmapped = 0
        readstats = {t: {'filtered': 0, 'aligned': 0}
                     for t in self.target_names}

        if filtered_cs:
            filtered_cols = ['query_name', 'filter_reason', 'filter_cs']
        else:
            filtered_cols = ['query_name', 'filter_reason']

        if to_csv:
            filtered = {t: os.path.join(csv_dir, t.replace(' ', '_') +
                                        '_filtered.csv')
                        for t in self.target_names}
            aligned = {t: os.path.join(csv_dir, t.replace(' ', '_') +
                                       '_aligned.csv')
                       for t in self.target_names}
            filenames = list(filtered.values()) + list(aligned.values())
            if (not overwrite_csv) and any(map(os.path.isfile, filenames)):
                raise IOError(f"existing file with name in: {filenames}")
        else:
            filtered = {t: [] for t in self.target_names}
            aligned = {t: [] for t in self.target_names}

        # Iterate samfile in contextlib so can handle files if using `to_csv`.
        with contextlib.ExitStack() as stack:
            # Define callback to delete CSV files on error. See here:
            # https://docs.python.org/3/library/contextlib.html#replacing-any-use-of-try-finally-and-flag-variables
            @stack.callback
            def delete_files_on_err():
                if to_csv:
                    for fname in filenames:
                        if os.path.isfile(fname):
                            os.remove(fname)

            if to_csv:
                # add files to stack so they are closed at end:
                # https://stackoverflow.com/a/19412700
                filtered_files = {t: stack.enter_context(open(f, 'w'))
                                  for t, f in filtered.items()}
                for f in filtered_files.values():
                    f.write(','.join(filtered_cols) + '\n')
                aligned_files = {t: stack.enter_context(open(f, 'w'))
                                 for t, f in aligned.items()}
                for t, f in aligned_files.items():
                    f.write(','.join(self._parse_returnvals(t)) + '\n')

            # iterate over samfile and process each alignment
            for aligned_seg in pysam.AlignmentFile(samfile):

                if aligned_seg.is_unmapped:
                    unmapped += 1
                    continue

                if aligned_seg.is_supplementary and skip_sups:
                    continue

                if aligned_seg.is_secondary and primary_only:
                    continue

                a = Alignment(aligned_seg, introns_to_deletions=True,
                              target_seqs=self.target_seqs)
                tname = a.target_name

                is_filtered, parse_tup = self._parse_single_Alignment(
                                                                a,
                                                                tname,
                                                                filtered_cs)

                if is_filtered:
                    readstats[tname]['filtered'] += 1
                    if to_csv:
                        filtered_files[tname].write(','.join(map(str,
                                                                 parse_tup)))
                        filtered_files[tname].write('\n')
                    else:
                        filtered[tname].append(parse_tup)

                else:
                    readstats[tname]['aligned'] += 1
                    if to_csv:
                        aligned_files[tname].write(','.join(map(str,
                                                                parse_tup)))
                        aligned_files[tname].write('\n')
                    else:
                        aligned[tname].append(parse_tup)

            # Done iterating over `samfile`, get values ready to return
            readstats = (pd.DataFrame(readstats)
                         .reset_index()
                         .melt(id_vars='index', value_name='count')
                         .assign(category=lambda x: (x['index'] + ' ' +
                                                     x['variable']))
                         [['category', 'count']]
                         .append({'category': 'unmapped', 'count': unmapped},
                                 ignore_index=True)
                         )
            if not to_csv:
                filtered = {t: pd.DataFrame(tlist, columns=filtered_cols)
                            for t, tlist in filtered.items()}
                aligned = {t: pd.DataFrame(tlist,
                                           columns=self._parse_returnvals(t))
                           for t, tlist in aligned.items()}
            else:
                for f in (list(aligned_files.values()) +
                          list(filtered_files.values())):
                    f.flush()

            stack.pop_all()  # no callback to delete CSV files if reached here

        return readstats, aligned, filtered

    def _parse_returnvals(self, targetname, featurename=None):
        """Values to return when parsing alignment.

        Parameters
        ----------
        targetname : str
        featurename : str or None

        Returns
        -------
        list
            If `featurename` is not `None`, gets list of all values
            to return for this feature and targets from `feature_parse_specs`.
            If `featurename` is `None`, gets list of all values for all
            features in target, prefixing with the feature name.

        """
        if featurename is None:
            returnlist = ['query_name', 'query_clip5', 'query_clip3']
            for featurename in self.features_to_parse(targetname, 'name'):
                for val in self._parse_returnvals(targetname, featurename):
                    returnlist.append(f"{featurename}_{val}")
            return returnlist
        else:
            return self._feature_parse_specs[targetname][featurename]['return']

    def _parse_single_Alignment(self, a, targetname, filtered_cs=False):
        """Parse a single alignment.

        Parameters
        ----------
        a : :class:`alignparse.cs_tag.Alignment`
        targetname : str
        filtered_cs : bool

        Returns
        --------
        2-tuple
            Tuple `(is_filtered, parse_tup)` where `is_filtered` is bool
            indicating if alignment fails `feature_parse_specs` filters.
            If it fails and `filtered_cs` is `False`, `parse_tup` is 2-tuple
            `(query_name, filter_reason)`. If it fails and `filtered_cs` is
            `True`, `parse_tup` is tuple `(query_name, filter_reason,
            filtered_cs)`. If it passes, `parse_tup` is list giving values
            specified by :meth:`_parse_returnvals` for `targetname`.

        """
        query_name = a.query_name
        parse_tup = [query_name]

        # get and filter on query clipping
        for clip in ['query_clip5', 'query_clip3']:
            clipval = getattr(a, clip)
            clipmax = self._feature_parse_specs[targetname][clip]
            if (clipmax is None) or clipval <= clipmax:
                parse_tup.append(clipval)
            else:
                return True, (query_name, clip)

        # get and filter on features
        target_parse_filters = self._parse_filters[targetname]
        for feature in self.features_to_parse(targetname):
            feat_info = a.extract_cs(feature.start, feature.end)
            if feat_info is None:
                # feature is fully clipped, determine which end
                if a.target_clip5 >= feature.end:
                    clip5 = feature.length
                    clip3 = 0
                    cs = ''
                elif a.target_lastpos <= feature.start:
                    clip5 = 0
                    clip3 = feature.length
                    cs = ''
                else:
                    raise ValueError(
                            f"Should not get here:\ntarget = {targetname}:\n"
                            f"feature = {feature}\n"
                            f"target_clip5 = {a.target_clip5}\n"
                            f"lastpos = {a.target_lastpos}\n"
                            )
            else:
                cs, clip5, clip3 = feat_info

            # apply filters
            featurename = feature.name
            filter_vals = {'clip5': clip5,
                           'clip3': clip3,
                           'mutation_nt_count': cs_to_nt_mutation_count(cs),
                           'mutation_op_count': cs_to_op_mutation_count(cs),
                           }
            for key, valmax in target_parse_filters[featurename].items():
                if filter_vals[key] > valmax:
                    if filtered_cs:
                        return True, (query_name, f"{featurename} {key}", cs)
                    else:
                        return True, (query_name, f"{featurename} {key}")

            # get return values
            for return_name in self._parse_returnvals(targetname, featurename):
                if return_name == 'clip5':
                    parse_tup.append(clip5)
                elif return_name == 'clip3':
                    parse_tup.append(clip3)
                elif return_name == 'mutations':
                    parse_tup.append(cs_to_mutation_str(cs, clip5))
                elif return_name == 'cs':
                    parse_tup.append(cs)
                elif return_name == 'sequence':
                    clippedseq = feature.seq[clip5: feature.length - clip3]
                    parse_tup.append(cs_to_sequence(cs, clippedseq))
                elif return_name == 'accuracy':
                    parse_tup.append(a.get_accuracy(feature.start,
                                                    feature.end))
                else:
                    allowednames = [name[1:] for name in self._return_suffixes]
                    raise ValueError(f"invalid `return_name` {return_name}, "
                                     f"should be one of {allowednames}")

        # if here, alignment not filtered: return it
        return False, parse_tup

    def _parse_alignment_cs(self, samfile, *, multi_align='primary',
                            skip_sups=True):
        """Parse alignment feature ``cs`` strings for aligned queries.

        Note
        ----
        This method returns the same information that can be better
        obtained via :meth:`Targets.parse_alignment` by setting
        to return 'cs', 'clip5', 'clip3' for every feature. It is
        currently retained only for debugging / testing purposes,
        and may eventually be removed.

        Parameters
        ----------
        See parameter definitions in :meth:`Targets.parse_alignment`.

        Returns
        -------
        dict
            Keyed be each target name. If there are no alignments for
            that target, value is `None`. Otherwise, value is a pandas
            DataFrame with a row for each feature in each query. The
            columns in this data frame are:

                - 'query_name' : name of query in `samfile`

                - 'query_clip5' : length at 5' end of query not in alignment

                - 'query_clip3' : length at 3' end of query not in alignment

                - for each feature listed for that target in
                  :meth:`Target.feature_parse_specs`, columns with name of the
                  feature and the following suffixes:

                    - '_cs' : the ``cs`` string for the aligned portion
                      of the target.

                    - '_clip5' : number of nucleotides clipped from 5' end
                      of feature in alignment.

                    - '_clip3' : number of nucleotides clipped from 3' end
                      of feature in alignment.

                  If the feature is not aligned at all, then the '_cs' suffix
                  column is an empty str, and either '_clip5' or '_clip3' will
                  be the whole length of feature depending on if feature is
                  upstream or downstream of aligned region.

            The returned dict also has a key 'unmapped' which gives
            the number of unmapped reads.

        """
        suffixes = self._return_suffixes[3:]
        d = {target: None for target in self.target_names}
        if 'unmapped' in self.target_names:
            raise ValueError('cannot have a target named "unmapped"')
        else:
            d['unmapped'] = 0

        if multi_align == 'primary':
            primary_only = True
        else:
            raise ValueError(f"invalid `multi_align` {multi_align}")

        for aligned_seg in pysam.AlignmentFile(samfile):
            if aligned_seg.is_unmapped:
                d['unmapped'] += 1
            else:
                if primary_only and aligned_seg.is_secondary:
                    continue
                if skip_sups and aligned_seg.is_supplementary:
                    continue
                a = Alignment(aligned_seg, introns_to_deletions=True,
                              target_seqs=self.target_seqs)
                tname = a.target_name
                if d[tname] is None:
                    d[tname] = {col: [] for col in self._reserved_cols}
                    for feature in self.features_to_parse(tname):
                        fname = feature.name
                        for suffix in suffixes:
                            d[tname][fname + suffix] = []

                d[tname]['query_name'].append(a.query_name)
                d[tname]['query_clip5'].append(a.query_clip5)
                d[tname]['query_clip3'].append(a.query_clip3)

                for feature in self.features_to_parse(tname):
                    feat_info = a.extract_cs(feature.start, feature.end)
                    fname = feature.name
                    if feat_info is None:
                        d[tname][fname + '_cs'].append('')
                        if a.target_clip5 >= feature.end:
                            d[tname][fname + '_clip5'].append(feature.length)
                            d[tname][fname + '_clip3'].append(0)
                        elif a.target_lastpos <= feature.start:
                            d[tname][fname + '_clip5'].append(0)
                            d[tname][fname + '_clip3'].append(feature.length)
                        else:
                            raise ValueError(
                                f"Should never get here for target {tname}:\n"
                                f"feature = {feature}\n"
                                f"target_clip5 = {a.target_clip5}\n"
                                f"lastpos = {a.target_lastpos}\n"
                                )
                    else:
                        d[tname][fname + '_cs'].append(feat_info[0])
                        d[tname][fname + '_clip5'].append(feat_info[1])
                        d[tname][fname + '_clip3'].append(feat_info[2])

        for target in d.keys():
            if target != 'unmapped':
                if d[target] is not None:
                    d[target] = pd.DataFrame.from_dict(d[target])

        return d


if __name__ == '__main__':
    import doctest
    doctest.testmod()
