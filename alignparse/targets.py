"""
=======
targets
=======

Defines :class:`Targets`, which holds :class:`Target` objects that define
alignment targets. Each :class:`Target` has some :class:`Feature` regions.

"""


import copy
import re
import tempfile

import Bio.SeqIO

import dna_features_viewer

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt

import pandas as pd

import pysam

import yaml

from alignparse.constants import CBPALETTE
from alignparse.cs_tag import Alignment


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
        if not hasattr(seqrecord, 'seq'):
            raise ValueError(f"`seqrecord` does not define a seq")
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
        Name parsed from `seqrecord`.

        """
        if not hasattr(seqrecord, 'name'):
            raise ValueError(f"`seqrecord` does not define a name")
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
        Indicates how :meth:`Targets.parse_alignment` parses the alignments.
        Should specify a dict or YAML file giving a dict. Keyed by name of
        of each target, values are target-level dicts keyed by feature names.
        These feature-level dicts can have two keys:

          - 'filter': a dict keyed by any of 'clip_count', 'mutation_nt_count',
            and 'mutation_op_count' which give the maximal clipping, number
            of nucleotide mutations, and number of ``cs`` tag mutation
            operations allowed for the feature. If 'filter' itself or any of
            the keys are missing, the value is set to zero. If the value
            is `None` ('null' in YAML notation), then no filter is applied.

          - 'return': a str or list of strings indicating what we return
            for this feature. Can be 'sequence' or 'mutations'. If 'returns'
            is absent or the value is `None` ('null' in YAML notation),
            nothing is returned for this feature.

        In addition, target-level dicts should have keys 'query_clip5' and
        'query_clip3' which give the max amount that can be clipped from
        each end of the query prior to the alignment. Use a value of
        `None` ('null' in YAML notation) to have no filter on this clipping.

    allow_extra_features : bool
        Can targets have features not in `feature_parse_specs`?
    seqsfileformat : {'genbank'}
        Format of `seqsfile`.

    Attributes
    ----------
    targets : list
        List of all :class:`Target` objects.
    target_names : list
        List of names of all targets.

    """

    def __repr__(self):
        """Get string representation."""
        return f"{self.__class__.__name__}(targets={self.targets})"

    def __init__(self, *, seqsfile, feature_parse_specs,
                 allow_extra_features=False, seqsfileformat='genbank'):
        """See main class docstring."""
        if isinstance(seqsfile, str):
            seqsfile = [seqsfile]

        seqrecords = []
        for f in seqsfile:
            seqrecords += list(Bio.SeqIO.parse(f, format=seqsfileformat))

        # name of columns with alignment clipping
        self._clip_cols = ['query_clip5', 'query_clip3']

        # read feature_parse_specs
        if isinstance(feature_parse_specs, str):
            with open(feature_parse_specs) as f:
                self._feature_parse_specs = yaml.safe_load(f)
        else:
            self._feature_parse_specs = copy.deepcopy(feature_parse_specs)

        # reserved columns for parsing, cannot be name of a feature
        self._reserved_cols = ['query_name'] + self._clip_cols

        # reserved suffixes, cannot be end of target name
        self._parse_alignment_cs_suffixes = ['_cs', '_clip5', '_clip3']
        self._parse_alignment_suffixes = ['_mutations', '_sequence']

        self.targets = []
        self.target_names = []
        self._target_dict = {}
        for seqrecord in seqrecords:
            targetname = Target.get_name(seqrecord)
            if targetname not in self._feature_parse_specs:
                raise ValueError(f"target {targetname} not in "
                                 '`feature_parse_specs`')
            target = Target(seqrecord=seqrecord,
                            req_features=(set(self._feature_parse_specs
                                              [targetname].keys()) -
                                          set(self._clip_cols)),
                            allow_extra_features=allow_extra_features,
                            )
            assert target.name == targetname
            if target.name in self._target_dict:
                raise ValueError(f"duplicate target name of {target.name}")
            self.target_names.append(target.name)
            self.targets.append(target)
            self._target_dict[target.name] = target
            # we cannot have feature names that match reserved names
            for feature in target.features:
                if feature.name in self._reserved_cols:
                    raise ValueError(f"cannot have a feature {feature.name}")
                if re.search('|'.join(s + '$' for s in
                                      self._parse_alignment_cs_suffixes),
                             feature.name):
                    raise ValueError(
                            'feature name cannot end in :\n' +
                            '\n'.join(self._parse_alignment_cs_suffixes))
                if re.search('|'.join(s + '$' for s in
                                      self._parse_alignment_suffixes),
                             feature.name):
                    raise ValueError('feature name cannot end in :\n' +
                                     '\n'.join(self._parse_alignment_suffixes))

        extra_targets = set(self._feature_parse_specs) - set(self.target_names)
        if extra_targets:
            raise ValueError('`feature_parse_specs` includes non-existent '
                             f"targets {extra_targets}")

        # check and set defaults in feature_parse_specs
        self._features_to_parse = {}  # features to parse for each target
        for targetname, targetspecs in self._feature_parse_specs.items():
            target = self.get_target(targetname)
            if targetname not in self.target_names:
                raise ValueError('`feature_parse_specs` includes non-existent '
                                 f"target {targetname}")
            extrafeatures = (set(targetspecs) -
                             set(self._clip_cols)
                             .union(target.feature_names)
                             )
            if extrafeatures:
                raise ValueError(f"`feature_parse_specs` for {targetname} has "
                                 f"specs for unknown features {extrafeatures}")
            if set(self._clip_cols) - set(targetspecs):
                raise ValueError(f"`feature_parse_specs` for {targetname} "
                                 f"lacks {self._clip_cols}")
            self._features_to_parse[targetname] = []
            for fname, fdict in targetspecs.items():
                if fname in self._clip_cols:
                    continue
                feature = target.get_feature(fname)
                self._features_to_parse[targetname].append(feature)
                if set(fdict.keys()) - {'return', 'filter'}:
                    raise ValueError(f"`feature_parse_specs` for {targetname} "
                                     f"{fname} has extra keys: only 'return' "
                                     "and 'filter' are allowed.")
                if 'return' not in fdict:
                    fdict['return'] = []
                else:
                    if isinstance(fdict['return'], str):
                        fdict['return'] = [fdict['return']]
                    for returnval in fdict['return']:
                        if returnval not in {'sequence', 'mutations'}:
                            raise ValueError(
                                    f"`feature_parse_specs` for {targetname} "
                                    f"{fname} has invalid return type "
                                    f"{returnval}")
                if 'filter' not in fdict:
                    fdict['filter'] = {}
                filterkeys = ['clip_count', 'mutation_nt_count',
                              'mutation_op_count']
                if set(fdict['filter'].keys()) - set(filterkeys):
                    raise ValueError(f"`feature_parse_specs` for {targetname} "
                                     f"{fname} has invalid filter type. Only "
                                     f"{filterkeys} are allowed.")
                for filterkey in filterkeys:
                    if filterkey not in fdict['filter']:
                        fdict['filter'][filterkey] = 0

    def feature_parse_specs(self, returntype):
        """Get the feature parsing specs.

        Parameter
        ---------
        returntype : {'dict', 'yaml'}
            Return a Python `dict` or a YAML string representation.

        Returns
        -------
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

    def plot(self, *, sharex=True, ax_width=5, ax_height=3, **kwargs):
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
        ``**kwargs``
            Keyword arguments passed to :meth:`Target.image`.

        Returns
        -------
        matplotlib.pyplot.figure
            Figure showing all targets.

        """
        fig, axes = plt.subplots(nrows=len(self.targets),
                                 ncols=1,
                                 sharex=True,
                                 squeeze=False,
                                 gridspec_kw={'hspace': 0.3},
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

    def parse_alignment(self, samfile, *, multi_align='primary'):
        """**Docs in progress.**"""
        raise RuntimeError('not yet implemented')

    def parse_alignment_cs(self, samfile, *, multi_align='primary'):
        """Parse alignment feature ``cs`` strings for aligned queries.

        Note
        ----
        The ``cs`` tags are in the short format returned by ``minimap2``;
        see here for details: https://lh3.github.io/minimap2/minimap2.html

        When an insertion occurs between two features, it is assigned to the
        end of the first feature.

        Parameters
        ----------
        samfile : str
            SAM file with ``minimap2`` alignments with ``cs`` tag, typically
            created by :meth:`Targets.align`.
        multi_align : {'primary'}
            How to handle multiple alignments. Currently only option is
            'primary', which indicates that we only retain primary alignments
            and ignore all secondary alignment.

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
        d = {target: None for target in self.target_names}
        if 'unmapped' in self.target_names:
            raise ValueError('cannot have a target named "unmapped"')
        else:
            d['unmapped'] = 0

        for a in pysam.AlignmentFile(samfile):
            if a.is_unmapped:
                d['unmapped'] += 1
            else:
                aligned_seg = Alignment(a)
                tname = aligned_seg.target_name
                features = self._features_to_parse[tname]
                if d[tname] is None:
                    d[tname] = {col: [] for col in self._reserved_cols}
                    for feature in features:
                        fname = feature.name
                        for suffix in self._parse_alignment_cs_suffixes:
                            d[tname][fname + suffix] = []

                d[tname]['query_name'].append(aligned_seg.query_name)
                d[tname]['query_clip5'].append(aligned_seg.query_clip5)
                d[tname]['query_clip3'].append(aligned_seg.query_clip3)

                for feature in features:
                    feat_info = aligned_seg.extract_cs(feature.start,
                                                       feature.end)
                    fname = feature.name
                    if feat_info is None:
                        d[tname][fname + '_cs'].append('')
                        if aligned_seg.target_clip5 >= feature.end:
                            d[tname][fname + '_clip5'].append(feature.length)
                            d[tname][fname + '_clip3'].append(0)
                        elif aligned_seg.target_lastpos <= feature.start:
                            d[tname][fname + '_clip5'].append(0)
                            d[tname][fname + '_clip3'].append(feature.length)
                        else:
                            raise ValueError(
                                f"Should never get here for target {tname}:\n"
                                f"feature = {feature}\n"
                                f"target_clip5 = {aligned_seg.target_clip5}\n"
                                f"lastpos = {aligned_seg.target_lastpos}\n"
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
