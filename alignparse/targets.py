"""
=======
targets
=======

Defines :class:`Targets`, which holds :class:`Target` objects that define
alignment targets. Each :class:`Target` has some :class:`Feature` regions.

"""


import copy
import tempfile

import Bio.SeqIO

import dna_features_viewer

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt

import pandas as pd

import pysam

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
    feature_parse_specs : dict
        Keyed by name of each target in `seqsfile`, values target-level dicts
        keyed by feature names with values feature-level dicts keyed by feature
        names indicating what :meth:`Targets.parse_alignment` filters for and
        and returns for each feature. These feature-level dicts can be keyed
        by any of the following:

           - 'sequence' : sequence of feature, clipping indicated as deletion

           - 'mutations' : string of mutations, clipping indicated as deletion

           - 'mutation_count' : number mutated nucleotides excluding clipping

           - 'clip_count' : number of clipped nucleotides, total both termini

        The values for 'sequence' and 'mutations' are ignored; the values for
        'mutation_count' and 'clip_count' give the max allowable count before
        alignment is filtered (filtered if value > this).

        In addition, target-level dicts can optionally have the following keys
        that give the maximum amount that can be clipped from target or query
        sequence prior to the aligned region:

            - 'query_clip5' : amount query extends 5' of alignment

            - 'query_clip3' : amount query extends 3' of alignment

            - 'target_clip5' : amount target extends 5' of alignment;
              equivalent information can be represented by 'clip_count'
              for the 5' feature(s)

            - 'target_clip3' : amount target extends 3' of alignment;
              equivalent information can be represented by 'clip_count'
              for the 3' feature(s)

        See :meth:`Targets.parse_alignment` for some additional details.

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
    feature_parse_specs : dict
        A copy of the parameter of the same name.

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

        self.feature_parse_specs = copy.deepcopy(feature_parse_specs)

        # name of columns with alignment clipping
        self._clip_cols = ['query_clip5', 'query_clip3',
                           'target_clip5', 'target_clip3']

        # reserved columns for parsing, cannot be name of a feature
        self._reserved_cols = ['query_name'] + self._clip_cols

        self.targets = []
        self.target_names = []
        self._target_dict = {}
        for seqrecord in seqrecords:
            targetname = Target.get_name(seqrecord)
            if targetname not in self.feature_parse_specs:
                raise ValueError(f"target {targetname} not in "
                                 '`feature_parse_specs`')
            target = Target(seqrecord=seqrecord,
                            req_features=(set(self.feature_parse_specs
                                              [targetname].keys()) -
                                          set(self._reserved_cols)),
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

        extra_targets = set(self.feature_parse_specs) - set(self.target_names)
        if extra_targets:
            raise ValueError('`feature_parse_specs` includes non-existent '
                             f"targets {extra_targets}")

        for targetname, targetspecs in self.feature_parse_specs.items():
            if targetname not in self.target_names:
                raise ValueError('`feature_parse_specs` includes non-existent '
                                 f"target {targetname}")
            extrafeatures = (set(targetspecs) -
                             set(self._clip_cols)
                             .union(self.get_target(targetname).feature_names)
                             )
            if extrafeatures:
                raise ValueError(f"`feature_parse_specs` for {targetname} has "
                                 f"specs for unknown features {extrafeatures}")

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

                - 'target_clip5' : length at 5' end of target not in alignment

                - 'target_clip3' : length at 3' end of target not in alignment

                - a column with the name of each feature in the target giving
                  the ``cs`` string for that feature's alignment. If feature's
                  alignment is clipped (incomplete), this is indicated by
                  adding '<clipN>' (where 'N' is the amount of clipping) to the
                  ``cs`` string. For instance: '<clip7>:5*cg:3<clip2>'
                  indicates 7 nucleotides clipped at 5' end, 2 nucleotides at
                  3' end, and a ``cs`` string of ':5*cg:3' for aligned portion.

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
                aligned_target = self.get_target(tname)
                features = aligned_target.features
                if d[tname] is None:
                    d[tname] = {col: [] for col in self._reserved_cols +
                                aligned_target.feature_names}

                d[tname]['query_name'].append(aligned_seg.query_name)
                d[tname]['query_clip5'].append(aligned_seg.query_clip5)
                d[tname]['query_clip3'].append(aligned_seg.query_clip3)
                d[tname]['target_clip5'].append(aligned_seg.target_clip5)
                d[tname]['target_clip3'].append(aligned_target.length -
                                                aligned_seg.target_lastpos)

                for feature in features:
                    feat_info = aligned_seg.extract_cs(feature.start,
                                                       feature.end)
                    if feat_info is None:
                        d[tname][feature.name].append(feat_info)
                    else:
                        feat_cs, clip5, clip3 = feat_info
                        if clip5 != 0:
                            feat_cs = f"<clip{clip5}>{feat_cs}"

                        if clip3 != 0:
                            feat_cs = f"{feat_cs}<clip{clip3}>"

                        d[tname][feature.name].append(feat_cs)

        for target in d.keys():
            if target != 'unmapped':
                if d[target] is not None:
                    d[target] = pd.DataFrame.from_dict(d[target])

        return d


if __name__ == '__main__':
    import doctest
    doctest.testmod()
