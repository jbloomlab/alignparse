"""
=======
targets
=======

Defines :class:`Targets`, which holds :class:`Target` objects that define
alignment targets. Each :class:`Target` has some :class:`Feature` regions.

"""


import tempfile

import Bio.SeqIO

import dna_features_viewer

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot as plt

from alignparse.constants import CBPALETTE


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

    """

    def __repr__(self):
        """Get string representation."""
        return (f"{self.__class__.__name__}(name={self.name}, seq{self.seq}, "
                f"features={self.features})")

    def __init__(self, *, seqrecord, req_features=frozenset(),
                 opt_features=frozenset(), allow_extra_features=False):
        """See main class docstring."""
        for attr in ['name', 'seq']:
            if not hasattr(seqrecord, attr):
                raise ValueError(f"`seqrecord` does not define a {attr}")
            setattr(self, attr, str(getattr(seqrecord, attr)))

        self.length = len(self.seq)

        allow_features = set(req_features) | set(opt_features)

        self._features_dict = {}
        self.features = []
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
            self._features_dict[feature_name] = feature

        missing_features = set(req_features) - set(self._features_dict)
        if missing_features:
            raise ValueError(f"{self.name} lacks features: {missing_features}")

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

    def __repr__(self):
        """Get string representation."""
        return f"{self.__class__.__name__}(targets={self.targets})"

    def __init__(self, *, seqsfile, req_features=frozenset(),
                 opt_features=frozenset(), allow_extra_features=False,
                 seqsfileformat='genbank'):
        """See main class docstring."""
        if isinstance(seqsfile, str):
            seqsfile = [seqsfile]

        seqrecords = []
        for f in seqsfile:
            seqrecords += list(Bio.SeqIO.parse(f, format=seqsfileformat))

        self.targets = []
        self._target_dict = {}
        for seqrecord in seqrecords:
            target = Target(seqrecord=seqrecord,
                            req_features=req_features,
                            opt_features=opt_features,
                            allow_extra_features=allow_extra_features,
                            )
            if target.name in self._target_dict:
                raise ValueError(f"duplicate target name of {target.name}")
            self.targets.append(target)
            self._target_dict[target.name] = target

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
        """Parse alignment of query to targets.

        Note
        ----
        **Link to describe ``cs`` tag format.**

        **Indels that occur at feature junctions are assigned as
        part of first feature? Look to see how ``minimap2`` handles
        indels in homopolymers.**

        Parameters
        ----------
        samfile : str
            SAM file with ``minimap2`` alignments with ``cs`` tag, typically
            created by :meth:`Targets.align`.
        multi_align : {'primary'}
            How to handle multiple alignments. **Needs some thought**.
            Maybe 'primary' = return only primary, etc...
        **unaligned? How do we handle unaligned queries?**
            Maybe either ignored or added to data frame with 'target'
            as `None`.

        Returns
        -------
        pandas.DataFrame
            Gives alignment for each feature in target. Columns are:
                - 'query' : alignment query name
                - 'target' : name of target to which query aligns.
                - 'feature' : ``cs`` tag equivalent for that feature, or
                  `None` if that feature not in that target.
                - **potentially another column related to `multi_align`?**

        """
        # I would recommend writing a separate function that somehow does
        # the splitting of the ``cs`` tags, perhaps in its own module. You
        # can then have some simple doctests for that. This could probably
        # go in its own module called `cs_tag` or something like that.
        # For reading the SAM file, I recommend `pysam`.
        # https://pysam.readthedocs.io/en/latest/usage.html#opening-a-file
        raise RuntimeError('not yet implemented')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
