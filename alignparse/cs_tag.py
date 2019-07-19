"""
=======
cs_tag
=======

Parse SAM entries with ``cs`` tag from ``minimap2`` output.

"""


import re


class Alignment:
    """Process a SAM alignment with a ``cs`` tag to extract features.

    Parameters
    ----------
    sam_alignment : pysam.AlignedSegment
        Aligned segment from `pysam <https://pysam.readthedocs.io>`_.

    Attributes
    ----------
    query_name : str
        Name of query in alignment.
    target_name : str
        Name of alignment target.
    cs : str
        The ``cs`` tag.
    query_clip5 : int
        Length at 5' end of query not in alignment.
    query_clip3 : int
        Length at 3' end of query not in alignment.
    target_clip5 : int
        Length at 5' end of target not in alignment.
    target_lastpos : int
        Last position of alignment in target (exclusive).

    """
    def __init__(self, sam_alignment):
        """See main class docstring."""
        self.query_name = sam_alignment.query_name
        self.target_name = same_alignment.reference_name
        self.query_clip5 = sam_alignment.query_alignment_start
        self.query_clip3 = (sam_alignment.query_length -
                            sam_alignment.query_alignment_end)
        self.target_clip5 = sam_alignment.reference_start
        self.target_lastpos = sam_alignment.reference_end

        if sam_alignment.is_reverse:
            raise ValueError(f"alignment {self.name} is reverse orientation")

        if sam_alignment.has_tag('cs'):
            self.cs = sam_alignment.get_tag('cs')
        else:
            raise ValueError(f"alignment {self.name} has no `cs` tag")


    def extract_cs(self, start, end, *,
                   max_clip5=0, max_clip3=0):
        """Extract ``cs`` tag corresponding to feature in target.

        Parameters
        ----------
        start : int
            Start of feature in target in 0, 1, ... numbering.
        end : int
            End of feature in target (not inclusive of this site).
        max_clip5 : int
            If alignment does not fully cover 5' end of feature,
            add gaps to returned ``cs`` string for uncovered region
            up to this length.
        max_clip3 : int
            Like `max_clip5` but for 3' end of feature.

        Returns
        -------
        str or `None`
            If region of target covered by feature is aligned, return
            str with ``cs`` tag for this portion of alignment. If that region
            is not covered (even adding any padding from `max_clip5` /
            `max_clip3`), return `None.

        """
        # make sure we have indexing correct (0- or 1-based)
        raise RuntimeError('not yet implemented')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
