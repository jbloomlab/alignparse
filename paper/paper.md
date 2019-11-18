---
title: 'alignparse: A Python package for parsing complex features out of sequence alignments'
tags:
  - Python
  - PacBio
  - deep mutational scanning
  - single-cell virus sequencing
  - genomics
  - sequencing
authors:
  - name: Katharine H.D. Crawford
    orcid: 0000-0002-6223-4019
    affiliation: "1, 2"
  - name: Jesse D. Bloom
    orcid: 0000-0003-1267-3408
    affiliation: "1, 2, 3"
affiliations:
  - name: Basic Sciences and Computational Biology, Fred Hutchinson Cancer Research Center, Seattle, Washington, USA
    index: 1 
  - name: Department of Genome Sciences and Medical Scientist Training Program, University of Washington, Seattle, Washington, USA
    index: 2
  - name: Howard Hughes Medical Institute, Seattle, Washington, USA
    index: 3
date: 18 November 2019
bibliography: paper.bib
---

# Summary & Purpose

Advances in sequencing technology to generate large numbers of long high-accuracy reads.
For instance, the new PacBio Sequel platform can generate XX high-quality circular consensus sequences in a single Y [CITATIONS].
Numerous programs exist for aligning these reads for purposes for genome assembly.
However, these long reads can also be used for other purposes such as sequencing long PCR amplicons that contain various sequence features of interest.
For instance, PacBio circular consensus sequences have been used to sequence influenza viruses in single cells [CITATION], or link barcodes to gene mutants in deep mutational scanning [CITATION].
For such applications, the alignment of the sequence to the template may be fairly trivial, but it is not trivial to then parse specific features of interest (such mutations, unique molecular identifiers, cell barcodes, and flanking sequences) from these alignments.

Here we describe ``alignparse``, a Python package for parsing complex sets of features from long sequences that map to known amplicons.
It provides researchers a flexible tool for aligning long-read sequencing data to user-specified target sequences and for extracting short sequences---corresponding to user-specified features---from these reads for further analyses.
**It's more important to say that it allows the user to provide a complex customized template in Genbank format that provides an arbitrary number of sequence features. It then aligns, and filters alignments based on whether the user-specified features are present with the desired ientity (which is also user-specified). It then gets out the sequences or mutations as requested (basically, what it does)...**

**Possibly then go into examples, then discuss a bit about API**

**Then describe a bit about how it does it.**
[alignparse](https://jbloomlab.github.io/alignparse/) uses [minimap2](https://github.com/lh3/minimap2) for aligning sequences and takes adavantage of the informative [cs tag](https://lh3.github.io/minimap2/minimap2.html#10) returned from [minimap2](https://github.com/lh3/minimap2) alignments to parse features from the long-read sequencing data [@Li:2018]. This parsing is all carried out using the [alignparse.targets](https://jbloomlab.github.io/alignparse/alignparse.targets.html) submodule. 

Further analysis of these features is facilitated by the [alignparse.consensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html) module that provides tools for grouping reads by shared barcodes, determining consensus sequences for barcoded reads, and further processing mutation information for downstream analyses. Such downstream analyses can be highly customized by the user as the output from ``alignparse`` is in an intuitive data frame format that can be easily used as input for additional analyses. Thus, ``alignparse`` provides a flexible and useful tool designed specifically for the initial processing of long-read sequencing data from the deep sequencing of known amplicons. 

As long-read sequencing continues to be applied to areas outside of genome assembly, our analysis tools also need to adapt. ``alignparse`` provides an important tool specifically for analyzing data that comes from the application of long-read sequencing technology to the deep sequencing of pre-defined amplicons that contain complex sets of sequence features to be parsed out for further analyses. **Maybe make this more succinct**

# Audience

``alignparse`` is designed to be used by genomics researchers analyzing long-read sequencing data that maps to pre-defined amplicons. Researchers using long-read sequencing to sequence viral genomes or libraries of variants of single genes for deep mutational scanning will find this tool particularly useful for extracting information from their sequencing reads. Nonetheless, any scientists with long-read sequencing data where each read fully maps to a pre-defined amplicon and from which they would like to extract information about pre-defined sequence features, will find this tool useful for the analysis of such features.

# Uses & Examples 

The flexibility built into ``alignparse`` allows researchers from many areas of genomics to use this tool to analyze long-read sequencing data from known amplicons. Below are two known use cases:

## Sequencing deep mutational scanning libraries

Long-read sequencing of mutant libraries has allowed researchers conducing deep mutational scanning experiments to link all of the mutations in a variant with a unique molecular barcode [@Matreyek:2018]. This has facilitated the study of the effects of genetic variants on gene function by improving the data generated by deep mutational scanning experiments. In analyzing their data, the Shendure lab wrote [custom Python scripts](https://github.com/shendurelab/AssemblyByPacBio) that are rather specific for their use case and construct design. The ``alignparse`` package can be used to parse long read sequences for analyzing deep mutational scanning libraries like those in @Matreyek:2018, with much more flexibility than previously available scripts. In particular, the ability to easily define any number of sequence features is very useful for customizing sequence analysis.  **I would not go into exactly how Shendure lab did this, but rather just describe use case and how it's applied.**

Currently, ``alignparse`` is being used to analyze deep mutational scanning libraries in the Bloom lab, including the RecA deep mutational scanning library in [this](https://jbloomlab.github.io/alignparse/recA_DMS.html) example.

## Single-cell viral sequencing

Given their small size, several viral genomes can be sequenced in full using long-read sequencing technology. Indeed, for influenza virus, each of its eight gene segments can be sequenced as a single long-read sequencing amplicon. Using long-read sequencing in conjunction with single cell technologies (such as 10x Chromium), it is possible to barcode single flu virions and sequence the entire genome of a virus that enters a single cell. In order to analyze these viral sequences, the viral genomic sequences and associated mutations, must be parsed in conjunction with cell-specific barcodes and additional sequence features introduced to ensure accurate variant calling. Such experiments and analyses have been carried out for influenza virus, resulting in expanding our understanding of influenza infection heterogeneity [@Russell:2019]. These analyses were initially carried out with project-specific python scripts. However, as seen in the [Single-cell virus sequencing](https://jbloomlab.github.io/alignparse/flu_virus_seq_example.html) example, this data can be readily analyzed using ``alignparse``. 


The flexibility of the ``alignparse`` package to be used to analyze long-read sequencing data from deep mutational scanning to single-cell virology experiments helps show the broad utility of this package. With the increasing availability of long-read sequencing technologies, we imagine the use cases for this package will grow to include other instances of deep sequencing long amplicons. ``alignparse`` will also hopefully increase the reproducibility of analyzing such experiments by reducing the current reliance on project-specific scripts and providing researchers an actively-maintained, publically available tool for the initial parsing of amplicon-based long-read sequencing data.


# Acknowledgements

We would like to thank members of the Bloom lab for helpful discussions and beta testing.

# References

