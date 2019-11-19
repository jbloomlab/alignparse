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
    affiliation: "1, 3"
affiliations:
  - name: Basic Sciences and Computational Biology, Fred Hutchinson Cancer Research Center, Seattle, Washington, USA
    index: 1 
  - name: Department of Genome Sciences and Medical Scientist Training Program, University of Washington, Seattle, Washington, USA
    index: 2
  - name: Howard Hughes Medical Institute, Seattle, Washington, USA
    index: 3
date: 19 November 2019
bibliography: paper.bib
---

# Summary & Purpose

Advances in sequencing technology have made it possible to generate large numbers of long, high-accuracy sequencing reads.
For instance, the new PacBio Sequel platform can generate hundreds of thousands of high-quality circular consensus sequences in a single run [@Rhoads:2015, @Hebert:2018].
Numerous programs exist for aligning these reads for purposes of genome assembly [@Chaisson:2012, @Li:2018].
However, these long reads can also be used for other purposes, such as sequencing PCR amplicons that contain various sequence features of interest.
For instance, PacBio circular consensus sequences have been used to identify the mutations in influenza viruses in single cells [@Russell:2019], or to link barcodes to gene mutants in deep mutational scanning [@Matreyek:2018].
For such applications, the alignment of the sequence to the template may be fairly trivial, but it is not trivial to then parse specific features of interest (such as mutations, unique molecular identifiers, cell barcodes, and flanking sequences) from these alignments.

Here we describe [alignparse](https://jbloomlab.github.io/alignparse/), a Python package for parsing complex sets of features from long sequences that map to known amplicons.
This package provides flexible tools for aligning long-read sequencing data to user-specified target sequences and then extracting specific subsequences corresponding to user-defined features for further analyses. 
Specifically, it allows the user to provide a complex "target" sequence in Genbank format that contains an arbitrary number of user-defined sub-sequence features. 
It then aligns the long sequencing reads to this target and filters alignments based on whether the user-specified features are present with the desired identity (which is also user-specified). 
Finally, it parses out the sequences, mutations, and/or accuracy of these features as specified by the user.
The flexibility of this package therefore fulfills the need for a tool to extract and analyze complex sets of features in large numbers of long sequencing reads.

# Uses & Examples 

Below are two example use cases of [alignparse](https://jbloomlab.github.io/alignparse/) from our research:

## Sequencing deep mutational scanning libraries

In deep mutational scanning experiments, researchers use mutant libraries to assay the effects of thousands of individual mutations to a gene-of-interest in one experiment [@Fowler:2014]. 
To examine the effects of mutations that appear together, researchers can link all mutations in a single variant to a unique molecular barcode [@Hiatt:2010]. 
Long-read PacBio sequencing can be used to link mutations and barcodes [@Matreyek:2018].   

The ``alignparse`` package provides a standard tool for parsing barcodes and their linked mutations from the long-read sequencing of deep mutational scanning libraries. 
It also allows for the parsing of any number of additional sequence features necessary for validating the quality of deep mutational scanning libraries. 
The [RecA deep mutational scanning library](https://jbloomlab.github.io/alignparse/recA_DMS.html) example demonstrates this use. Currently, ``alignparse`` is also being used to analyze additional deep mutational scanning libraries in the Bloom lab.

## Single-cell viral sequencing

Given their small size, several viral genomes can be sequenced in full using long-read sequencing technology. 
Using long-read sequencing in conjunction with single cell technologies (such as 10x Chromium), it is possible to barcode single virions and sequence the full-length genome of a virus that enters a single cell [@Russell:2019]. 
To analyze these viral sequences, the viral genomic sequences and associated mutations must be parsed in conjunction with cell-specific barcodes and additional sequence features introduced to ensure accurate variant calling. 
``alignparse`` provides a highly-customizable tool for parsing this data from single-cell viral sequencing data sets.
The [Single-cell virus sequencing](https://jbloomlab.github.io/alignparse/flu_virus_seq_example.html) example shows how these data can be readily analyzed using ``alignparse``. 


# How ``alignparse`` works

The ``alignparse`` package is designed to parse features--provided in a user-defined Genbank file--from alignments of long sequencing reads aligned to a user-specified reference. 
Inputs to ``alignparse`` include a user-defined Genbank file containing the sequence of one or more reference amplicons with defined features; a YAML file containing parsing specifications (such as filter cutoffs and desired outputs) for each feature; and the long-read sequencing data processed into ``fastq`` format.
These inputs are used to define a [Targets](https://jbloomlab.github.io/alignparse/alignparse.targets.html#alignparse.targets.Targets) object. 
The main functionality of ``alignparse`` uses a ``Targets`` object to create sequence alignments and parse out sequence features as defined in the input Genbank and YAML files. 

``alignparse`` aligns sequencing reads to user-defined reference sequences using [minimap2](https://github.com/lh3/minimap2) and the [alignparse.minimap2](https://jbloomlab.github.io/alignparse/alignparse.minimap2.html) submodule provides alignment specifications opitimzed for the two use cases described above. 
Alignments to be parsed by ``alignparse`` must contain short [cs tags](https://lh3.github.io/minimap2/minimap2.html#10).
These tags record descriptive information about the alignment in a succint format and are necessary for parsing features from ``Targets`` objects. Following parsing, features (and user-defined information about them) are output for each sequencing read in an intuitve data frame format. 

Downstream analyses of parsed features are facilitated by the [alignparse.consensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html) submodule.
This submodule provides tools for grouping reads by shared barcodes, determining consensus sequences for barcoded reads, and further processing mutation information for downstream analyses.
Since the main outputs from this package are in intuitive data frame formats, downstream analyses of reads parsed using ``alignparse`` can be highly customized by the user.
Thus, ``alignparse`` provides a flexible and useful tool designed specifically for parsing any number of complex features from long-read sequencing data of pre-defined amplicons.

# Acknowledgements

We would like to thank members of the Bloom lab for helpful discussions and beta testing. We also thank ___ **funding sources?** **ADD THE HIV DMS AND GENERAL VIRUS DMS R01s, PLUS SAYING THAT I AM AN INVESTIGATOR OF THE HOWARD HUGHES MEDICAL INSTITUTE.**

# References

