---
title: 'alignparse: A Python package for parsing complex features from high-throughput long-read sequencing'
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
For instance, the new PacBio Sequel platform can generate hundreds of thousands of high-quality circular consensus sequences in a single run [@Rhoads:2015; @Hebert:2018].
Good programs exist for aligning these reads for genome assembly [@Chaisson:2012; @Li:2018].
However, these long reads can also be used for other purposes, such as sequencing PCR amplicons that contain various features of interest.
For instance, PacBio circular consensus sequences have been used to identify the mutations in influenza viruses in single cells [@Russell:2019], or to link barcodes to gene mutants in deep mutational scanning [@Matreyek:2018].
For such applications, the alignment of the sequences to the targets may be fairly trivial, but it is not trivial to then parse specific features of interest (such as mutations, unique molecular identifiers, cell barcodes, and flanking sequences) from these alignments.

Here we describe [alignparse](https://jbloomlab.github.io/alignparse/), a Python package for parsing complex sets of features from long sequences that map to known targets.
Specifically, it allows the user to provide complex target sequences in Genbank Flat File format that contain an arbitrary number of user-defined sub-sequence features [@Sayers:2019]. 
It then aligns the sequencing reads to these targets and filters alignments based on whether the user-specified features are present with the desired identities (which can be set to different thresholds for different features). 
Finally, it parses out the sequences, mutations, and/or accuracy (sequence quality) of these features as specified by the user.
The flexibility of this package therefore fulfills the need for a tool to extract and analyze complex sets of features in large numbers of long sequencing reads.

# Uses & Examples 

Below are two example use cases of [alignparse](https://jbloomlab.github.io/alignparse/) from our research.
Code, data, and example output are included in the [alignparse](https://jbloomlab.github.io/alignparse/) documentation.

## Sequencing deep mutational scanning libraries

In deep mutational scanning experiments, researchers use mutant libraries to assay the effects of tens of thousands of individual mutations to a gene-of-interest in one experiment [@Fowler:2014]. 
One way to make deep mutational scanning of long gene variants work efficiently with short-read Illumina sequencing is to link the mutations in each variant to a unique molecular barcode [@Hiatt:2010]. 
This barcode linking can be done by long-read PacBio sequencing of the variant library [@Matreyek:2018], but it is then necessary to parse the resulting long reads to associate the barcode with the mutations in the variant.   

The [alignparse](https://jbloomlab.github.io/alignparse/) package provides a standard tool for parsing barcodes and linked mutations from the long-read sequencing data. 
It also allows for the parsing of additional sequence features necessary for validating the quality of deep mutational scanning libraries, such as the presence of terminal sequences or other identifying tags. 
The [RecA deep mutational scanning library example](https://jbloomlab.github.io/alignparse/recA_DMS.html) demonstrates this use. 

## Single-cell viral sequencing

Some viral genomes are sufficiently small to be sequenced in their entirety using long-read sequencing technology. 
Recent work has shown that such long-read sequencing of viral genes can be combined with standard single-cell transcriptomic technologies (such as 10x Chromium) to simultaneously sequence the infecting virus and characterize the transcriptome in single infected cells [@Russell:2019]. 
Such experiments require parsing the long-read viral sequences to identify viral mutations as well as cell barcodes, unique molecular identifiers, and other flanking sequences. 
The [single-cell virus sequencing example](https://jbloomlab.github.io/alignparse/flu_virus_seq_example.html) shows how such parsing can readily be performed using [alignparse](https://jbloomlab.github.io/alignparse/).

# How alignparse works

[alignparse](https://jbloomlab.github.io/alignparse/) takes the following inputs:

 1. One or more user-defined Genbank files containing the sequence of one or more alignment targets with an arbitrary number of user-defined features. These Genbank files can be readily generated using sequence editing programs, such as [ApE](https://jorgensen.biology.utah.edu/wayned/ape/) or [Benchling](https://www.benchling.com/).
 
 2. A YAML file containing parsing specifications for each feature. These specifications include filters indicating the maximal allowed mutations in each feature, as well as information on what output should be parsed for each feature (e.g., its sequence, its mutations, or simply if it is present).

3. A FASTQ file containing the long-read sequencing data. This file can be gzipped. There is no need to decompress gzipped FASTQ files first.

These inputs are used to define a [Targets](https://jbloomlab.github.io/alignparse/alignparse.targets.html#alignparse.targets.Targets) object. 
[alignparse](https://jbloomlab.github.io/alignparse/) then uses this [Targets](https://jbloomlab.github.io/alignparse/alignparse.targets.html#alignparse.targets.Targets) object to create sequence alignments and parse sequence features defined in the input Genbank and YAML files. 

[alignparse](https://jbloomlab.github.io/alignparse/) aligns sequencing reads to the targets using [minimap2](https://github.com/lh3/minimap2).
The [alignparse.minimap2](https://jbloomlab.github.io/alignparse/alignparse.minimap2.html) submodule provides alignment specifications optimized for the two example use cases described above. 
[alignparse](https://jbloomlab.github.io/alignparse/) uses the [cs tags](https://lh3.github.io/minimap2/minimap2.html#10) generated by [minimap2](https://github.com/lh3/minimap2) to extract the relevant features from the alignments into intuitive data frames or CSV files. 

We expect most users to align sequences and parse features in a single step using the [alignparse.targets.Targets.align_and_parse](https://jbloomlab.github.io/alignparse/alignparse.targets.html#alignparse.targets.Targets.align_and_parse) function.
However, these aligning and parsing steps can be carried out separately, as seen in the [Lassa virus glycoprotein](https://jbloomlab.github.io/alignparse/lasv_pilot.html) example.
Indeed, the [alignparse.targets.Targets.parse_alignment](https://jbloomlab.github.io/alignparse/alignparse.targets.html#alignparse.targets.Targets.parse_alignment) function should be able to parse features from any alignment file (in [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf)) as long as the alignments have [cs tags](https://lh3.github.io/minimap2/minimap2.html#10) and a corresponding [Targets](https://jbloomlab.github.io/alignparse/alignparse.targets.html#alignparse.targets.Targets) object has been defined that identifies the targets to which the query sequences were aligned and specifies the features to parse and filters to use. 

Downstream analyses of parsed features are facilitated by the [alignparse.consensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html) submodule.
This submodule provides tools for grouping reads by shared barcodes, determining consensus sequences for barcoded reads, and further processing mutation information for downstream analyses.
Since the main outputs from [alignparse](https://jbloomlab.github.io/alignparse/) are in intuitive data frame formats, downstream analyses can be highly customized by the user.
Thus, [alignparse](https://jbloomlab.github.io/alignparse/) provides a flexible and useful tool for parsing complex sets of features from high-throughput long-read sequencing of pre-defined targets.

# Code Availability 

The [alignparse](https://jbloomlab.github.io/alignparse/) source code is on GitHub at https://github.com/jbloomlab/alignparse and the documentation is at https://jbloomlab.github.io/alignparse.

# Acknowledgements

We would like to thank members of the Bloom lab for helpful discussions and beta testing.
This work was supported by the following grants from NIAID of the NIH: R01 AI141707 and R01 AI140891.
JDB is an Investigator of the Howard Hughes Medical Institute.

# References

