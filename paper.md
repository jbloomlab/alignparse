---
title: 'alignparse: A Python package for aligning and parsing long-read sequencing data'
tags:
 - Python
 - PacBio
 - genomics
 - sequencing
 - deep mutational scanning
authors:
 - name: Katharine H.D. Crawford
   orcid:
 - name: Jesse D. Bloom
   orcid: 


# Summary

In the past decade, technology for genome sequencing has been advancing at an unprecented rate. (citation) Along with becoming cheaper, faster, and easier to sequence short segments of DNA, it has recently become feasible to sequence long reads of DNA (> 1000 bp) at scale. This has led to an increase in the amount of long-read sequencing data that is being produced and an expansion of using long-read sequencing into areas other than large-scale genomics, such as viral sequencing and deep mutatinoal scanning. 

With this expansion of long-read sequencing, also comes an increased need for tools to analzye such data. Several tools exist to map long sequencing reads to reference genomes, but most of these tools were designed for eukaryotic genomics and do not have the flexibility needed to adapt to the ever increasing uses of long-read sequencing.

`alignparse` provides researchers a flexible tool for aligning their long-read sequencing data to a user-specified target sequence and to extract short sequences - corresponding to user-specified features - from their sequencing reads for further analysis. This flexibility allows for the indepth analysis of long-read sequencing results to be more widely applied. `alignparse` uses `minimap2` for aligning sequences and takes adavantage of the informative `cs` tag returned from `minimap2` alignments in order to parse features from the long-read sequencing results and characterize their mutations and accuracy. 

In our examples, we apply `alignparse` to analyze sequencing of bacterial and viral constructs. These examples show the utility of this package to analyze both full-length viral genomes and laboratory designed mutant libraries of single genes. 

`alignparse` is currently being actively used to analyze deep mutational scanning libraries of several proteins, but we expect its utility to span across fields. It is designed to be used by genomics researchers with long-read sequencing data that they need to align to known reference sequences and from which they need to extract important sequence features for further analysis.

# Acknowledgements

# References