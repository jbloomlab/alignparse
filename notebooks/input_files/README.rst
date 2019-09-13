Files for example notebooks
===========================

RecA data
---------
Data from Danny Lawrence's deep mutational scanning of RecA:

  - `recA_amplicon.gb <recA_amplicon.gb>`_: PacBio amplicon of RecA

  - `recA_feature_parse_specs.yaml <recA_feature_parse_specs.yaml>`_: YAML file spcifying how to parse features from `recA_amplicon.gb <recA_amplicon.gb>`_

  - `lib-1_ccs.fastq <lib-1_ccs.fastq>`_ and `lib-2_ccs.fastq <lib-2_ccs.fastq>`_: FASTQ files from running PacBio ``ccs``, shortened to include just a handful of entries.

  - `lib-1_report.txt <lib-1_report.txt>`_ and `lib-2_report.txt <lib-2_report.txt>`_: ``ccs`` report files corresponding the the FASTQ files. Manually edited to also have stats for just a handful of entries.

flu WSN data
-------------
Snippets of data from the single-cell virus sequencing study [Russell et al (2019)](https://jvi.asm.org/content/93/14/e00500-19). Note that the data in that study was **not** actually analyzed using ``alignparse``:

 - `flu_WSN_amplicons.gb <flu_WSN_amplicons.gb>`_: PacBio amplicons of flu genes
 
 - `flu_WSN_feature_parse_specs.yaml <flu_WSN_feature_parse_specs.yaml>`_: YAML file specifying how to parse features from PacBio amplicons.

 - `flu_WSN_ccs.fastq <flu_WSN_ccs.fastq>`_: FASTQ file with a snippet of data from one of the PacBio ``ccs`` runs.

VEP data
---------
Data from Kate Dusenbury's initial pilot PacBio sequencing run of VEP constructs. 

    - `LASV_G1959_WT.gb <LASV_G1959_WT.gb>`_: PacBio amplicon of LASV GP from strain G1959, wildtype sequence.
    - `LASV_G1959_OPT.gb <LASV_G1959_OPT.gb>`_: PacBio amplicon of LSAV GP fro strain G1959, codon optimized sequence.
    - `vep_pilot_report.txt <vep_pilot_report>`_: Report file from the VEP pilot PacBio sequencing. The ZMW numbers are scaled to be consistent with the size (250 CCSs) of the accompanying ``ccs`` file. 
    - `vep_pilot_ccs.fastq <vep_pilot_ccs.fastq>`_: FASTQ file containing the first 250 CCSs from the VEP pilot PacBio sequencing run. This should contain ~50 reads matching the LASV G1959 sequences.
