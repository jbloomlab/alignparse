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

LASV data
---------
Data from Kate Dusenbury's initial pilot PacBio sequencing run of viral entry protein constructs. The data analyzed here is from the LASV glycoprotein constructs only.

    - `LASV_Josiah_WT.gb <LASV_Josiah_WT.gb>`_: PacBio amplicon of LASV GP from strain Josiah, wildtype sequence.
    - `LASV_Josiah_OPT.gb <LASV_Josiah_OPT.gb>`_: PacBio amplicon of LASV GP from strain Josiah, codon optimized sequence.
    - `lasv_pilot_report.txt <lasv_pilot_report>`_: Report file from the VEP pilot PacBio sequencing. The ZMW numbers are scaled to be consistent with the size (250 CCSs) of the accompanying ``ccs`` file. 
    - `lasv_pilot_ccs.fastq <lasv_pilot_ccs.fastq>`_: FASTQ file containing 250 CCSs mapping to LASV Josiah from the pilot PacBio sequencing run. All of these reads should match the LASV Josiah sequences, but ~50 of them should be filtered and ~200 of them should be aligned.