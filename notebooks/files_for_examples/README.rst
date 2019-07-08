Files for example notebooks
===========================

RecA data
---------
Data from Danny Lawrence's deep mutational scanning of RecA:

  - `recA_amplicon.gb <recA_amplicon.gb>`_: PacBio amplicon of RecA

  - `lib-1_ccs.fastq <lib-1_ccs.fastq>`_ and `lib-2_ccs.fastq <lib-2_ccs.fastq>`_: FASTQ files from running PacBio ``ccs``, shortened to include just a handful of entries.

  - `lib-1_report.txt <lib-1_report.txt>`_ and `lib-2_report.txt <lib-2_report.txt>`_: ``ccs`` report files corresponding the the FASTQ files. Manually edited to also have stats for just a handful of entries.

VEP data
---------
Data from Kate Dusenbury's initial pilot PacBio sequencing run of VEP constructs. 

    - `LASV_G1959_WT.gb <LASV_G1959_WT.gb>`_: PacBio amplicon of LASV GP from strain G1959, wildtype sequence.
    - `LASV_G1959_OPT.gb <LASV_G1959_OPT.gb>`_: PacBio amplicon of LSAV GP fro strain G1959, codon optimized sequence.
    - `vep_pilot_2019-06-05_report.txt <vep_pilot_2019-06-05_report>`_: Report file for all CCSs from the VEP pilot PacBio sequencing. This contains info on all CCSs. 
    - `vep_pilot_test_reads.fastq <vep_pilot_test_reads.fastq>`_: FASTQ file containing the first 1000 CCSs from the VEP pilot PacBio sequencing run. This should contain ~100 reads matching the LASV G1959 WT sequence.
    - `LASV_G1959_targets.fasta <LASV_G1959_targets.fasta>`_: FASTA file with the sequences for the LASV G1959 targets - wildtype and codon optimized sequences.
    - `g1959_test_alignment.paf <g1959_test_alignment.paf>`_: PAF alignment output from aligning the G1959 targets with the `vep_pilot_test_reads`. Used minimap2 with `--cs=long` argument for alignment.
