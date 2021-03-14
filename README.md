# Project Description

From https://bf528.readthedocs.io/en/latest/content/projects/project_2_rnaseq_1/project_2_rnaseq_1.html :
“High-throughput sequencing technologies have helped researchers expand their understanding of genome-wide gene expression through their ability to detect potentially any mRNA molecule in a sample, rather than measuring only molecules with specific sequences as with microarrays. This agnostic approach to molecular profiling allows us to ask different and somewhat more complex questions of a sequencing dataset than a microarray dataset, (e.g. novel splice patterns, lincRNAs, coding mutations, etc), but the most basic information, mRNA abundance, is analogous to microarray expression measurements and extremely useful. In this project, you will download, QC, process, and analyze sequencing data that was generated to better understand how neonatal mice are able to regenerate their heart tissue but lose this ability later in development.”

In this project, we will analyze 8 samples from O’Meara et al 2015, and attempt to replicate the some of the findings shown in Figure 1 and 2.

This repository contains the scripts used to process and analyze the data.

O’Meara, C. C., Wamstad, J. A., Gladstone, R. A., Fomovsky, G. M., Butty, V. L., Shrikumar, A., Gannon, J. B., Boyer, L. A., & Lee, R. T. (2015). Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circulation Research, 116(5), 804–815. https://doi.org/10.1161/CIRCRESAHA.116.304269
s

# Contributors  
Data Curator: Allison Nau  
Programmer: Mae Rose Gott  
Analyst: Sheila Yee  
Biologist: Abhishek Thakar  


# Repository Contents
run_extract.qsub : script to unpack SRA file into FASTQ files. To submit as a qsub, login to SCC, navigate to the appropriate folder, and enter the command “qsub -P bf528 run_extract.qsub”.

20210224 proj2 data curator notes.txt: outlines the commands used to initially process the sample and set up the GIT repository. This is NOT a shell script, merely an outline of the command that were used in command line.

library_size.py: Python script to calculate average library size for the 8 samples looked at. (Data curator).
Will print the minimum, maximum, and average library size.
Expects info files to be in a subdirectory "prep_reads_info" in the current working directory.
Library size of P0_1 is loaded in by default WITHOUT an info file.


