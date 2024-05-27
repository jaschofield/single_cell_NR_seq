# single_cell_NR_seq
Code associated with Schofield and Hahn 2024 publication.

## Raw data
Raw data can be obtained from the GEO database under accession GSE247795. The series contains raw single cell NR-seq FASTQ datasets and partially processed files containing assigned cell barcodes (CellRanger outputs).

## Original code description
1. Raw FASTQ data should first be processed using the TimeLapse-seq mutation calling pipeline, generating a large .csv file containing mutation counts for each read ("counts.csv" file). The annotation file for pipeline execution is the Cellranger_PolyA.gtf provided in this repository to be consistent with the CellRanger analysis pipeline.
2. CellRanger is run using the script "Cellranger_run.sh" and the PolyA_10X_genome provided in this repository, generating a list of identified cell barcodes as well as a .bam file containing cell identity information.
3. The TimeLapse-seq analysis output counts.csv file is joined with the CellRanger .bam output file using the "Cellranger_TL_join.R" script.
4. Potential cell doublet droplets are filtered using the "Doublet_removal.R" script.
5. Fano calculation and bootstrapping analysis is performed using the "fano_bootstrapping.R" script, which can be executed on a Slurm scheduler using the "fano_boot_execute.sh" script.
6. Fon calculation is run using "Fon_calculation.R" script, using intermediate output of "Cellranger_TL_join.R" script.
7. Cell cycle analysis is run using "cell_cycle.R" script.
8. Fraction new estimates and half life calculation are generated using the "bakR_run.R" script.

## Outside code description
The TimeLapse-seq mutation calling pipeline can be found here: https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/

The estimation of transcript half-lives uses the bakR package found here: https://simonlabcode.github.io/bakR/

Cell cycle data (from Spellman et al. 1998) can be downloaded here: https://www.molbiolcell.org/doi/suppl/10.1091/mbc.9.12.3273
