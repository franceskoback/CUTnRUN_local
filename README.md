# CUTnRUN Local Version
## Bioinformatics Pipeline For CUT&RUN Analysis

This repo will contain the code needed to run the CUT&RUN pipeline. Currently in development.

For all of these steps, you will need to run the python script for the step, and then run the associated bash script that python script generates before moving onto the next step in the pipeline.

To start with fastqs, use scripts 1 and 2. If you already have sorted bedgraphs, skip to the "STARTING FROM 3:" Portion below. Else 
1. Run python script one
2. On first run, you'll need to make a bwa index, so you'll need to install the bwa software:
3. Once you have bwa downloaded (typing bwa in command window returns bwa help page) run the following command: **bwa index /wynton/group/gladstone/users/franceskoback/Projects/CUTnRUN/data/index/mm10.fa (where mm10.fa is the reference genome you have downloaded. See below step 23 to see an example of how to download this) **
4. 


These are run on each set of paired fastqs until you get a list of sorted bedgraphs, then proceed with the steps below.

### Steps to Install STARTING FROM 3: ##
1. Git clone this repository 
2. make results and data directories within this folder 
3. copy your data (.sorted.bedGraphs) into this repository or know where it is (cp -r /path/to/data ./data within this repo folder)
4. Run 3rd python script 
5. nano that bash script output and put #!/usr/bin/env bash on top of bash script if submitting as job 
6. Rename 03_SEACRcall.sh to three_SEACRcall.sh if submitting as a job-- DONE change this to be three, four, etc instead of starting with numbers
7. qsub -cwd -pe smp 12 -l mem_free=12G -l scratch=1000G -l h_rt=50:00:00 -m bea -M frances.koback@gladstone.ucsf.edu three_SEACRcall.sh
8. run four_Annotation.py as shown below 
9. make sure you have all packages downloaded in R (ChIPseeker,TxDb.Mmusculus.UCSC.mm10.knownGene, clusterProfiler,ReactomePA,tidyverse,ggupset, ggimage)
10. change working directory in four_Annotation.R : setwd("/wynton/group/gladstone/users/franceskoback/CUTnRUN/results/three_calledpeaks")
11. module load CBI r/4.1.3 , then run "Rscript 04_Annotation.R" this will generate plots in the directory above. 
12. To move them to a new folder called "four_ChIPSeeker", run bash script four_moveplots.sh 
13. On initial run, to make a Homer folder in your software directory (or wherever you want to have it downloaded), navigate to that software folder and type wet http://homer.ucsd.edu/homer/configureHomer.pl
14. perl configureHomer.pl -install 
15. PATH=$PATH:/wynton/home/srivastava/franceskoback/software/Homer/.//bin/
16. perl /wynton/home/srivastava/franceskoback/software/Homer/.//configureHomer.pl -install mm10
18. run 5th python script 
19. bash five_cut_n_run_homer_motifs_v1beds.sh
20. on wynton, do module load CBI bedops 
21. mkdir six_MemeMotifs in results folder 
22. go to three_calledpeaks and make relaxed and stringent folders then do mv *.relaxed.bed ./relaxed and mv *.stringent.bed ./stringent-- DONE make this automatically happen after running third step 
23. go to your data folder and do wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz to download the mm10.fa genome  
24. navigate to 6th script in src folder and nano to change the path names at the top of the Rscript to match your paths-- including the path to the mm10 genome above  
25. Make sure you're using R > 4.0 ( module load r/4.1.3 ) 
26. Rscript 06_MemeMotifs.R (DONE: make this six_MemeMotifs.R)




## Usage example: ##
- python3 **one_cut_n_run_pairedReads_filter_align.py** "cells_FLAG_S4" "/Users/fkoback/Documents/Projects/Arun/CUTnRUN/CnRAP/cutNrun_fastq_jan2022/cells_FLAG_S4_R1_001.fastq" "/Users/fkoback/Documents/Projects/Arun/CUTnRUN/CnRAP/cutNrun_fastq_jan2022/cells_FLAG_S4_R2_001.fastq" 8 "/Users/fkoback/Documents/Projects/Arun/CUTnRUN/CnRAP/results"
- python3 **two_cut_n_run_bamToBed_normalize_SEACRPrepv1.py** "/Users/fkoback/Documents/Projects/Arun/CUTnRUN/CnRAP/results" "/Users/fkoback/Documents/Projects/Arun/CUTnRUN/CnRAP/results_2_normalizedbeds" "/Users/fkoback/Documents/Projects/Arun/CUTnRUN/CnRAP/cutNrun_fastq_jan2022/mm10.chrom.sizes.txt"
- python **three_SEACR.py** "/Users/fkoback/software/SEACR/SEACR_1.3.sh" "/Users/fkoback/Documents/Projects/CUTnRUN/data"  "/Users/fkoback/Documents/Projects/CUTnRUN/three_calledpeaks" "y"
- **or** python **three_SEACR.py** "/wynton/home/srivastava/franceskoback/software/SEACR/SEACR_1.3.sh" "/wynton/group/gladstone/users/franceskoback/CUTnRUN/data"  "/wynton/group/gladstone/users/franceskoback/CUTnRUN/results/3_calledpeaks" "y"
- python **four_Annotation.py** "/Users/fkoback/Documents/Projects/CUTnRUN/results/three_calledpeaks"
- or python **four_Annotation.py** "/wynton/group/gladstone/users/franceskoback/CUTnRUN/results/three_calledpeaks"

- run **four_Annotation.R**

- python **five_HomerMotifs.py** /Users/fkoback/Documents/Projects/CUTnRUN/results_3_calledpeaks/relaxed_beds relaxed.bed /Users/fkoback/Documents/Projects/CUTnRUN/05_homer_peak_motifs
- or **python five_HomerMotifs.py** /wynton/group/gladstone/users/franceskoback/CUTnRUN/results/three_calledpeaks "relaxed.bed" /wynton/group/gladstone/users/franceskoback/CUTnRUN/results/five_HomerMotifs

- Run **six_MemeMotifs.R**

## Steps of Analysis (and to-dos for development)
  
  - **Script 1**
      - Trimming with trimmomatic ([CnRAP](https://star-protocols.cell.com/protocols/944#key-resources-table)) or TrimGalore([nf-core](https://nf-co.re/cutandrun))-- TBD. Wynton has TrimGalore module available). Modifying this workflow to make it use TrimGalore. The original script uses kseq to trim, but that is unneccessary with TrimGalore (Trimmomatic fails to trim reads containing 6 bp or less, which is why Kseq is used in conjunction with Trimmomatic. 
      - Alignment with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
      - Filtering (removing unmapped reads and gathering sort, index, and alignment statistics): [samtools](http://www.htslib.org/)
    
  - **Script 2**
      - Bams to Bedgraphs 
      - Optional normalization to spike-in genome


  - **Script 3**
      - Call [SEACR](https://github.com/FredHutch/SEACR)
      - **Argument 1**: path to SEACR
      - **Argument 2**: path to data folder containing bedgraphs
      - **Argument 3**: path to folder where you want to store your results (and bash script this python script generates)
      - **Argument 4**: "y" or "n"-- "y" if there is an igg control 

  - **Script 4 Annotation**
      - **four_Annotation.py Argument 1**: path to data containing bed file outputs from SEACR (*relaxed.bed and *stringent.bed)
      - **four_Annotation.R** run in RStudio or R to generate plots 

  - **Script 5: Homer Motifs**
      - **five_Motifs.py** make sure you have genome downloaded (ie perl /opt/anaconda3/envs/python385/share/homer-4.10-0/.//configureHomer.pl -install mm10
      - **Note**: make sure beds are in correct Homer format. The script "check_bed_format.py" makes sure the bed files are tab-delimited and allows you to change them to be tab-delimited if they are space delimited instead

  - **Script 6: Meme Motifs**
      - **Download Bedops** https://bedops.readthedocs.io/en/latest/content/installation.html
      - **Note**: change path names in R script to match your environment 

Polishing in progress to make it reproducible in other compute environments!  

