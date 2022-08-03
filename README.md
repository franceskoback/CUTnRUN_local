# CUTnRUN Pipeline
## Bioinformatics Pipeline For CUT&RUN Analysis

This repo will contain the code needed to run the CUT&RUN (cleavage under targets and release using nuclease) pipeline, a more efficient alternative to the standard ChIP-Seq method. 

For all of these steps, you will need to run the python script for the step, and then run the associated bash script that python script generates before moving onto the next step in the pipeline.
### Software Requirements: ###
1. [TrimGalore](https://github.com/FelixKrueger/TrimGalore) wrapper to apply adapter and quality trimming to fastq files -- wynton has this already 
2. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) for aligning. Other pipelines use BWA, but Bowtie2 seems to perform better and does not require Stampy after alignment. 
3. [samtools](http://www.htslib.org/download/)
4. [bamtools](https://github.com/pezmaster31/bamtools) install with **conda install -c bioconda bamtools**  **for these, make sure your conda folder is in your path, ie export PATH="/your/path/to/miniconda3:$PATH"**
5. [deeptools](https://deeptools.readthedocs.io/en/develop/) **conda install -c conda-forge -c bioconda deeptools** 
6. [Bwa](https://github.com/lh3/bwa) follow the instructions in that link to download bwa and make your reference index by running the following command: **bwa index /path/to/genome.fa** 
7. Make sure you have all required packages installed in R (**ChIPseeker,TxDb.Mmusculus.UCSC.mm10.knownGene, clusterProfiler,ReactomePA,tidyverse,ggupset, ggimage**)
8. [Homer](http://homer.ucsd.edu/homer/) To download this, make a Homer folder in your software directory (or in whatever location you want to have Homer downloaded). Then navigate to that folder and type **wet http://homer.ucsd.edu/homer/configureHomer.pl** , **perl configureHomer.pl -install, PATH=$PATH:/your/path/to/Homer/.//bin/** ,then once you have Homer added to your path, install the relevant genome, ie  **perl /wynton/home/srivastava/franceskoback/software/Homer/.//configureHomer.pl -install mm10**
9. [Bedops](https://bedops.readthedocs.io/en/latest/content/installation.html)

To start with fastqs, use scripts 1 and 2. If you already have sorted bedgraphs, skip to the "STARTING FROM SORTED BEDGRAPHS:" Portion below. Else 
### STARTING FROM PAIRED FASTQS: ### 
1. Git clone this repository 
2. make results and data directories within this folder (mkdir results, mkdir data)
3. copy your paired fastq data into this repository (..R1_001.fastq.gz and R2_001.fastq.gz or something similar) or cp -r /path/to/data ./data within this repo folder)
5. Run python script one and associated bash script following the format in the "Usage example" portion of the readme below. This will output a bash script in your results/one_aligned folder that you will then run to do the first step in this pipeline, ie **bash one_cut_n_run_cells_BRD4_dia_trim_align.sh** **NOTE: Sometimes this alignment will output data that contains lines of the form: chr1_GL456221_random, chr4_JH584293_random, etc. In order to fix this issue, and get rid of any lines that did not align to  chromomes properly, in your results/one_aligned folder:  samtools view -o out.bam in.bam `seq 1 21 | sed 's/^/chr/'`**  where 1 and 21 correspond to the mouse genome. 
8. Run python script two following the format in the "Usage example" portion of the readme below. This will output a bash script in your results/two_normalize folder that you will then run to do the first step in this pipeline, ie **bash two_cut_n_run_bamToBed_normalize.sh** 
9. **Repeat steps 1-8 for every set of paired fastq data you wish to run. Then, once you have a bunch of "...Aligned_Filtered.mappedSorted.normalized.bed" files in your results/two_normalize folder, move on to the portion below. Alternatively, if you already have sorted bedgraphs, skip to the step below** 

### STARTING FROM SORTED BEDGRAPHS/ Pipeline Continued : ### 
11. Run python script three following the format in the "Usage example" portion of the readme below. This will output a bash script in your three_calledpeaks folder that you will then run to do the first step in this pipeline, ie **bash three_SEACRcall.sh** 
13.The SEACR output data structure is: 
```
<chr>	<start>	<end>	<total signal>	<max signal>	<max signal region>
```
where chr is a number, without "chr" appended. To add this "chr" to make the files compatible with further analysis, **run four_Annotation.py as shown in the "Usage example" below**-- **this will modify your bed files so be sure about this step before you run! Ie do not run this twice!!** 

14. Change working directory in four_Annotation.R: setwd("/your/path/to/results/three_calledpeaks") four_Annotation.R 
15. Make sure you're using R> 4.0 module load r/4.1.3 then run **Rscript four_Annotation.R** 
16. The four_Annotation.py python script will also output a bash script in your results/four__ChIPSeeker folder that will move all plots generated in the previous step into folders to clean up your directory. So once you've run four_Annotation.R, go to results/four__ChIPSeeker and run **bash four_moveplots.sh**
17. Run 5th python script as shown in the usage example below
18. Run 6th R script as shown in the usage example elow 


These are run on each set of paired fastqs until you get a list of sorted bedgraphs, then proceed with the steps below.


## Usage example: ##
- python **one_Align_wynton.py** "cells_BRD4_dia" "/local/path/to/cells_BRD4_dia_S3_R1_001.fastq.gz" "/local/path/to/cells_BRD4_dia_S3_R2_001.fastq.gz" 8 "/local/path/to//results/one_aligned" mm10
- **bash one_cut_n_run_cells_BRD4_dia_trim_align.sh** from within your one/aligned folder 
- python **two_Normalize.py** "/local/path/to/results/one_aligned" "/local/path/to/results/two_normalize" 
- **bash two_cut_n_run_bamToBed_normalize.sh** from within your results/two_normalize folder
- python **three_SEACR.py** "/path/to/SEACR/SEACR_1.3.sh" "/local/path/to/results/two_normalize" "/local/path/to/three_calledpeaks" "n"
- **bash three_SEACRcall.sh** from within your results/three_calledpeaks folder
- python **four_Annotation.py** "/local/path/to/results/three_calledpeaks" "/local/path/to/results/four_ChIPSeeker"
- module load r/4.1.3 then **Change working directory in four_Annotation.R: setwd("/your/path/to/results/three_calledpeaks")** run **Rscript four_Annotation.R**
- **bash four_moveplots.sh** from within your results/four_ChIPSeekr folder
- **python five_HomerMotifs.py** "/local/path/to/results/three_calledpeaks "relaxed.bed" "/local/path/to/results/five_HomerMotifs"
- **Change path names in R script to match your environment** run **Rscript six_MemeMotifs.R** 

## Steps of Analysis Explained
  
  - **Script 1**
      - Trimming with TrimGalore([nf-core](https://nf-co.re/cutandrun))-- TBD. Wynton has TrimGalore module available). The [CnRAP](https://star-protocols.cell.com/protocols/944#key-resources-table) pipeline uses trimmomatic and  uses kseq to trim, but kseq is unneccessary with TrimGalore (Trimmomatic fails to trim reads containing 6 bp or less, which is why Kseq is used in conjunction with Trimmomatic.) 
      - Alignment with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
      - Filtering (removing unmapped reads and gathering sort, index, and alignment statistics): [samtools](http://www.htslib.org/)
    
  - **Script 2**
      - Converts bams to bedgraphs 
      - Optional normalization to spike-in genome


  - **Script 3**
      - Takes your normalized bedgraphs, sorts them, and makes coverage files before running [SEACR](https://github.com/FredHutch/SEACR) to call peaks. 
      - **Argument 1**: path to SEACR
      - **Argument 2**: path to data folder containing bedgraphs
      - **Argument 3**: path to folder where you want to store your results (and bash script this python script generates)
      - **Argument 4**: "y" or "n"-- "y" if there is an igg control 
      - **Note**: If not using igg control, this script has a built-in manual threshold value of 0.1. To modify this threshlold, simply enter the python script and modify the value manually. 

  - **Script 4 Annotation**
      - **four_Annotation.py Argument 1**: path to data containing bed file outputs from SEACR (*relaxed.bed and *stringent.bed)
      - **four_Annotation.R** run in RStudio or R to generate plots 

  - **Script 5: Homer Motifs**
      - **five_Motifs.py** make sure you have genome downloaded (ie perl /opt/anaconda3/envs/python385/share/homer-4.10-0/.//configureHomer.pl -install mm10
      - **Note**: make sure beds are in correct Homer format. The script "check_bed_format.py" makes sure the bed files are tab-delimited and allows you to change them to be tab-delimited if they are space delimited instead

  - **Script 6: Meme Motifs**
      - **Note**: change path names in R script to match your environment 
      
This pipeline was generated for use locally. To see the version designed for use in a HPC like wynton, navigate to [CUTnRUN](https://github.com/franceskoback/CUTnRUN) 
