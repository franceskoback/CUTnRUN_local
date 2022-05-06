# import packages required
import os
import sys

# read in input arguments that are required
sample_id		= sys.argv[1] # what to name the samples - easier to define than to subset
read1_fq_gz		= sys.argv[2] # read1.fq.gz
read2_fq_gz		= sys.argv[3] # read2.fq.gz
num_cores		= sys.argv[4] # how many cores to use for the analysis
aligned_folder	= sys.argv[5] # where to dump the aligned reads
index_touse = sys.argv[6]

# folders/programs to be used for initializing the script

index		=	index_touse # need to put all the bwa index files in the working directory!! Or bwa aln will not run!! because you did bwa index -p mm10, this needs to be "mm10" here
#read_length_file		=	"<path>/zzz_read_length/length" # maybe shouldn't have commented out but it points to text file with a single number so i just made this a variable defined below instead because i thought that was easier 

read_length=76

multimap_limit			=	str( 20 )

# pull out how long the reads are from the file- should be 42 though
#with open(read_length_file, 'r') as file:
#	read_length = str( int( file.read().strip() ) )


# create output directories if don't exist
if not os.path.exists(aligned_folder):
	os.makedirs(aligned_folder)

trim_folder	=	"" +aligned_folder+ "/trimmed_reads/"
if not os.path.exists(trim_folder):
	os.makedirs(trim_folder)

intermediates_folder	=	"" +aligned_folder+ "/intermediate_files/"
if not os.path.exists(intermediates_folder):
	os.makedirs(intermediates_folder)

unapired_trim_folder	=	"" +trim_folder+ "unpaired/"
if not os.path.exists(unapired_trim_folder):
	os.makedirs(unapired_trim_folder)

logs_folder = "" +trim_folder+ "logs/"
if not os.path.exists(logs_folder):
	os.makedirs(logs_folder)


# start actually doing work now
os.chdir( aligned_folder )
# start saving output script
script_name = "one_cut_n_run_" +sample_id+ "_trim_align.sh"
output_script = open( script_name, 'w' )

# load necessary modules 
output_command = "echo \"loading necessary modules" + "\""
output_script.write(output_command)
output_script.write("\n")

output_command = "module load CBI trimgalore samtools "
output_script.write(output_command)
output_script.write("\n")

# make output command write to file
output_command = "echo \"TrimGalore Trimming : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")

output_command = "echo \"TrimGalore Trimming : " +sample_id+ "\""


# trim the files with Trim Galore


output_command = "trim_galore --illumina --paired --fastqc -o trim_galore/ " +read1_fq_gz+ " " +read2_fq_gz
output_script.write(output_command)
output_script.write("\n\n")



# make output command write to file
output_command = "echo \"Alignment using Bowtie2: " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")

# align the reads with bowtie2
# bwa aln -q10 -t8 hg18 reads_1.fastq > 1.sai
# bowtie2 -p 8 -x mm10 -1 out_1.fq  -2 out_2.fq | samtools view -@ 10 -Sb -o bowtie2.bam 1>1.txt 2>2.txt 
# bowtie2 -x mm10 -1 cells_BRD4_dia_S3_R1_001.fastq.gz -2 cells_BRD4_dia_S3_R2_001.fastq.gz|samtools view -bS - > out.bam
output_command = "bowtie2 -p" +num_cores+ " -x" +index+ " -1 " + read1_fq_gz + " -2 " +read2_fq_gz+ " | samtools view -@ 10 -Sb -o " +aligned_folder+ "/" +sample_id+  "_Aligned_Filtered.bam 1>1.txt 2>2.txt" 
output_script.write(output_command)
output_script.write("\n")


# make output command write to file
output_command = "echo \"Removing Unmapped reads : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")
# remove unmapped reads
output_command = "bamtools filter -isMapped true -in " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.bam > " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mapped.bam 2> " +logs_folder+ "" +sample_id+ "_bamtools_removeUnmapped.err"
output_script.write(output_command)
output_script.write("\n")
# make output command write to file
output_command = "echo \"Sort, Index, Alignment stats: " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")
# samtools sort -l 9 -O bam -@ " +num_cores+ " -o " +output_dir+ "" +sample_name+ ".sorted.bam " +output_dir+ "" +sample_name+ "Aligned.sortedByCoord.out.bam"
output_command = "samtools sort -l 9 -O bam -n -@ " +num_cores+ " -o " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mapped.bam"
output_script.write(output_command)
output_script.write("\n")
# fixmates
output_command = "samtools fixmate " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# sort fixed mates
output_command = "samtools sort -l 9 -O bam -@ " +num_cores+ " -o " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# index
output_command = "samtools index -b " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bam"
output_script.write(output_command)
output_script.write("\n")
# genome coverage map
output_command = "bamCoverage -b " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bam -o " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bw 2> " +logs_folder+ "" +sample_id+ "_bamCoverage.err"
output_script.write(output_command)
output_script.write("\n")
# sort fixed mates in way to convert to beds! - sort by names since thats whats required
output_command = "samtools sort -l 9 -O bam -n -@ " +num_cores+ " -n -o " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSortedBed.bam " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# move intermediate files
output_command = "mv " +sample_id+ "_Aligned.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_Aligned_Filtered.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_Aligned_Filtered.mappedSortedFixed.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_Aligned_Filtered.mapped.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_Aligned_r1.sai " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_Aligned_r2.sai " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
# get the mapping stats since bwa doesnt provide
output_command = "bamtools stats -in " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedSorted.bam > " +aligned_folder+ "/" +sample_id+ "_Aligned_Filtered.mappedStats.txt"
output_script.write(output_command)
output_script.write("\n\n")

# close the script file
output_script.close()

# make script executable
os.system("chmod +x " +script_name)


