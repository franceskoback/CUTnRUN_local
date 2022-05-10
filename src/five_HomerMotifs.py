# import packages required
import os
from os import listdir
import sys

# read in input arguments that are required
peaks_dir				= sys.argv[1] # where peak files are located /results/three_calledpeaks
output_base_dir			= sys.argv[2] # where to put results, ie /results/five_HomerMotifs
#homer_path= sys.argv[4] # path to Homer ie /path/to/Homer

# Move to the folder with the saved bam and stats files
path = str(peaks_dir)+"/relaxed"
print(path)
os.chdir(path)


# figure out which are the nuclear files, which are the regular files
peak_beds = []
for file_name in listdir():
	if file_name.endswith("relaxed.bed"):
		peak_beds.append(peaks_dir+ "/relaxed/" +file_name)
  
path = str(peaks_dir)+"/stringent"
os.chdir(path)
for file_name in listdir():
    if file_name.endswith("stringent.bed"):
        peak_beds.append(peaks_dir+ "/stringent/" +file_name)
        
# sort alphabetically
peak_beds.sort()

# folders/programs to be used for initializing the script
output_directory		=	output_base_dir

# create output directories if don't exists
if not os.path.exists(output_directory):
	os.makedirs(output_directory)

# Move to the output folder
os.chdir( output_directory )

# save commands to the output script
script_name = "five_cut_n_run_homer_motifs_v1beds.sh"
output_script = open( script_name, 'w' )


# for each file re-arrange to make in the format for Homer to call peaks
for sample_counter in range(len(peak_beds)):
	# update the user as to whats going on
	output_command = "echo \"Running Motif Discovery on " +peak_beds[sample_counter].split("/")[-1].split("_R1")[0]+ "\""
	output_script.write(output_command)
	output_script.write("\n")
	output_command = "cd " + output_base_dir
	output_script.write(output_command)
	output_script.write("\n")
	# call Homer motif discovery
	output_command ="findMotifsGenome.pl " +peak_beds[sample_counter]+ " mm10 ./" +peak_beds[sample_counter].split("/")[8][:-35]+ "_given -size given -mask -p 20 -S 50"
	output_script.write(output_command)
	output_script.write("\n")
	output_command = "findMotifsGenome.pl " +peak_beds[sample_counter]+ " mm10 ./" +peak_beds[sample_counter].split("/")[8][:-35]+ " -size 50,100,200 -mask -p 20 -S 50"
	output_script.write(output_command)
	output_script.write("\n\n\n")

# close the script file
output_script.close()

# make script executable
os.system("chmod +x " +script_name)

