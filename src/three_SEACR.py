# Import Required Packages 
import os
from os import listdir
import sys
import numpy as np
import pandas as pd
from tabulate import tabulate

# Read in Input Arguments
seacr_path = sys.argv[1] # path/to/SEACR_1.3.sh
bedgraphs_folder = sys.argv[2] # Path to bedgraphs 
output_folder = sys.argv[3] # Output folder to store results and resulting bash script this script generates 
igg_control = sys.argv[4] # "y" or "n" if there is an igg control -- todo: "n" will become numberic threshold to return top n fraction of peaks 


threshold = 0.1
# create output directories if don't exist
if not os.path.exists(output_folder):
	os.makedirs(output_folder)

# Change Directories to Bedgraphs Folder
os.chdir(bedgraphs_folder)

if igg_control=="y":
    igg_bedgraphs = []
    for file_name in listdir():
        if file_name.startswith("igg"):
            igg_bedgraphs.append(bedgraphs_folder+ "/" +file_name)
    igg_bedgraphs.sort()
    
## TO DO : make this script work if there is no igg control
# (Use a numeric threshold n between 0 and 1 to return the top n fraction of peaks based on total signal within peaks) 

# Grouping the Nuclei Data
nuclei_bedgraphs = []
for file_name in listdir():
	if file_name.startswith("nuclei"):
		nuclei_bedgraphs.append(bedgraphs_folder+ "/" +file_name)
# sort alphabetically
nuclei_bedgraphs.sort()

# Grouping the Cell Data
cell_bedgraphs = []
for file_name in listdir():
	if file_name.startswith("cells"):
		cell_bedgraphs.append(bedgraphs_folder+ "/" +file_name)
# sort alphabetically
cell_bedgraphs.sort()

# Grouping the Unfixed Cells Data
unfixed_bedgraphs = []
for file_name in listdir():
	if file_name.startswith("unfixedCells"):
		unfixed_bedgraphs.append(bedgraphs_folder+ "/" +file_name)
# sort alphabetically
unfixed_bedgraphs.sort()

# Move to the output folder
os.chdir( output_folder )

# save commands to the output script
script_name = "three_SEACRcall.sh"
output_script = open( script_name, 'w' )

# move to the output directory
output_command = "cd " +output_folder
output_script.write(output_command)
output_script.write("\n\n")

for count in range(len(nuclei_bedgraphs)):
    output_command = "echo \"Calling SEACR on Nuclei Data  - Non and Relaxed\""
    # Non is no normalization of control to target data
    output_script.write(output_command)
    output_script.write("\n")
    prefix = nuclei_bedgraphs[count].rsplit('/', 1)[-1].rsplit('.')[0]
    current_command= nuclei_bedgraphs[count]
    if igg_control=="y":
        output_command = "bash " +seacr_path+ " " +str(nuclei_bedgraphs[count])+ " " +str(igg_bedgraphs[0])+ " non relaxed " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    else: 
        output_command = "bash " +seacr_path+ " " +str(nuclei_bedgraphs[count])+ " " +str(threshold)+ " non relaxed " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
        
    

output_script.write("\n\n")

for count in range(len(cell_bedgraphs)):
    output_command = "echo \"Calling SEACR on Cell Data  - Non and Relaxed\""
    # Non is no normalization of control to target data
    output_script.write(output_command)
    output_script.write("\n")
    prefix = cell_bedgraphs[count].rsplit('/', 1)[-1].rsplit('.')[0]
    current_command= cell_bedgraphs[count]
    if igg_control=="y":
        output_command = "bash " +seacr_path+ " " +str(cell_bedgraphs[count])+ " " +str(igg_bedgraphs[0])+ " non relaxed " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    else: 
        output_command = "bash " +seacr_path+ " " +str(cell_bedgraphs[count])+ " " +str(threshold)+ " non relaxed " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    
output_script.write("\n\n")

for count in range(len(unfixed_bedgraphs)):
    output_command = "echo \"Calling SEACR on Unfixed Data  - Non and Relaxed\""
    # Non is no normalization of control to target data
    output_script.write(output_command)
    output_script.write("\n")
    prefix = unfixed_bedgraphs[count].rsplit('/', 1)[-1].rsplit('.')[0]
    current_command= cell_bedgraphs[count]
    if igg_control=="y":
        output_command = "bash " +seacr_path+ " " +str(unfixed_bedgraphs[count])+ " " +str(igg_bedgraphs[0])+ " non relaxed " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    else: 
        output_command = "bash " +seacr_path+ " " +str(unfixed_bedgraphs[count])+ " " +str(threshold)+ " non relaxed " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")

output_script.write("\n\n")    
    
for count in range(len(nuclei_bedgraphs)):
    output_command = "echo \"Calling SEACR on Nuclei Data  - Non and Stringent\""
    # Non is no normalization of control to target data
    output_script.write(output_command)
    output_script.write("\n")
    prefix = nuclei_bedgraphs[count].rsplit('/', 1)[-1].rsplit('.')[0]
    current_command= nuclei_bedgraphs[count]
    if igg_control=="y":
        output_command = "bash " +seacr_path+ " " +str(nuclei_bedgraphs[count])+ " " +str(igg_bedgraphs[0])+ " non stringent " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    else: 
        output_command = "bash " +seacr_path+ " " +str(nuclei_bedgraphs[count])+ " " +str(threshold)+ " non stringent " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")

output_script.write("\n\n")

for count in range(len(cell_bedgraphs)):
    output_command = "echo \"Calling SEACR on Cell Data  - Non and Stringent\""
    # Non is no normalization of control to target data
    output_script.write(output_command)
    output_script.write("\n")
    prefix = cell_bedgraphs[count].rsplit('/', 1)[-1].rsplit('.')[0]
    current_command= cell_bedgraphs[count]
    if igg_control=="y":
        output_command = "bash " +seacr_path+ " " +str(cell_bedgraphs[count])+ " " +str(igg_bedgraphs[0])+ " non stringent " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    else: 
        output_command = "bash " +seacr_path+ " " +str(cell_bedgraphs[count])+ " " +str(threshold)+ " non stringent " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    
output_script.write("\n\n")

for count in range(len(unfixed_bedgraphs)):
    output_command = "echo \"Calling SEACR on Unfixed Data  - Non and Stringent\""
    # Non is no normalization of control to target data
    output_script.write(output_command)
    output_script.write("\n")
    prefix = unfixed_bedgraphs[count].rsplit('/', 1)[-1].rsplit('.')[0]
    current_command= cell_bedgraphs[count]
    if igg_control=="y":
        output_command = "bash " +seacr_path+ " " +str(unfixed_bedgraphs[count])+ " " +str(igg_bedgraphs[0])+ " non stringent " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")
    else: 
        output_command = "bash " +seacr_path+ " " +str(unfixed_bedgraphs[count])+ " " +str(threshold)+ " non stringent " + str(prefix)
        output_script.write(output_command)
        output_script.write("\n")

output_script.write("\n\n")

output_command="mkdir relaxed"
output_script.write(output_command)
output_script.write("\n")
output_command="mkdir stringent"
output_script.write(output_command)
output_script.write("\n")
output_command="mv *.relaxed.bed ./relaxed"
output_script.write(output_command)
output_script.write("\n")
output_command="mv *.stringent.bed ./stringent"
output_script.write(output_command)
output_script.write("\n")

