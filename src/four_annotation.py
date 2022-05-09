# import packages required
import os
from os import listdir
import sys
import numpy as np
import pandas as pd
from tabulate import tabulate

# Import Required Packages
normalized_beds_folder	= sys.argv[1]
output_folder = sys.argv[2]


# create output directories if don't exist
if not os.path.exists(output_folder):
	os.makedirs(output_folder)

def add_chr_to_bed_file(file_path):
    df = pd.read_csv(file_path,sep='\t',header=None)
    df.iloc[:, 0] = 'chr' + df.iloc[:, 0].astype(str)
    df.to_csv(file_path, sep='\t',header=None,index=False)

# iterate through all file
path = os.path.join(normalized_beds_folder,"/relaxed")
os.chdir( path )
for file in os.listdir():
    # Check whether file is in text format or not
    if file.endswith("relaxed.bed") | file.endswith("stringent.bed"):
        file_path = f"{normalized_beds_folder}/{file}"
  
        # call read text file function
        add_chr_to_bed_file(file_path)

# iterate through all file
path = os.path.join(normalized_beds_folder,"/stringent")
os.chdir( path )
for file in os.listdir():
    # Check whether file is in text format or not
    if file.endswith("relaxed.bed") | file.endswith("stringent.bed"):
        file_path = f"{normalized_beds_folder}/{file}"
  
        # call read text file function
        add_chr_to_bed_file(file_path)
	
	
# Move to the output folder
os.chdir( output_folder )

# save commands to the output script
script_name = "four_moveplots.sh"
output_script = open( script_name, 'w' )
command= "cd " + normalized_beds_folder
output_script.write(command)
output_script.write("\n")

command= "mv *.txt " + output_folder
output_script.write(command)
output_script.write("\n")

command= "mv *.png " + output_folder
output_script.write(command)
output_script.write("\n")


