##############
# Import packs
##############

import os
import pandas as pd
import argparse
import shutil

##############

from find_structure_list import find_structure_list
from validate_sequence import validate_sequence

#parser creation
parser = argparse.ArgumentParser(prog = "PDBminer",
                                 usage = '%(prog)s [-h] help',
                                 description='PDBminer test software')

# Arguments 
# Currently there will only be one argument, the input file. This may expand. 
    
parser.add_argument('Input',
                    metavar = 'input_file',
                    type=str,
                    help='The input file containing the Hugo_name, Uniprot ID, isoform and Mutations.')

args = parser.parse_args()

input_data = args.Input

path = os.getcwd()
#this is problematic when running on the server, because it makes it very 
#difficult to point to another directory. How do I fix this?


if not os.path.exists("structures"):
    os.makedirs("structures")

if not os.path.exists("output_files"):
    os.makedirs("output_files")

################################################################################################################################
#Functions

def import_csv(csv_full_path):
    
    input_dataframe = pd.read_csv(csv_full_path)
    
    input_dataframe['mutation_positions'] = input_dataframe['Mutations'].str.split(';').apply(lambda x: [int(y[1:-1]) for y in x])
    
    #if all input_dataframe.ClusterIDs are "N/A"
    if len(set(input_dataframe.ClusterID)) == 1 and input_dataframe.ClusterID[0] == "N/A":
        
        #Set all clusters as 999  
        #To accompas later code that require ClusterID to be an interger.
        #200 is chosen because it is unlikely to have 999 clusters. 
        input_dataframe.ClusterID = 999
        
    return input_dataframe

def clean_up(path):
    
    #move the output files into the output folder, if they were created.
    if os.path.exists(path + "/inspection_structures.csv"):
        
        shutil.move(path+"/inspection_structures.csv", path+"/output_files")
    
    if os.path.exists(path + "/potential_structures.csv"):
        
        shutil.move(path+"/potential_structures.csv", path+"/output_files")

    if os.path.exists(path + "/missing_ID.txt"):
        shutil.move(path+"/missing_ID.txt", path+"/output_files")
    
    #move all the downloaded stuctures into the structure directory. 
    #should consider deleting them
    
    for dir in os.listdir(path):
        if len(dir) == 2:
            shutil.move(path+"/"+dir, path+"/structures")    
    
    return

################################################################################################################################
#Running the script

input_dataframe = import_csv(path+"/"+input_data)
found_structures = find_structure_list(input_dataframe)
structure_lists = validate_sequence(found_structures, input_dataframe, path)
clean_up(path)
