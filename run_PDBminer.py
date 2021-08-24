
"""
This file contains the code required to run the PDBminer program

"""

##############
# Import packs
##############

import os
import pandas as pd
import argparse

##############

from find_structure_list import find_structure_list
from validate_sequence import validate_sequence

#from evaluate_mutations import evaluate_mutations
#from choose_model import choose_model

#############

#python3 run_PDBminer.py input_test_file.csv


#parser creation
parser = argparse.ArgumentParser(prog = "PDBminer",
                                 usage = '%(prog)s [options] ',
                                 description='PDBminer version 0.02 (29-July-2021)')

# Arguments 
# Currently there will only be one argument, the input file. This may expand. 
    
parser.add_argument('Input',
                    metavar = 'input_file',
                    type=str,
                    help='The input file containing the Hugo_name, Uniprot ID, isoform and Mutations.')

args = parser.parse_args()

input_data = args.Input

path = os.getcwd()  

print("####################################################################")

#upload csv file
input_dataframe = pd.read_csv(path+"/"+input_data)

#run find_structure_list
print("Step 1: input file is loaded and lists of 3D structures are being generated...")
find_structure_list(input_dataframe, path)

print("Step 1 has finalized and the list of PDB sturctures per uniprot ID are avialable in the directory")

#might be nice to input a progress bar here

print("Step 2: Checking for the domains containing the mutations...")
validate_sequence(path, input_dataframe)

print("Step 2 has finalized and the .csv file containig the PDB file, and mutational sites for each Uniprot ID is avialable in the directory")
