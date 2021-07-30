
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

#Answer = input("Do you have a new input file? (Y/N):  ")

#if Answer == "Y":
#    input_data = input("Please write the name of the input file, if you do not have one, use /input_test_file.csv")
#elif Answer == "N":
#    input_data = "/input_test_file.csv" #why "/"
#else: print("Your input was neither Y or N")    

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

#path = '/Users/krde/Documents/Research/TCGA_3D/structure_lists'
#names = import_names(path)
#test = import_csvs(path, names[0], names[1])
#mutations = seq_pos_dict(input_df)

print("This what PDBminer is capable of doing at the moment. Program has finished.")
#print("Step 3: Control if PDBs contain the correct AA mutations..")
#evaluate_mutations()

#print("Step 4: ... choosing a model")
#choose_model()
