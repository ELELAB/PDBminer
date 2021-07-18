#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file contains the code required to run the PDBminer program

###################################################################

This program is used to build 3D structures   

####################################################################

"""

##############
# Import packs
##############
import os
import pandas as pd

##############

from find_structure_list import find_structure_list
from validate_sequence import validate_sequence

#from evaluate_mutations import evaluate_mutations
#from choose_model import choose_model

#############

print("####################################################################")
print("Version 0.01 of PDBminer. This text was written 13-July-2021.")
print("####################################################################")
print("To run this program, please make sure that this python file package is placed in the same directory as you input file and that the input file format is described in the read me file")
print("####################################################################")

path = os.getcwd()

Answer = input("Do you have a new input file? (Y/N):  ")

if Answer == "Y":
    input_data = input("Please write the name of the input file, if you do not have one, use /input_test_file.csv")
elif Answer == "N":
    input_data = "/input_test_file.csv"
else: print("Your input was neither Y or N")    

print("####################################################################")

#upload csv file
input_dataframe = pd.read_csv(path+input_data)

#run find_structure_list
print("Step 1: input file is loaded and lists of 3D structures are being generated...")
find_structure_list(input_dataframe, path)
print("Step 1 has finalized and the list of PDB sturctures per uniprot ID are avialable in the directory")

#might be nice to input a progress bar here

print("Step 2: .... validating the sequence")
validate_sequence(path, input_dataframe)

#print("Step 3: .... evaluating mutations the sequence")
#evaluate_mutations()

#print("Step 4: ... choosing a model")
#choose_model()
