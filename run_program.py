##############
# Import packs
##############

import os
import pandas as pd
import argparse

##############

from find_structure_list import find_structure_list
from validate_sequence import validate_sequence

#parser creation
parser = argparse.ArgumentParser(prog = "PDBminer",
                                 usage = '%(prog)s [options] ',
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

################################################################################################################################

def import_csv(csv_full_path):
    
    input_dataframe = pd.read_csv(csv_full_path)
    
    input_dataframe['mutation_positions'] = input_dataframe['Mutations'].str.split(';').apply(lambda x: [int(y[1:-1]) for y in x])
    
    return input_dataframe

#csv_full_path = "/Users/krde/Documents/Research/TCGA_3D/input_test_file.csv"

print("####################################################################")

#upload csv file
input_dataframe = import_csv(path+"/"+input_data)

found_structures = find_structure_list(input_dataframe)
stucture_lists = validate_sequence(found_structures, input_dataframe, path)

