"""
DESRIPTION TO BE UPDATED 
"""

import pandas as pd
import numpy as np
import json
import logging as log
import requests
import os
import re
import glob

def import_names(Path):
    #create an empty name list for each uniprot ID with corresponding pdbs
    names_list = []
    
    #ensure correct placement
    os.chdir(Path)
    
    #create a list of all files in directory
    for files in glob.glob("*.csv"):
        names_list.append(files)
    
    uniprot_id = []
    
    for i in range(len(names_list)):
        name = names_list[i][:-4]
        uniprot_id.append(name)
    
    return names_list, uniprot_id
    

def import_csvs(Path, names_list, uniprot_id):
    """Returns a list, where each item is a pandas dataframe containing the list of pdb sturctures"""
    
    #create an empty list for all the dataframes to be imported into.
    dataframes_list = []
    
    #ensure correct placement
    os.chdir(Path)
    
    #import all files in directory into a list of dataframes 
    for csvs in range(len(names_list)):
        temp_df = pd.read_csv(names_list[csvs])
        temp_df = temp_df.rename(columns={'Unnamed: 0': uniprot_id[csvs]})
        dataframes_list.append(temp_df)
        
        #update - consider adding names
    
    return dataframes_list

#Create a dictionary with the uniprot ID and mutational sites.

def seq_pos_dict(input_dataframe):
    """returns a dictionary with the uniprot ID as key and all the positions of the mutations as values"""
    
    #make a list of all uniprotIDs
    IDs = list(input_dataframe.Uniprot)
    
    #make an empty list
    mut = []
    
    #iterate over all the mutations, to extract the position of the mutation
    for i in range(len(input_dataframe)):
        
        mut_pos = re.split('(\d+)',input_dataframe.Mutations[i])
        mut_pos = [int(word) for sublist in map(str.split, mut_pos) for word in sublist if word.isdigit()]
        
        mut.append(mut_pos)
    
    #collect the information in a dictionary
    dictionary = {IDs[i]: mut[i] for i in range(len(IDs))}
        
    return dictionary

 ################################################################################
 
def find_ranges(segments):
    #output from get_uniprot_segments
    
    segment_array = np.array(list(segments.values()))
    
    #reshape to fit ranges
    segment_array = segment_array.reshape([-1,2])
    
    #create empty list for removal:
    to_be_removed = []
    
    for i in range(len(segment_array)):
        if None in segment_array[i]:
            to_be_removed.append(i)
        else:
            segment_array[i] == segment_array[i]
    
    tbr = np.array(to_be_removed)
    
    if len(tbr) > 0:
        segment_array = np.delete(segment_array, tbr, 0)
    else: 
        segment_array = segment_array
    
    return segment_array

def make_searchable_ranges(residue_range):
    #the output from find_segments
    
    a_list = []
    
    for jj in range(len(residue_range)):
        
        ranges_ = np.arange(residue_range[jj][0], residue_range[jj][1]+1)
        
        a_list.append(ranges_)
    
    ranges_ = np.hstack(a_list)
    
    return ranges_

def find_mutational_sites(mutation_location, range_of_model):
    
    presence_of_mutations = []
    
    for mutation in range(len(mutation_location)):
        
        if mutation_location[mutation] in range_of_model:
            r = 1
        
        else:
            r = 0
            
        presence_of_mutations.append(r)
    
    if 0 in presence_of_mutations:
        mut_result = False 
    
    else: 
        mut_result = True
    
    return mut_result

##################################################################################
#Slimfast Functions

# Get the module logger
logger = log.getLogger(__name__)

def get_uniprot_segments(pdb_id, uniprot_id):
    """Get which portions (segments) of a protein, identified by its
    UniProt ID, are covered in a given PDB file, given its PDB ID. SLiMfast code.
    """

    # URL where the mapping between PDB entities in a PDB entry
    # and their corresponding UniProt data is stored on PDBe
    mapping_url = \
        "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{:s}"

    # Get the data from PDBe
    response = requests.get(mapping_url.format(pdb_id))
    
    # Convert to text and read the resulting dictionary through JSON
    response_text = json.loads(response.text)

    # Get the data about the mappings
    data = \
        response_text[pdb_id.lower()]["UniProt"][uniprot_id]["mappings"]
    
    # Create an empty dictionary to store the mapping between
    # the PDB ID and the chains corresponding to the protein
    # identified by the UniProt ID
    mappings = {}

    # For each mapping found
    for m in data:

        # Get the chain ID
        pdb_chain = m["chain_id"]

        # Get the start and end of the chain as stored in the PDB
        pdb_start = m["start"]["author_residue_number"]
        pdb_end = m["end"]["author_residue_number"]

        # Get the start and end of the chain as they would be
        # in UniProt numbering
        uniprot_start = m["unp_start"]
        uniprot_end = m["unp_end"]

        if pdb_start == "null" or pdb_end == "null":
            warnstr = \
                f"Mapping for UniProt ID {uniprot_id} residues " \
                f"{unp_start}-{unp_end} not available in PDB " \
                f"{pdb_id} chain {pdb_chain}. Therefore, this " \
                f"segment will not be considered."
            logger.warning(warnstr)
            continue

        # Create a tuple of tuples to store those ranges
        ranges = ((pdb_start, pdb_end), (uniprot_start, uniprot_end))

        # If the chain ID was already encountered (i.e. the PDB
        # chain corresponds to multiple discontinuous protein
        # segments)
        if pdb_chain in mappings.keys():
            # Add the ranges to the mapping
            mappings[pdb_chain].append(ranges)

        # If it is the first time that the chain ID is encountered
        else:
            # Create a new list and add the mapping as first element
            mappings[pdb_chain] = [ranges]

    # Return the mappings
    return mappings

##############################################################

def validate_sequence(path, input_df):
    
    #create an empty list to hold the found pdb id's
    pdb_list = []
    
    path = path + "/structure_lists"
    
    #import all the Uniprot ID's to which pdb IDs were found
    names = import_names(path)
    
    #Import the csv files of each uniprot id's PDBs as a list of dataframes
    csv_ = import_csvs(path, names[0], names[1])
    
    #Find all the mutational positions.
    mutations = seq_pos_dict(input_df)
    
    for i in range(len(names[1])):
        
        print(names[1][i])
        
        #create a list of all the pdb's realted to one uniprot ID
        r_list = list(csv_[i][names[1][i]])
        
        for j in range(len(r_list)):
            
            #search for the segments covered by the model in PDB
            r = get_uniprot_segments(r_list[j], names[1][i])
            
            #Clean up the segments, removing all the sequence elements of the PDB, for which 
            #there is not a model
            res_range = find_ranges(r)
            
            #Convert the segment descriptors [start, stop] into searchable ranges
            searchable_range = make_searchable_ranges(res_range)
            
            #Seach through the ranges for the relevant mutational site 
            placement = find_mutational_sites(mutations[names[1][i]], searchable_range)
            
            #if all the mutational sites realated to a specific uniprot ID is found in a singular PDB model, 
            #the model is kept, if not the search contunues. 
            if placement == True:
                pdb_list.append(r_list[j]) 
                break
            else:
                continue
    
    pdb_uniprot_df = pd.DataFrame(list(zip(names[1], pdb_list)), columns =['Uniprot_ID', 'PDB'])   
    pdb_uniprot_df.to_csv("pdb_uniprot_dataframe.csv")      
    
    #something if they are not found - start over, look for all -1 ? 
            
    return(pdb_uniprot_df)
