
"""
This file contains the code required to answer
the question: 

Does the found structure contain the sequence where the mutations are?

"""

##############
# Import packs
##############

import pandas as pd
import numpy as np
import json
import logging as log
import requests
import os
import re
import glob

##############

def import_names(Path):
    """Function that takes the path and returns two lists, 
    names list =  A list of all he PDBs related to a single Uniprot ID
    Uniprot_id = A list with the uniprot IDs to which there are solved structures. 
    """
    
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
    """Takes the path, and the output from import_names and 
    Returns a list, where each item is a pandas dataframe 
    containing the list of pdb sturctures"""
    
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
    """Takes the input file as input and returns a dictionary with 
    the uniprot ID as key and all the positions of the mutations as 
    values"""
    
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
    """The function takes the segments found in the quality control function
       get_uniprot_segments() and removes sequences that does not match.
    
    Input:
    --------------------------------------------------------------------------
     {'Chain A': [((pdb_start, pdb_end), (Uniprot_start, Uniprot_end))], 
      'Chain B': [(pdb_start, pdb_end), (Uniprot_start, Uniprot_end)]}
    
    Example: 
        {'A': [((19, 351), (39, 371))], 
         'B': [((None, 395), (377, 415))]}
        
    Chain A:
        
    PDB sequence       19 ##################################### 351
    Uniprot sequence          39 ##################################### 371
    
    Chain B
    
    PDB sequence         ? XXXXXXXXXX######### 395
    Uniprot sequence          377 ################### 415    
    
    
    Aim and Methodology of function:
    --------------------------------------------------------------------------
    The aim of find_ranges(segemnts) is to:
        i) remove the segments that does not match, such as chain B above
        ii) to correct the numbering to match the Uniprot numbering
            
    Output:
    --------------------------------------------------------------------------
    The output is an array of the corrected segment. 
    
    """
    
    #Convert the values of the input dictionary to an array    
    segment_array = np.array(list(segments.values()))
        
    #create empty list for removal:
    to_be_removed = []
    
    #control if any ranges are not covered
    for i in range(len(segment_array)):
        if None in segment_array[i]:
            to_be_removed.append(i)
        else:
            segment_array[i] == segment_array[i]
    
    #create a list of the segments the PDB does not cover if any. 
    #e.g. a chain not modelled. 
    tbr = np.array(to_be_removed)
    
    #if the list of uniprot segments that are not covered by the PDB
    #is longer than 0 some segments are to be removed. 
    if len(tbr) > 0:
        segment_array = np.delete(segment_array, tbr, 0)
    else: 
        segment_array = segment_array

    #For the remaining segments, a control of the length is to be conducted

    #First the segment is reshaped to access both pdb and uniprot numbering
    #seperately.    
    segment_array = segment_array.reshape([-1,2])
    
    #Length of the segment is controlled to be identical and is either 
    #discarded or the uniprot numbering is kept.
    
    #Collect all even numbers from the segment_array 
    array_number =  list(range(len(segment_array)))
    even_numbers = [num for num in array_number if num % 2 == 0]
    
    #create an empty list
    collected_segments = []
    
    #assess the length of each segment
    for j in even_numbers:
        if segment_array[j][0]-segment_array[j][1] == segment_array[j+1][0]-segment_array[j+1][1]:
            #choose the uniprot segment
            segment_ = segment_array[j+1]
            
            collected_segments.append(segment_)
    
    if len(collected_segments) == 0:
        segment_array = segment_array
    else:
        segment_array = collected_segments
    
    return segment_array

def make_searchable_ranges(residue_range):
    """The function takes the output range from find_ranges() and
    convert this to an array of AA positions that can be searched."""
    
    #create an empty list
    a_list = []
    
    #create a for loop to convert each residue range into an array
    for jj in range(len(residue_range)):
        
        ranges_ = np.arange(residue_range[jj][0], residue_range[jj][1]+1)
        
        #append each range
        a_list.append(ranges_)
    
    ranges_ = np.hstack(a_list)
    
    return ranges_

def find_mutational_sites(mutation_location, range_of_model):
    """A function that takes each mutational location and the range of 
    the PDB model and controls if the PDB model covers the mutational site. 
    The output is twofold:
        
    mut_result = True or False, if a specific mutation is present in pdb 
    mutation_location = The mutational location associated with each true or false
    
    """
    
    #empty list of collect the presence of mutations iteratively
    presence_of_mutations = []

    #for loop over the mutations to decide on a binary scale if the 
    # mutational location exist in the PDB file. 
    for mutation in range(len(mutation_location)):
        
        if mutation_location[mutation] in range_of_model:
            r = mutation_location[mutation]
            
            presence_of_mutations.append(r)

    
    return presence_of_mutations


def find_sequence(name,csv,mutations):
    
    """
    A function that takes a Uniprot ID, the corresponding CSV file and the mutations and search 
    trough the PDBs in the CSV file to find the matching values. 
    
    """
    
    #create a list of pdb ids
    r_list = list(csv[name])
    
    #create empty lists for collection of information    
    pdb_id_success = []
    ranges_list = []
    mutational_sites = []
    wt_list = []
    method = []
    res = []
        
    for j in range(len(r_list)):
            
        #search for the segments covered by the model in PDB
        try:
            r = get_uniprot_segments(r_list[j], name)
            
            #Clean up the segments, removing all the sequence elements of the PDB, for which 
            #there is not a model
            res_range = find_ranges(r)
                
            #create a dictionary of each PDBid and the range.
            if len(res_range) > 0:
                
                pdb_id_success.append(r_list[j])
                ranges_list.append(res_range)
                wt_list.append(csv.loc[csv[name]==r_list[j],"model_mutations"].iloc[0])
                method.append(csv.loc[csv[name]==r_list[j],"experimental_method"].iloc[0])
                res.append(csv.loc[csv[name]==r_list[j],"resolution"].iloc[0])
                
        except:
            print("failed to analyze " + r_list[j])
                
    for jj in range(len(pdb_id_success)):
            
        #Convert the segment descriptors [start, stop] into searchable ranges
        searchable_range = make_searchable_ranges(ranges_list[jj])
            
        #Seach through the ranges for the relevant mutational site 
        placement = find_mutational_sites(mutations[name], searchable_range)
            
        mutational_sites.append(placement)        

    return pdb_id_success, wt_list, mutational_sites, method, res

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

    # Get the data about the mappings - #p53, pos (csv_[2].iloc[90]) complains here for some reason.
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
    
    """
    A wrapper function to validate the sequence of each of the PDBs
    found in "find_structure.py"
    
    The input is:
        
        path = the path where the input file is.
        input_df = the input file
        
    The output is:
        
        a csv file per uniprot_id in the following format:
            
        uniprot_id | PDB_id | WT_status | Found_sites         
        
    """
    
    path = path + "/structure_lists"
    
    #import all the Uniprot ID's to which pdb IDs were found
    names = import_names(path)
    
    #Import the csv files of each uniprot id's PDBs as a list of dataframes
    csv_ = import_csvs(path, names[0], names[1])
    
    #Find all the mutational positions.
    mutations = seq_pos_dict(input_df)
    
    for i in range(len(names[1])):
        print(names[1][i])
        
        #creation of a tuple containing three lists:
            #i) all the PDBids that contain the domains where the mutations are
            #ii) if the PDB is a model of a mutated or WT protein
            #iii) the found mutational sites, not all mutations may be found.
        pdb_data = find_sequence(names[1][i],csv_[i], mutations)
        
        uniprot_name_list = list([names[1][i]]*len(pdb_data[0]))
        
        d = {'uniprot_id': uniprot_name_list, 'PDB_id': pdb_data[0], 'WT_status': pdb_data[1], "Exp_method": pdb_data[3], "Resolution": pdb_data[4], "Found_sites": pdb_data[2]}
        
        df = pd.DataFrame(data=d)

        df.to_csv(names[1][i]+"_PDB_Uniprot_df.csv")
            
    return
