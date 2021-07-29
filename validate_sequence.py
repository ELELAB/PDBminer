
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
    """Takes the input file as inpit and returns a dictionary with 
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
    """The function takes the segments found in get_uniprot_segments() 
    (SLIM) which is ranges where the uniprot ID is covered by a given PDB 
    file. find_ranges() removes the segments that are not covered by the 
    PDB and output an array of the AA range the specific PDB covers"""
    
    #output from get_uniprot_segments
    
    segment_array = np.array(list(segments.values()))
    
    #reshape to fit ranges
    segment_array = segment_array.reshape([-1,2])
    
    #create empty list for removal:
    to_be_removed = []
    
    #control if any ranges are not covered by the PDB
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
            r = 1
        
        else:
            r = 0
            
        presence_of_mutations.append(r)
    
    if 0 in presence_of_mutations:
        mut_result = False 
    
    else: 
        mut_result = True
    
    return mut_result, mutation_location

def find_mutational_sites_2nd(mutation_location, range_of_model, iteration):
    """A function that is not yet utilized - the aim is to have a function
    that can find mutational sites, if a specific model does cover two out 
    of three mutations, perhaps this solved PDB is better than a homolog.
    """
    
    #create an empty list to store the presence of mutations, presence = 1, missing = 0
    presence_of_mutations = []
    
    sites = []
    
    #a for loop to investigate if the mutational sites (locations) are present in the model
    for mutation in range(len(mutation_location)):
        
        #when the mutational sites are present or not present, this is stored in a binary list as 1 or 0. 
        if mutation_location[mutation] in range_of_model:
            r = 1
            #when a site has been found, the information is stored in case not all mutations can be found.
            sites.append(mutation_location[mutation])
        
        else:
            r = 0
            
        presence_of_mutations.append(r)

    #for loop to ensure, that even if not all mutations are found, the best model 
    #that satisfy most of the criteria is chosen  
    if sum(presence_of_mutations) >= len(presence_of_mutations)-iteration:
        mut_result = True

    else:
        mut_result = False

    return mut_result, sites


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
    
    """
    A wrapper function to validate the sequence of each of the PDBs
    found in "find_structure.py"
    
    The input is:
        
        path = the path where the input file is.
        input_df = the input file
        
    The output is:
        
        pdb_uniprot_df = a csv file in the following format:
            
        Uniprot | PDB | found_mutational_sites         
        
    """
    
    
    #create an empty list to hold the found pdb id's
    pdb_list = []
    
    #create an empty list to hold the found pdb id's
    mutational_sites = []
    
    #create an empty list for the Uniprot IDs where the pdbs does not cover the mutational sites
    # uniprot_missing_mut_site = []
    
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
            if placement[0] == True:
                pdb_list.append(r_list[j]) 
                mutational_sites.append(placement[1])
                break
            
            else:
                continue
                
#   For later;          
#            if placement[0] == False:
#                placement_ = find_mutational_sites_2nd(mutations[names[1][i]], searchable_range, j)
#                
#                if placement_[0] == True:
#                    pdb_list.append(r_list[j]) 
#                    mutational_sites.append(placement_[1])
#                    break
                    
#                else:
#                    continue
            
#            if placement[0] == False:
#                uniprot_missing_mut_site.append(r_list[j])
                 
    #write into a file
#    textfile = open("missing_mutational_sites.txt", "w")
#    for element in uniprot_missing_mut_site:
#        textfile.write(element + "\n")
#        textfile.close()
    
    pdb_uniprot_df = pd.DataFrame(list(zip(names[1], pdb_list, mutational_sites)), 
                                  columns =['Uniprot', 'PDB', 'found_mutational_sites'])         
    
    pdb_uniprot_df.to_csv("pdb_uniprot_df.csv")
            
    return(pdb_uniprot_df)  

