
"""
This file contains the code required to answer
the question: 

Does the found structure contain the sequence where the mutations are?
And what is the AA at that position?

last update: 24th aug 2021

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
from Bio.PDB import *
from Bio import SeqIO

###############################################################################
#Preperation functions start
###############################################################################

def import_names(Path):
    """Function that takes the path and returns basic information regarding the 
    files found in the path.
    
    Parameters
    --------------
    Path:       A string indicating the placement of the .csv files to import.
    
    Returns
    -------
    names_list: A list of strings, all he PDBs related to a single Uniprot ID
    uniprot_id: A list of strings, with the uniprot IDs to which there are 
                solved structures.
    
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
    containing the list of pdb sturctures
    
    Parameters
    ------------
    Path:              A string indicating the placement of the .csvs found in 
                       import_names(Path)
    names_list:        A list of strings of all the .csvs found in 
                       import_names(Path)
    Uniprot_id:        A list of strings of all the uniprot_ids covered by 
                       the names_list.
    
    Returns
    -----------
    dataframes_list:   A list of pandas dataframes. dataframes_list[0] concern
                       uniport_id[0]
                       
                       The dataframe contains; (Example)
                           
    PDB_id | experimental_method | resolution | deposition_date | model_mutations
    ------------------------------------------------------------------------------
    2XN6   | X-RAY DIFFRACTION   |     1.29   |   2010-07-30    |     Mutated
    
    
    """
    
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
    values
    
    Parameters
    -----------
    input_dataframe:  The pandas dataframe used as input for the entire pipeline.
        
    Returns
    -----------
    dictionary:       A dictionary with the uniprot ids as keys and the 
                      mutational sites as values. 
    """
    
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

###############################################################################
#Preperation functions end
###############################################################################

###############################################################################
#Biological & Individual Calutation functions start
###############################################################################

###############################################################################
#Slimfast Function

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
 
def find_ranges(segments):
    """The function takes the segments found in the quality control function
       get_uniprot_segments() and removes sequences that does not match.
    
    Parameters
    ------------
    segments: A dictionary with chains as keys and sequence numbers as values.
    
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
            
    Resturns
    --------------------------------------------------------------------------
    The output is an array of the corrected segment. This is a quite stringent 
    code, that removes everything, that contains differencies. For the future, 
    this might need to be more dynamic. 
    
    chain_list: A list of integers, [0, 1, 2] in the length of the number of 
                chains.
    
    collected segments: A list in the length of chain_list containing 
                        np.arrays for each segment of a chain; example:
                            
                        [[array([ 94, 291]), array([322, 356])],
                         [array([ 94, 291]), array([322, 356])]])
    """
    
    #Convert the values of the input dictionary to an array    
    segment_array = np.array(list(segments.values()))
    
    ###########################    
    #    Control of chains    # 
    ###########################
    
    #create empty list for removal:
    to_be_removed = []
    
    #control if any chains are not covered
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
    
    ###################################    
    # Control of Segment within chain #
    ###################################    
    
    #For the remaining segments, a control of the length is to be conducted
    chain_list = []
    collected_segments = []
    
    #Length of the segment is controlled to be identical and is either 
    #discarded or the uniprot numbering is kept.
        
    for chain in range(len(segment_array)):
        
        statement = []
        _segment_ = []
        
        for segment in range(len(segment_array[chain])):
            
            n_uniprot = 1
            n_pdb = 0
            
            if segment_array[chain][segment][n_uniprot][1]-segment_array[chain][segment][n_uniprot][0] == segment_array[chain][segment][n_pdb][1]-segment_array[chain][segment][n_pdb][0]:
                
                segmentcontent = True
                
                statement.append(segmentcontent)
            
            if False not in statement:
                
                chain_list.append(chain)
                
                segment_ = segment_array[chain][segment][1]
                
                _segment_.append(segment_)
            
        collected_segments.append(_segment_)    

    return chain_list, collected_segments            

#        if segment_array[chain][0][1][1]-segment_array[chain][0][1][0] == segment_array[chain][0][0][1]-segment_array[chain][0][0][0] and segment_array[chain][1][1][1]-segment_array[chain][1][1][0] == segment_array[chain][1][0][1]-segment_array[chain][1][0][0]:
            # for the true statements, this should be saved
            #save ii to capture chain information    


def make_searchable_ranges(residue_range):
    """The function takes the output range from find_ranges() and
    convert this to an array of AA positions that can be searched.
    
    Parameters
    --------------
    residue_range:     A list containing a number of arrays for each
                       segment in a chain of the protein; eg. 
                       
                       residue_range = [[array([ 94, 291]), array([322, 356])]]
                       
                       Where residue_range[0] indicate the first segment, which
                       is a list consisting of two np.arrays.
    
    Returns
    -------
    ranges_:           A singular array of intergers indicating the numerical 
                       sequence of a particular chain.
    
    """
    
    #create an empty list
    a_list = []
    
    #create a for loop to convert each residue range into an array
    for jj in range(len(residue_range)):
        
        ranges_ = np.arange(residue_range[jj][0], residue_range[jj][1]+1)
        
        #append each range
        a_list.append(ranges_)
    
    ranges_ = np.hstack(a_list)
    
    return ranges_



def find_mutational_sites(mutation_location, range_of_chain):
    """A function that takes each mutational location and the range of 
    the PDB model and controls if the PDB model covers the mutational site. 
    This exist to narrow down the searchpool before importing sequences.
    
    Parameters
    --------------
    mutation_location:       An array of intergers indicating the mutational 
                             location, e.g. np.array([278, 337, 344])
    range_of_chain:          An array of intergers indicating the numerical 
                             sequence as calculated in the function 
                             "make_searchable_ranges()"
    
    Returns
    -------
    presence_of_mutations:   List of intergers indicating the 
                             mutation_location items present in range_of_model.
    
    """
    
    #empty list of collect the presence of mutations iteratively
    presence_of_mutations = []

    #for loop over the mutations to decide on a binary scale if the 
    # mutational location exist in the PDB file. 
    for mutation in range(len(mutation_location)):
        
        if mutation_location[mutation] in range_of_chain:
            r = mutation_location[mutation]
            
            presence_of_mutations.append(r)
    
    return presence_of_mutations


def get_sequence(pdb_id, path):
    
    """ 
    A function that takes a pdb id and the path to where it sould be downloaded 
    and find the sequence and chains.
    
    Parameters
    --------------
    pdb_id:     A string corresponding to a pdb id, e.g. "3Q01"
    path:       A string indicating the path where the pdb file 
                should be downloaded to.
    
    Returns
    -------
    A tuple containing;
    
    chains:     A list containing strings, the name of the chains, 
                e.g. ['A', 'B']
    sequences:  A list containing srings, the Amino Acid sequence, 
                e.g., [Seq('SSSVP....'), Seq('FRLHG...')] 
    
    """
    
    pdb = PDBList()
    
    sequences = []
    chains = []
            
    print(pdb_id)
    pdb.retrieve_pdb_file(pdb_id)
        
    pdb_ = pdb_id[1:3]
    pdb_ = pdb_.lower()
    os.chdir(pdb_)
        
    for record in SeqIO.parse(pdb_id+".cif", "cif-seqres"):
        print('> chain_'+record.annotations['chain']+'\n'+record.seq)
        chains.append(record.annotations['chain'])
        sequences.append(record.seq)
        
    os.chdir(path)
        
    return chains, sequences

def AA_finder(sequence, range_, mutation_list):
    """
    A function that finds the specific Amino Acid placed on the mutational 
    site.

    Parameters
    ----------
    sequence :      A sting of Amino Acids
    range_ :        A range of numbers corresponding to the sequence positions
    mutation_list : List of mutational sites, integers, to seach for

    Returns
    -------
    AA_list :       List of Positions and the corresponding Amino Acid.

    """
    
    #empty list to collect the Amino Acids in the correct position
    AA_list = []
    
    #for loop to look for all the mutations
    for mutation in range(len(mutation_list)):
        
        #find the numerical index of the sequential placement
        placement_index = np.where(range_ == mutation_list[mutation])
        
        #convert to an integer
        mut_placement_ = int(placement_index[0])
        
        #find the placement in the Amino Acid sequence
        AA = sequence[mut_placement_]
        
        #combine AA and placement
        
        Mut = str(mutation_list[mutation]) + AA
        
        #append AA to the list
        AA_list.append(Mut)
    
    return AA_list

def find_numbering(pdb_id, ranges_, mutations, path):
   
    """
    A function that takes the numerical input and pdb_id and convert
    this to specific amino acids on the mutational sites

    Parameters
    ----------
    pdb_id:    A string, The string of the pdb_id in question, e.g. '3Q01'
    ranges_:   A list of arrays, The consolidated sequence range for each chain
    mutations: An array of integers, the mutational sites
    path:      A string, the path to which the pdb should be downloaded 

    Returns
    -------
    A tuple containing a variable number of elements following the form:
        
        ('chain','Pos_AA', 
         'A', ['278P', '337R', '344R'], 
         'B', ['278P', '337R', '344R'])
    
    for a PDB model with chains A and B

    """
    
    #Convert the segment descriptors [start, stop] into searchable ranges
    #and find the 

    AA_tuple = ('chain','Pos_AA')
    
    sequence = get_sequence(pdb_id, path)  
    #output chains, sequences
    
    for i in range(len(ranges_)):
        chain = sequence[0][i]
        
        searchable_range = make_searchable_ranges(ranges_[i]) 
        
        #Seach through the ranges for the relevant mutational site 
        placement = find_mutational_sites(mutations, searchable_range)
        
        #find the amino acid at the mutational site
        amino_acid_at_site = AA_finder(sequence[1][i], searchable_range, placement)
        
        inf = (chain, amino_acid_at_site)
    
        AA_tuple += (inf)
        
    return AA_tuple

###############################################################################
#Biological & Individual Calutation functions end
###############################################################################

###############################################################################
#Wrapper functions start
###############################################################################

def collect_range_information(name,csv):
    
    """
    A function that takes a Uniprot ID, the corresponding CSV file and the mutations and search 
    trough the PDBs in the CSV file to find the matching values. 
    
    
    Parameters
    ---------------
    name:           A string indicating a uniprot id, e.g. 'P05543'
    csv:            A dataframe following the below structure:
    
    Uniprot_id | experimental_method | resolution | deposition_date | model_mutations
    ----------------------------------------------------------------------------------
    2XN6       | X-RAY DIFFRACTION   |     1.29   |   2010-07-30    |     Mutated
   
    
    Returns
    -------------
    pdb_id_success: A list of stings, each pdb ID which lives up to 
                    the quality contol parameters.
    wt_list:        A list of strings, each indicating if the pdb model in 
                    pdb_id_success is wildtype or mutated.
    method:         list of strings, each indicating which experimental model 
                    was used to solve the structure in pdb_id_success.
    chain_list:     A list of intergers indicating the chains in the pdb_model
    ranges_list:    A list of arrays indicating the ranges covered by the 
                    model.
    
    """
    
    #create a list of pdb ids
    r_list = list(csv[name])
    
    #create empty lists for collection of information    
    pdb_id_success = []
    chain_list = []
    ranges_list = []
    wt_list = []
    method = []
    res = []
        
    for j in range(len(r_list)):
        print(r_list[j] + " is being analyzed")
        #search for the segments covered by the model in PDB
        try:
            r = get_uniprot_segments(r_list[j], name)
            
            #Clean up the segments, removing all the sequence elements of the PDB, for which 
            #there is not a model
            res_range = find_ranges(r)
                
            #create a dictionary of each PDBid and the range.
            #count number of chains
            if len(res_range[0]) > 0:
                
                print(r_list[j] + ' saved for further investigation')
                #keep the pdb id thatwas approved
                pdb_id_success.append(r_list[j])
                
                #capture the relevant experiemental data regarding the pdb_id
                wt_list.append(csv.loc[csv[name]==r_list[j],"model_mutations"].iloc[0])
                method.append(csv.loc[csv[name]==r_list[j],"experimental_method"].iloc[0])
                res.append(csv.loc[csv[name]==r_list[j],"resolution"].iloc[0])
                
                #capture the approved range
                chain_list.append(res_range[0])
                ranges_list.append(res_range[1])
                
        except:
            print(r_list[j] + " does not pass the Uniprot quality control")
    
    return pdb_id_success, wt_list, method, chain_list, ranges_list, res


def validate_sequence(path, input_df):
    
    """
    A wrapper function to validate the sequence of each of the PDBs
    found in "find_structure.py"
    
    The input is:
        
        path = the path where the input file is.
        input_df = the input file
        
    The output is:
        
        a csv file per uniprot_id in the following format:
            
        uniprot_id | PDB_id | WT_status | Exp_method | Resolution | Found_sites | Range         
        
    """
    amino_acids = []
    path = path + "/structure_lists"
    
    #import all the Uniprot ID's to which pdb IDs were found
    names = import_names(path)
    
    #Import the csv files of each uniprot id's PDBs as a list of dataframes
    csv_ = import_csvs(path, names[0], names[1])
    
    #Find all the mutational positions.
    mutations = seq_pos_dict(input_df)
    
    for uniprot_id in range(len(names[1])):
        
        mut = mutations[names[1][uniprot_id]]
        
        print(names[1][uniprot_id])
        
        #creation of a tuple containing three lists:
            #i) all the PDBids that contain the domains where the mutations are
            #ii) if the PDB is a model of a mutated or WT protein
            #iii) the found mutational sites, not all mutations may be found.
        pdb_data = collect_range_information(names[1][uniprot_id],csv_[uniprot_id])
        
        for pdb in range(len(pdb_data[0])):
        
            #find_numbering(pdb_id, ranges_, mutations, path)
            AA = find_numbering(pdb_data[0][pdb], pdb_data[4][pdb], mut, path)
            amino_acids.append(AA)
        
            d = {'PDB_id': pdb_data[0], 'WT_status': pdb_data[1], 
                 "Exp_method": pdb_data[2], "Resolution": pdb_data[5], "Found_sites": amino_acids}
        
            df = pd.DataFrame(data=d)

            df.to_csv(names[1][uniprot_id]+"_PDB_Uniprot_df.csv")
            
    return

###############################################################################
#Wrapper functions end
###############################################################################
