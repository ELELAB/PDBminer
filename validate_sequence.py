##############
# Import packs
##############

import pandas as pd
import numpy as np
import os
from Bio.PDB import *
from Bio import SeqIO
import shutil


#############

from SLiM_codes import get_uniprot_segments

#############

def identify_mutational_positions(segments, mutation_list, pdb_id, path):
    
    #This function is too strict in its quality control and should be reworked somehow. 
    #Idea is to introduce alignment rather than segmentation. 
    
    """The function takes the segments found in the quality control function
       get_uniprot_segments(), removes sequences that does not match, identifty mutational sites
       and find the related amino acids.
    
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
    
    
    mutation_list: A list of intergers of the mutational sites, derived from the 
                    input datafile. 
                    
    pdb_id:       The PDB id to which the segments belong.
            
    Resturns
    --------------------------------------------------------------------------
    The output is a dictionary of the amino acids at the mutational sites if 
    covered by the structure. 
    
    chain_dict format is either a dictionary with each chain as key; e.g.:
    
    {'A': ['P', 'X', 'X'], 'B': ['P', 'X', 'X']}
    
    Values are either the Amino acid at the mutational site or "X" indicating that the 
    model does not cover that particular mutation. 
    
    The output can also be "0" which indicates that the model is not sufficient 
    for the tool to use, either because the uniprot sequence and the PDB sequence
    does not match (quality control) or because it does not cover any of the 
    mutational sites the function looks for.
        
    """
    
    ##########################
    # Step 1: Clean segments #
    ##########################
    
    
    #looping through the segments dict
    for k,v in list(segments.items()):
        
        #Convert the values of the input dictionary to an array 
        segment_array = np.array(v)

        #Remove the chains of the model where the pdb file does not
        #cover the begining or the end of the segment, e.g. due to
        #disorderd regions or other model construction issues.        
        if None in segment_array.flatten():
            segments.pop(k)
            continue
            
        #For the remaining segments, a control of the length is 
        #carried out - the length of the uniprot chain and the pdb
        #chain have to match exactly. 
    
        length_of_segment = segment_array[:,:,1] - segment_array[:,:,0]
    
        if not np.all((length_of_segment[:,0] - length_of_segment[:,1])) == 0:
            segments.pop(k)
            
            continue
    
        segments[k] = segments[k][0][1]
    
    #Rename purely for human comprehension/for my own understanding 
    clean_segments = segments
    del segments

    ############################################################
    # Step 2: search through the segments for mutational sites #
    ############################################################
    
    #loop trough the clean segments converting them to searchable numerical ranges
    for k,v in list(clean_segments.items()):
        
        #convert the values to an array
        segment_array = np.array(v)
        
        #Expand array to cover the entire range of amino acids
        segment_range = np.arange(segment_array[0], segment_array[1]+1)
        
        #Introducing an empty list for collection of mutational placements in the sequence
        placement_list = []
        
        #loop over each mutation, and find out where it is placed in the string:
        
        #e.g. mutation_list = [50, 57, 101], mutation_list[mutation] = 50
        # segment_range = [45, 46, 47, 48, 49, 50, 51, 52 ...]
        
        #placement of mutation 1, at AA 50 is: #5 in segment_range, hence, we need to look at 
        # AA 5 in sting of Amino acids; e.g., seq = [SSRTYGG], then we are for G, that is 
        # Amino acid substituted at mutation_list[mutation] = 50
        
        for mutation in range(len(mutation_list)):
            
            if mutation_list[mutation] in segment_range:
                
                #find the numerical index of the sequential placement
                placement_index = np.where(segment_range == mutation_list[mutation])
                
                #convert to an integer
                placement_index = int(placement_index[0])
                
                #append to list
                placement_list.append(placement_index)
            
            #If the mutation is not covered by the model, the "X" is used to indicate this
            else: 
                placement_list.append("X")
            
            #reforming the segment dictionary for each amino acid
            clean_segments[k] = placement_list   
       
        # This should not be there???? how can this situation happen?
        if clean_segments[k] == []:
            clean_segments.pop(k)
            continue #?
    
    #renaming for comprehension
    mutational_site_dict = clean_segments
    del clean_segments
    # eg: mutational_site_dict = {'A': [239, 298, 305]}, or {'B': ['X', 276, 'X']}, or {}
    
    ### retrieve amino acid sequence if the dictionary is not empty 
    
    #######################################
    # Step 3: Find the Amino Acid at site #
    #######################################
    
    #removal of empty dictionaries
    if mutational_site_dict == {}: 
        
        return "Mutational Domains Not found"     
    
    #Retrieval of structures
    else: 
        
        pdb = PDBList()

        sequences = []
        chains = []

        pdb.retrieve_pdb_file(pdb_id)

        pdb_name = pdb_id[1:3]
        pdb_name = pdb_name.lower()
        os.chdir(pdb_name) #this is not a good solution, but I end up creating a russian doll of structures if I dont. 

        for record in SeqIO.parse(pdb_id+".cif", "cif-seqres"):
            chains.append(record.annotations['chain'])
            sequences.append(record.seq)
        
        os.chdir(path) # this is needed as long as os.chdir(pdb_name) is in use
        
	#Create an empty dictionary to capture the amino acids
        chain_dict = {}
        
        #loop over the doctionary containing the mutational sites. 
        for k,v in list(mutational_site_dict.items()):

            pl = [] #empty list to capture the amino acids or "X"

            for chain in range(len(chains)):

                if k == chains[chain]: #ensureing that found chain segments are retained. 

                    for position in v:                

                        try: 
                            amino_acid = sequences[chain][position]
                            pl.append(amino_acid)

                        except:
                            position = 'X'
                            pl.append(position)
                            

            #If the chains does not match. 
            
            #example; a structrue was retrieved from in PDB based on the p53 Uniprot ID
            #The chains of the structure model that satisfy the numerical sequence of 
            #the mutational site, is chain E and F, but chain E and F is each the 
            #DNA strands that binds to the tetramer A, B, C, D. Hence, the chains does 
            #not match. 
            
            if pl == []:
                return "The chains does not match"
            
            #Capture amino acids unless it is just "X" for all indicating that the mutational sites where not found.
            elif len(set(pl)) > 1 or pl[0] != "X":
                chain_dict[k] = pl
    
       # shutil.move(path+"/"+pdb_name, path+"/structures") # now I move all the structures to a dir, but they should proberbly just be deleted.             
            
    return chain_dict


def validate_sequence(sturcture_df, input_dataframe, path):
    
    """
    A function that takes the dataframe created in find_structure_list
    and retrieves information regarding the mutational sites and out-
    put the same dataframe with additional columns, mutation_coverage and
    Research_mutations. 
    
    Parameters: 
    ----------------------------
    structure_df: dataframe created in find_structure_list
    
    input_dataframe: the input file
    
    Returns:
    ----------------------------
    A csv file & df for inspection of all the stuctures discarded
    
    A csv file & df containing the potential models. 
    
    """
    
    #creation of empty lists 
    mutation_coverage = []
    searching_mutations = []
    
    for i in range(len(input_dataframe)):
    
        for j in range(len(sturcture_df)):

            if sturcture_df.iloc[j]["Uniprot_ID"] == input_dataframe.Uniprot[i] and sturcture_df.iloc[j]["ClusterID"] == input_dataframe.ClusterID[i]: 

                try: 
                    segments = get_uniprot_segments(sturcture_df.iloc[j].name, input_dataframe.Uniprot[i])
                    mutation_list = input_dataframe["mutation_positions"][i]

                    outcome = identify_mutational_positions(segments, mutation_list, sturcture_df.iloc[j].name, path)
                    
                    if outcome == {}: #removing empty dictionaries
                        mutation_coverage.append("segments not found")
                        searching_mutations.append(input_dataframe["Mutations"][i])
                         
                    else: 
                        mutation_coverage.append(outcome)
                        searching_mutations.append(input_dataframe["Mutations"][i])
                   
                except:
                    "KeyError:" + input_dataframe.Uniprot[i]
                    mutation_coverage.append("segments not found")
                    searching_mutations.append(input_dataframe["Mutations"][i])
                

    sturcture_df.insert(5, "mutation_coverage", mutation_coverage, True) 
    sturcture_df.insert(6, "Research_mutations", searching_mutations, True)
    
    if len(set(sturcture_df.ClusterID)) == 1 and sturcture_df.ClusterID[0] == 999:
        
        #Return all clusters set at 999 to N/A  
        input_dataframe.ClusterID = "N/A"
    
    #Notice 
    # 1) The information on isoforms are not used here!
    # 2) PTM information should also be captured in some capacity!
    
    #Seperate out potential choices for structures and a list of structures to insepect
    binary = []

    #The last quality check that each coverage is added as a dictionary
    for i in range(len(sturcture_df)):
        if type(sturcture_df.mutation_coverage[i]) == dict and sturcture_df.mutation_coverage[i] != {}:

            binary.append(0)

        else:
            
            binary.append(1)

    sturcture_df.insert(0, "type_of_inf", binary, True)
    
    if len(set(sturcture_df['type_of_inf'])) > 1:
        
        df1, df2 = [x for _, x in sturcture_df.groupby(sturcture_df['type_of_inf'] > 0)]

        potential_structures = df1.drop("type_of_inf", axis=1)
        potential_structures.to_csv(path+'/potential_structures.csv')

        inspection_structures = df2.drop("type_of_inf", axis=1)
        inspection_structures.to_csv(path+'/inspection_structures.csv')
        
    else:
        
        if sturcture_df['type_of_inf'].iloc[0] == 0:
            potential_structures = sturcture_df.drop("type_of_inf", axis=1)
            potential_structures.to_csv(path+'/potential_structures.csv')
            
        else: 
            inspection_structures = sturcture_df.drop("type_of_inf", axis=1)
            inspection_structures.to_csv(path+'/inspection_structures.csv')

    return #potential_structures, inspection_structures #how do I deal with this when not both are there. 
