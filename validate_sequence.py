#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WORK UNDER CONSTRUCTION
"""


#Create a dictionary with the uniprot ID and mutational sites.

def seq_pos_dict(input_dataframe):
    """returns a dictionary with the uniprot ID as key and all the positions of the mutations as values"""
    
    #make a list of all uniprotIDs
    IDs = list(input_dataframe.Uniprot)
    
    #iterate over all the mutations, to extract the position of the mutation
    for i in range(len(input_dataframe)):
        
        mut_pos = re.split('(\d+)',input_dataframe.Mutations[i])
        mut_pos = [int(word) for sublist in map(str.split, mut_pos) for word in sublist if word.isdigit()]
        
        mut.append(mut_pos)
    
    #collect the information in a dictionary
    dictionary = {IDs[i]: mut[i] for i in range(len(IDs))}
        
    return dictionary


def imoport_all_pdb_lists(Path, input_dataframe):
    
    IDs = list(input_dataframe.Uniprot) 
    dataframes_list = []
    
    for i in range(len(IDs)):
        try:
            temp_df = pd.read_csv(Path+"/"+IDs[i])
            dataframes_list.append(temp_df)
        
        except:
            print("a PDB file does not exist for %s" % IDs[i])
            
    return dataframes_list

Path = '/Users/krde/Documents/Research/TCGA 3D'
pdb_lists = imoport_all_pdb_lists(Path, input_df)

def imoport_all_pdb_lists(Path, input_dataframe):
    
    IDs = list(input_dataframe.Uniprot) 
    dataframes_list = []
    
    for i in range(len(IDs)):
        try:
            temp_df = pd.read_csv(Path+"/"+IDs[i])
            dataframes_list.append(temp_df)
        
        except:
            print("a PDB file does not exist for %s" % IDs[i])
            
    return dataframes_list