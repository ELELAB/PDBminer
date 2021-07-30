
"""
This file contains the code required to answer
the question: 

Does the found structure contain the mutations?

"""

##############
# Import packs
##############
import pandas as pd
import os
from Bio.PDB import *
from Bio import SeqIO

##############


def import_csv(Path, pdb_uniprot_df):
    """This function takes the path of and csv file created in the 
    validate_sequence.py script, the pdb_uniprot_df and returns it as 
    a dataframe."""

    #ensure correct placement
    os.chdir(Path)
    
    pdb_df = pd.read_csv(pdb_uniprot_df)
    
    pdb_df.drop('Unnamed: 0', axis='columns', inplace=True)
    
    return pdb_df

def merge_information(start_df, pdb_df):
    """This function takes pdb_uniprot_df (pdb_df) imported in import_csv()
    and merge the information with the input file to return a singular 
    dataframe containing all the relevant information"""
    
    #start by removing the uniprots with no related PDBs 
    new_df = start_df[start_df["Uniprot"].isin(list(pdb_df["Uniprot"]))]
    
    #join together based on the Uniprot id
    df = new_df.join(pdb_df.set_index("Uniprot"),on="Uniprot")
    
    return df

def get_sequence(df, path):
    
    """This function takes the dataframe created in merge_information() 
    as input.
    
    The aim is to import the sequence of amino acids from the PDB file and
    use the mutational position to check if the mutation is present in the 
    PDB file. However, the sequence I import below is not the same in length
    as the segment ranges I have used to find the domains. Therefore 
    there a discrepency that I don't have an answer as to why it is there. 
    
    Planed actions: 
        
        - Revisit the way I retrive the segments from slimfast, perhaps
            I have not really understood and therefore I have a problem. 
        - PDB versus PDBe 
    
    The output is supposed to be an additional column to the merged df
    which includes information about the specific mutation, if present 
    or not. 
    
    """
    
    pdb = PDBList()
    
    sequences = []
    chains = []
    
    for i in range(len(df)):
        
        _pdb = df["PDB"][i]
        print(_pdb)
        pdb.retrieve_pdb_file(df["PDB"][i])
        
        pdb_ = _pdb[1:3]
        pdb_ = pdb_.lower()
        os.chdir(pdb_)
        
        for record in SeqIO.parse(_pdb+".cif", "cif-seqres"):
            print('> chain_'+record.annotations['chain']+'\n'+record.seq)
            chains.append(record.annotations['chain'])
            sequences.append(record.seq)
        
        os.chdir(path)
    
    d = {'Chain':chains,'Sequence':sequences}
    seq_df = pd.DataFrame(d)
    
    return seq_df
