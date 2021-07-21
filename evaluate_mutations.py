#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os
from Bio.PDB import *
from Bio import SeqIO


def import_csv(Path, pdb_uniprot_df):
    """Returns..."""

    #ensure correct placement
    os.chdir(Path)
    
    pdb_df = pd.read_csv(pdb_uniprot_df)
    
    pdb_df.drop('Unnamed: 0', axis='columns', inplace=True)
    
    return pdb_df

def merge_information(start_df, pdb_df):
    """Returns..."""
    
    #start by removing the uniprots with no related PDBs 
    new_df = start_df[start_df["Uniprot"].isin(list(pdb_df["Uniprot"]))]
    
    #join together based on the Uniprot id
    df = new_df.join(pdb_df.set_index("Uniprot"),on="Uniprot")
    
    return df

def get_sequence(df, path):
    
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


######## there is something here that does not add up biologically. 
# unfortunately I think the problem is in the way i have constructed the ranges. 

