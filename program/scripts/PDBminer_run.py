#!/usr/bin/env python

# PDBminer_run: classes and functions for the snakefile.
# Copyright (C) 2022, Kristine Degn
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import pandas as pd
import shutil

from PDBminer_functions import find_structure_list, combine_structure_dfs, align, collect_complex_info, cleanup_all, filter_all

def run_list(full_path):
    """
    Run list runs PDBminer for a unique uniprot id. It relies on the fundtions
    in the script PDBminer_functions. 
    
    Functionality: 
        
        1)  reads the input dataframe for the specific Uniprot id
        2)  find all structures related to the uniprot id 
        3)  combine the input 1) and the structures 2)
        4)  For each sequence of structure (PDBid) an alignment to the fasta of 
            the specified isoform. Here the differencies between the PDB and
            fasta sequence is identified, the area of the sequence covered
            by the PDB is annotated and the amino acids at the mutational sites
            are found. NB! this does not account for quality of the structure. 
        5)  The PDB files are then analyzed in terms of other present proteins, 
            indicating a complex, ligands and other molecules present. 
        6)  All these informations is reported in all_{uniprot_id}_structural_df.csv
        7)  A cleanup of the path removing structures and if structures.
        8)  The structural_df if cleaned only keeping the PDB files that at 
            least cover one of the specifed mutations. This is reported as
            clean_{uniprot_id}_structural_df.csv. If there is no structures 
            that cover the mutations, a txt file is reported indicating 
            that an alphafold structure may be the next path. 
            
    If these steps cannot take place for one reason or another, an issue log 
    will be reported.

    Parameters
    ----------
    uniprot_list : The list of unique uniprot ids as output of import_csv

    Returns
    -------
    None. The function creates files and directories running the scripts
    but does not have any output variables itself. 
    
    """

    full_path = str(full_path)
    
    path = full_path.split("/results")[0]
    uniprot_id = full_path.split("/")[-2]

    os.chdir(f"{path}/results/{uniprot_id}")

    input_dataframe = pd.read_csv(f"{uniprot_id}_input.csv", index_col=0)
    print(input_dataframe)
    #only one of each 
    
    #prior issue with one uniprot id, still to be fixed this is a patch
    found_structures = find_structure_list(input_dataframe)
    print("structures found")
    
    if len(found_structures) != 0:
        
        combined_structure = combine_structure_dfs(found_structures, input_dataframe)
        print("files combined")
        combined_structure = collect_complex_info(combined_structure)
        print("complexes found")
        structural_df = align(combined_structure, os.getcwd())
        print("files downloaded and aligned")
        
        if len(structural_df) != 0:
            #prep and export all file
            all_df = cleanup_all(structural_df)
            print("filed cleaned")
            
            if type(input_dataframe.mutations[0]) != str: 
                all_df = all_df.drop(columns=('mutations'))
            if set(input_dataframe.cluster_id) == {999}:
                all_df = all_df.drop(columns=('cluster_id'))
            
            all_df.to_csv(f"{path}/results/{uniprot_id}/{uniprot_id}_all.csv")
            
            #prep and export the filered file
            if input_dataframe['mutations'].isnull().values.any() == False:
                if input_dataframe['mutations'][0] not in ["N/A", "NA"]:
                    filtered_df = filter_all(structural_df, input_dataframe)   
            
                if len(filtered_df) == 1:            
                    if len(filtered_df[0]) > 0:
                        if set(input_dataframe.cluster_id) == {999}:
                            filtered_df[0] = filtered_df[0].drop(columns=('cluster_id'))
                        filtered_df[0].to_csv(f"{path}/results/{uniprot_id}/{uniprot_id}_filtered.csv")
                else:
                    for cluster_df in filtered_df:
                        if len(cluster_df) > 0:
                            cluster_df.to_csv(f"{path}/results/{uniprot_id}/{uniprot_id}_cluster{cluster_df.cluster_id[0]}_filtered.csv")
                
            #remove the alphafold model
            os.chdir(f"{path}/results/{uniprot_id}")
            os.system("find . -maxdepth 1 -name '*.pdb*' -type f -delete")
            
            if os.path.exists("structure"):
                shutil.rmtree("structure")
            if os.path.exists("obsolete"):
                shutil.rmtree("obsolete")
                
            if os.path.getsize("missing_ID.txt") == 0:
                os.remove("missing_ID.txt")            
    
    with open(f"{uniprot_id}_done.txt", "w") as textfile: 
        textfile.write(f"Analysis complete for {uniprot_id}")
          
    os.chdir(f"{path}/results/")
  
    return

run_list(snakemake.input)
