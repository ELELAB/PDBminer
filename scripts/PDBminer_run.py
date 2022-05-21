import os
import pandas as pd
import shutil

from PDBminer_functions import find_structure_list, combine_structure_dfs, align, collect_complex_info, cleanup_df

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
    
    #prior issue with one uniprot id, still to be fixed this is a patch
    found_structures = find_structure_list(input_dataframe)
    
    if len(found_structures) != 0:
        
        combined_structure = combine_structure_dfs(found_structures, input_dataframe)

        structural_df = align(combined_structure, os.getcwd())
        
        if len(structural_df) != 0:

            structural_df = collect_complex_info(structural_df)
    
            structural_df.to_csv(f"all_{uniprot_id}_structural_df.csv")
            shutil.rmtree("structure")
            if os.path.getsize("missing_ID.txt") == 0:
                os.remove("missing_ID.txt")
    
            #regular clean up
            clean_df = cleanup_df(structural_df)
    
            if len(clean_df) > 0:
                clean_df.to_csv(f"clean_{uniprot_id}_structural_df.csv")
            
            else: 
                with open("alphafold.txt", "w") as textfile: 
                   textfile.write("use the alphafold model, API to come")
        else:
            with open("issue_log.txt", "w") as textfile: 
                textfile.write("structural df issue")
    
    else:
        with open("issue_log.txt", "w") as textfile: 
            textfile.write("No structures found")
    
    with open(f"{uniprot_id}_done.txt", "w") as textfile: 
        textfile.write(f"Analysis complete for {uniprot_id}")
            
    os.chdir(f"{path}/results/")
  
    return 

run_list(snakemake.input)
