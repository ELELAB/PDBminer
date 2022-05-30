import pandas as pd
import os
import numpy as np

configfile: 'config.yaml'

path = config['path']

#importing the file
df = pd.read_csv(path+"/"+config['input_file']) 
if len(set(df.ClusterID)) == 1 and df.ClusterID[0] == "N/A":
    df.ClusterID = 999

#creating the uniprot list
uniprot_list = list(set(df.Uniprot))

#creating the results directory
if os.path.exists(path+"/results") == False: 
    os.mkdir(path+"/results")                    

target = expand("{path}/results/{uniprot_id}/{uniprot_id}_done.txt", uniprot_id=uniprot_list, path=path)                   

rule all:
    input:
        target     

#works as intended
rule prepare_files_and_dirs:
    input:
        "{path}"
    output:
        "{path}/results/{uniprot_id}/{uniprot_id}_input.csv"
    run:
        os.chdir(f"{input}/results")
        
        for uniprot_id in uniprot_list:
            if os.path.exists(uniprot_id) == False:
                os.mkdir(uniprot_id)
                os.chdir(uniprot_id)
            else:
                os.chdir(uniprot_id)
            
            uniprot_dataframe = df[df.Uniprot == uniprot_id]        
            uniprot_dataframe = uniprot_dataframe.reset_index(drop=True)
            uniprot_dataframe.to_csv(f"{uniprot_id}_input.csv")
            os.chdir(f"{input}/results")

rule run_PDBminer:
    input:  
        "{path}/results/{uniprot_id}/{uniprot_id}_input.csv"
    output:    
        "{path}/results/{uniprot_id}/{uniprot_id}_done.txt"
    script:
        "scripts/PDBminer_run.py"
