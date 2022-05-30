##############################################
###		PDBminer		   ###
##############################################

#This read me describes how to run this interim version 
#of PDBminer, last update 21May2022 

#NB! Only structures reported in 
#clean_{unipot_id}_structural_df.csv should be 
#Considered for further analysis. The next step is
#to check the quality of the structure. This step 
#May be implemented later. 

##############################################
# Environment
##############################################

#Everytime:
module load conda/4.9.2/modulefile

#First time:
conda env create -f environment_python.yml 
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0

#all subsequent times
conda activate PDBminer

##############################################
# Input Files & Scripts
##############################################

#config.yaml
        #input file
        #input path

#scripts
#PDBminer comprises a two python scripts both present in the directory "scripts", 

        # PDBminer_run.py
        # PDBminer_functions.py

#Input file:

#The input file have to be a .csv file format and include the following comma
#separated values: Hugo_name, Uniprot, Uniprot-isoform,Mutations, ClusterID.

#         Hugo_name       String, e.g., AKT1
#         Uniprot         Uniprot ID e.g., P31749 
#         Uniprot-isoform Integer
#         Mutations       strings separated by a semicolon, ;, e.g. E17K or L755S;K128M, deletion with *
#         ClusterID       Integer or N/A.

##############################################
# Running
##############################################

tsp -N 4 snakemake --cores 4   

##############################################
# Output
##############################################

# Each Uniprot ID will have its own directory 
# containing: 

#	1) the input {unipot_id}_input.csv
#	2) a common file indicating that the calculations 
#	   have finished: {unipot_id}_done.txt

#Moreover, there will be varying other files: 

#	3) all_{unipot_id}_structural_df.csv
#	   An output file with all PDBs associated 
#	   with the uniprot_id

#	4) clean_{unipot_id}_structural_df.csv
#	   An output file with the PDBs associated 
#	   with the uniprot_id that covers at least 
#	   one mutation.

#	5) missing_id.txt
#	   If there are no PDBids associated with the 
#	   uniprot id. This should become 
#	   an alphafold api.

#	6) alphafold.txt
#	   If there are PDBs associated but none of them 
#	   covers a single mutation. This should become 
#	   an alphafold api.

#	7) issue.log

#content of clean_{unipot_id}_structural_df.csv and 
#all_{unipot_id}_structural_df.csv:

#Output Columns and explanations
#================================================================
#Hugo_name                      Gene names from the input file
#Uniprot_id                     ID 
#Uniprot_isoform                Number identifying which isoform to which the alignment was done. 
#Cluster_ID                     Number of cluster from the input file. 
#PDB_id                         Identifier of the PDB file.
#Research_mutations             The mutations specified in the input file
#experimental_method            By which the PDB was generated.
#resolution                     Estimation of PDB quality
#deposition_date                Timing of file placement in PDB
#chains                         Letter describing the chains covering the Uniprot ID 
#coverage                       Range where the PDB file and uniprot ID are aligned. “;” separated for multiple chains.
#AA_in_PDB                      The amino acids at the mutational spots described in the input file.
#mutations_in_pdb               All mutations found in the pdb file compared to the uniprot sequence if none: []
#complex_protein                Binary, either “NA” or “Protein in complex” indicates if the PDB file contains a protein complex.
#complex_protein_details        Details regarding the protein complex indicating the Uniprot ID of the other protein and the chains.
#complex_nucleotide             Binary, indicates if the protein is bound to a nucleotide string such as DNA.
#complex_nucleotide_details     Details regarding the DNA or RNA binding. 
#complex_ligand                 Binary indicating if metal or small molecules are present in the pdb file.
#complex_ligand_details         Details describing which “other” things are in the file.
