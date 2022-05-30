# PDBminer

## Introduction to the Program 
PDBminer is a snakemake pipeline that takes a singular input file containing information about a protein and its mutations and outputs an overview of the best possible structural models in the 
Protein data bank to cover the protein and its mutations.  

PDBminer is currently only applicable for solved structures and does therefore not find alphafold 
models.  

## Dependencies

It is recommended to create an environment to run PDBminer, all specifications for this is 
Described in environment_python.yml.

#First time:
conda env create -f environment_python.yml 
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0

#all subsequent times
conda activate PDBminer


Within the envrionment_python.yml it is the following packs are described: 

* python=3.8.8
* pandas=1.2.4
* requests=2.25.1

## Setup
PDBminer Requires a Input file in a csv file format, containing the following:

```
Hugo_name | Uniprot | Uniprot-isoform | Mutations
----------------------------------------------------------
TP53      |  P04637 |         2       | P278L;R337C;L344P
MAT1A     |  Q00266 |         1       | P30N;W300H
SSTR3     |  P05543 |         1       | T11S;C191S;R330L
SAMD4A    |  Q9UPU9 |         3       | L10R;I80A

        
```
This should be specified in the config.yaml file. 

This file is present in this repository under the name "inputfile.csv" and can be used as an example.

run_program.py is the main script that will guide the program and currently contains crude user prompting. Run from the terminal using: 

```
$ snakemake --cores X

#An example of the output can be found in results.

```
The output:
Each Uniprot ID will have its own directory containing: 

* {unipot_id}_input.csv, the input
* a common file indicating that the calculations have finished: {unipot_id}_done.txt

Moreover, there will be varying other files: 

* all_{unipot_id}_structural_df.csv, An output file with all PDBs associated with the uniprot_id
* clean_{unipot_id}_structural_df.csv, An output file with the PDBs associated with the uniprot_id that covers at least one mutation.
* missing_id.txt, If there are no PDBids associated with the uniprot id. This should become an alphafold api.
* alphafold.txt, If there are PDBs associated but none of them covers a single mutation. This should become an alphafold api.
* issue.log

#content of clean_{unipot_id}_structural_df.csv and 
#all_{unipot_id}_structural_df.csv:

```
Output Columns and explanations
================================================================
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

```
