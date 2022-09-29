# PDBminer

## Introduction to the Program 
PDBminer is a snakemake pipeline that uses a singular file as input. The file must contain information about a protein of interest and its mutations. PDBminer outputs an overview of the best possible structural models in the Protein Data Bank.

PDBminer find and annotate both structures from the Protein Data Bank and Alphafold.

## Dependencies

It is recommended to create a virtual environment to run PDBminer. The environment can be created using environment_python.yml as described below:

### First time:

```
conda env create -f environment_python.yml 
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0
conda install -c conda-forge biopandas=0.4.1
```

### All subsequent times

```
conda activate PDBminer
```


PDBminer is dependent on the following packs (as described in envrionment_python.yml):

* python=3.8.8
* pandas=1.2.4
* requests=2.25.1

## Setup
PDBminer requires an input file in a comma-separated values file format, with four mandatory columns; "Hugo_name", "Uniprot", "Uniprot-isoform", and "Mutations". It is optional to include a fifth column titled "ClusterID". The cluster ID should be an integer and can be used to annotate mutations forming a spacial cluster. See "inputfile.csv" as an example. 

```
Hugo_name | Uniprot | Uniprot-isoform | Mutations         | ClusterID
-----------------------------------------------------------------------
TP53      |  P04637 |         2       | P278L;R337C;L344P | 1
MAT1A     |  Q00266 |         1       | P30N;W300H        | 1
SSTR3     |  P05543 |         1       | T11S;C191S;R330L  | 1
SAMD4A    |  Q9UPU9 |         3       | L10R;I80A         | 1

        
```
The name of the input file should be specified in the config.yaml file. 

## Running the Program
Once the config.yaml has been updated and the dependencies installed, the program can be run from the terminal:
```
$ snakemake --cores X
```

## The Output
For each Uniprot ID in the input file, a directory is created. After a successful run, this directory will contain the following: 

* {unipot_id}_input.csv, the input
* a common file indicating that the calculations have finished: {unipot_id}_done.txt

Moreover, there will be varying other files: 

* all_{unipot_id}_structural_df.csv, An output file with all PDBs and alphafold structure associated with the uniprot_id
* clean_{unipot_id}_structural_df.csv, An output file with the PDBs and alphafold structure associated with the uniprot_id that covers at least one mutation.
* missing_id.txt, If there are no structures from the protein data bank or alphafold structures available for the uniprot_id

See an example of the output of the example input file in the directory results.

#content of clean_{unipot_id}_structural_df.csv and 
#all_{unipot_id}_structural_df.csv:

```
Output Columns and explanations
================================================================
#Hugo_name                      Gene names from the input file
#Uniprot_id                     ID 
#Uniprot_isoform                Number identifying which isoform to which the alignment was done. 
#Cluster_ID                     Number of cluster from the input file. 
#structure_id                   Identifier of the PDB file or Alphafold Model.
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
