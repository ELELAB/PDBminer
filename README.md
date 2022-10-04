# PDBminer

## Introduction to the Program 
PDBminer is a snakemake pipeline that uses a single inputfile. The file must contain information about a protein of interest and its mutations. 
PDBminer outputs a ranked overview of the possible structural models in the Protein Data Bank and the most current version of the Alphafold2 model, if any.

## Dependencies

It is recommended to create a virtual environment to run PDBminer. The environment can be created using environment_python.yml as described below:

### First time:

```
conda env create -f program/environment_python.yml 
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0
conda install -c conda-forge biopandas=0.4.1
```

### All subsequent times

```
conda activate PDBminer
```

PDBminer is dependent on the following packs (as described in program/envrionment_python.yml):

* python=3.8.8
* pandas=1.2.4
* requests=2.25.1

## Setup
PDBminer requires an input file in a comma-separated values file format, with four mandatory columns; "hugo_name", "uniprot", "uniprot_isoform", and "mutations". It is optional to include a fifth column titled "cluster_id". 
The cluster id should be an integer and can be used to annotate mutations forming a spacial cluster. See "inputfile.csv" as an example. 

```
hugo_name | uniprot | uniprot_isoform | mutations         | cluster_id
-----------------------------------------------------------------------
TP53      |  P04637 |         2       | P278L;R337C;L344P | 1
MAT1A     |  Q00266 |         1       | P30N;W300H        | 1
SSTR3     |  P05543 |         1       | T11S;C191S;R330L  | 1
SAMD4A    |  Q9UPU9 |         3       | L10R;I80A         | 1

        
```
The name of the input file should be specified in the commandline. 

## Running the Program

```
$ python PDBminer -i [input file name] -n [cores]
```
See example directory.

## The Output
The output is found in the directory "results".
For each uniprot id in the input file, a directory is created. After a successful run, this directory will contain the following: 

* {unipot_id}_input.csv, the input
* {unipot_id}_done.txt, indicating that the job has finished for this protein.

Moreover, there will be varying other files: 

* {unipot_id}_all.csv, An output file with all PDBs and alphafold structure associated with the uniprot_id regardless of mutational coverage.
* {unipot_id}_filtered.csv, An output file with the PDBs and alphafold structure associated with the uniprot_id that covers at least one mutation.
* missing_id.txt, If there are no structures from the protein data bank or alphafold structures available for the uniprot_id.

See an example of the output of the example input file in the directory results.

#content of {unipot_id}_clean.csv and {unipot_id}_all.csv:

```
Output Columns and explanations
================================================================
#hugo_name                      Gene names from the input file
#uniprot_id                     ID 
#uniprot_isoform                Number identifying which isoform to which the alignment was done. 
#cluster_id                     Number of cluster from the input file, or an assigned cluster "999" if not used. 
#structure_id                   Identifier of the PDB file or Alphafold Model.
#Research_mutations             The mutations specified in the input file
#experimental_method            By which the PDB was generated.
#resolution                     Estimation of PDB quality
#deposition_date                Timing of file placement in PDB or model generation in the Alphafold Database
#chains                         Letter describing the chains covering the Uniprot ID 
#coverage                       For structures from the Protein Data bank, this indicates the range where the PDB file and uniprot ID are aligned. The alignment is done based on the reported sequence in the PDB and does not take low quality areas or breaks in the model into account. A quality control should ALWAYS be conducted on PDB structures. For the alphafold models the coverage is based on pLDDT scores > 70, which means that the coverage of a single chain can be split for different domains, such as [(9,22),(30,60)] indicating that the positions within these intervals are of high quality. “;” separated for multiple chains.
#mutations_in_pdb               All mutations found in the pdb file compared to the uniprot sequence if none: [] indicating WT structure.
#complex_protein                Binary, either “NA” or “Protein in complex” indicates if the PDB file contains a protein complex.
#complex_protein_details        Details regarding the protein complex indicating the Uniprot ID of the other protein and the chains.
#complex_nucleotide             Binary, indicates if the protein is bound to a nucleotide string such as DNA.
#complex_nucleotide_details     Details regarding the DNA or RNA binding. 
#complex_ligand                 Binary indicating if metal or small molecules are present in the pdb file.
#complex_ligand_details         Details describing which “other” things are in the file.

```
