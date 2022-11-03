# PDBminer

## Introduction to the Program 
PDBminer is a snakemake pipeline that generates a ranked overview of the available structural models in the 
Protein Data Bank and the most current version of the Alphafold2 model, if any.

## Dependencies

It is recommended to create a virtual environment to run PDBminer. The environment can be created using 
environment_python.yml as described below:

### First time:

```
conda env create -f program/environment_python.yml 
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0
conda install -c conda-forge biopandas=0.4.1
conda install -c conda-forge matplotlib=3.2.2
conda install -c anaconda seaborn=0.12.0
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
There are two ways of running PDBminer. Either by using and input file listing one or more proteins, or by using the command line
to find the available structures for a single protein. 

### Using an input file

The input file should be in a comma-separated values file format, with two mandatory columns; "hugo_name", and "uniprot".
Additionally the optional columns are "uniprot_isoform", "mutations" and "cluster_id". 

Example of input_file.csv:

```
hugo_name | uniprot | uniprot_isoform | mutations         | cluster_id
-----------------------------------------------------------------------
TP53      |  P04637 |         2       | P278L;R337C;L344P | 1
MAT1A     |  Q00266 |         1       | P30N;W300H        | 1
SSTR3     |  P05543 |         1       | T11S;C191S;R330L  | 1
SAMD4A    |  Q9UPU9 |         3       | L10R;I80A         | 1

        
```
The name of the input file should be specified in the commandline: 

### Running the Program with an input file

```
$ python PDBminer -i [input file name] -n [cores]
```
### Using the commandline directly

When PDBminer is run on only a single protein it may sometimes be beneficial to run it directly in the 
commandline. To do so, a input file does not need to be constructured and the content can be specified 
with flags. Agian is the hugo_name and uniprot options mandatory while the rest is optional. 

```
$ python PDBminer -g [hugo_name] -u [uniprot_id] -s [uniprot_isoform] -m [mutations] -c [cluster_id] -n [cores]

$ python PDBminer -g TP53 -u P04637 -m "P278L;R337C;L344P" -n 1
``` 

NOTICE: when isoform is not specified 1 is assumed.

## The Output
A log.txt file is created for each run. 

The output is found in the directory "results".
For each uniprot id in the input file, a directory is created. 
After a successful run, this directory can contain the following: 

* missing_id.txt, If there are no structures from the protein data bank or alphafold structures available for the uniprot_id.

However, after the advent of alphafold, this is not a very likely output. It is more likely to have 
a file called {uniprot_id}_all.csv:

* {unipot_id}_all.csv, An output file with all PDBs and alphafold structure associated with the uniprot_id regardless of mutational coverage.

If mutations are included in the input, a filtered version of all will also be available if the mutations
are covered by any structure.
* {unipot_id}_filtered.csv, An output file with the PDBs and alphafold structure associated with the uniprot_id that covers at least one mutation.

Notice that multiple filtered files are available when multiple clusters are parsed. 
* {uniprot_id}_cluster{cluster_id}_filtered.csv

See examples of the in- and  output of the example directories.

#content of {unipot_id}_clean.csv and {unipot_id}_all.csv:

```
Output Columns and explanations
================================================================
#structure_rank			Index, the lower the value the seemingly better the model.
#hugo_name			Gene name from the input.
#uniprot_id			Uniprot id from the input.	
#uniprot_isoform		Uniprot isoform from the input or 1 if none given.
#mutations			The mutations from the input if any. 
#cluster_id			The cluster from the input if any.
#structure_id			Identifier of the PDB file or Alphafold Model.
#deposition_date		Timing of file placement in PDB or model generation in the Alphafold Database
#experimental_method		By which the PDB was generated. 
#resolution			Estimation of PDB quality (for X-ray structures)
#chains				Letter describing the chains covering the Uniprot ID 
#coverage			For structures from the Protein Data bank, this indicates the range where the 
#				PDB file and uniprot ID are aligned. The alignment is done based on the reported 
#				sequence in the PDB and does not take low quality areas or breaks in the model into account. 
#				A quality control should ALWAYS be conducted on PDB structures. For the alphafold models 
#				the coverage is based on pLDDT scores > 70, which means that the coverage of a single chain 
#				can be split for different domains, such as [(9,22),(30,60)] indicating that the positions 
#				within these intervals are of high quality. “;” separated for multiple chains.			 
#mutations_in_pdb		All mutations found in the pdb file compared to the uniprot sequence if none: 
#				[] indicating WT structure.
#missing_residues		A string of all the missing residues on the different chains, e.g. 'chain c: M1;L2;W3;W4;E5;E6;V7;E8'
#complex_protein                Binary, either “NA” or “Protein in complex” indicates if the PDB file contains a protein complex.
#complex_protein_details        Details regarding the protein complex indicating the Uniprot ID of the other protein and the chains.
#complex_nucleotide             Binary, indicates if the protein is bound to a nucleotide string such as DNA.
#complex_nucleotide_details     Details regarding the DNA or RNA binding. 
#complex_ligand                 Binary indicating if metal or small molecules are present in the pdb file.
#complex_ligand_details         Details describing which “other” things are in the file.

```
# Plotting 

## PDBminer2coverage
PDBminer2coverage is a plotting tool creating an overview of the protein sequence on the x-axis and the different models covering the sequence on the y-axis. The area the model covers is colored grey. If any positions are mutated, the position will be colored in with a transparent blue hue across all entries on the y-axis. 
PDBminer2coverage takes the results directory and input file as required input, per default these are input_file.csv and the current working directory/results. Hence, the PDBminer2coverage module can be run in the same place as PDBminer without any arguments.
The output is one or more plots. If there are both a filtered- and an all-output file, both will be plotted. If there are multiple clusters, these will be plotted separately. In cases where the protein sequence is longer than 500 amino acids, the plot will be split into multiple output files, termed "chunks". Additionally, it is possible to narrow down the plotting area with the flag -s --sequence. 

PDBminer2coverage -s 1-20,50-95 

Would, for example, only plot the sequence 1-20 and 50-95 in the same plot.  

Additionally you can also set a limit on the x-axis, indicating how many positions you want plotted using -t. 

PDBminer2coverage -t 50 

would, for example only plot 50 amino acids per chunk. Default is 500. 
