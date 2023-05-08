# PDBminer

Repository associated to the Preprint:
```
PDBminer to Find and Annotate Protein Structures for Computational Analysis
Kristine Degn, Ludovica Beltrame, Matteo Tiberti, Elena Papaleo
bioRxiv 2023.05.06.539447; doi: https://doi.org/10.1101/2023.05.06.539447
```

## Introduction to the Program 
PDBminer is a snakemake pipeline that generates a ranked overview of the available structural models in the 
Protein Data Bank and the most current version of the AlphaFold2 model, if any.

* Installation
* Setup
* Running PDBminer
* Understanding the Output
* Plotting the Output

<a name="-dependencies"></a>
## Dependencies

It is recommended to create a virtual environment to run PDBminer. The environment can be installed 
via conda using the environment.yml or via pip with requirements.txt.

### Conda 
#### First time:

```
git clone https://github.com/ELELAB/PDBminer.git
conda env create -f environment.yml
conda activate PDBminer
```

#### All subsequent times

```
conda activate PDBminer
```

### Pip
#### First time:
```
git clone https://github.com/ELELAB/PDBminer.git
python3 -m venv PDBminer
source PDBminer/bin/activate
python3 -m pip --default-timeout=1000 install -r requirements.txt
```

#### All subsequent times

```
source PDBminer/bin/activate
```

<a name="-setup"></a>
## Setup

When running PDBminer, you must specify the location of the program using the -f flag, however, if you wish to avoid this 
you can change the default placement in the PDBminer script. 
The default is default = "/usr/local/envs/PDBminer/PDBminer/program/snakefile" when you clone, which may not 
fit your setup. If you choose to make this change, disregard -f in the following and the examples. 

## Running PDBminer the first time
There are two ways of running PDBminer. Either by using and input file, or by using the command line
to find the available structures for a single protein. 
In the directory examples are three examples and their commands in the do.sh file. Consider testing the installation and use
by running one or more of these. Notice that the -f flag is used in the do.sh files. 

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
The name of the input file should be specified in the command line: 

### Running the Program with an input file

```
$ python PDBminer -i [input file name] -n [cores]
```
### Using the command line directly

When PDBminer is run on only a single protein it may sometimes be beneficial to run it directly in the 
commandline. To do so, a input file does not need to be constructured and the content can be specified 
with flags. Again is the hugo_name and uniprot options mandatory while the rest is optional. 

```
$ python PDBminer -g [hugo_name] -u [uniprot_id] -s [uniprot_isoform] -m [mutations] -c [cluster_id] -n [cores] -f [path to snakefile]

$ python PDBminer -g TP53 -u P04637 -m "P278L;R337C;L344P" -n 1 -f program/snakefile
``` 

NOTICE: when isoform is not specified 1 is assumed.

<a name="-the-output"></a>
## The Output
A log.txt file is created for each run. 

The output is found in the directory "results".
For each uniprot id in the input file, a directory is created. 
After a successful run, this directory can contain the following: 

* missing_id.txt, If there are no structures from the protein data bank or alphafold structures available for the uniprot_id.

However, after the advent of alphafold, this is not a very likely output. It is more likely to have 
a file called {uniprot_id}_all.csv:

* {unipot_id}_all.csv, An output file with all PDBs and AlphaFold structure associated with the uniprot_id regardless of mutational coverage.

If mutations are included in the input, a filtered version of all will also be available if the mutations
are covered by any structure.
* {unipot_id}_filtered.csv, An output file with the PDBs and alphafold structure associated with the uniprot_id that covers at least one mutation.

Notice that multiple filtered files are available when multiple clusters are parsed. 
* {uniprot_id}_cluster{cluster_id}_filtered.csv

See examples of the in- and  output of the example directories.

#content of {unipot_id}_clean.csv and {unipot_id}_all.csv:

## Output Columns and Explanations:

* structure_rank: Index, the lower the value the seemingly better the model.
* hugo_name: Gene name from the input.
* uniprot_id: Uniprot id from the input.
* uniprot_isoform: Uniprot isoform from the input or 1 if none given.
* mutations: A list of input mutations, only visible in the filtered output.
* cluster_id: If clusters are specified, the cluster will be visible in this column.  
* structure_id: Identifier of the PDB file or Alphafold Model.
* deposition_date: Timing of file placement in PDB or model generation in the Alphafold Database.
* experimental_method: By which the PDB was generated. For AlphaFold model "predicted" is used. 
* resolution: Estimation of PDB quality (for X-ray structures).
* complex_protein: A column defining if a complex or fusion is present in the PDB file.
* complex_protein_details: Details regarding the protein complex indicating the Uniprot ID of the other protein and the chains.
* complex_nucleotide: Binary, indicates if the protein is bound to a nucleotide string such as DNA.
* complex_nucleotide_details: Details regarding the DNA or RNA binding.
* complex_ligand: Binary indicating if metal or small molecules are present in the pdb file.
* complex_ligand_details: Details describing which “other” things are in the file.
* chains: Letter or letters describing the chains covering the Uniprot ID 
* coverage: Coverage is a range of numbers indicating the area the model covers the uniprot sequence using the uniprot numbering of the sequence. One chain can have multiple sequences, for PDB structures indicating missing residues and for AlphaFold structures when the pLDDT score is below 70.  
* mutations_in_pdb: When structures contain amino acid sequence that differ from the sequence specified with the isoform, they are annotated as mutations and can be found here. 
* warnings: This coloumn contains information regarding dissimilarities or discrepancies that cannot be explained by fusions or other known alterations. This includes expression tags and amino acids added to terminals for structure solving purposes.

For all columns ";" seperate data on the annotated chains and "NA" indicates that no relevant data is present.

<a name="plotting"></a>
# Plotting 

## PDBminer2coverage
PDBminer2coverage is a plotting tool creating an overview of the protein sequence on the x-axis and the different models 
covering the sequence on the y-axis. The area the model covers is colored grey. If any positions are mutated, the position 
will be colored in with a transparent blue hue across all entries on the y-axis. 
PDBminer2coverage takes the results directory and input file as required input, per default these are input_file.csv and 
the current working directory/results. Hence, the PDBminer2coverage module can be run in the same place as PDBminer without any arguments.
The output is one or more plots. If there are both a filtered- and an all-output file, both will be plotted. 
If there are multiple clusters, these will be plotted separately. In cases where the protein sequence is longer than 500 amino acids, 
the plot will be split into multiple output files, termed "chunks". Additionally, it is possible to narrow down the 
plotting area with the flag -s --sequence. 

```
$ PDBminer2coverage -s 1-20,50-95 
```

Would, for example, only plot the sequence 1-20 and 50-95 in the same plot.  

Additionally you can also set a limit on the x-axis, indicating how many positions you want plotted using -t. 

```
$ PDBminer2coverage -t 50 
```
would, for example only plot 50 amino acids per chunk. Default is 500. 

Flags:
```
-r: choosing the results path if not default.
-i: The input file.
-u: uniprot id can be added if only one of the proteins in a multi protein run should be visualized.
-s: sequence, dash and comma separated values such as 1-10 or 1-20,30-40 for the sequence of the protein to plot. 
-t: threshold, a integer of the sequence length to be placed in each individual plot, default is 500, why each plot is maximum 500 amino acid broad. 
-c: color coverage, the plot is white and the coverage of each structure is colored in. Default is grey, but any color can be used, e.g '#64b2b5'. 
-m: mutation color, is the overlay color on the sites of mutation. It is blue per default but any color can be used, e.g '#183233'
-d: color for mutations in the PDBfile, meaning the inherent mutations within the file. These can be colored different from the WT amino acids. fault is light grey.
```
All options, example:

```
$ PDBminer2coverage -r PDBminer_run/results -i PDBminer_run/input_file.csv -u P00000 -s 30-120 -t 100 -c '#64b2b5' -m '#183233' -d '#a3cacc'
```

## PDBminer2network
PDBminer2network is a plotting tool creating an overview of the protein complexes within the protein data bank for the protein of interest.
PDBminer2network takes the current working directory/results directory as default input, just as the PDBminer2coverage module.
The output is a network with the Uniprot ID of your protein at the center, branching out to "protein complexes" and/or "fusion products" depending
on the content of the output file. From here each bound or fused protein is the first node and the second is the PDBid with the complex or fusion. 
The output is one or more plots. If there are both a filtered- and an all-output file, both will be plotted. If there are multiple clusters, 
they will be plotted seperately.

```
$ PDBminer2network -h
```

Would, for example, write out the help information.

The plot can be adjusted using the following flags: 

```
-r: choosing the results path if not default. 
-i: The input file. 
-u: uniprot id can be added if only one of the proteins in a multi protein run should be visualized.
-c: node color for center of graph
-p: node color for proteins (fused and bound)
-s: node color for structures (PDBid)
-t: color for the nodes with "protein complex" and "fusion product".
```
The edges of the graph are black. Any color can be used e.g. '#183233’.

All options, example:

```
$ PDBminer2network -r PDBminer_run/results -i PDBminer_run/input_file.csv -u P00000 -c '#64b2b5' -p '#183233' -s '#1fc6cc' -t '#a3cacc'
```
