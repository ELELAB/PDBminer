# PDBminer
29-July-2021

## Introduction to the Program 
PDBminer is a collection of scripts that takes a singular input file containing information about a protein and its mutations and outputs the best possible structural model in the PDB file format covering the protein and its mutations. 

To identify the best model, PDBminer works on the premise that the hierarchy is as follows: 
* The best option is a solved structure with good resolution.  
* Secondly, a solved structure that has undergone mutagenesis
* Thirdly, a structure created using homology modeling
* Lastly, a de novo model. 

Hence, the program runs iteratively aiming to find a suitable model through the hierarchy. 

PDBminer is currently only applicable for solved structures.  

## Dependencies

This project is writen in Python 3.8

The following packs are required: 

* Pandas 1.3.1
* Biopython 1.79
* glob 
* numpy

## Setup
Requires a Input file in a csv file format, containing the following:

```
Hugo_name | Uniprot | Uniprot-isoform | Mutations
-------------------------------------------------
          |         |                 |           
```

code is written as:

run_program.py is the main script that will guide the program and currently contains crude user prompting. Run from the terminal using: 

```
$ .....

```
The output:
* A directory called "structure_lists" containing a csv file for each uniprot ID. (to be expanded upon when the directory does contain more)
* A directory called "sturctures" containing the final output structures in a PDB file format for further analysis.  
