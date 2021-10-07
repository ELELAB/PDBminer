# PDBminer

THIS READ ME IS OBSOLETE

## Introduction to the Program 
PDBminer is a collection of scripts that takes a singular input file containing information about a protein and its mutations and outputs tree files, a list containing the proteins (uniprot IDs) that does not have any solved models associated with it, a .csv file of all the best solved structures for that particular protein and a .csv file of all the structures associated with the uniprot ID that was removed and the reason for their removal. 

## Dependencies

This project is writen in Python 3.8

The following packs are required: 

* Pandas 1.3.1
* Biopython 1.79
* numpy

It is recommended to create a virtuel environment, virtualenv -p /usr/bin/python3.7 [name of env]

## Input file
Requires a Input file in a csv file format, containing the following:

Hugo_name	String, e.g., AKT1
Uniprot		Uniprot ID e.g., P31749 
Uniprot-isoform	Integer
Mutations	strings separated by a semicolon, ;, e.g. E17K or L755S;K128M, deletion with *
ClusterID	Integer or N/A.

```
Hugo_name | Uniprot | Uniprot-isoform | Mutations         | ClusterID
-----------------------------------------------------------------------
TP53      |  P04637 |         2       | P278L;R337C;L344P | 1
MAT1A     |  Q00266 |         1       | P30N;W300H        | 1
SSTR3     |  P05543 |         1       | T11S;C191S;R330L  | 2
SAMD4A    |  Q9UPU9 |         3       | L10R;I80A         | 5
        
```
## Running the program

```
$ python3 [script placement]/run_PDBminer.py [input file] 

# e.g. 

$ python3 /Desktop/PDBminer_program/run_PDBminer.py /Desktop/input_files/fast_input_tcga_3d.csv


```
## The Output

The program creates two directories

* structures: Where all the downloaded PDB structures reside. Can be deleted or kept after finalization. 
* output_files: Here one to three files are placed:

        1) missing_IDs.txt: a textfile containing a list of Uniprot IDs where there are no PDB structures.
        2) potential_structures.csv: A csv file that contains the PDB structures that have passed each quality control step
        3) inspection_structures.csv: A csv file that contains the PDB structures that have not passed each quality control step and why. 

Example of the potential_structures.csv:
```
PDB_ID | Uniprot_ID | ClusterID | experimental_method  | resolution | deposition_date | mutation_coverage     |  Research_mutations  | model_mutated
------------------------------------------------------------------------------------------------------------------------------------------------------
6IPU   | P68431     | 1         | X-RAY DIFFRACTION    | 1.99       | 2018-11-04      | {'A': ['D', 'K', 'E', |  D82N;K80N;E74Q;F79L;| Not-Mutated
       |            |           |                      |            |                 | 'F', 'E', 'R'], 'E':  |  E74K;R73P           |
       |            |           |                      |            |                 | ['D', 'K', 'E', 'F',  |                      |
       |            |           |                      |            |                 | 'E', 'R']}            |                      |
------------------------------------------------------------------------------------------------------------------------------------------------------
6KE9   | P68431     | 1         | X-RAY DIFFRACTION    | 2.22       | 2019-07-04      | {'E': ['D', 'K', 'E', |  D82N;K80N;E74Q;F79L;| Not-Mutated
       |            |           |                      |            |                 | 'F', 'E', 'R']}       |  E74K;R73P           |
------------------------------------------------------------------------------------------------------------------------------------------------------
6JXD   | P68431     | 1         | X-RAY DIFFRACTION    | 2.25       | 2019-04-23      | {'A': ['D', 'K', 'E', |  D82N;K80N;E74Q;F79L;| Not-Mutated
       |            |           |                      |            |                 | 'F', 'E', 'R'], 'E':  |  E74K;R73P           |
       |            |           |                      |            |                 | ['D', 'K', 'E', 'F',  |                      |
       |            |           |                      |            |                 | 'E', 'R']}            |                      |

```
