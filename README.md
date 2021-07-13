# PDBminer
scripts to find PBD structures for cancer driver proteins

This collection of files are ready for an inital review. 

run_program.py is the main script that will guide the program and currently contains crude user prompting. 



The input is a .csv file that should contain the following information:

Hugo_name | Uniprot | Uniprot-isoform | Mutations



The output is currently (13-July-2021) a directory called "structure_lists" which contains two primary things:
  1:  A .csv file for each uniprot ID with all the avialable pdb structures including specifications about the experimental method, 
      resolution and desposition date
  2: A .txt file containing a list of all the uniprot IDs to which no pdb file exist. 
