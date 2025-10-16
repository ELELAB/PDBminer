# AFTER HAVING SUCCESSFULLY RUN PDBMINER YOU CAN FILTER ON COMPLEXES:

#activate the python environment (module load python)

#SIMPLE USE, just specify input and output
python ../find_PDBminer_complexes.py -i P78352_all.csv -o P78352_filtered.csv

#INTERFACE USE WITH DEFINED DISTANCE:
python ../find_PDBminer_complexes.py -i P78352_all.csv -o P78352_filtered_interface.csv --binding_interface
#this downloads the relevant PDBs, removed here.
