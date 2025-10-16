# AFTER HAVING SUCCESSFULLY RUN PDBMINER YOU CAN FILTER ON COMPLEXES:

#activate the python environment (module load python)

#SIMPLE USE, just specify input and output
python ../find_PDBminer_complexes.py -i P04637_all.csv -o P04637_filtered.csv

#DOMAIN USE, define your domain of interest
python ../find_PDBminer_complexes.py -i P04637_all.csv -o P04637_filtered_DBD.csv -s 91 -e 289

#INTERFACE USE:
python ../find_PDBminer_complexes.py -i P04637_all.csv -o P04637_filtered_interface.csv --binding_interface

#INTERFACE USE WITH DEFINED DISTANCE:
python ../find_PDBminer_complexes.py -i P04637_all.csv -o P04637_filtered_DBD_interface.csv --binding_interface -d 5 -s 91 -e 289

#NB! all the pdb files are downloaded when running the interface, they are removed here after the run.
