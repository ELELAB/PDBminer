#This read me describes how to run PDBminer, 
#on the Bioinfo server after cloning the repro
#from github.

#Everytime:
module load conda/4.9.2/modulefile

#First time:
conda env create -f program/environment_python.yml
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0
conda install -c conda-forge biopandas=0.4.1

#all subsequent times
conda activate PDBminer

# prep files as described on Github
# Running it
tsp -N 4 python PDBminer -i input_file.csv -n 4

#once PDBminer is added to our path, we can 
#omit "python". 
