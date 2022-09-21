#This read me describes how to run PDBminer, 
#on the Bioinfo server after cloning the repro
#from github.

#Everytime:
module load conda/4.9.2/modulefile

#First time:
conda env create -f environment_python.yml
conda activate PDBminer
conda install -c conda-forge biopython=1.78
conda install -c bioconda -c conda-forge snakemake=7.7.0

#all subsequent times
conda activate PDBminer

# prep files as described on Github
# Running it
tsp -N 4 snakemake --cores 4
