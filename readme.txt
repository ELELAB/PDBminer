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
# Running it: 

#### option 1/2:
#	Create a input file. This is recommended if you are
#	running multiple proteins. 
#	notice that only hugo_name and uniprot columns 
#	are required. 

tsp -N 4 python PDBminer -i {input_file} -n {cores}

tsp -N 4 python PDBminer -i input_file.csv -n 4

#### option 2/2:
#	write the protein of interest in the commandline.
#	this can be done when only one protein is of
#	interest. 

tsp -N {cores} python PDBminer -g {hugo_name} -u {uniprot_id} -s {uniprot_isoform} -m {mutations} -c {cluster_id} -n {cores}

#only -g, -u and -n is required. example: 

python PDBminer -g RNH -u P0A7Y7 -n 1

python PDBminer -g RNH -u P0A7Y7 -m "D134N;D134K" -n 1

#notice that only 1 core is used when only one protein is used, 
#this is because parallelization via snakemake is done per
#protein. A generated input_file.csv will be available after 
#the run as well. 
