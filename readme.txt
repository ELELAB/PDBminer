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
conda install -c conda-forge matplotlib=3.2.2
conda install -c anaconda seaborn=0.12.0
conda install -c anaconda networkx=2.8.4

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

#Additionally, a plotting modules can be run

python PDBminer2coverage

#Plotting of the sequence of the protein and the relative mapping of the 
#structural coverage.

#PDBminer2coverage takes up to tree arguments and requires the inputfile, 
#either supplied or generated and the results folder. 

#If you run PDBminer2coverage in the same directory you ran PDBminer in, you
#do not need to add any flags. Alternatively you can:

python PDBminer2coverage -r PDBminer_run/results/ -i PDBminer_run/input_file.csv

#If you wish only to plot a section of the protein use the flag -s, example:

python PDBminer2coverage -s 1-20,60-140 

#Output: One or more plots illustrating the coverage of the found structures.

python PDBminer2network

#plotting of the protein complexes in the output file. 

#PDBminer2network requires the inputfile and result folder.
#If you run PDBminer2network in the same directory you ran PDBminer in, you
#do not need to add any flags. Alternatively you can:

python PDBminer2network -r PDBminer_run/results/ -i PDBminer_run/input_file.csv


