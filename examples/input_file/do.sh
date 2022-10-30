module load conda/4.9.2/modulefile
conda activate PDBminer

#require input file to run

../../PDBminer -i input_file.csv -n 1  -f ../../program/snakefile
../../PDBminer2coverage
