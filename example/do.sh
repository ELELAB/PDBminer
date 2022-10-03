module load conda/4.9.2/modulefile
conda activate PDBminer
cp -r program/ example/
../PDBminer -i input_file.csv -n 1
