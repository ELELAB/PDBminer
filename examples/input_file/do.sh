module load conda/4.9.2/modulefile
conda activate PDBminer
cp -r ../../program/ .
../../PDBminer -i input_file.csv -n 1
rm -r program

