module load conda/4.9.2/modulefile
conda activate PDBminer
cp -r ../../program/ .
../../PDBminer -g SCP2 -u P22307 -n 1
