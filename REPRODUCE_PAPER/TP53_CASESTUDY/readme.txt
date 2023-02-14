#activate environment
../../PDBminer -g TP53 -u P04637 -f ../../program/snakefile
../../PDBminer2coverage -i input_file.csv -u P04637 -c '#238A8DFF' -d '#482677FF'
../../PDBminer2network -i input_file.csv -u P04637 -c '#1F968BFF' -p '#238A8DFF' -s '#B8DE29FF' -t '#39568CFF'


