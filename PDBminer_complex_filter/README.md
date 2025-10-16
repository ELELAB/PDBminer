# Filtering PDBminer output to identify Complexes in the Protein Data Bank

## Description

## Requirements
- Python >= 3.10
- pandas
- os
- MDAnalysis
- ast
- argparse
- PDBminer output file generated generated in at least 2024

## Running the script

- Step 1: Activate the envrionment `module load python`
- Step 2: Run the script: `python find_PDBminer_complexes.py`

For step 2, remember to add your flags.

### Flags 

- `-i`: The input file, should be the PDBminer output file either {uniprot}_all.csv or {uniprot}_all.json
- `-o`: Name of the output file, per default {uniprot}_filtered.csv
- `-s`: Start residue if you are interested in a specific domain, int. 
- `-e`: End residue if you are interested in a specific domain, int. 
- `--binding_interface`: If you want the interface residues for each self_chain (the target protein chains) binding to any other chain in the complex to be reported. 
- `-d`: If you use `--binding_interface` you can adjust the distance of contact, it is 5 per default, set accordingly, int

### Example
`python find_PDBminer_complexes.py -i uniprot_all.csv -o uniprot_filtered.csv --binding_interface -d 5 -s start_residue -e end_residue`

### Output
The output file contains the following columns and an example
- `structure_id` - 7YGI
- `self_chains` - A;B
- `coverage` - {'A': '[(92, 289)]', 'B': '[(92, 289)]'}
- `mutations` - NA
- `complex_details`, "P53_HUMAN, P04637, chain_A;AZUR_PSEAE, P00282, chain_C"
- `binding_partners`, A residues binding C
- `residues` - [551, 555, 569, 584, 606, 617, 625, 63]
- `complex_type` - protein complex

Where the 'complex_details' here indicate that p53 (chain A of PDB) is in complex with the Azurin protein (chain C of PDB) and binds on 
residues 551, 555, 569 .. on chain A. 
