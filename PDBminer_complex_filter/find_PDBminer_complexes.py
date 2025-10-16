#!/usr/bin/env python3

# Copyright (C) 2024, Kristine Degn
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import os
import MDAnalysis as mda
import ast
import argparse
import requests

parser = argparse.ArgumentParser(prog = "PDBminer complexes")

parser.add_argument("-i", '--inputfile',
                    metavar = 'input file',
                    type=str,
                    help='the output from PDBminer in csv or json format')

parser.add_argument("-o", '--outputfile',
                    metavar = 'output file',
                    type=str,
                    help='name of the output')

parser.add_argument("-d", "--distance",
                    metavar = "distance in Å",
                    type=float,
                    default = 5.0,
                    help="The local interaction distance")

parser.add_argument("-s", "--start",
                    metavar = "start_residue",
                    type=int,
                    default = 1,
                    help="The minimum residue of domain")

parser.add_argument("-e", "--end",
                    metavar = "end_residue",
                    type=int,
                    default = 50000,
                    help="The maximum residue of domain")

parser.add_argument("-b", "--binding_interface",
                    action='store_true',
                    help="if the interface should be identified")

args = parser.parse_args()

def import_data(file_name):
    if file_name.endswith('.json'):
        df = pd.read_json(file_name)
    elif file_name.endswith('.csv'):
        df = pd.read_csv(file_name)
    else:
        raise ValueError("Unsupported file format. Supported formats: .json, .csv")
    return df

def download_pdb_file(structure_id):
    pdb_file_path = f"{structure_id}.pdb"
    if os.path.isfile(pdb_file_path):
        print(f"The file {pdb_file_path} already exists.")
    else:
        response = requests.get(f"https://files.rcsb.org/view/{structure_id}.pdb")
        if response.status_code == 200:
            download_command = f"wget https://files.rcsb.org/view/{structure_id}.pdb"
            os.system(download_command)
            print(f"File {pdb_file_path} downloaded successfully.")
        else:
            print(f"{structure_id} PDB cannot be downloaded.")

def get_residues_protein_complex(self_chain, binding_chain, distance, universe):

    interacting_residues_self = set()

    for atom in universe.select_atoms(f'segid {self_chain} and around {distance} segid {binding_chain}'):
        interacting_residues_self.add(atom.resid)

    interacting_residues_binder = set()

    for atom in universe.select_atoms(f'segid {binding_chain} and around {distance} segid {self_chain}'):
        interacting_residues_binder.add(atom.resid)

    return interacting_residues_self,  interacting_residues_binder

def create_universe_collect_data(pdb, self_chains, complex_partners, distance):
    download_pdb_file(pdb)

    binder_dictionary = {}

    if os.path.isfile(f"{pdb}.pdb"):
        try:
            u = mda.Universe(pdb+".pdb")
            for chain in self_chains:
                for binding_chain in complex_partners:
                    #retuns a set of numbers
                    interacting_residues_self, interacting_residues_binder = get_residues_protein_complex(chain, binding_chain, distance, u)
                    interacting_residues_self = list(interacting_residues_self)
                    interacting_residues_self.sort()
                    binder_dictionary[f'{chain} residues binding {binding_chain}'] = interacting_residues_self
                    interacting_residues_binder = list(interacting_residues_binder)
                    interacting_residues_binder.sort()
                    binder_dictionary[f'{binding_chain} residues binding {chain}'] = interacting_residues_binder
        except Exception as e:
            print(f"Error loading PDB into MDAnalysis: {e}")
    return binder_dictionary

if __name__ == '__main__':
    file_name = args.inputfile
    df = import_data(file_name)

    output_df = pd.DataFrame(columns=["structure_id", "self_chains", "coverage",
                                      "mutations", "complex_details", "binding_partners",
                                      "residues", "complex_type"])
    for structure in df.structure_id:
        info = df[df.structure_id == structure].iloc[0]
        if not pd.isna(info.complex_protein):
            self_chains = info.chains.split(";")
            complex_protein_details = ast.literal_eval(info.complex_protein_details) if type(info.complex_protein_details) == str else info.complex_protein_details
            complex_partners = complex_protein_details[0].split(';')
            coverage_dict = ast.literal_eval(info.coverage) if type(info.coverage) == str else info.coverage

            for key, value_str in coverage_dict.items():
                if value_str != 'NA' and ';' not in value_str:
                    try:
                        value_list = ast.literal_eval(value_str)
                        for start, end in value_list:
                            if args.start <= start <= args.end  or args.start <= end <= args.end:
                                binding_chain_dict = {}
                                for chain in self_chains:
                                    binding_chains = []
                                    for entry in complex_partners:
                                        try:
                                            relevant_chain = entry.split(', ')[2].split('_')[1]
                                            if relevant_chain != chain:  # prevent self-binding
                                                binding_chains.append(relevant_chain)
                                        except IndexError:
                                            continue
                                    binding_chain_dict[chain] = binding_chains

                                if args.binding_interface:
                                    binder_dictionary = {}
                                    for self_chain in self_chains:
                                        if self_chain in binding_chain_dict:
                                            bd = create_universe_collect_data(structure, [self_chain], binding_chain_dict[self_chain], args.distance)
                                            binder_dictionary.update(bd)
                                    if binder_dictionary != {}:
                                        for key, value in binder_dictionary.items():
                                            new_row = {
                                                "structure_id": structure,
                                                "self_chains": info.chains,
                                                "coverage": info.coverage,
                                                "mutations": info.mutations_in_pdb,
                                                "complex_details": complex_protein_details[0],
                                                "binding_partners": key,
                                                "residues": value,
                                                "complex_type": info.complex_protein
                                            }
                                            output_df = pd.concat([output_df, pd.DataFrame([new_row])], ignore_index=True)

                                else:
                                    new_row = {
                                        "structure_id": structure,
                                        "self_chains": info.chains,
                                        "coverage": info.coverage,
                                        "mutations": info.mutations_in_pdb,
                                        "complex_details": complex_protein_details[0],
                                        "binding_partners": "NA",
                                        "residues": "NA",
                                        "complex_type": info.complex_protein
                                    }
                                    output_df = pd.concat([output_df, pd.DataFrame([new_row])], ignore_index=True)

                    except (ValueError, SyntaxError):
                        print(f"Error parsing value for Coverage {key}")

    subset_columns = ['structure_id', 'self_chains', 'complex_details', 'binding_partners']
    output_df.drop_duplicates(subset=subset_columns, keep="first", inplace=True)
    if args.binding_interface:
        output_df = output_df[output_df.residues.apply(lambda x: isinstance(x, list) and len(x) > 0)]
    output_df = output_df.reset_index(drop=True)
    output_df.to_csv(args.outputfile)

