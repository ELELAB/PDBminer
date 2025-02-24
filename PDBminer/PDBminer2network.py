#!/usr/bin/env python3

# Copyright (C) 2023 & 2025, Kristine Degn
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
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
import argparse
import os
import logging

parser = argparse.ArgumentParser(prog = "PDBminer2network",
                                 usage = "PDBminer2network [-h] help, all flags are optional: [-r] results path, [-i] input file, [-u] uniprot id, can be added if the user wish to visualize only one protein network, [-c] color center, the center of the plot is colored dark blue by default, but any color can be used, e.g '#64b2b5'. [-p] protein color is grey by default but any color can be used, e.g '#183233’. [-s] PDBid color is blue by default but any color can be used, e.g '#1fc6cc'")

parser.add_argument("-r", '--result',
                    metavar = 'result path',
                    type=str,
                    default = "./results/",
                    help='The path to the result path, default ./results')

parser.add_argument("-i", "--inputfile",
                    metavar = 'input file',
                    type=str,
                    default = "input_file.csv",
                    help="input file created in the PDBminer run or original file")

parser.add_argument("-u", "--uniprot",
                    metavar = 'uniprot',
                    default = "none",
                    type=str,
                    help="Uniprot id when multiple proteins are run")

parser.add_argument("-c", "--color_center",
                    metavar = "color for center of graph",
                    default = "darkblue",
                    type=str,
                    help="plot color for center")

parser.add_argument("-p", "--color_proteins",
                    metavar = "color for proteins",
                    default = "grey",
                    type=str,
                    help="plot color for proteins")

parser.add_argument("-s", "--color_structures",
                    metavar = "color for structures",
                    default = "lightblue",
                    type=str,
                    help="plot color for structures")

parser.add_argument("-t", "--type_color",
                    metavar = "color for the protein and fusion clades",
                    default = "blue",
                    type=str,
                    help="plot color for structures")

parser.add_argument('-v', '--verbose', 
                    action='store_true', 
                    help='Enable verbose output')

def plot_protein_network(path, uniprot, file_name, output_name,
                         center_color, protein_color, structure_color, 
                         type_color):
    from_list = []
    to_list = []
    
    p = file_name.dropna(subset=['complex_protein'])
    p = p[p.complex_protein != "NA"]
    if len(p) == 0:
        return logging.debug(f"plot for {uniprot}/{output_name} omitted, no protein complex found")
        
    #regular protein complex
    protein_complex=p[p.complex_protein.str.startswith('protein complex')]

    #fusion products
    fusion=p[p.complex_protein.str.startswith('fusion product')]

    if len(protein_complex) > 0:
        from_list.append(uniprot)
        to_list.append("Protein Complex")
        protein_complex = protein_complex[['structure_id', 'complex_protein_details']]
        protein_complex = protein_complex.reset_index(drop = True)
        for i, info in enumerate(protein_complex.values):
            structure = info[0]
            i = info[1]
            if type(i) == str:
                i = i.strip('[]')
            else:
                i = i[0]
            i = i.split(";")
            for item in i:
                if uniprot in item:
                    i.remove(item)
            if len(i) > 0:
                for protein in i:                
                    from_list.append("Protein Complex")
                    protein_specifics = protein.split(", ")[1]
                    to_list.append(protein_specifics)
                    from_list.append(protein_specifics)
                    to_list.append(structure)
    
        if len(fusion) > 0:
            from_list.append(uniprot)
            to_list.append("Fusion Product")
            fusion = fusion[['structure_id', 'complex_protein_details']]
            fusion = fusion.reset_index(drop = True)
            for i, info in enumerate(fusion.values):
                structure = info[0]
                i = info[1]
                if type(i) == str:
                    i = i.strip('[]')
                else:
                    i = i[0]
                i = i.split(";")
                for item in i:
                    if uniprot in item:
                        i.remove(item)
                if len(i) > 0:
                    for protein in i:                
                        from_list.append("Fusion Product")
                        protein_specifics = protein.split(", ")[1]
                        to_list.append(protein_specifics)
                        from_list.append(protein_specifics)
                        to_list.append(structure)

        f = pd.DataFrame({'from':from_list, 'to':to_list})

        # Build your graph
        G=nx.from_pandas_edgelist(f, 'from', 'to')

        rcParams['figure.figsize'] = 20, 20
        pos = nx.spring_layout(G, scale=20, k=3/np.sqrt(G.order()))
    
        colors = []
        size = []
        for node in G.nodes():
            if len(node) == 4:
                colors.append(structure_color)
                size.append(700)
            elif node == "Protein Complex":
                colors.append(type_color)
                size.append(6000)
            elif node == "Fusion Product":
                colors.append(type_color)
                size.append(6000)
            elif node == uniprot:
                colors.append(center_color)
                size.append(10000)
            else:
                colors.append(protein_color)
                size.append(1500)
            
        nx.draw(G, pos, 
                with_labels=True, 
                node_size = size,
                node_color=colors)
    
        plt.savefig(f"{path}/{uniprot}/{output_name}_network.pdf", format="pdf", bbox_inches="tight")
        plt.clf()
        out = f"plot for {uniprot}/{output_name} done"
    else:
        out = f"plot for {uniprot}/{output_name} omitted, no protein complex found"
    return out


def run(df, path, center_color, protein_color, structure_color, type_color):
    """
    A run file to generate the plots based on the original input file.

    Parameters
    ----------
    df : Dataframe of the inputfile.
    path : path to output
    center_color : color for center of graph (node)
    protein_color : color for proteins (nodes)
    structure_color: color for structures - PDBids (nodes)
    type_color: color for the protein and fusion clades (nodes)


    Returns
    -------
    None.

    """
    
    for i in range(len(df)): 
        logging.debug(f"{df.iloc[i].uniprot}")
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.json"):
            logging.debug("filtered file identified")
            filtered = pd.read_json(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.json")
            logging.debug("filtered is loaded")
            o = plot_protein_network(path, df.iloc[i].uniprot, filtered, "filtered",
                                     center_color, protein_color, structure_color, type_color)
            logging.debug(f"{o}")
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.csv"):
            logging.debug("filtered file identified")
            filtered = pd.read_csv(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.csv", na_values=[], keep_default_na=False)
            o = plot_protein_network(path, df.iloc[i].uniprot, filtered, "filtered",
                                     center_color, protein_color, structure_color, type_color)
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.json"):
            logging.debug("filtered cluster file identified")
            filtered = pd.read_json(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.json")
            o = plot_protein_network(path, df.iloc[i].uniprot, filtered,"filtered",
                                     center_color, protein_color, structure_color, type_color)
            logging.debug(f"{o}")
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.csv"):
            logging.debug("filtered cluster file identified")
            filtered = pd.read_csv(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.csv", na_values=[], keep_default_na=False)
            o = plot_protein_network(path, df.iloc[i].uniprot, filtered,"filtered",
                                     center_color, protein_color, structure_color, type_color)
            logging.debug(f"{o}")
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.json"):
            logging.debug("All file identified")
            full = pd.read_json(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.json")
            o = plot_protein_network(path, df.iloc[i].uniprot, full,"all",
                                     center_color, protein_color, structure_color, type_color)
            logging.debug(f"{o}")
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.csv"):
            logging.debug("All file identified")
            full = pd.read_csv(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.csv", na_values=[], keep_default_na=False)
            o = plot_protein_network(path, df.iloc[i].uniprot, full,"all",
                                     center_color, protein_color, structure_color, type_color)
            logging.debug(f"{o}")        
    return

def main():
    args = parser.parse_args()

    path = args.result

    log_file = f'{path}/log_network.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.addHandler(file_handler)

    df = pd.read_csv(args.inputfile)

    uniprot_id = args.uniprot
    if uniprot_id != "none":
        df = df[df.uniprot == uniprot_id]
        df = df.reset_index()
        
    center_color = args.color_center
    protein_color = args.color_proteins
    structure_color = args.color_structures
    type_color = args.type_color
    
    run(df, path, center_color, protein_color, structure_color, type_color)
    file_handler.close()    
    if os.path.exists(log_file) and os.path.getsize(log_file) == 0:
        os.remove(log_file)
        
if __name__ == '__main__':
    main()
