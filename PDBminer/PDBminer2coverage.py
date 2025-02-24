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


#coverage and mutation coverage
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import requests
from Bio import SeqIO
from io import StringIO
import os
import seaborn as sns
import math
from matplotlib.colors import ListedColormap
import argparse
import re
from sys import exit
from requests.exceptions import ConnectionError
import ast
import logging

parser = argparse.ArgumentParser(prog = "PDBminer2coverage",
                                 usage = "PDBminer2coverage [-h] help, all flags are optional: [-r] results path, [-i] input file, [-u] uniprot id, can be added if the user wish to visualize only one coverage, [-s] sequence, dash and comma separated values such as 1-10 or 1-20,30-40 for the sequence of the protein to plot, [-t] threshold, a integer of the sequence length to be placed in each individual plot, default is 500, why each plot is maximum 500 amino acid broad, then you can specify the coloring based on b_factors: [-bb] a threshold value for the best b-factors, < 15 per default, [-bg] a threshold value for good b-factors, < 25 per default, [-bp] a threshold value for poor b-factors, > 40 per default")

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

parser.add_argument("-s", "--sequence",
                    metavar = "sequence",
                    type=str,
                    help="The area to be plotted, e.g., 1-200 or 1-20,50-66, comma seperated")

parser.add_argument("-t", "--threshold",
                    metavar = "theshold",
                    default = 500,
                    type=int,
                    help="Max number of residues plotted on the y-axis.")

parser.add_argument("-b", "--b_factor_best",
                    metavar = "b_factor_best",
                    default = 15,
                    type=int,
                    help="Threshold for the best b-factors")

parser.add_argument("-g", "--b_factor_good",
                    metavar = "b_factor_good",
                    default = 25,
                    type=int,
                    help="Threshold for good b-factors")

parser.add_argument("-p", "--b_factor_poor",
                    metavar = "b_factor_poor",
                    default = 40,
                    type=int,
                    help="Threshold for poor b-factors")

parser.add_argument('-v', '--verbose', 
                    action='store_true', 
                    help='Enable verbose output')

def get_uniprot_sequence(uniprot_id, isoform):
    """
    A function that takes the uniprot accession number and the isoform 
    and retrieves the fasta file with the sequence. 

    Parameters
    ----------
    uniprot_id : Uniprot accession number, String.
    isoform : Integer of isoform.

    Returns
    -------
    uniprot_sequence : A one letter amino acid sequence of the protein. String.
    uniprot_numbering : A list of numerical values describing the residue position.

    """
    logging.debug(f"FUNCTION: get_uniprot_sequence({uniprot_id}, {isoform})")
    
    if isoform==1:
    
        try: 
            response = requests.post(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
        except ConnectionError as e:
            logging.error(f"EXITING: Uniprot database API rejected for {uniprot_id}.")
            exit(1) 
    
    else:
        try: 
            response = requests.post(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}-{isoform}.fasta")
        except ConnectionError as e:
            logging.error(f"EXITING: Uniprot database API rejected for {uniprot_id}.")
            exit(1) 
    
    if response.status_code == 200:
        sequence_data=''.join(response.text)
        Seq=StringIO(sequence_data)
        pSeq=list(SeqIO.parse(Seq,'fasta'))
        uniprot_sequence = str(pSeq[0].seq)
        uniprot_numbering = list(range(1,len(uniprot_sequence)+1,1)) 

    else:
        logging.error(f"EXITING: The canonical sequence could not be retrieved, ensure that isoform {isoform} exist, {uniprot_id}.")
        exit(1) 
    
    uniprot_sequence = list(uniprot_sequence)
    columns = []
    for i in range(len(uniprot_numbering)):
        columns.append(f"{uniprot_sequence[i]}{uniprot_numbering[i]}")
    return uniprot_sequence, columns


def get_a_range(a, chain, uniprot_numbering):
    """
    A function that generates an array the length of the protein
    with the b-factor values as entries where there are coverage.

    Parameters
    ----------
    a : One row from the result file from PDBminer
    chain : The relevant Chain for the coverage
    uniprot_numbering : The canonical numbering of the uniprot sequence.

    Returns
    -------
    output : An array with 0 if no structure coverage and the b-factor if
    applicable.

    """
    UN = np.array(uniprot_numbering)
    output = np.zeros(UN.shape)  
    
    if type(a.coverage) == str:
        coverage = ast.literal_eval(a.coverage)
        coverage = coverage[chain]
    else:
        coverage = a.coverage[chain]
    if coverage == "NA" or coverage == "[]":
        return output
    #problems when working with csv files compared to json files. Everything is
    #converted to a string.
    if type(a.b_factor) == str:    
        if "\n" in a.b_factor:
            items_to_remove = ["\n", " ", "[", "]"]
            a.b_factor = ''.join(char for char in a.b_factor if char not in items_to_remove)
        if 'array' in a.b_factor:
            a.b_factor = a.b_factor.replace('array', '')
        bfactor = ast.literal_eval(a.b_factor)
        bfactor = bfactor[chain]
    else:
        bfactor = a.b_factor[chain]    
    for r in (("[", ""), ("]", ""), ("'", ""), ("(", ""), (")", "")): 
        coverage = coverage.replace(*r)
    coverage = coverage.split(";")
    for segment in coverage:
        start, end = map(int, segment.split(', '))  # Parse the coverage range
        output[start - 1:end] = bfactor[0:end - start + 1]
    return output

def convert_bfactor(a_range, bfactortheshold_best, bfactortheshold_good, 
                    bfactortheshold_poor):
    """
    A function that converts the b-fcator values to category values.
    The categroies are determined based on the variables. 

    Parameters
    ----------
    a_range : The full length array of the protein and b-factor values where 
    the structure coveres the protein.
    bfactortheshold_best : A integer defining a value where lower b-factors are 
    considered the best quality.
    bfactortheshold_good : A integer defining a value where lower b-factors are 
    considered good quality.
    bfactortheshold_poor : A integer defining a value where higher b-factors are 
    considered poor quality.

    Returns
    -------
    scored_values : An array shaped like a_range with the categorized b-factors.

    """
    scored_values = []
    for score in a_range:
        if score == 0.0:
            scored_values.append(0)
        elif score <= bfactortheshold_best:
            scored_values.append(4)
        elif score <= bfactortheshold_good:
            scored_values.append(3)
        elif score <= bfactortheshold_poor:
            scored_values.append(2)
        elif score > bfactortheshold_poor:
            scored_values.append(1)     
    scored_values = np.array(scored_values)
    return scored_values

def get_pLLDT_score(score):
    if score > 90:
        return 4
    elif score > 70:
        return 3
    elif score > 50:
        return 2
    else:
        return 1
    
    
def get_mutations(mutations_in_pdb, a_range, chain): 
    """
    A function that annotates the mutations in the pdb file.

    Parameters
    ----------
    mutations_in_pdb : A dictionary
    a_range : The b-factor categorized coverage array
    chain : Chain in structure.

    Returns
    -------
    a_range : A b-factor categorized coverage array where the mutations
    are marked with a novel category.

    """
    if mutations_in_pdb == "NA":
        return a_range
    if type(mutations_in_pdb) == str:
        mutations_in_pdb = ast.literal_eval(mutations_in_pdb)
        if mutations_in_pdb[chain] == "NA":
            return a_range        
    for mutation in mutations_in_pdb[chain]:
        mut_pos = int(mutation[1:-1])
        a_range[mut_pos-1] = 5
    return a_range       



def plot_coverage(path, isoform, uniprot_id, mutations, cluster, file, name, 
                  sequence, threshold, bfactortheshold_best, 
                  bfactortheshold_good, bfactortheshold_poor):
    """
    The plotting function, taking the output files from PDBminer and coverting 
    them into a plot.

    Parameters
    ----------
    path : The path to the relevant Uniprot ID folder in the results.
    isoform : From the input file of a particular protein.
    uniprot_id : From the input file of a particular protein.
    mutations : From the input file of a particular protein.
    cluster : From the input file of a particular protein.
    file : The PDBminer output file, either filtered or all
    name : Filtered or all
    sequence : User input of area to be plotted.

    Returns
    -------
    One or more plots of the coverage. The number of plots depends on the
    length of the protein (max 500 AA per plot), number of clusters and 
    availablity of a filtered and all file.

    """
    uniprot_numbering, columns = get_uniprot_sequence(uniprot_id, isoform)
    vectorized_get_pLLDT_score = np.vectorize(get_pLLDT_score)
    logging.debug("    columns calculated imported from uniprot sequence")
    f = file[['structure_id', 'chains', 'coverage', 'mutations_in_pdb', 'b_factor']]
    index = []
    ranges = []
    logging.debug("    converting ranges on chains to b-factor dependent values")
    for i in range(len(f)):
        a = f.iloc[i]
        logging.debug(f"        {a.structure_id}")
        chains = a.chains.split(";")
        for chain in chains:
            #creating an array the length of the uniprot sequence,
            #where the pLDDT and B-factors are converted into scores
            #that will illustrate the coverage in colors.
            if a.structure_id.startswith('AF-'):
                if type(a.b_factor) == str:
                    bfactor = ast.literal_eval(a.b_factor)
                else:
                    bfactor = a.b_factor    
                a_range = np.array(bfactor[chain])
                a_range = vectorized_get_pLLDT_score(a_range)
            else:
                a_range = get_a_range(a, chain, uniprot_numbering)
                a_range = convert_bfactor(a_range, 
                                              bfactortheshold_best, 
                                              bfactortheshold_good, 
                                              bfactortheshold_poor)
            # we need to color in mutations
            # if any a_range is assigned 5.
            a_range = get_mutations(a.mutations_in_pdb, a_range, chain)
            if set(a_range) != {0.0}: 
                index.append(f"{a.structure_id}_{chain}")                
                ranges.append(a_range)
    logging.debug("    combine to dataframe")
    heatmap = np.stack(ranges)
    df = pd.DataFrame(data=heatmap, index=index, columns=columns)
    
    if np.array_equal(sequence, 0) == False:
        if len(sequence) > 2:
            evens = [x for x in range(len(sequence)) if x%2 == 0]
            df_list = []
            for i in evens:
                df_list.append(df.iloc[: , sequence[i]:sequence[i+1]])
            df = pd.concat(df_list, axis=1, join="inner")
        else: 
            df = df.iloc[: , sequence[0]:sequence[1]]
    logging.debug("   plotting")
    length_threshold = threshold
    if len(df.columns) > length_threshold:
        logging.debug(f"    sequence > {length_threshold} AA, multiple plots")
        chunks = len(uniprot_numbering)/length_threshold
        chunks = math.ceil(chunks)            
        for i in range(chunks): 
            logging.debug(f"    plotting chunk {i+1} of {chunks}")
            df_chunk = df.iloc[:,i*length_threshold:(i*length_threshold)+length_threshold]
            fig, ax = plt.subplots(figsize=((len(df_chunk.columns)/6),(len(df_chunk.index)/4)))
            ax = sns.heatmap(df_chunk, vmin=0, vmax=5, cmap=ListedColormap(['white', 'orange', 'yellow', 'lightblue', 'darkblue', 'red']), 
                             linewidths=1, linecolor='black', cbar = False, xticklabels=True, yticklabels=True)            
            ax.xaxis.tick_top()
            plt.xticks(rotation=90)
            if "cluster_id" not in file.columns: 
                plt.savefig(f"{path}/{uniprot_id}/{name}_coverage_plot_chunk_{i}.pdf", format="pdf", bbox_inches="tight")
            else:
                plt.savefig(f"{path}/{uniprot_id}/{name}_coverage_plot_cluster_{cluster}_chunk_{i}.pdf", format="pdf", bbox_inches="tight")
    else:
        logging.debug(f"    sequence < {length_threshold} AA, one plot")
        fig, ax = plt.subplots(figsize=((len(df.columns)/6),(len(df.index)/4)))
        ax = sns.heatmap(df, vmin=0, vmax=5, cmap=ListedColormap(['white', 'orange', 'yellow', 'lightblue', 'darkblue', 'red']), 
                         linewidths=1, linecolor='black', cbar = False, xticklabels=True, yticklabels=True)            
        ax.xaxis.tick_top()
        plt.xticks(rotation=90)
        if "cluster_id" not in file.columns: 
            plt.savefig(f"{path}/{uniprot_id}/{name}_coverage_plot.pdf", format="pdf", bbox_inches="tight")
        else:
            plt.savefig(f"{path}/{uniprot_id}/{name}_coverage_plot_cluster_{cluster}.pdf", format="pdf", bbox_inches="tight")
    return

def run(df, path, threshold, sequence, bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor):
    """
    A run file to generate the plots based on the original input file.

    Parameters
    ----------
    df : Dataframe of the inputfile.
    path : The path to the output files.
    threshold : Max number of residues plotted on the y-axis
    sequence : User input of area to be plotted.
    bfactortheshold_best : A integer defining a value where lower b-factors are
                           considered the best quality.
    bfactortheshold_good : A integer defining a value where lower b-factors are
                           considered good quality.
    bfactortheshold_poor : A integer defining a value where higher b-factors are
                           considered poor quality.

    Returns
    -------
    None.

    """
    for i in range(len(df)): 
        logging.debug(df.iloc[i].uniprot)
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.json"):
            logging.debug("filtered file identified")
            filtered = pd.read_json(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.json")
            plot_coverage(path, df.iloc[i].uniprot_isoform, df.iloc[i].uniprot, 
                          df.iloc[i].mutations, df.iloc[i].cluster_id, filtered, 
                          "filtered", sequence, threshold, 
                          bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
            logging.debug("plot done")
        elif os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.csv"):
            logging.debug("filtered file identified")
            filtered = pd.read_csv(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_filtered.csv", na_values=[], keep_default_na=False)
            plot_coverage(path, df.iloc[i].uniprot_isoform, df.iloc[i].uniprot, 
                          df.iloc[i].mutations, df.iloc[i].cluster_id, filtered, 
                          "filtered", sequence, threshold, 
                          bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
            logging.debug("plot done")
        elif os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.json"):
            logging.debug("filtered cluster file identified")
            filtered = pd.read_json(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.json")
            plot_coverage(path, df.iloc[i].uniprot_isoform, df.iloc[i].uniprot, 
                          df.iloc[i].mutations, df.iloc[i].cluster_id, filtered, 
                          "filtered", sequence, threshold,
                          bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
            logging.debug("plot done")
        elif os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.csv"):
            logging.debug("filtered cluster file identified")
            filtered = pd.read_csv(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_cluster{df.iloc[i].cluster_id}_filtered.csv", na_values=[], keep_default_na=False)
            plot_coverage(path, df.iloc[i].uniprot_isoform, df.iloc[i].uniprot, 
                          df.iloc[i].mutations, df.iloc[i].cluster_id, filtered, 
                          "filtered", sequence, threshold,
                          bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
            logging.debug("plot done")
        if os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.json"):
            logging.debug("All file identified")
            full = pd.read_json(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.json")
            plot_coverage(path, df.iloc[i].uniprot_isoform, df.iloc[i].uniprot, 
                          df.iloc[i].mutations, df.iloc[i].cluster_id, full, 
                          "all", sequence, threshold,
                          bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
            logging.debug("plot done")
        elif os.path.isfile(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.csv"):
            logging.debug("All file identified")
            full = pd.read_csv(f"{path}/{df.iloc[i].uniprot}/{df.iloc[i].uniprot}_all.csv", na_values=[], keep_default_na=False)
            plot_coverage(path, df.iloc[i].uniprot_isoform, df.iloc[i].uniprot, 
                          df.iloc[i].mutations, df.iloc[i].cluster_id, full, 
                          "all", sequence, threshold,
                          bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
            logging.debug("plot done")
        
    return

def main():
    args = parser.parse_args()
    path = args.result

    log_file = f'{path}/log_coverage.txt'
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

    threshold = args.threshold

    bfactortheshold_best = args.b_factor_best
    bfactortheshold_good = args.b_factor_good
    bfactortheshold_poor = args.b_factor_poor

    if args.sequence: 
        sequence = []
        sequence_pairs = args.sequence.split(",")
        for pair in sequence_pairs:
            pair = pair.split("-")
            sequence.append(int(pair[0])-1)
            sequence.append(int(pair[1]))
    else:
        sequence = 0

    logging.debug(args)
    run(df, path, threshold, sequence, bfactortheshold_best, bfactortheshold_good, bfactortheshold_poor)
    file_handler.close()    
    if os.path.exists(log_file) and os.path.getsize(log_file) == 0:
        os.remove(log_file)
        
if __name__ == '__main__':
    main()
