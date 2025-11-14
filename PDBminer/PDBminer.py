#!/usr/bin/env python3

# Copyright (C) 2023, 2024 & 2025, Kristine Degn
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

import os
import argparse
from sys import exit
import requests
from requests.exceptions import ConnectionError
from datetime import datetime
import logging
import re

import pandas as pd
import numpy as np

import json
import shutil
import glob
from multiprocessing import Pool
import itertools

from Bio.PDB import *
from Bio import Align
from Bio import SeqIO
from io import StringIO


parser = argparse.ArgumentParser(prog = "PDBminer",
                                 usage = 'PDBminer [-h] help [-n] cores and then either [-i] input file, or [-u] uniprot id')

parser.add_argument("-i", '--inputfile',
                    metavar = 'input file',
                    type=str,
                    help='The name of the input_file if any')

parser.add_argument("-g", '--hugo_name',
                    metavar = 'hugo name',
                    type=str,
                    help='The gene name in hugo formating')

parser.add_argument("-u", '--uniprot_id',
                    metavar = 'uniprot_id',
                    type=str,
                    help='The uniprot id')

parser.add_argument("-s", '--uniprot_isoform',
                    metavar = 'uniprot_isoform',
                    type=int,
                    help='The uniprot isoform, an interger')

parser.add_argument("-m", '--mutations',
                    metavar = 'mutations',
                    type=str,
                    help='string with each mutation seperated by a ";" with "" around the expression, e.i. "E120K;I390P"')

parser.add_argument("-c", '--cluster_id',
                    metavar = 'cluster id',
                    type=int,
                    help='The mutational cluster, an interger')

parser.add_argument("-n", "--cores",
                    metavar = "cores for run",
                    type=int,
                    default = 1,
                    help="the number of cores to allocate the run, if nothing is choosen: 1 core")

parser.add_argument("-f", "--format",
                    metavar = "output format, json or csv",
                    type=str,
                    default = "json",
                    help="the output format of the run, default json")

parser.add_argument("-p", '--peptide_length',
                    metavar = 'peptide_length',
                    type=int,
                    default = 5,
                    help='The minimum length of a peptide to be aligned. It is recomended not to go below 5 due to the risk of false alignment of a very small peptide')

parser.add_argument("-rm", "--remove_mutations",
                    action='store_true',
                    help="Option to add a processed, additional output file, where PDBs with mutations have been removed")

parser.add_argument("-ri", "--remove_interactors",
                    action='store_true',
                    help="Option to add a processed, additional output file, where PDBs with protein and DNA/RNA interactors are removed.")

parser.add_argument("-a", "--file_save",
                    metavar = "file_save_strategy",
                    type=str,
                    default = "none",
                    choices=["none", "retain"],
                    help="the file save strategy of the run (mmcif file format), default is none, no files are saved. Alternative is 'retain', which keeps all downloded files. Choose none or remove")

parser.add_argument('-v', '--verbose', 
                    action='store_true', 
                    help='Enable verbose output')

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(filename='log.txt', level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s')
else:
    logging.basicConfig(filename='log.txt', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')


#####################

# INPUT HANDLING 

#####################

def fetch_uniprot_data(uniprot_id):
    """
    This function takes the uniprot assccesion number and checks that it
    does indeed exist.

    Parameters
    ----------
    uniprot_id : uniprot accession number.

    Returns
    -------
    Nothing, exist program if faulty. .

    """
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        response = requests.get(url)
        if response.status_code == 400:
            logging.error(f"ERROR: The UniProt accession number {uniprot_id} does not exist. Exiting...")
            exit(1)
        elif response.status_code == 200:
            return response.json()
        else:
            logging.error(f"ERROR: Failed to fetch data. Status code: {response.status_code}")
            exit(1)

    except requests.exceptions.RequestException as e:
        logging.error(f"ERROR: Network problem occurred: {e}")
        exit(1)

def uniprot_accession_number_type(df):
    """
    This function takes the input to PDBminer read as a pandas dataframe 
    and check if the uniprot accession numbers are correctly formatted.

    Parameters
    ----------
    df : The input dataframe generated from the input CSV.

    Returns
    -------
    problems : If the uniprot accession number is not a string, the format
    is wrong and will return a list of string(s) with the offending input.

    """
    problems = []
    
    #Uniprot Accession number check:
    problematic_uniprots = df.uniprot[~df.uniprot.apply(lambda x: isinstance(x, str))]
    for uniprot_value in problematic_uniprots:
        problems.append(f"there is a problem with the uniprot input: {uniprot_value}")
    
    # If there are issues with the uniprot input, the program should end, 
    # since this is the only mandatory element of the input.
    if problems != []:
        return problems
    
def retrieve_hugo_name(uniprot_accession_number):
    """
    This functon takes a uniprot accession number and retrieves a gene name
    from the uniprot database. If Uniprot DB is down, the program will
    end with an error.

    Parameters
    ----------
    uniprot_accession_number : Accession number. 

    Returns
    -------
    hugo_name : The listed gene name for the uniprot accession number.

    """
    
    #try to get the hugoname from the uniprot database. 
    # If the uniprot datbase is down, the program will exit.
    try: 
        response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_accession_number}.json")
    except ConnectionError as e:
        logging.error(f"EXITING: Could not connect to Uniprot database via API for {uniprot_accession_number} when finding the Hugo Name.")
        exit(1)
    
    # If the response is successfull, hugo name is collected.
    if response.status_code == 200:
        uniprot_file_content = json.loads(response.text)
        if 'genes'in uniprot_file_content and 'geneName'in uniprot_file_content['genes'][0]:
            hugo_name = uniprot_file_content['genes'][0]['geneName']['value']
        else:
            hugo_name = "NA"
            logging.warning(f"Uniprot has records in this id, but no Hugo Name could be retrived with UniProt for {uniprot_accession_number}. Gene name assigned NA.")
    # Otherwise the field is left empty
    else: 
        hugo_name = "NA"
        logging.warning(f"No Uniprot record in this id, no Hugo Name could be retrived with UniProt for {uniprot_accession_number}. Gene name assigned NA.")
        
    return hugo_name
    
    
def check_input_file(df):

    problems = []
    
    # The rest of the columns are optional and may be populated in this function.
    for new_column_label in ["hugo_name", "mutations", "cluster_id"]:
        if not new_column_label in df.columns:
            df[new_column_label] = 'NA'
    
    if not 'uniprot_isoform' in df.columns:
        df['uniprot_isoform'] = 1

    for i, r in df.iterrows():
        #UNIPROT
        data = fetch_uniprot_data(r['uniprot'])

        # HUGO NAME
        if r['hugo_name'] in ["N/A", "NA"] or pd.isna(r['hugo_name']):
            hugo_name_specific = retrieve_hugo_name(r['uniprot'])
            df.at[i, "hugo_name"] = hugo_name_specific

        if not isinstance(r['hugo_name'], str):
            problems.append(f"there is a problem with the hugo name input on the line of {r['uniprot']}")

        # ISOFORM
        try:
            df.at[i, "uniprot_isoform"] = int(r['uniprot_isoform'])
        except ValueError:
            problems.append(f"there is a problem with the isoform input on the line of {r['uniprot']}")
        
        #MUTATIONS
        if r['mutations'] in ["N/A", "NA"] or pd.isna(r['mutations']):
            df.at[i, "mutations"] = "NA"
        elif not isinstance(r['mutations'], str):
            problems.append(f"there is a problem with the mutation input on the line of {r['uniprot']}")
        
        #CLUSTER ID
        if r['cluster_id'] in ["N/A", "NA"] or pd.isna(r['cluster_id']):
            df.at[i, "cluster_id"] = "NA"
        else:
            try:
                int(r['cluster_id'])
            except ValueError:
                problems.append(f"there is a problem with the cluster_id input on the line of {r['uniprot']}")
        
    df = df[['hugo_name', 'uniprot', 'uniprot_isoform', 'mutations', 'cluster_id']]
    
    return df, problems

################################

# LOCATING RELEVANT STRUCTURES 

################################


def get_alphafold_basics(uniprot_id, uniprot_isoform=None):
    logging.debug(f"FUNCTION: get_alphafold_basics({uniprot_id})")
    """
    Function that takes a uniprot id and retrieve data from the alphafold
    database, and prepare the basic information of that model in alignment 
    with the PDB data.

    Parameters
    ----------
    uniprot_id : The uniprot accession number. 

    Returns
    -------
    A tuple of information fitting as a line in the structure_df captured 
    in get_structure_df. 

    """
    
    if uniprot_isoform is not None:
        acc = f"{uniprot_id}-{uniprot_isoform}"
    else:
        acc = uniprot_id
    
    try:
        response = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
    except ConnectionError:
        logging.error(f"WARNING: Could not connect to AlphaFold database API for {acc}.")
        return
    
    # If the response is successfull, data on the model is collected.
    if response.status_code == 200:
        result = response.json()[0]
        deposition_date = result['modelCreatedDate'] 
        Alphafold_ID = result['pdbUrl'].split('/')[-1][:-4]

        return Alphafold_ID, uniprot_id, deposition_date, "PREDICTED", "NA", 0
    
    # If the response was unsuccessfull, e.g. if there are no alphafold structure 
    # in the AlphaFold database, a warning is written to the log.txt file.
    else:
        logging.warning(f"The Alphafold Database returned an error for the request of {uniprot_id}.")
        return

    
def get_pdbs(uniprot_id):
    logging.debug(f"FUNCTION: get_pdbs({uniprot_id})")
    """
    Function is taken from SLiMfast, documentation and comments there.
    Credit: Valentina Sora
    
    """
    #try to get the list of structures associated with a uniprot accession number. 
    # If the Uniprot is down, the program will exit.
    
    try: 
        response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.txt")
    except ConnectionError as e:
        logging.error(f"WARNING: Could not connect to Uniprot database API for {uniprot_id}.")
    
    # If the response is successfull, data on the models are collected.
    if response.status_code == 200:

        pdbs = []        
        for line in response.text.split("\n"):
            if line.startswith("DR   PDB;"):
                db, pdb_id, exp, res, chain_res = \
                    [item.strip(" ") for item \
                     in line.rstrip(".\n").split(";")]
                pdbs.append(pdb_id)    
        return pdbs
    
    # If the response was unsuccessful, e.g. if there are no pdb structures 
    # associated in UniPort, a warning is written to the log.txt file.
    else:
        logging.warning(f"The Uniprot Database returned an error for the request of {uniprot_id}.")

def get_structure_metadata(pdb_id, uniprot):
    logging.debug(f"FUNCTION: get_structure_metadata({pdb_id})")
    
    """
    Function that takes each pdb_id and retrive metadata from the PDBe.
    The metadata consist of desposition date to the PDB, the experimental
    metod used to solve the structure and the resolution if reported. 

    Parameters
    ----------
    pdb_id : four letter code defining a structure reported in the protein 
             data bank.

    Returns
    -------
    A tuple of metadata: deposition_date, experimental_method, resolution

    """
    #try to get the basic metadata on the pdb id from PDBe. 
    # If the PDBe Database is down, the program will exit.
    try: 
        response = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}")
    except ConnectionError as e:
        logging.error(f"WARNING: Could not connect to PDBe database via API for {pdb_id}, {uniprot}.")   
    
    # If the response was unsuccessfull, e.g. if there are no metadata 
    # associated to a structure, a warning is written to the log.txt file.
    if response.status_code != 200:
        logging.warning(f"The PDBe Database returned an error for the request of summary data for {pdb_id}, {uniprot}.")
        return
    
    # If the response is successfull, data on the metadata on the models are collected.
    response_text = json.loads(response.text)
    metadata_dictionary = response_text[pdb_id.lower()]
    metadata_dictionary = metadata_dictionary[0]

    #Change the date format
    deposition_date = f"{metadata_dictionary['deposition_date'][:4]}-{metadata_dictionary['deposition_date'][4:6]}-{metadata_dictionary['deposition_date'][6:]}"
    #Find the experimental method
    experimental_method = str(metadata_dictionary['experimental_method'][0]).upper()
    
    #Retrieve information regarding resolution
    
    #exlude NMR structures (all NMR types)
    if "NMR" in experimental_method:
        
        resolution = "NA"
        #Would be nice to have ResProx here, but there is not an API.
    
    #include all others
    else:
        #try to get the detailed metadata on the pdb id from PDBe. 
        # If the PDBe Database is down, the program will exit.
        try: 
            response_experiment = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{pdb_id}")
        except ConnectionError as e:
            logging.error(f"WARNING: Could not connect to PDB database via API (experiments) for {pdb_id}, {uniprot}.")
        
        # If the response was unsuccessfull, e.g. if there are no experimental metadata 
        # associated to a structure, a warning is written to the log.txt file.
        if response_experiment.status_code != 200:
            logging.warning(f"The PDBe Database returned an error for the request of eperimental data of {pdb_id}, {uniprot}.")
            return

        response_text_exp = json.loads(response_experiment.text)
        dictionary_exp = response_text_exp[pdb_id.lower()]
        dictionary_exp = dictionary_exp[0]
        resolution = dictionary_exp['resolution']
    
    return deposition_date, experimental_method, resolution

def get_PDBredo(pdb, uniprot_id):
    logging.debug(f"FUNCTION: get_PDBredo({pdb}, {uniprot_id})")
    """
    
    A function that takes the PDB id and checks the availability and r-free of
    the structure in the PDBREDO data bank. 

    Parameters
    ----------
    pdb : four letter PDB code.

    Returns
    -------
    str: YES/NO for availability in PDB-REDO database.
    rfree_improve: a string detailing the PDBredo r-free value or "NA".
    """
    
    try: 
        response = requests.get(f"https://pdb-redo.eu/db/{pdb}/data.json")
    except ConnectionError as e:
        logging.error(f"WARNING: Could not connect to PDB-REDO database via API for {pdb}, {uniprot_id}. ")
        return "NO", "NA"     
    
    if response.status_code == 200: 
        response_data = response.json()
        r_free_pdbredo = response_data['properties']['RFFIN']
        return "YES", r_free_pdbredo
    else:
        return "NO", "NA"     
    
def get_structure_df(uniprot_id): 
    """
    This function takes a single uniprot ID and outputs a 
    dataframe containing a sorted list of PDB ids and their metadata. 

    Parameters
    ----------
    uniprot_id : The uniprot accession number. 

    Returns
    -------
    structure_df :  A pandas dataframe containing all the names and details
                    regarding the solved structures related to the 
                    uniprot id in terms of pdb files.  

    """
    
    logging.debug(f"FUNCTION: get_structure_df({uniprot_id})")
    
    pdb = []
    deposition_date = []
    experimental_method = [] 
    resolution = []
    
    # Attempt to get data from 3D-Beacons
    # As of October 2023 the 3D becons API does not always return a code 200
    # even for entries on the website. So far this has been the case for 
    # a set amount of the test structures. It is seemingly not at random. 
    try: 
        response = requests.get(f"https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/{uniprot_id}.json?provider=pdbe")
    except ConnectionError as e:
        logging.error(f"WARNING: Could not connect to 3D-Beacons database via API for {uniprot_id}.")

    if response.status_code == 200:
        # Process data from 3D-Beacons
        response_text = json.loads(response.text)
        structures = response_text['structures']
           
        for structure in structures:
            pdb.append(structure['summary']['model_identifier'].upper())
            deposition_date.append(structure['summary']['created'])
            experimental_method.append(structure['summary']['experimental_method'])
            resolution.append(structure['summary']['resolution'])
                
    else:
        # Try to get data from Uniprot
        pdbs = get_pdbs(uniprot_id)
        
        if len(pdbs) != 0:
            logging.warning(f"Uniprot returned {len(pdbs)} structures. NOTICE that structures deposited in the PDB within 8 weeks may not be included in this list, {uniprot_id}.")
            
            for pdb_id in pdbs:
                metadata = get_structure_metadata(pdb_id, uniprot_id)

                if metadata is None:
                    logging.warning(f"PDBe database API retured an error for {pdb_id}, {uniprot_id}.")

                else:       
                    pdb.append(pdb_id)
                    deposition_date.append(metadata[0])
                    experimental_method.append(metadata[1]) 
                    resolution.append(metadata[2])
        else:
            logging.warning(f"Uniprot did not return any structures for {uniprot_id}")
            
    structure_df = pd.DataFrame({"pdb": pdb, 
                                 "uniprot_id": [uniprot_id]*len(pdb), 
                                 "deposition_date": deposition_date, 
                                 "experimental_method": experimental_method, 
                                 "resolution": resolution}) 
        
    rank_dict = {'X-RAY DIFFRACTION': 1,'ELECTRON MICROSCOPY': 2,'ELECTRON CRYSTALLOGRAPHY': 3,'SOLUTION NMR': 4, 'SOLID-STATE NMR': 5}
    
    structure_df['method_priority'] = structure_df['experimental_method'].map(rank_dict).fillna(6)
    
    AF_model = get_alphafold_basics(uniprot_id, args.uniprot_isoform)
            
    if AF_model is not None:
        structure_df.loc[len(structure_df)] = tuple(AF_model)
    else:
        logging.warning(f"AlphaFold database returned an error for {uniprot_id}. This may indicate that there are no structure for {uniprot_id} in the Alphafold Database.")
                
    structure_df.sort_values(["method_priority", "resolution", "deposition_date"], ascending=[True, True, False], inplace=True)
    structure_df = structure_df.drop(['method_priority'], axis=1)
    
    structure_df['PDBREDOdb'] = structure_df.apply(lambda row: get_PDBredo(row['pdb'], row['uniprot_id']), axis=1)

    if structure_df.empty:
        logging.error(f"No structures found for UniProt ID {uniprot_id}. Skipping.")
    
    try:
        structure_df[['PDBREDOdb', 'PDBREDOdb_details']] = structure_df['PDBREDOdb'].apply(pd.Series)
    except ValueError as e:
        logging.warning(f"Skipping PDBREDO expansion for {uniprot_id}: {e}")
                                                            
    return structure_df   

def find_structure_list(input_dataframe):  
    logging.debug("FUNCTION: find_structure_list(input_dataframe)")
    """
    Takes the input file and outputs
    a directory with a csv file for each uniprot id input and a txt file 
    including all the uniprot ids that does not have any solved structures.
    
    parameters
    ------------
    input_dataframe             The input df, as described in the readme file.  
    
    Returns          
    --------------
    log.txt:             Containing the uniprot id string which there 
                                are no solved structurs.
    
    found_structure_list:       A pandas datafrane where each solved structure
                                and a number of describtors are detailed. 
    """

    df_collector = []
    
    #take all uniprot id's from the input file
    all_uniprot_ids = list(input_dataframe.uniprot)
    all_uniprot_ids = sorted(set(all_uniprot_ids), key=all_uniprot_ids.index)
           
    for row in range(len(all_uniprot_ids)):
        
        logging.debug(all_uniprot_ids[row])
        
        structure_info = get_structure_df(all_uniprot_ids[row]) 
    
        if type(structure_info) != str: 
            df_collector.append(structure_info)
        
        else:
            logging.warning(f"No structures found in any resource for {structure_info}, {all_uniprot_ids[row]}. ")
            
    if len(df_collector) > 0:
        found_structure_list = pd.concat(df_collector) 
        
    else:
        found_structure_list = []
        
    return found_structure_list


def combine_structure_dfs(input_dataframe, found_structures):
    logging.debug("FUNCTION: combine_structure_dfs(found_structures, input_dataframe)")
    """
    This function takes the found structures and the input dataframe
    and combine these for continues computation. 

    Parameters
    ----------
    found_structures : A pandas dataframe created in step 2. 
    input_dataframe :  The original input dataframe with mutational information.

    Returns
    -------
    final_df :  A dataframe with nine columns including hugo_name, uniprot_id,
                uniprot_isoform, mutations, cluster_id, structure_id, deposition_date
                experimental_method, resolution

    """
    
    input_dataframe = input_dataframe.rename(columns={"uniprot": "uniprot_id"}) 
    # Merge the input_dataframe with found_structures based on uniprot_id
    merged_df = input_dataframe.merge(found_structures, on=['uniprot_id'])

    # Create the final DataFrame with the required columns
    final_df = pd.DataFrame({
        "hugo_name": merged_df['hugo_name'],
        "uniprot_id": merged_df['uniprot_id'],
        "uniprot_isoform": merged_df['uniprot_isoform'],
        "mutations": merged_df['mutations'],
        "cluster_id": merged_df['cluster_id'],
        "structure_id": merged_df['pdb'],
        "deposition_date": merged_df['deposition_date'],
        "experimental_method": merged_df['experimental_method'],
        "resolution": merged_df['resolution'],
        "PDBREDOdb": merged_df['PDBREDOdb'],
        "PDBREDOdb_rfree": merged_df['PDBREDOdb_details']
    })

    return final_df


################################

# ALIGNING RELEVANT STRUCTURES 
# TO THE ISOFORM DETERMINED 
# FASTA

################################

def to_ranges(iterable):
    logging.debug("FUNCTION: to_ranges()")
    """
    This function converts a list of integers into a list of ranges.

    This function groups consecutive integers into ranges.
    For example, [1,2,3,4,5,7,8,9,13,14,15] becomes [(1, 5), (7, 9), (13, 15)].

    Parameters
    ----------
    iterable : list of integers
        Input list of integers to be converted into ranges.

    Returns
    -------
    ranges : list of tuples
        List of tuples representing ranges.

    Example
    -------
    >>> to_ranges([1, 2, 3, 4, 5, 7, 8, 9, 13, 14, 15])
    [(1, 5), (7, 9), (13, 15)]
    """
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
                                        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]
        

def alignment(uniprot_sequence, AA_pdb, peptide_min_length, b_factors, uniprot_numbering):
    """
    A function that aligns the uniprot and AA sequence from the PDB.

    Parameters
    ----------
    uniprot_sequence : A canonical string of amino acid letters
    AA_pdb : A list of amino acid one letter codes covered by the PDB 

    Returns
    -------
    uniprot_aligned : A string of letters. 
    pdb_aligned : A string of letters including "-" for areas not covered by 
    the pdb.

    """
    pdb_sequence = ''.join(AA_pdb)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -5
    aligner.open_gap_score = -20
    aligner.extend_gap_score = -10
    
    # X is code for any amino acid and will be used here to indicate a gap, 
    # which represent the missing residues in the structure, as  
    # Align.PairwiseAligner() does not take "-" as input. 
    uniprot_sequence = uniprot_sequence.replace("-", "X" )
    pdb_sequence = pdb_sequence.replace("-", "X" )
    
    alignments = aligner.align(uniprot_sequence, pdb_sequence)
    
    if alignments.score < peptide_min_length: 
        return None
    
    a = format(alignments[0])
    a = a.split("\n")    
    
    uniprot_seq_len = len(uniprot_sequence)
    
    #as the uniprot sequence as target and pdb is query, we get the
    #values for the start and end positions of the alignment
    target_strings = [string for string in a if 'target' in string]
    query_strings = [string for string in a if 'query' in string]
    alignment_strings = [string for string in a if not string.startswith(('target', 'query')) and string.startswith(' ')]
    
    #get where the areas are locally aligned:
    start_value_uniprot = int(re.findall(r'\d+', target_strings[0])[0])
    start_value_pdb = int(re.findall(r'\d+', query_strings[0])[0])
    end_value_uniprot = int(re.findall(r'\d+', target_strings[-1])[-1])
    end_value_pdb = int(re.findall(r'\d+', query_strings[-1])[-1])
    
    #combine to one dataframe:
    dfs = []
    for i in range(len(target_strings)):
        data = {'pdb_seq': list(query_strings[i]), 'alignment': list(alignment_strings[i]), 'uniprot_seq': list(target_strings[i])}
        df = pd.DataFrame(data)
        df = df[20:] #removes "target, "query and these values"
        dfs.append(df)
    
    df_all = pd.concat(dfs)
    #mask numbers
    mask_numbers = df_all.apply(lambda x: x.astype(str).str.isdigit())
    #mask empty rows
    mask_empty = df_all.apply(lambda x: x.astype(str).str.strip().eq('')).all(axis=1)
    # Use the combined mask to filter out rows with numbers or empty rows
    df_all = df_all[~(mask_numbers.any(axis=1) | mask_empty)]
    # Reset the index if needed
    df_all.reset_index(drop=True, inplace=True)
    
    #get the uniprot seqpos + possible insertions
    uniport_seqpos = [f'{item1}{item2}' for item1, item2 in zip(list(uniprot_sequence[start_value_uniprot:end_value_uniprot]), list(uniprot_numbering[start_value_uniprot:end_value_uniprot]))]
    dash_positions = [i for i, char in enumerate(list(df_all['uniprot_seq'])) if char == '-']
    
    for position in dash_positions:
        uniport_seqpos.insert(position, '-')

    #this might be a problem for instertions
    df_all['uniprot_sequence'] = uniport_seqpos
    df_all[['uniprot_letter', 'uniprot_pos']] = df_all['uniprot_sequence'].str.extract(r'([A-Za-z]+)(\d+)')
    
    #get b-factors
    dash_positions = [i for i, char in enumerate(list(df_all['pdb_seq'])) if char == '-']
    relevant_b_factors = b_factors[start_value_pdb:end_value_pdb]

    for position in dash_positions:
        relevant_b_factors.insert(position, '-')
    
    df_all['b_factor']=relevant_b_factors
    
    df_all.replace("X", "-", inplace=True)

    df_align = df_all[['uniprot_seq', 'uniprot_pos', 'pdb_seq', 'b_factor']]
    
    df_align = df_align.reset_index(drop=True)

    return df_align

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
    if isoform is None:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}-{isoform}.fasta"

    try:
        response = requests.get(url)  
    except requests.exceptions.RequestException as e:
        logging.error(f"EXITING: connection error when trying to connect to Uniprot database API, {uniprot_id}, isoform {isoform}.")
        raise RuntimeError(
            f"Connection error when trying to connect to Uniprot for "
            f"{uniprot_id}, isoform {isoform}: {e}")
   
    if response.status_code == 200:
        sequence_data=''.join(response.text)
        Seq=StringIO(sequence_data)
        pSeq=list(SeqIO.parse(Seq,'fasta'))
        uniprot_sequence = str(pSeq[0].seq)
        uniprot_numbering = list(range(1,len(uniprot_sequence)+1,1)) 
    else:
        logging.error(f"EXITING: The sequence could not be retrieved, ensure that isoform {isoform} exist, {uniprot_id}.")
        raise RuntimeError(
            f"The sequence could not be retrieved, ensure that isoform "
            f"{isoform} exists for {uniprot_id} (status {response.status_code}).")

    return uniprot_sequence, uniprot_numbering
    
def get_mutations(mut_pos, df):
    """
    A function finding the amino acids in the structure at the sites of 
    mutation in the input file.

    Parameters
    ----------
    mut_pos : A list of intergers representing mutational sites.
    df : A pandas dataframe with the aligned sequence.

    Returns
    -------
    muts : A list of strings with the sites and amino acids, e.g. ['R249H','..'] 

    """
    logging.debug("FUNCTION: get_mutations(mut_pos, df)")
    mut_pos = list(map(int, mut_pos))
    
    muts = []
    
    filtered_df = df[df['uniprot_pos'].isin(mut_pos)]

    for _, row in filtered_df.iterrows():
        mutation = row['uniprot_seq'] + str(row['uniprot_pos']) + row['pdb_seq']
        muts.append(mutation)

    return muts

def get_all_discrepancies(df, residue_names): 
    """
    A function that takes the aligned dataframe and identify problems if any. 

    Parameters
    ----------
    df : The aligened uniprot and pdb structure sequence. 

    Returns
    -------
    df : An updated input dataframe.
    all_mutations : The mutations in the PDB compared to the uniprot sequence, list.
    warnings : Qualitative warinings regarding the structure, list of strings. 

    """
    logging.debug("FUNCTION: get_all_discrepancies(df, pdb_residues)")

    warnings = []
        
    df = df[df.uniprot_seq != "-"]
    df = df[df.pdb_seq != "-"]        

    if len(df[df["uniprot_seq"] != df["pdb_seq"]]) != 0:
        all_mutations = df[df["uniprot_seq"] != df["pdb_seq"]].apply(lambda row: f"{row['uniprot_seq']}{row['uniprot_pos']}{row['pdb_seq']}", axis=1).tolist()
    else:
        all_mutations = []
        
    pdb_residues = "".join(residue_names)
    # check if these are part of a tag
    if "HHH" in pdb_residues[:20]:
        warnings.append("The structure likely contains an expression tag")
        
    if "HHH" in pdb_residues[-20:]:
        warnings.append("The structure likely contains an expression tag")
    
    df = df.reset_index(drop=True)

    return df, all_mutations, warnings

def get_coverage(df):
    """
    A function that findes the area the protein structure covers the canonical 
    sequence from uniprot allowing amino acid substitutions.

    Parameters
    ----------
    df :  The aligned uniprot and pdb structure sequence. 

    Returns
    -------
    ranges_covered : A string with the coverage ranges. 

    """
    logging.debug("FUNCTION: get_coverage(df)")
    
    #list(df.uniprot_pos) example: [1, 2, 3, 4, 5, 7, 8, 9, 13, 14, 15]
    
    df['uniprot_pos'] = df['uniprot_pos'].astype(int)
    
    ranges_covered = list(to_ranges(list(df.uniprot_pos)))
    
    #ranges covered example: [(1, 5), (7, 9), (13, 15)]
    ranges_covered = str(ranges_covered)
    #example: "[(1, 5), (7, 9), (13, 15)]"
    ranges_covered = ranges_covered.replace("), (", ");(" )
    #example: "[(1, 5);(7, 9);(13, 15)]"
    return ranges_covered


def processing_normal_alignment(uniprot_sequence, residue_names, uniprot_numbering, b_factors, peptide_min_length):    
    """
    Function that takes the sequence from uniprot and pdb and create an aligned dataframe
    and finds mutations.

    Parameters
    ----------
    uniprot_sequence : A list of amino acids from uniprot.
    residue_names : A list of amino acids from the structure.
    uniprot_numbering : A list of numbers pertaining to the uniprot sequence.
    b_factors : A list of CA b-factors pertaining to the pdb sequence.

    Returns
    -------
    df_align : A dataframe with the uniprot sequence, uniprot positions, the 
    aligned pdb seuqnde and b-factor.
    mutations_in_all : The mutations in the pdb file.
    warnings : A list of potential warnings. 

    """
    logging.debug("FUNCTION: processing_normal_alignment(uniprot_sequence, residue_names, uniprot_numbering, b_factors)")
    
    #Align the PDB sequence to the uniprot sequence.
    df_align = alignment(uniprot_sequence, residue_names, peptide_min_length, b_factors, uniprot_numbering)   
    
    if type(df_align) != pd.core.frame.DataFrame:
        return "NA", "NA", "NA"
    
    #Find mutations and create a list of warnings. 
    df_align, mutations_in_all, warnings = get_all_discrepancies(df_align, residue_names)

    return df_align, mutations_in_all, warnings


def align_uniprot_pdb(pdb_id, uniprot_sequence, uniprot_numbering, mut_pos, 
                      path, self_chains, uniprot_id, complex_status, 
                      complex_details, peptide_min_length):
    """
    This function takes the pdbs and align to the uniprot sequence per chain. 
    Fusion products is handled in a different way than single protein chains. 

    Parameters
    ----------
    pdb_id : Four letter protein data bank ID.
    uniprot_sequence : A string of the amino acid sequence.
    uniprot_numbering : A list of each number corresponding to the letters in
    the amino acid sequence.
    mut_pos : The query mutations in the input file.
    path : Where to download the structures to.
    self_chains : A list of the chains associated with the uniprot accession number.
    uniprot_id : The uniprot accession number.
    complex_status : The complex status, such as fusion product.
    complex_details : The deatils of the complex or fusion, including the fused chain,
    and the uniprot accession number associated.

    Returns
    -------
    output_array: An array of mixed data for each pdb, indluding dictionaries
    for the b-fcator, mutations in the sequence and the r-free value.

    """    
    logging.debug(f"FUNCTION: align_uniprot_pdb({pdb_id}, uniprot_sequence, uniprot_numbering, mut_pos, path, {self_chains}, {uniprot_id}, complex_status, complex_details)")
    
    allowed_amino_acids = ['D','V','A','M','N','F','G','L','I',
                           'T','P','C','S','R','Y','Q','K','W',
                           'E','H', 'X']

    #download the mmcif of the pdb_id
    download_mmcif(pdb_id, path)
    
    #import the mmcif into a dataframe per chain.
    chain_dfs, chain_names, r_free = import_mmcif(pdb_id, self_chains, path)
    
    #The function returns empty if the struture could not be downloaded. 
    if chain_dfs == "NA":
        return ['']
    
    # Empty lists for popultion by function
    muts = {}
    coverage = {}
    mutation_list = {}
    warning_column = {}
    b_factor_dictionary = {}
    
    for i, chain in enumerate(chain_names): 
        #For each self chain remove residues named "X" 
        logging.debug(f"CHAIN: {chain}")
        chain_df = chain_dfs[i]
        chain_df = chain_df[chain_df['residue_name'].isin(allowed_amino_acids)]        
        chain_df = chain_df.replace("X", "-")
        
        if len(chain_df) < peptide_min_length:
            #peptide_min_length recommendended to be 5
            #remove very small fragments
            output_array = np.array(['', '', '', '', '', ''], dtype=object)
            return output_array
        
        #insert "-" for missing residues to aid the alignment. 
        lst = list(chain_df.seq_num)
        missing_residues = sorted(set(range(lst[0], lst[-1])) - set(lst))
        
        if missing_residues != []:       
            insertion_df = pd.DataFrame({"residue_name":["-"]*len(missing_residues), "seq_num":missing_residues, "b_factor":["-"]*len(missing_residues)})        
            chain_df = pd.concat([chain_df, insertion_df]).sort_values(by="seq_num")
            
        #If there are no residue name, there may still be a b-fcator. This tend
        #to be the case with "X" residues and with additions to the structure 
        #that is not an amino acid.
        # Overwrite b_factor with "-" where residue_name is "-"
        chain_df.loc[chain_df["residue_name"] == "-", "b_factor"] = "-"
        
        uniprot_numbering_copy = uniprot_numbering.copy() 
        df_align, mutations_in_all, warnings = processing_normal_alignment(uniprot_sequence, 
                                                                           chain_df.residue_name.to_list(),
                                                                           uniprot_numbering_copy, 
                                                                           chain_df.b_factor.to_list(),
                                                                           peptide_min_length)
        #input mutations
        if type(mut_pos) == float or type(mut_pos) == np.float64:
            mut_pos = "NA"
            
        if mut_pos not in ["NA", "N/A"]:
            
            mutation_positions = [x[1:-1] for x in mut_pos.split(';')]
            mutational_sites = get_mutations(mutation_positions, df_align)
            muts[chain] = mutational_sites       
        
        
        warnings = ";".join(warnings)
        if len(warnings) > 0 and warnings[0] == ";":
            warning_column[chain] = warnings
                
        if mutations_in_all == []:
            mutation_list[chain] = "NA"
        else:
            mutation_list[chain] = mutations_in_all

        #capture the area the PDB covers according to the alignment. 
        if isinstance(df_align, pd.DataFrame):
            ranges_covered = get_coverage(df_align)
            coverage[chain] = ranges_covered
            b_factor_dictionary[chain] = np.array((df_align[df_align.b_factor != "-"].b_factor).astype(float))
        else:
            coverage[chain] = "NA"
            b_factor_dictionary[chain] = "NA"

    if warning_column == {}:
        warning_column = "NA"
        
    chains_string = ';'.join([str(elem) for elem in chain_names])
    output_array = np.array([chains_string, coverage, muts, mutation_list, b_factor_dictionary, warning_column, r_free], dtype=object)

    return output_array

def download_mmcif(structure_identifier, path):
    """
    A function that downloads a mmcif file from a known source.

    Parameters
    ----------
    structure_identifier : string of pdb id or alphafold identifier

    Returns
    -------
    None.
    
    This functions downloads the mmCIF file of the structure to be analyzed.

    """
    #check if structures directory exist
    if os.path.isdir(f"{path}/structures") == False:
        os.mkdir(f"{path}/structures")
    
    #check if file already exists
    if os.path.isfile(f"structures/{structure_identifier}.cif") == False:
        
        if structure_identifier.startswith("AF-"):
            #download the alphafold file - we know the structure exist. If ther is a download issue, 
            #this is handled in the calling function. 
            try: 
                response = requests.get(f"https://alphafold.ebi.ac.uk/files/{structure_identifier}.cif")
            except ConnectionError as e:
                logging.error(f"WARNING: Could not connect to AlphaFold database via API for {structure_identifier}, {path}.")
            
            if response.status_code == 200:
                with open(f'{path}/structures/{structure_identifier}.cif', 'w') as fh:
                    fh.write(response.text)
            
        else:
            #download the pdb file - we know the structure exist. If ther is a download issue, 
            #this is handled in the calling function.             
            try: 
                response = requests.get(f"https://files.rcsb.org/download/{structure_identifier}.cif")
            except ConnectionError as e:
                logging.error(f"WARNING: Could not connect to PDB via API for {structure_identifier}, {path}.")
            
            if response.status_code == 200:
                with open(f'{path}/structures/{structure_identifier}.cif', 'w') as fh:
                    fh.write(response.text)
    
    return

def three_to_one(amino_acid):
    amino_acid = amino_acid.upper()
    amino_acid_mapping = {
        'ASP': 'D', 'VAL': 'V', 'ALA': 'A', 'MET': 'M', 'ASN': 'N',
        'PHE': 'F', 'GLY': 'G', 'LEU': 'L', 'ILE': 'I', 'THR': 'T',
        'PRO': 'P', 'CYS': 'C', 'SER': 'S', 'ARG': 'R', 'TYR': 'Y',
        'GLN': 'Q', 'LYS': 'K', 'TRP': 'W', 'GLU': 'E', 'HIS': 'H'}
    return amino_acid_mapping.get(amino_acid, None)

def import_mmcif(structure_identifier, self_chains, path):
    """
    A function that imports the downloaded mmcif file and handles the 
    lack of a file.

    Parameters
    ----------
    structure_identifier : string of pdb id or alphafold identifier
    self_chains : chains reated to the protein of interest

    Returns
    -------
    chain_dfs : A list of dataframes one for each self chain, containing; 
                residue_name: 1 letter amino acid
                seq_num: residue number in structure.
                b-factor: the b-factor of the CA.
    chain_names : A list of chains found. 

    """
    #import the mmcif dictionary
    
    if os.path.isfile(f"{path}/structures/{structure_identifier}.cif") == False:
        return "NA","NA","NA"
        
    mmcif_dict = MMCIF2Dict.MMCIF2Dict(f"{path}/structures/{structure_identifier}.cif")
    
    #create a dataframe 
    df = pd.DataFrame({"ATOM": mmcif_dict['_atom_site.label_atom_id'], "seq_num": mmcif_dict['_atom_site.auth_seq_id'],
                       "residue": mmcif_dict['_atom_site.auth_comp_id'], "chain": mmcif_dict['_atom_site.auth_asym_id'],
                       "occupancy": mmcif_dict['_atom_site.occupancy'], "b_factor": mmcif_dict['_atom_site.B_iso_or_equiv'],
                       "model_num": mmcif_dict['_atom_site.pdbx_PDB_model_num']})
       
    #reduce to only the carbon alpha
    df = df[df.ATOM == "CA"]
    
    #ensure we choose the one with the highest occupancy
    df = df.sort_values(by='occupancy', ascending=False)
    df = df.drop_duplicates(subset=['ATOM', 'seq_num', 'residue', 'chain'], keep="first")
    df.seq_num = df.seq_num.astype(int)
    df = df.sort_values(by="seq_num")
    
    #prepare putput dataframes
    chain_dfs = []
    chain_names = []
    chains = df.groupby("chain")
    
    for chain_name, chain_df in chains:
        if chain_name in self_chains:
            chain_names.append(chain_name)
            chain_df['residue_name'] = [three_to_one(aa) for aa in chain_df["residue"]]
            chain_df = chain_df[["residue_name", "seq_num", "b_factor"]]
            df.b_factor = df.b_factor.astype(float)
            chain_df = chain_df[chain_df.residue_name != '']
            chain_dfs.append(chain_df)

    try:
        r_free = float(mmcif_dict['_refine.ls_R_factor_R_free'][0])
    except ValueError:
        r_free = "NA" #In cases where r-free has been entered as: ?/NA/X/Yes/No/none and more.
    except KeyError: # relevant for AlphaFold structures
        r_free = "NA"

    return chain_dfs, chain_names, r_free

def align_alphafold(alphafold_id, mutations, path):
    logging.debug("FUNCTION: align_alphafold(alphafold_id, mutations, path)")
    """
    This function takes an alphafold ID and the mutational positions. 
    The aim is to find the high quality areas of the protein, and set these in 
    relation to the mutations. 

    Parameters
    ----------
    alphafold_id : String of the ID
    mutations : List of mutations

    Returns
    -------
    coverage : [(x, y)] per chain high quality areas (pDDLT > 70)
    AA_in_PDB : If the high quality portions of the AF structure covers 
    mutations. 

    """
    download_mmcif(alphafold_id, path)
    
    chain_dfs, chain_names, r_free = import_mmcif(alphafold_id, ["A"], path)
    
    #since there only is one chain
    AF_df = chain_dfs[0]
    
    AF_df["seq_num"] = AF_df["seq_num"].astype(int)
    AF_df["b_factor"] = AF_df["b_factor"].astype(float)

    confidence_categories = ["high" if i > 70 else "low" for i in AF_df.b_factor]
    AF_df["category"] = confidence_categories

    confident_seq = AF_df[AF_df.category == "high"].seq_num.values
    
    if len(confident_seq) > 0:
        # Create coverage string (PDBminer output style)
        f = [int(confident_seq[0])]
        for i in range(len(confident_seq) - 1):
            if int(confident_seq[i]) + 1 != int(confident_seq[i + 1]):
                f.extend([int(confident_seq[i]), int(confident_seq[i + 1])])
        f.append(int(confident_seq[-1]))
        f = np.array(f)
    
        coverage = [(f[i], f[i + 1]) for i in range(0, len(f), 2)]
        coverage = str(coverage).replace("), (", ");(")
        coverage = dict({"A":coverage})
    
        AA_in_PDB = []
        
        if type(mutations) == float or type(mutations) == np.float64: 
            mutations="NA"

        if mutations not in ["NA", "N/A"]:
            mutation_positions = [x[1:-1] for x in mutations.split(';')]
            for mutant in mutation_positions:
                mutation = f"{AF_df.residue_name[AF_df.seq_num == int(mutant)].values[0]}{int(mutant)}{AF_df.residue_name[AF_df.seq_num == int(mutant)].values[0]}"
            
            AA_in_PDB.append(mutation)
    
        AA_in_PDB = ",".join(AA_in_PDB)
        AA_in_PDB = dict({"A": AA_in_PDB})

    else: 
        coverage = "NA"
        AA_in_PDB = "NA"
    
    b_factor_dict = dict({"A": list(AF_df.b_factor)})
    
    return coverage, AA_in_PDB, b_factor_dict, r_free

def align(combined_structure, path, peptide_min_length):
    logging.debug("FUNCTION: align(combined_structure, path, peptide_min_length)")
    """
    This functions takes the pandas dataframe containing all the structures, 
    their metadata and information regarding mutations of interest, import 
    the relevant fastafiles and conduct alignment, which captures information
    regarding the structure in terms of mutations. All is outputted as 
    additional columns in the input file.

    Parameters
    ----------
    combined_structure : Pandas dataframe created in step 3.
    path : string, directory of interest.

    Returns
    -------
    combined_structure : Updated input. 
    """
    
    combined_structure["chains"] = " "
    combined_structure["coverage"] = " "
    combined_structure["AA_in_PDB"] = " "
    combined_structure["mutations_in_pdb"] = " "
    combined_structure["b_factor"] = " "
    combined_structure["warnings"] = " "
    combined_structure["r_free"] = " "
    
    iso_raw = combined_structure['uniprot_isoform'].iloc[0]
    if pd.isna(iso_raw) or str(iso_raw).strip().upper() in ("", "NA", "N/A"):
        iso_norm = None
    else:
        iso_norm = int(iso_raw)

    uniprot_sequence, uniprot_numbering = get_uniprot_sequence(combined_structure.uniprot_id.iloc[0], iso_norm)     

    for i in range(len(combined_structure)):
                           
        if combined_structure['structure_id'][i].startswith("AF-"):
            alignment_info = align_alphafold(combined_structure['structure_id'][i], combined_structure['mutations'][i], path)
            
            combined_structure.at[i, 'chains'] = "A"  
            combined_structure.at[i, 'coverage'] = alignment_info[0] 
            combined_structure.at[i, 'AA_in_PDB'] = alignment_info[1] 
            combined_structure.at[i, 'mutations_in_pdb'] = "NA"
            combined_structure.at[i, 'b_factor'] = alignment_info[2]
            combined_structure.at[i, 'warnings'] = "NA"
            combined_structure.at[i, 'r_free'] = alignment_info[3]
            
        else:    
        #uniprot_sequence, uniprot_numbering = get_uniprot_sequence(combined_structure.uniprot_id[i], int(combined_structure['uniprot_isoform'][i]))     
            alignment_info = align_uniprot_pdb(combined_structure.structure_id[i],
                                           uniprot_sequence, 
                                           uniprot_numbering,
                                           combined_structure['mutations'][i],
                                           path,
                                           combined_structure.self_chains[i],
                                           combined_structure.uniprot_id[i],
                                           combined_structure.complex_protein[i],
                                           combined_structure.complex_protein_details[i],
                                           peptide_min_length)            
            if alignment_info[0] != '':
                combined_structure.at[i, 'chains'] = alignment_info[0]  
                combined_structure.at[i, 'coverage'] = alignment_info[1] 
                combined_structure.at[i, 'AA_in_PDB'] = alignment_info[2] 
                combined_structure.at[i, 'mutations_in_pdb'] = alignment_info[3]
                combined_structure.at[i, 'b_factor'] = alignment_info[4]
                combined_structure.at[i, 'warnings'] = alignment_info[5]
                combined_structure.at[i, 'r_free'] = alignment_info[6]
            #drop missing values. 
            else:
                combined_structure = combined_structure.drop([i])
    
    combined_structure = combined_structure.reset_index(drop=True)  
                        
    return combined_structure

################################

# FINDING INFORMATION ON NON-
# SELF ENTRIES IN THE PDB FILE.

################################

def get_complex_information(pdb_id, uniprot):
    logging.debug(f"FUNCTION: get_complex_information(pdb_id, uniprot), {pdb_id}, {uniprot}")
    """
    This function takes a PDB id and analyzes its content to estabilish if 
    there is any other elements within the file such as a ligand. 
    

    Parameters
    ----------
    pdb_id : Four letter code, string. 

    Returns
    -------
    output_array: A np.array contoning of 
                1) protein_complex_list: binary  
                2) protein_info; description of complex if any
                3) self_chain; a double check list of chains annotated as self
                4) nucleotide_complex_list: binary 
                5) nuleotide_info: description of complex if any 
                6) ligand_complex_list: binart 
                7) ligand_info: description of complex if any. Include metal.

    """    
    #finding protein complexes and their related uniprot_ids 
    #this step also serves as a quality control of uniprot id's and chains.
    
    try: 
        response = requests.get(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{pdb_id}")
    except ConnectionError as e:
        logging.error(f"WARNING: Could not connect to PDBe database via API for {pdb_id}, {uniprot}.")
    
    if response.status_code == 200:
        response_text = json.loads(response.text)
        protein_segment_dictionary = response_text[pdb_id.lower()]
        
        self_chain = []
    
        if len(protein_segment_dictionary['UniProt']) <= 1:
            protein_complex_list = "NA"
            protein_info = "NA"    
            
            #check that uniprot found the correct PDB:
            if list(protein_segment_dictionary['UniProt'].keys())[0] != uniprot:
                output_array = np.array(['NA', 'NA', 'NA', 'NA', 
                                         'NA', 'NA', 'NA'], dtype=object)
                return output_array
            
            else:               
                #allow isoforms
                for i in next(v for k,v in protein_segment_dictionary['UniProt'].items() if uniprot in k)['mappings']:
                    self_chain.append(i['chain_id'])
                self_chain = list(set(self_chain))
                
        else: 
            if uniprot not in list(protein_segment_dictionary['UniProt'].keys()):
                output_array = np.array(['NA', 'NA', 'NA', 'NA', 
                                         'NA', 'NA', 'NA'], dtype=object)
                return output_array

            else:            
                #allow isoforms
                for i in next(v for k,v in protein_segment_dictionary['UniProt'].items() if uniprot in k)['mappings']:
                    self_chain.append(i['chain_id'])
                self_chain = list(set(self_chain))
                
                info = []
                fusion_test = []
                for i in protein_segment_dictionary['UniProt']:
                    prot_info = f"{protein_segment_dictionary['UniProt'][i]['identifier']}, {i}, chain_{protein_segment_dictionary['UniProt'][i]['mappings'][0]['chain_id']}"
                    info.append(prot_info)
                    
                for item in enumerate(info):
                    fusion_test.append(item[1][-1])
                    if uniprot in item[1]:
                        value = item[0]
                
                if len(set(fusion_test)) == 1:
                    protein_complex_list = 'fusion product'
                elif len(fusion_test) > len(set(fusion_test)):
                    if fusion_test.count(fusion_test[value]) > 1:
                        protein_complex_list = 'fusion product in protein complex'
                    else:
                        protein_complex_list = 'protein complex with fusion product'
                else:
                    protein_complex_list = 'protein complex'
            
                info = ';'.join(info)
                protein_info = [info] 

    else:
        protein_complex_list = "NA"
        protein_info = "NA"
    
    #finding complexes with other ligands by identifying other molecules
    try: 
        response = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}")
    except ConnectionError as e:
        logging.error(f"WARNING: Could not connect to PDBe database via API for molecules,{pdb_id}, {uniprot}.")
        
    if response.status_code == 200: 
        response_text = json.loads(response.text)
        molecule_dictionary = response_text[pdb_id.lower()]

        nucleotide_info = []
        ligand_info = []
        
        for i in range(len(molecule_dictionary)):
            
            if "nucleotide" in molecule_dictionary[i]['molecule_type']:
                n1 = molecule_dictionary[i]['molecule_name'][0]
                n2 = molecule_dictionary[i]['in_chains']
                n2 = ",".join(n2)
                n = f"{n1}, chain {n2}"
                
                nucleotide_info.append(n)
            
            elif molecule_dictionary[i]['molecule_type'] != 'polypeptide(L)':
                if molecule_dictionary[i]['molecule_type'] != 'water': 
                    l1 = molecule_dictionary[i]['molecule_name'][0]
                    l2 = molecule_dictionary[i]['in_chains']
                    l2 = ",".join(l2)
                    l = f"{l1}, chain {l2}"
                    
                    ligand_info.append(l)
       
        if len(nucleotide_info) > 0:
            nucleotide_complex_list = 'nucleotide complex'
            if len(nucleotide_info) > 1: 
                nucleotide_info = ';'.join(nucleotide_info)
        
        else: 
            nucleotide_complex_list = "NA"
            nucleotide_info = "NA"
               
        if len(ligand_info) > 0:
            ligand_complex_list = "Other ligands"
            if len(ligand_info) > 1: 
                ligand_info = ';'.join(ligand_info)
        
        else: 
            ligand_complex_list = "NA"
            ligand_info = "NA"
    
    else:
        nucleotide_complex_list = "NA"
        nucleotide_info = "NA"
        ligand_complex_list = "NA"
        ligand_info = "NA"
        
    output_array = np.array([protein_complex_list, protein_info, self_chain, nucleotide_complex_list, 
                             nucleotide_info, ligand_complex_list, ligand_info], dtype=object)
    
    return output_array

def collect_complex_info(structural_df):
    logging.debug("FUNCTION: collect_complex_info(structural_df)")
    """
    A function that parses all pdbids though the get_complex_information
    function and capture and merge with the input file. 

    Parameters
    ----------
    structural_df: The pandas dataframe containing the original input file,
                   PDBs and subsequent additional columns. 
    
    Returns
    -------
    df_combined: A format of the structural_df including complex information. 

    """
    
    uniprot_id = structural_df['uniprot_id'][0]
    
    df = pd.DataFrame(columns=['structure_id','complex_protein',
                           'complex_protein_details', 'self_chains', 'complex_nucleotide',
                           'complex_nucleotide_details','complex_ligand', 
                           'complex_ligand_details'])

    for pdb in set(structural_df['structure_id']):
        
        if pdb.startswith("AF-"):
            list_of_values = [pdb,'NA',
                              'NA',['A'],
                              'NA','NA',
                              'NA','NA']
            
        else:
            complex_info = get_complex_information(pdb, uniprot_id)
    
            list_of_values = [pdb, complex_info[0], 
                              complex_info[1], complex_info[2], 
                              complex_info[3], complex_info[4], 
                              complex_info[5], complex_info[6]]
            
        if list_of_values[3] != 'NA':
        
            df.loc[pdb] = np.array(list_of_values, dtype="object") 
    
    df_combined = pd.merge(structural_df, df, how='inner', on = 'structure_id')
    
    return df_combined

################################

# FILTER AND EXPORT BASED ON 
# INPUT

################################

def cleanup_all(structural_df,input_dataframe, output_format):
    logging.debug("FUNCTION: cleanup_all(structural_df,input_dataframe, output_format)")
    if structural_df.empty:
        uniprot_id = input_dataframe.iloc[0]['uniprot']
        logging.warning(f"No structures for {uniprot_id} — writing header-only output.")
        output_path = f"results/{uniprot_id}"
        os.makedirs(output_path, exist_ok=True)

        if output_format == "json":
            structural_df.to_json(f"{output_path}/{uniprot_id}_all.json", orient="columns")
        else:
            structural_df.to_csv(f"{output_path}/{uniprot_id}_all.csv", index=False)
        return

    uniprot_id = structural_df.iloc[0]['uniprot_id']
    #rearrange columns
    structural_df = structural_df[['hugo_name', 'uniprot_id', 'uniprot_isoform', 
                                   'mutations', 'cluster_id', 'structure_id', 
                                   'deposition_date', 'experimental_method', 
                                   'resolution', 'r_free','PDBREDOdb', 
                                   'PDBREDOdb_rfree', 'complex_protein',
                                   'complex_protein_details', 'complex_nucleotide',
                                   'complex_nucleotide_details', 'complex_ligand',
                                   'complex_ligand_details', 'chains', 'coverage', 
                                   'mutations_in_pdb','AA_in_PDB','b_factor', 'warnings']]
    
    df = structural_df.reset_index(drop=True)
    
    #if there are no coverage in any chain, remove the structure: 
    for index, row in df.iterrows():
        all_values_empty = all(value == '[]' for value in row['coverage'].values())
        if all_values_empty:
            df = df.drop([index])

    df = df.reset_index(drop=True)
    
    # Iterate through mutations_in_pdb
    for i,v in enumerate(df.mutations_in_pdb):
        if isinstance(v, dict):
            if all(value == 'NA' for value in v.values()):
                df.loc[i,'mutations_in_pdb'] = "NA"
    
    # Iterate through warnings
    for i in range(len(df.warnings)):
        if isinstance(df.warnings[i], dict):
            for key, value in df.warnings[i].items():
                if all(c in ['', ' ', ';'] for c in value):
                    df.warnings[i][key] = 'NA'
            if all(value == 'NA' for value in df.warnings[i].values()):
                df.loc[i,'warnings'] = "NA"
    
    key = df.map(str).agg('|'.join, axis=1)
    df = df.loc[~key.duplicated(keep='first')].reset_index(drop=True)

    df.index.name = 'structure_rank'
    
    #### find if there are mutations in the file containing non numbers.
    if type(df.mutations[0]) != float or type(df.mutations[0]) != np.float64:
        filtered_df = df.copy()
        filtered_df = filtered_df[filtered_df['AA_in_PDB'].apply(lambda d: any(d.values()))]
        
        if type(df.cluster_id[0]) == float or type(df.cluster_id[0]) == np.float64:
            df.drop('cluster_id', inplace=True, axis=1)  
        
        if len(filtered_df) > 0:
            if output_format == "json":
                filtered_df.to_json(f"results/{uniprot_id}/{uniprot_id}_filtered.json", orient="columns")
            else:
                filtered_df.to_csv(f"results/{uniprot_id}/{uniprot_id}_filtered.csv")
        
    if 'cluster_id' in df.columns:
        grouped = df.groupby(['cluster_id'])
        for cluster_id, cluster_data in grouped:
            cluster_data = cluster_data[cluster_data['AA_in_PDB'].apply(lambda d: any(d.values()))]
            if len(cluster_data) > 0:
                if output_format == "json":
                    cluster_data.to_json(f"results/{uniprot_id}/{uniprot_id}_cluster_{cluster_id}.json", orient="columns")
                else: 
                    cluster_data.to_csv(f"results/{uniprot_id}/{uniprot_id}_cluster_{cluster_id}.csv")
        
    if 'mutations' in df.columns:
        df.drop('mutations', inplace=True, axis=1)
        df.drop('AA_in_PDB', inplace=True, axis=1)
    
    if 'cluster_id' in df.columns:
        df.drop('cluster_id', inplace=True, axis=1)
    
    if output_format == "json":
        df.to_json(f"results/{uniprot_id}/{uniprot_id}_all.json", orient="columns")
    else:
        df.to_csv(f"results/{uniprot_id}/{uniprot_id}_all.csv")
        
    #post processing the file:
    filtering_applied = False

    #remove non-self interactors. Proteins in homodimers and similar is not listed as a complex, 
    #but rather with multiple chains covering the sequence. 
    #PDBminer does not know about the direction and configuration - e.g. if 
    #multiple chains are a true dimer or just the result of the asymmetic unit.
    if args.remove_interactors:
        df = df[(df.complex_nucleotide == "NA") & (df.complex_protein == "NA")]
        filtering_applied = True  
    
    #remove PDB structures harbouring mutations:
    if args.remove_mutations:
        df = df[df.mutations_in_pdb == "NA"]
        filtering_applied = True  
    
    #export file, if there has been any processing:
    if filtering_applied:
        if output_format == "json":
            df.to_json(f"results/{uniprot_id}/{uniprot_id}_processed.json", orient="columns")
        else:
            df.to_csv(f"results/{uniprot_id}/{uniprot_id}_processed.csv")
        
    return 

################################

# RUNNING THE TOOL

################################

def process_uniprot(uniprot_id, df, output_format, peptide_min_length, file_save_strategy):
    uniprot_dir = f"results/{uniprot_id}"
    os.makedirs(uniprot_dir, exist_ok=True)
    
    # Create a copy of the DataFrame for this process
    uniprot_dataframe = df[df['uniprot'] == uniprot_id].reset_index(drop=True)

    found_structures = find_structure_list(uniprot_dataframe)

    if len(found_structures) != 0:
        combined_structure = combine_structure_dfs(uniprot_dataframe, found_structures)
        combined_structure = collect_complex_info(combined_structure)
        structural_df = align(combined_structure, uniprot_dir, peptide_min_length)
    else:
        structural_df = pd.DataFrame(columns=[
        'hugo_name','uniprot_id','uniprot_isoform','mutations','cluster_id',
        'structure_id','deposition_date','experimental_method','resolution',
        'r_free','PDBREDOdb','PDBREDOdb_rfree','complex_protein',
        'complex_protein_details','complex_nucleotide','complex_nucleotide_details',
        'complex_ligand','complex_ligand_details','chains','coverage',
        'mutations_in_pdb','AA_in_PDB','b_factor','warnings'])
    
    cleanup_all(structural_df, uniprot_dataframe, output_format)
    
    if file_save_strategy == "none" and os.path.isdir(f"{uniprot_dir}/structures"):
        shutil.rmtree(f"{uniprot_dir}/structures")

def run_list(input_file, cores, output_format, peptide_min_length, file_save_strategy):
    df = pd.read_csv(input_file)
    uniprot_list = df['uniprot'].unique()

    os.makedirs("results", exist_ok=True)

    with Pool(processes=cores) as pool:
        #use starmap to add multiple arguments.
        try:    
            pool.starmap(process_uniprot, [(uniprot_id, df, output_format, peptide_min_length, file_save_strategy) for uniprot_id in uniprot_list])
        except Exception as e:
            logging.error(f"Fatal error in worker: {e}")
            exit(1)

def main():
    start = datetime.now()
    logging.info("PDBminer is starting")
    
    if args.inputfile:
        input_file = args.inputfile
        logging.info("input file mode is choosen and input file has been identified.")

        df = pd.read_csv(input_file, na_values=[], keep_default_na=False)
        
        uniprot_problems = uniprot_accession_number_type(df)
        if uniprot_problems is not None:
            logging.error(f"ERROR: The input file is not correctly formated. {uniprot_problems}. Exiting...")
            exit(1)
        
        df,problems = check_input_file(df)

        if problems == []:
            df.to_csv(input_file, index=False)

        elif len(problems) > 0:
            logging.error(f"ERROR: The input file is not correctly formated. {problems}. Exiting...")
            exit(1)

    elif args.uniprot_id:
        uniprot_id = args.uniprot_id
        data = fetch_uniprot_data(uniprot_id)
        #will exit if problem, data only occur if the response is 200.

        if args.hugo_name:
            hugo_name = args.hugo_name
        else:  
            hugo_name = retrieve_hugo_name(uniprot_id)
        
        df = pd.DataFrame({'hugo_name':[hugo_name],'uniprot':[uniprot_id]})
        
        if args.uniprot_isoform:
            df['uniprot_isoform'] = [args.uniprot_isoform]
        else:
            df['uniprot_isoform'] = [None]

        if args.mutations:
            mutations = [args.mutations]
            df['mutations'] = mutations
        else: df['mutations'] = ["NA"]
        
        if args.cluster_id:
            cluster = [args.cluster_id]
            df['cluster_id'] = cluster
        
        else: df['cluster_id'] = ["NA"]
        df.to_csv('input_file.csv')
        input_file = 'input_file.csv'
        logging.info("input file created based on hugo name and uniprot id, file name: input_file.csv")

    else:
        logging.error("ERROR: input file or input hugo name and uniprot id is missing. Exiting...")
        exit(1)
        
    logging.info("Pipeline is starting")
    
    run_list(input_file, args.cores, args.format, args.peptide_length, args.file_save)
    
    finished = datetime.now() 
    
    logging.info(f"total runtime: {finished-start}")
    logging.info("PDBminer has finished.")

if __name__ == '__main__':
    main()

