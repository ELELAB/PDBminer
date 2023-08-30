#!/usr/bin/env python

# PDBminer_functions: functions for PDBminer_run.py
# Copyright (C) 2023, Kristine Degn
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

#============================================================================#
#Importing relevant packages
#============================================================================#
import os
import pandas as pd
import requests
import json
import itertools
import numpy as np
from Bio.PDB import *
from Bio import pairwise2
from Bio import SeqIO
from io import StringIO
from Bio.SeqUtils import seq1
from requests.exceptions import ConnectionError

def get_alphafold_basics(uniprot_id):
    print(f"FUNCTION: get_alphafold_basics({uniprot_id})")
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
    
    #try to get the metadata on the alphafold ID. 
    # If the AlphaFold Database is down, the program will exit.
    try: 
        response = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}")
    except ConnectionError as e:
        with open("log.txt", "a") as textfile:
            textfile.write(f"EXITING: Could not connect to AlphaFold database API for {uniprot_id}.\n")
        exit(1)
    
    # If the response is successfull, data on the model is collected.
    if response.status_code == 200:
        result = response.json()[0]
        deposition_date = result['modelCreatedDate'] 
        Alphafold_ID = result['pdbUrl'].split('/')[-1][:-4]

        return Alphafold_ID, uniprot_id, deposition_date, "PREDICTED", "NA", 0
    
    # If the response was unsuccessfull, e.g. if there are no alphafold structure 
    # in the AlphaFold database, a warning is written to the log.txt file.
    else:
        with open("log.txt", "a") as textfile:
            textfile.write(f"WARNING: The Alphafold Database returned an error for the request of {uniprot_id}.\n")
        return

    
def get_pdbs(uniprot_id):
    print(f"FUNCTION: get_pdbs({uniprot_id})")
    """
    Function is taken from SLiMfast, documentation and comments there.
    Credit: Valentina Sora
    
    """
    #try to get the list of structures associated with a uniprot accession number. 
    # If the Uniprot is down, the program will exit.
    
    try: 
        response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.txt")
    except ConnectionError as e:
        with open("log.txt", "a") as textfile: #"a" to append.
            textfile.write(f"EXITING: Could not connect to Uniprot database API for {uniprot_id}. \n")
        exit(1) 
    
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
    
    # If the response was unsuccessfull, e.g. if there are no pdb structures 
    # associated in UniPort, a warning is written to the log.txt file.
    else:
        with open("log.txt", "a") as textfile: #"a" to append.
            textfile.write(f"WARNING: The Uniprot Database returned an error for the request of {uniprot_id}.\n")

def get_structure_metadata(pdb_id):
    print(f"FUNCTION: get_structure_metadata({pdb_id})")
    
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
        with open("log.txt", "a") as textfile:
            textfile.write(f"EXITING: Could not connect to PDBe database via API for {pdb_id}. \n")
        exit(1) 
    
    
    # If the response was unsuccessfull, e.g. if there are no metadata 
    # associated to a structure, a warning is written to the log.txt file.
    if response.status_code != 200:
        with open("log.txt", "a") as textfile: #"a" to append.
            textfile.write(f"WARNING: The PDBe Database returned an error for the request of summary data for {pdb_id}.\n")
        return
    
    # If the response is successfull, data on the metadata on the models are collected.
    response_text = json.loads(response.text)
    dictionary = response_text[pdb_id.lower()]
    dictionary = dictionary[0]

    #Change the date format
    deposition_date = f"{dictionary['deposition_date'][:4]}-{dictionary['deposition_date'][4:6]}-{dictionary['deposition_date'][6:]}"
    #Find the experimental method
    experimental_method = str(dictionary['experimental_method'][0]).upper()
    
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
            with open("log.txt", "a") as textfile:
                textfile.write(f"EXITING: Could not connect to PDB database via API (experiments) for {pdb_id}. \n")
            exit(1) 
        
        # If the response was unsuccessfull, e.g. if there are no experimental metadata 
        # associated to a structure, a warning is written to the log.txt file.
        if response_experiment.status_code != 200:
            with open("log.txt", "a") as textfile: #"a" to append.
                textfile.write(f"WARNING: The PDBe Database returned an error for the request of eperimental data of {pdb_id}.\n")
            return

        response_text_exp = json.loads(response_experiment.text)
        dictionary_exp = response_text_exp[pdb_id.lower()]
        dictionary_exp = dictionary_exp[0]
        resolution = dictionary_exp['resolution']
    
    return deposition_date, experimental_method, resolution

def get_PDBredo(pdb):
    print(f"FUNCTION: get_PDBredo({pdb})")
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
        with open("log.txt", "a") as textfile:
            textfile.write(f"EXITING: Could not connect to PDB-REDO database via API for {pdb_id}. \n")
        exit(1) 
    
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
    
    print(f"FUNCTION: get_structure_df({uniprot_id})")
    
    try: 
        response = requests.get(f"https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/{uniprot_id}.json?provider=pdbe")
    except ConnectionError as e:
        with open("log.txt", "a") as textfile:
            textfile.write(f"EXITING: Could not connect to 3D-Beacons database via API for {uniprot_id}.\n")

        exit(1) 
    
    # As of August 2023 the 3D becons API is unstable, and routinely return an empty json file. 
    # to overcome this limitation, a warning is issues and the uniprot API is accessed.
    if response.status_code != 200:
        with open("log.txt", "a") as textfile:
            textfile.write(f"WARNING: 3D-Beacons did not return any PDBe structures for {uniprot_id}.\n")
        pdbs = get_pdbs(uniprot_id)
        
        if len(pdbs) != 0:
            with open("log.txt", "a") as textfile:
                textfile.write(f"WARNING: Uniprot returned {len(pdbs)} structures. NOTICE that structures deposited in the PDB within 8 weeks may not be included in this list.\n")
            
            pdb = []
            deposition_date = []
            experimental_method = [] 
            resolution = []
            
            for pdb_id in pdbs:
                metadata = get_structure_metadata(pdb_id)

                if metadata is None:
                    textfile.write(f"WARNING: PDBe database API retured an error for {pdb_id}.\n")
                    return uniprot_id
            
                pdb.append(pdb_id)
                deposition_date.append(metadata[0])
                experimental_method.append(metadata[1]) 
                resolution.append(metadata[2])
            
            structure_df = pd.DataFrame({"pdb": pdb, 
                                         "uniprot_id": [uniprot_id]*len(pdb), 
                                         "deposition_date": deposition_date, 
                                         "experimental_method": experimental_method, 
                                         "resolution": resolution}) 
            
            rank_dict = {'X-RAY DIFFRACTION': 1,'ELECTRON MICROSCOPY': 2,'ELECTRON CRYSTALLOGRAPHY': 3,'SOLUTION NMR': 4, 'SOLID-STATE NMR': 5}
            
            structure_df['method_priority'] = structure_df['experimental_method'].map(rank_dict).fillna(6)
            
            AF_model = get_alphafold_basics(uniprot_id)
            
            if AF_model is not None:
                structure_df.loc[len(structure_df)] = tuple(AF_model)
                        
            structure_df.sort_values(["method_priority", "resolution", "deposition_date"], ascending=[True, True, False], inplace=True)
            structure_df = structure_df.drop(['method_priority'], axis=1)
            
            structure_df['PDBREDOdb'] = structure_df.apply(lambda row: get_PDBredo(row['pdb']), axis=1)
            structure_df[['PDBREDOdb', 'PDBREDOdb_details']] = structure_df['PDBREDOdb'].apply(pd.Series)
            
            structure_df = structure_df.set_index('pdb')
            
        else:
            with open("log.txt", "a") as textfile:
                textfile.write(f"WARNING: Uniprot did not return any structures for {uniprot_id}.\n")
            AF_model = get_alphafold_basics(uniprot_id)
            if AF_model is None:
                textfile.write(f"WARNING: AlphaFold database returned an error for {uniprot_id}. This may indicate that there are no structure for {uniprot_id} in the Alphafold Database.\n")
                return uniprot_id
            
            structure_df = pd.DataFrame({'uniprot_id': [AF_model[1]], 
                                         'deposition_date': [AF_model[2]], 
                                         'experimental_method': [AF_model[3]], 
                                         'resolution': [AF_model[4]],
                                         'PDBREDOdb': "NO", 
                                         'PDBREDOdb_details': "NA"})

            structure_df.index = [AF_model[0]]           
        
        return structure_df

    else: 
        response_text = json.loads(response.text)
        structures = response_text['structures']
           
        pdb = []
        deposition_date = []
        experimental_method = [] 
        resolution = []
    
        for structure in structures:
            pdb.append(structure['summary']['model_identifier'].upper())
            deposition_date.append(structure['summary']['created'])
            experimental_method.append(structure['summary']['experimental_method'])
            resolution.append(structure['summary']['resolution'])

        structure_df = pd.DataFrame({"pdb": pdb, "uniprot_id": [uniprot_id]*len(pdb), "deposition_date": deposition_date, "experimental_method": experimental_method, "resolution": resolution})    
    
        rank_dict = {'X-RAY DIFFRACTION': 1,'ELECTRON MICROSCOPY': 2,'ELECTRON CRYSTALLOGRAPHY': 3,'SOLUTION NMR': 4, 'SOLID-STATE NMR': 5}
        
        structure_df['method_priority'] = structure_df['experimental_method'].map(rank_dict).fillna(6)
        
        AF_model = get_alphafold_basics(uniprot_id)
        
        if AF_model is not None:
            structure_df.loc[len(structure_df)] = tuple(AF_model)
        
        structure_df.sort_values(["method_priority", "resolution", "deposition_date"], ascending=[True, True, False], inplace=True)
        structure_df = structure_df.drop(['method_priority'], axis=1)
        
        structure_df['PDBREDOdb'] = structure_df.apply(lambda row: get_PDBredo(row['pdb']), axis=1)
        structure_df[['PDBREDOdb', 'PDBREDOdb_details']] = structure_df['PDBREDOdb'].apply(pd.Series)
        
        structure_df = structure_df.set_index('pdb')
                                                        
    return structure_df     

def find_structure_list(input_dataframe):  
    print("FUNCTION: find_structure_list(input_dataframe)")
    """
    Takes the input file and the path where it is placed and outputs
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
        
        print(all_uniprot_ids[row])
        
        structure_info = get_structure_df(all_uniprot_ids[row]) 
    
        if type(structure_info) != str: 
            df_collector.append(structure_info)
        
        else:
            with open("log.txt", "w") as textfile:
                textfile.write(f"No structures found in any resource for {structure_info}. \n")
            
    if len(df_collector) > 0:
        found_structure_list = pd.concat(df_collector) 
        
    else:
        found_structure_list = []
        
    return found_structure_list


def combine_structure_dfs(found_structures, input_dataframe):
    print("FUNCTION: combine_structure_dfs(found_structures, input_dataframe)")
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
    
    df_collector = []

    for i in range(len(input_dataframe)):
        sub_df = found_structures[found_structures.uniprot_id == input_dataframe.uniprot[i]]
        num_sub_df = len(sub_df)
        
        data = {
            "hugo_name": [input_dataframe.hugo_name[i]] * num_sub_df,
            "uniprot_id": [input_dataframe.uniprot[i]] * num_sub_df,
            "uniprot_isoform": [input_dataframe['uniprot_isoform'][i]] * num_sub_df,
            "mutations": [input_dataframe.mutations[i]] * num_sub_df,
            "cluster_id": [input_dataframe.cluster_id[i]] * num_sub_df,
            "structure_id": list(sub_df.index),
            "deposition_date": list(sub_df.deposition_date),
            "experimental_method": list(sub_df.experimental_method),
            "resolution": list(sub_df.resolution),
            "PDBREDOdb":list(sub_df.PDBREDOdb),
            "PDBREDOdb_rfree": list(sub_df.PDBREDOdb_details)
        }
        
        df_collector.append(pd.DataFrame(data))
    
    final_df = pd.concat(df_collector, ignore_index=True)
           
    return final_df

def to_ranges(iterable):
    print("FUNCTION: to_ranges()")
    """
    Function to make each mutational group iterable, called to make a range
    interable.

    Parameters
    ----------
    iterable : a range of numbers e.g. (1, 10)

    """
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
                                        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]
        
def alignment(uniprot_sequence, AA_pdb):
    """
    A function that takes the canonical sequence from uniprot and align the
    sequence from the pdb file to it. 
    
    1) Local alignment: Aim to find the area of the uniprot
       sequence the PDB covers.
    2) Match (identical amino acids) =  1 point
    3) Non-match (not identical) = -10 points
    4) Opening a gap: -20 points
    5) keeping a gap open: -10 
    
    reasoning: we may be OK with a mutation, but a missing residue
    insertion or deleting will be punished much harder.                
    The options have been chosen to force the best fit locally
    #and highly discourage gaps, as they do not make sense in 

    Parameters
    ----------
    uniprot_sequence : String of single letter amino acids. 
    AA_pdb : A string of single letter amino acids. 

    Returns
    -------
    uniprot_aligned : The updated string based on alginment. 
    pdb_aligned : The updated string based on alignment. 
    
    In both cases a gap is indicated with "-"
    """
    
    print("FUNCTION: alignment(uniprot_sequence, AA_pdb)")
   
    #this structral space. 
    pdb_sequence = ''.join(AA_pdb)
    alignments = pairwise2.align.localms(uniprot_sequence, pdb_sequence, 1, -10, -20, -10)
    
    if alignments != []: # and alignments[0][2] >= 10:
        uniprot_aligned = alignments[0][0]
        pdb_aligned = alignments[0][1]
        
    else:
        uniprot_aligned = []
        pdb_aligned = []
    
    return uniprot_aligned, pdb_aligned

def build_fusion_df(df, uniprot_id):
    print(f"FUNCTION: build_fusion_df(df, {uniprot_id})")
    """
    A function that takes the dataframe, df, with the protein alignment, 
    and aligns to the other uniprot_ids in the fusion product. This
    allow a discrimination between mutations and fusion. 

    Parameters
    ----------
    df : The dataframe, df_align
    uniprot_id : The fused uniprot accession number.

    Returns
    -------
    df : An updated input dataframe. 

    """
    #GET THE FUSED SEQUENCE:
    uniprot_sequence_fusion, uniprot_numbering_fusion = get_uniprot_sequence(uniprot_id, 1)

    #CREATE STRINGS AND LISTS FROM THE DATAFRAME
    pdb_sequence = "".join(list(df.pdb_seq))
    uniprot_aligned = "".join(list(df.uniprot_seq))
    uniprot_num = list(df.uniprot_pos)
    b_factor = list(df.b_factor)
    
    #ALIGN THE FUSION SEQUENCE TO THE PDB_SEQUENCE. 
    alignments_fusion = pairwise2.align.localms(uniprot_sequence_fusion, pdb_sequence, 1, -10, -20, -10)
    # The aligned fusion uniprot
    uniprot_fusion_aligned = alignments_fusion[0][0]
    # The aligned pdb with additional gaps
    pdb_aligned_fusion = alignments_fusion[0][1]
    #get the original uniprot to match      
    pdb_uniprot_align = pairwise2.align.localms(uniprot_aligned, pdb_aligned_fusion, 1, -10, -20, -10)
    pdb_uniprot_fusion = pdb_uniprot_align[0][0]
    #pdb_fusion_align = pairwise2.align.localms(uniprot_fusion_aligned, pdb_self_fusion, 1, -10, -20, -10)
    #pdb_fusion_align = pdb_fusion_align[0][0]
    
    for i, v in enumerate(list(pdb_uniprot_fusion)):
        if i >= len(uniprot_num) or (v == "-" and uniprot_num[i] != "-"):
            uniprot_num.insert(i, "-")
    
    for i, v in enumerate(list(pdb_aligned_fusion)):
        if i >= len(b_factor) or (v == "-" and b_factor[i] != "-"):
            b_factor.insert(i, "-")
    
    df_align = pd.DataFrame(data={"uniprot_seq": list(pdb_uniprot_fusion),
                                  "uniprot_pos": uniprot_num,
                                  "pdb_seq": list(pdb_aligned_fusion),
                                  "fusion_seq":  list(uniprot_fusion_aligned),
                                  "b_factor": b_factor})    
    
    df = df_align[df_align.uniprot_seq != "-"]
    df = df[df.pdb_seq != "-"]
    
    condition = lambda x: x['fusion_seq'] == x["pdb_seq"]
    
    if uniprot_aligned[0] == "-":
        #N-terminal fusion
        while len(df) > 0 and condition(df.iloc[0]):
            df = df.iloc[1:]
        
    elif uniprot_aligned[-1] == "-":
        #C-terminal fusion
        while len(df) > 0 and condition(df.iloc[-1]):
            df = df.iloc[:-1]
    
    df = df[["uniprot_seq","uniprot_pos",
             "pdb_seq","b_factor"]]
    
    return df

def alignment_fusion(uniprot_sequence, AA_pdb, uniprot_ids, uniprot_num, b_factor):     
    """
    A function that takes the uniprot and pdb sequence, numbers and b-factor and
    aligns as well as call a function to align towards other fused uniprots.

    Parameters
    ----------
    uniprot_sequence : A tring with amino acids
    AA_pdb : A sting of amino acids from the PDB.
    uniprot_ids : A list of fused uniprot_ids
    uniprot_num : A numberical list with the uniprot numbering.
    b_factor : A list of b-fcators related to the AA_pdb resdues.

    Returns
    -------
    A dataframe with the alignment. 

    """
    print("FUNCTION: alignment_fusion")
    
    uniprot_num_copy = uniprot_num.copy()  
    
    pdb_sequence = ''.join(AA_pdb)
    alignments = pairwise2.align.localms(uniprot_sequence, pdb_sequence, 1, -10, -20, -10)
    
    if alignments == [] or alignments[0][2] <= 10:
        uniprot_aligned = []
        pdb_aligned = []
        return []
        
    uniprot_aligned = alignments[0][0]
    pdb_aligned = alignments[0][1]

    for i, v in enumerate(list(uniprot_aligned)):
        if i >= len(uniprot_num_copy) or (v == "-" and uniprot_num_copy[i] != "-"):
            uniprot_num_copy.insert(i, "-")
    
    for i, v in enumerate(list(pdb_aligned)):
        if i >= len(b_factor) or (v == "-" and b_factor[i] != "-"):
            b_factor.insert(i, "-")
            
    # Create an aligned dataframe for the chain.
    df_aligned = pd.DataFrame(data={"uniprot_seq": list(uniprot_aligned),
                                  "uniprot_pos": uniprot_num_copy,
                                  "pdb_seq": list(pdb_aligned),
                                  "b_factor": b_factor})    
        
    for uniprot in uniprot_ids:
        df_aligned = build_fusion_df(df_aligned, uniprot)
        
    return df_aligned   

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
    print(f"FUNCTION: get_uniprot_sequence({uniprot_id}, {isoform})")
    
    if isoform==1:
    
        try: 
            response = requests.post(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
        except ConnectionError as e:
            with open("log.txt", "a") as textfile:
                textfile.write(f"EXITING: Uniprot database API rejected for {uniprot_id}. \n")
            exit(1) 
        
        if response.status_code == 200:
            sequence_data=''.join(response.text)
            Seq=StringIO(sequence_data)
            pSeq=list(SeqIO.parse(Seq,'fasta'))
            uniprot_sequence = str(pSeq[0].seq)
            uniprot_numbering = list(range(1,len(uniprot_sequence)+1,1)) 

        else:
            with open("log.txt", "a") as textfile:
                textfile.write(f"EXITING: The canonical sequence could not be retrieved. \n")
            exit(1) 
    
    else:
        
        try: 
            response = requests.post(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}-{isoform}.fasta")
        except ConnectionError as e:
            with open("log.txt", "a") as textfile:
                textfile.write(f"EXITING: Uniprot database API rejected for {uniprot_id}. \n")
            exit(1) 
            
        if response.status_code == 200:
            sequence_data=''.join(response.text)
            Seq=StringIO(sequence_data)
            pSeq=list(SeqIO.parse(Seq,'fasta'))
            uniprot_sequence = str(pSeq[0].seq)
            uniprot_numbering = list(range(1,len(uniprot_sequence)+1,1)) 
        else:
            with open("log.txt", "a") as textfile:
                textfile.write(f"EXITING: The canonical sequence could not be retrieved, ensure that isoform {isoform} exist. \n")
            exit(1) 
    
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
    print("FUNCTION: get_mutations(mut_pos, df)")
    mut_pos = list(map(int, mut_pos))
    
    muts = []
    
    filtered_df = df[df['uniprot_pos'].isin(mut_pos)]

    for _, row in filtered_df.iterrows():
        mutation = row['uniprot_seq'] + str(row['uniprot_pos']) + row['pdb_seq']
        muts.append(mutation)

    return muts

def get_all_discrepancies(df): 
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
    print("FUNCTION: get_all_discrepancies(df)")

    warnings = []
        
    df = df[df.uniprot_seq != "-"]
    df = df[df.pdb_seq != "-"]        

    if len(df[df["uniprot_seq"] != df["pdb_seq"]]) != 0:
        all_mutations = df[df["uniprot_seq"] != df["pdb_seq"]].apply(lambda row: f"{row['uniprot_seq']}{row['uniprot_pos']}{row['pdb_seq']}", axis=1).tolist()
    else:
        all_mutations = []
        
    # check if these are part of a tag
    if "HHH" in "".join(list(df.pdb_seq[:20])):
        warnings.append("The structure likely contains an expression tag")
        
    if "HHH" in "".join(list(df.pdb_seq[-20:])):
        warnings.append("The structure likely contains an expression tag")
    
    df = df.reset_index(drop=True)

    return df, all_mutations, warnings

def get_coverage(df):
    """
    A function that findes the area the protein structure covers the canonical 
    sequence from uniprot allowing amino acid substitutions.

    Parameters
    ----------
    df :  The aligened uniprot and pdb structure sequence. 

    Returns
    -------
    ranges_covered : A string with the coverage ranges. 

    """
    print("FUNCTION: get_coverage(df)")
    ranges_covered = list(to_ranges(list(df.uniprot_pos)))
    ranges_covered = str(ranges_covered)
    ranges_covered = ranges_covered.replace("), (", ");(" )
    return ranges_covered

def get_aligned_df(uniprot_al, pdb_al, uniprot_num, b_factor):
    """
    A function that builds the dataframe with the aligned sequences. 

    Parameters
    ----------
    uniprot_al : Aligned sequence uniprot string.
    pdb_al : Aligned pdb sequence string. 
    uniprot_num : The list of residue numbers in the uniprot sequence. 

    Returns
    -------
    df_align : A dataframe containing the relevant alignment. 

    """
    uniprot_al = list(uniprot_al)  
    pdb_al = list(pdb_al) 

    for i, v in enumerate(uniprot_al):
        if i >= len(uniprot_num) or (v == "-" and uniprot_num[i] != "-"):
            uniprot_num.insert(i, "-")
    
    for i, v in enumerate(pdb_al):
        if i >= len(b_factor) or (v == "-" and b_factor[i] != "-"):
            b_factor.insert(i, "-")

    # Create an aligned dataframe for the chain.
    df_align = pd.DataFrame(data={"uniprot_seq": uniprot_al,
                                  "uniprot_pos": uniprot_num,
                                  "pdb_seq": pdb_al,
                                  "b_factor": b_factor})    
    return df_align


def processing_normal_alignment(uniprot_sequence, residue_names, uniprot_numbering, b_factors):
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
    
    #Align the PDB sequence to the uniprot sequence.
    uniprot_aligned, pdb_aligned = alignment(uniprot_sequence, residue_names)    
    uniprot_numbering_copy = uniprot_numbering.copy()        
    #create a dataframe with the uniprot sequence, numbering, aligned pdb sequence and b-factor.
    df_align = get_aligned_df(uniprot_aligned, pdb_aligned, uniprot_numbering_copy, b_factors)
    
    #Find mutations and create a list of warnings. 
    df_align, mutations_in_all, warnings = get_all_discrepancies(df_align)

    return df_align, mutations_in_all, warnings

def processing_fusion_alignment(uniprot_sequence, residue_names, uniprot_numbering, b_factors, uniprot_ids):
    """
    Function that takes the sequence from uniprot and pdb and create an aligned dataframe
    and finds mutations.

    Parameters
    ----------
    uniprot_sequence : A list of amino acids from uniprot.
    residue_names : A list of amino acids from the structure.
    uniprot_numbering : A list of numbers pertaining to the uniprot sequence.
    b_factors : A list of CA b-factors pertaining to the pdb sequence.
    uniprot_ids: a list of the uniprots in the fusion.

    Returns
    -------
    df_align : A dataframe with the uniprot sequence, uniprot positions, the 
    aligned pdb seuqnde and b-factor.
    mutations_in_all : The mutations in the pdb file.
    warnings : A list of potential warnings. 

   """
    print("FUNCTION: processing_fusion_alignment")
    
    #Align the PDB sequence to the uniprot sequence.    
    uniprot_numbering_copy = uniprot_numbering.copy()
    df_align = alignment_fusion(uniprot_sequence, residue_names, uniprot_ids, uniprot_numbering_copy, b_factors) 
    if isinstance(df_align, pd.DataFrame):
        df_align, mutations_in_all, warnings = get_all_discrepancies(df_align)
        return df_align, mutations_in_all, warnings
    else:
        return df_align, "NA", "NA"   
        
def align_uniprot_pdb(pdb_id, uniprot_sequence, uniprot_numbering, mut_pos, path, self_chains, uniprot_id, complex_status, complex_details):
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
    print("===========================================")
    print(pdb_id)
    print("===========================================")
    
    print("FUNCTION: align_uniprot_pdb(pdb_id, uniprot_sequence, uniprot_numbering, mut_pos, path, self_chains, uniprot_id, complex_status, complex_details)")

    #download the mmcif of the pdb_id
    download_mmcif(pdb_id)
    
    #import the mmcif into a dataframe per chain.
    chain_dfs, chain_names, r_free = import_mmcif(pdb_id, self_chains)
    
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
        print(f"CHAIN: {chain}")
        chain_df = chain_dfs[i].replace("X", "-")
        #print(chain_df[:10])
        
        if len(chain_df) < 5:
            #remove very small fragments
            os.chdir(path)
            output_array = np.array(['', '', '', '', '', ''], dtype=object)
            return output_array
        
        #insert "-" for missing residues to aid the alignment. 
        lst = list(chain_df.seq_num)
        missing_residues = sorted(set(range(lst[0], lst[-1])) - set(lst))
        #print(missing_residues)
        
        if missing_residues != []:       
            insertion_df = pd.DataFrame({"residue_name":["-"]*len(missing_residues), "seq_num":missing_residues, "b_factor":["-"]*len(missing_residues)})        
            chain_df = pd.concat([chain_df, insertion_df]).sort_values(by="seq_num")
            
        #If there are no residue name, there may still be a b-fcator. This tend
        #to be the case with "X" residues and with additions to the structure 
        #that is not an amino acid.
        # Overwrite b_factor with "-" where residue_name is "-"
        chain_df.loc[chain_df["residue_name"] == "-", "b_factor"] = "-"
        
        #In cases of fusions, alignment can be a challenge. 
        if "fusion" in complex_status:

            #identify the other protein fused to the protein of interest.
            complex_list = complex_details[0].split(";")
            filtered_data = [item for item in complex_list if f'chain_{chain}' in item and uniprot_id not in item]
            
            if filtered_data != []:
                #popuate a list of fused uniprot_ids in the chain, that is not the
                #uniprot of interst.
                uniprot_ids = []
                for fusion in filtered_data:
                    uniprot_ids.append(fusion.split(", ")[1])
                uniprot_numbering_copy = uniprot_numbering.copy() 
                df_align, mutations_in_all, warnings = processing_fusion_alignment(uniprot_sequence, 
                                                                                   list(chain_df.residue_name), 
                                                                                   uniprot_numbering_copy, 
                                                                                   list(chain_df.b_factor),
                                                                                   uniprot_ids)                
            else:
                uniprot_numbering_copy = uniprot_numbering.copy() 
                df_align, mutations_in_all, warnings = processing_normal_alignment(uniprot_sequence, 
                                                                                   list(chain_df.residue_name), 
                                                                                   uniprot_numbering_copy, 
                                                                                   list(chain_df.b_factor))

        else:

            uniprot_numbering_copy = uniprot_numbering.copy() 
            df_align, mutations_in_all, warnings = processing_normal_alignment(uniprot_sequence, 
                                                                               list(chain_df.residue_name), 
                                                                               uniprot_numbering_copy, 
                                                                               list(chain_df.b_factor))
        
        #input mutations
        if type(mut_pos) == float or type(mut_pos) == np.float64:
            mut_pos = "NA"
        #if np.isnan(mut_pos):
        #    mut_pos="NA"
            
        if mut_pos not in ["NA", "N/A"]:
            
            mutation_positions = [x[1:-1] for x in mut_pos.split(';')]
            mutational_sites = get_mutations(mutation_positions, df_align)
            muts[chain] = mutational_sites       
        
        
        warnings = ";".join(warnings)
        if len(warnings) > 0 and warnings[0] == ";":
            warnings = warnings[1:]
            warning_column[chain] = warnings
                
        if mutations_in_all == []:
            mutation_list[chain] = "NA"
        else:
            mutation_list[chain] = mutations_in_all

        #capture the area the PDB covers accourding to the alignment. 
        if isinstance(df_align, pd.DataFrame):
            ranges_covered = get_coverage(df_align)
            coverage[chain] = ranges_covered
            b_factor_dictionary[chain] = np.array((chain_df[chain_df.b_factor != "-"].b_factor).astype(float))
        else:
            coverage[chain] = "NA"
            b_factor_dictionary[chain] = "NA"

    if warning_column == {}:
        warning_column = "NA"

    chains_string = ';'.join([str(elem) for elem in chain_names])
    output_array = np.array([chains_string, coverage, muts, mutation_list, b_factor_dictionary, warning_column, r_free], dtype=object)
    
    os.chdir(path)
    return output_array

def download_mmcif(structure_identifier):
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
    if os.path.isdir("structures") == False:
        os.mkdir("structures")
    
    #check if file already exists
    if os.path.isfile(f"structures/{structure_identifier}.cif") == False:
        
        if structure_identifier.startswith("AF-"):
            #download the alphafold file - we know the structure exist. If ther is a download issue, 
            #this is handled in the calling function. 
            os.system(f"wget https://alphafold.ebi.ac.uk/files/{structure_identifier}.cif -P structures")
            
        else:
            #download the pdb file - we know the structure exist. If ther is a download issue, 
            #this is handled in the calling function. 
            os.system(f"wget https://files.rcsb.org/download/{structure_identifier}.cif -P structures")
    
    return

def import_mmcif(structure_identifier, self_chains):
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
    
    if os.path.isfile(f"structures/{structure_identifier}.cif") == False:
        return "NA","NA","NA"
        
    mmcif_dict = MMCIF2Dict.MMCIF2Dict(f"structures/{structure_identifier}.cif")
    
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
            chain_df["residue_name"] = [seq1(i) for i in chain_df["residue"]]
            chain_df = chain_df[["residue_name", "seq_num", "b_factor"]]
            df.b_factor = df.b_factor.astype(float)
            chain_df = chain_df[chain_df.residue_name != '']
            chain_dfs.append(chain_df)

    try:
        r_free = mmcif_dict['_refine.ls_R_factor_R_free']
        try:
            r_free = float(r_free[0])
        except ValueError:
            r_free = "NA" #In cases where r-free has been entered as: ?/NA/X/Yes/No/none and more.
    except KeyError: # relevant for AlphaFold structures
        r_free = "NA"

    return chain_dfs, chain_names, r_free

def align_alphafold(alphafold_id, mutations):
    print("FUNCTION: align_alphafold(alphafold_id, mutations)")
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
    download_mmcif(alphafold_id)
    
    chain_dfs, chain_names, r_free = import_mmcif(alphafold_id, ["A"])
    
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

def align(combined_structure, path):
    print("FUNCTION: align(combined_structure, path)")
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
    
    uniprot_sequence, uniprot_numbering = get_uniprot_sequence(combined_structure.uniprot_id[0], int(combined_structure['uniprot_isoform'][0]))     

    for i in range(len(combined_structure)):
                           
        if combined_structure['structure_id'][i].startswith("AF-"):
            alignment_info = align_alphafold(combined_structure['structure_id'][i], combined_structure['mutations'][i])
            
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
                                           combined_structure.complex_protein_details[i])

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

def get_complex_information(pdb_id, uniprot):
    print(f"FUNCTION: get_complex_information(pdb_id, uniprot), {pdb_id}, {uniprot}")
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
        with open("log.txt", "a") as textfile:
            textfile.write(f"EXITING: Could not connect to PDBe database via API for {pdb_id}. \n")
        exit(1) 
    
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
        with open("log.txt", "a") as textfile:
            textfile.write(f"EXITING: Could not connect to PDBe database via API for molecules,{pdb_id}. \n")
        exit(1) 

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
            #n = f"[{molecule_dictionary[i]['molecule_type']}, {molecule_dictionary[i]['in_chains']}]"
            nucleotide_info.append(n)
        
        elif molecule_dictionary[i]['molecule_type'] != 'polypeptide(L)':
            if molecule_dictionary[i]['molecule_type'] != 'water': 
                l1 = molecule_dictionary[i]['molecule_name'][0]
                l2 = molecule_dictionary[i]['in_chains']
                l2 = ",".join(l2)
                l = f"{l1}, chain {l2}"
                #l = f"[{molecule_dictionary[i]['molecule_name'][0]}, {molecule_dictionary[i]['in_chains']}]"
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
            
    output_array = np.array([protein_complex_list, protein_info, self_chain, nucleotide_complex_list, 
                             nucleotide_info, ligand_complex_list, ligand_info], dtype=object)
    
    return output_array

def collect_complex_info(structural_df):
    print("FUNCTION: collect_complex_info(structural_df)")
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

def filter_all(structural_df, input_dataframe):
    print("FUNCTION: filter_all(structural_df, input_dataframe)")
    """
    This function cleans up sloppy coding from earlier in the pipeline
    which is needed for further investigation by removing structures that 
    does not satisfy the criteria.

    Parameters
    ----------
    structural_df : The pandas dataframe containing the original input file,
                   PDBs and subsequent additional columns. 

    Returns
    -------
    structural_df : The pandas dataframe containing the original input file,
                   PDBs that cover at least ONE mutation. 

    """
    final_dfs = []
    
    filtered_df = structural_df[structural_df['AA_in_PDB'].apply(lambda d: any(d.values()))]
    if not filtered_df.empty:
        filtered_df.index.name = 'structure_rank'
        final_dfs.append(filtered_df)
      
    if len(set(input_dataframe.cluster_id)) != 1:
        for cluster in list(input_dataframe.cluster_id):
            cluster_df = structural_df[structural_df.cluster_id == cluster]
            cluster_df = cluster_df.reset_index(drop=True)
            cluster_df.index.name = 'structure_rank'
            final_dfs.append(cluster_df)
    
    else:
        structural_df.index.name = 'structure_rank'
        final_dfs.append(structural_df)

    return final_dfs

def cleanup_all(structural_df,input_dataframe, path):
    print("FUNCTION: cleanup_all(structural_df,input_dataframe, path)")
    
    final_dfs = filter_all(structural_df, input_dataframe)
    for df in final_dfs:
        uniprot_id = df.iloc[0]['uniprot_id']
        #rearrange columns
        df = df[['hugo_name', 'uniprot_id', 'uniprot_isoform', 
                                       'mutations', 'cluster_id', 'structure_id', 
                                       'deposition_date', 'experimental_method', 
                                       'resolution', 'r_free','PDBREDOdb', 
                                       'PDBREDOdb_rfree', 'complex_protein',
                                       'complex_protein_details', 'complex_nucleotide',
                                       'complex_nucleotide_details', 'complex_ligand',
                                       'complex_ligand_details', 'chains', 'coverage', 
                                       'mutations_in_pdb','AA_in_PDB','b_factor', 'warnings']]
        
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
        
        df.index.name = 'structure_rank'
        df.to_json(f"{path}/results/{uniprot_id}/{uniprot_id}_all.json", orient="index")
        
        if type(df.mutations[0]) == float or type(df.mutations[0]) == np.float64:
            df.drop('mutations', inplace=True, axis=1)
            df.drop('AA_in_PDB', inplace=True, axis=1)
        
        if type(df.cluster_id[0]) == float or type(df.cluster_id[0]) == np.float64:
            df.drop('cluster_id', inplace=True, axis=1)    
        
        if 'mutations' in df.columns: #no cluster, but with mutations
            df.to_json(f"{path}/results/{uniprot_id}/{uniprot_id}_filtered.json", orient="index")
            
        if 'cluster_id' in df.columns:
            grouped = df.groupby(['cluster_id'])
            for cluster_id, cluster_data in grouped:
                cluster_data.to_json(f"{path}/results/{uniprot_id}/{uniprot_id}_cluster_{cluster_id}.json", orient="index")
                
    return 

