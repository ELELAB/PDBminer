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
from biopandas.pdb import PandasPdb
from Bio.SeqUtils import seq1

#============================================================================#
# Following the Documentation from PDBminer_run.py
# Step 1 - import is done in PDBminer_run
#============================================================================#

#============================================================================#
# 2)  find all structures related to the uniprot id 
#============================================================================#
# This step consist of four functions. 
#
#   get(pdbs)
#   get_structure_metadata(pdb_id)
#   get_structure_df(uniprot_id)
#   get_alphafold_basics(uniprot_id)
#   find_structure_list(input_dataframe)
#
#   The functions are called through find_structure_list
#   the aim is to take the input dataframe as input and output
#   a dataframe "found_structure_list" with all avialable PDB & newest AF 
#   structure including their metadata for further analysis. 

def get_pdbs(uniprot_id):
    """
    Function is taken from SLiMfast, documentation and comments there.
    Credit: Valentina Sora
    
    """
    uniprot_url = "https://www.uniprot.org/uniprot/{:s}.txt"
    response = requests.get(uniprot_url.format(uniprot_id))

    pdbs = []
    
    for line in response.text.split("\n"):
        
        if line.startswith("DR   PDB;"):
            
            db, pdb_id, exp, res, chain_res = \
                [item.strip(" ") for item \
                 in line.rstrip(".\n").split(";")]
            
            pdbs.append(pdb_id)

    return pdbs

def get_structure_metadata(pdb_id): 
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
    
    #Introduce the API to PDBe.
    mapping_url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{:s}"
    
    response = requests.get(mapping_url.format(pdb_id))
    response_text = json.loads(response.text)
    
    #exlucsion of PDBs without any metadata
    if response_text == {}:
        return {}
    
    #inclusion of PDBs with metadata
    else:
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
            response_experiment = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{pdb_id}")
            response_text_exp = json.loads(response_experiment.text)
            dictionary_exp = response_text_exp[pdb_id.lower()]
            
            dictionary_exp = dictionary_exp[0]
            
            resolution = dictionary_exp['resolution']
        
    return deposition_date, experimental_method, resolution

def get_alphafold_basics(uniprot_id):
    """
    Function that takes a uniprot id and retrieve data from the alphafold
    database, and prepare the basic information of that model in alignment 
    with the PDB data.

    Parameters
    ----------
    uniprot_id : A sting, e.g. 'P04637'

    Returns
    -------
    A tuple of information fitting as a line in the structure_df captured 
    in get_structure_df. 

    """
    
    result = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}")
    result = json.loads(result.text)
    
    if result == {}:
        return {}

    else:     
        result=result[0]
        deposition_date = result['modelCreatedDate'] 
        url = result['pdbUrl']
        local_filename = url.split('/')[-1] 
        Alphafold_ID = local_filename[:-4]
        
    experimental_method = "Predicted"
    resolution = "NA"
    rank = 0

    return Alphafold_ID, rank, uniprot_id, deposition_date, experimental_method, resolution


def get_structure_df(uniprot_id): 
    """
    This function takes a single uniprot ID and outputs a 
    dataframe containing a sorted list of PDB ids and their metadata. 

    Parameters
    ----------
    uniprot_id : A sting, e.g. 'P04637'

    Returns
    -------
    structure_df :  A pandas dataframe containing all the names and details
                    regarding the solved structures related to the 
                    uniprot id in terms of pdb files.  

    """
    #find all pdbs for a uniprot id
    pdb_list = get_pdbs(uniprot_id)
    
    #if there are no PDBs check for alphafold models
    if len(pdb_list) == 0:
        
        AF_model = get_alphafold_basics(uniprot_id)
        
        #if no PDB structures nor AF models exist,return and end the function.
        if AF_model == {}:
            return uniprot_id
        
        #if an alphagold model exist, use that as input for structure_df
        else:
            structure_df = pd.DataFrame({'uniprot_id': [AF_model[2]], 
                                         'deposition_date': [AF_model[3]], 
                                         'experimental_method': [AF_model[4]], 
                                         'resolution': [AF_model[5]]})
            structure_df.index = [AF_model[0]]           
            
            return structure_df
    
    #if there are PDBs these should be captured. 
    else: 
        #empty list for collection of PDBs 
        currated_pdb_list = []
        #create an empty df
        structure_df = pd.DataFrame()
    
        #make a for loop to populate the empty df with metadata for each pdb
        for pdb in pdb_list:
            
            #retrive all metadata        
            structure_metadata = get_structure_metadata(pdb)
          
            if type(structure_metadata) == tuple:
                
                data = [(structure_metadata[0]), (structure_metadata[1]), (structure_metadata[2])]
            
                structure_metadata_df = pd.DataFrame([data], columns = ['deposition_date', 
                                                                      'experimental_method', 
                                                                      'resolution'])
                
                structure_df = structure_df.append(structure_metadata_df)
                currated_pdb_list.append(pdb)
        
            else: 
                continue
                
        #Add pdb ID & uniprot ID as columns
        structure_df.index = currated_pdb_list
        structure_df.insert(0, "uniprot_id", [uniprot_id]*len(structure_df), True)
        
        #Split based on experimental method
        method_dfs = [x for _, x in structure_df.groupby(structure_df['experimental_method'])]
        
        #Sort within the method
        for method in range(len(method_dfs)):
            method_dfs[method] = method_dfs[method].sort_values(["resolution", "deposition_date"], ascending = [True, False])
            
            #add numerical values in each of the dataframes to ensure that the sorting is kept in place 
            if method_dfs[method].experimental_method[0] == 'X-RAY DIFFRACTION':
                method_dfs[method].insert(0, "Rank", sorted(range(1000,1000+len(method_dfs[method]))), True)
                
            elif method_dfs[method].experimental_method[0] == 'ELECTRON MICROSCOPY':
                method_dfs[method].insert(0, "Rank", sorted(range(2000,2000+len(method_dfs[method]))), True)
                
            elif "NMR" in method_dfs[method].experimental_method[0]:
                method_dfs[method].insert(0, "Rank", sorted(range(3000,3000+len(method_dfs[method]))), True)

        #concatenate all the dataframes
        structure_df = pd.concat(method_dfs)
        
        AF_model = get_alphafold_basics(uniprot_id)
        
        if AF_model != {}: 
            structure_df.loc[AF_model[0]] = list(AF_model[1:])
        
        #Sort based on rank & remove rank when sorted
        structure_df = structure_df.sort_values(by="Rank")
        structure_df = structure_df.drop(['Rank'], axis=1)
                                                            
        return structure_df 

def find_structure_list(input_dataframe):    
    """
    Takes the input file and the path where it is placed and outputs
    a directory with a csv file for each uniprot id input and a txt file 
    including all the uniprot ids that does not have any solved structures.
    
    parameters
    ------------
    input_dataframe             The input df, as described in the readme file.  
    
    Returns          
    --------------
    missing_ID.txt:             Containing the uniprot id string which there 
                                are no solved structurs.
    
    found_structure_list:       A pandas datafrane where each solved structure
                                and a number of describtors are detailed. 
    """

    df_collector = []
    
    #take all uniprot id's from the input file
    all_uniprot_ids = list(input_dataframe.uniprot)
    all_uniprot_ids = sorted(set(all_uniprot_ids), key=all_uniprot_ids.index)
    
    with open("missing_ID.txt", "w") as textfile:
        
        for row in range(len(all_uniprot_ids)):
            
            print(all_uniprot_ids[row])
            
            structure_info = get_structure_df(all_uniprot_ids[row]) 
        
            if type(structure_info) != str: 
                df_collector.append(structure_info)
            
            else:
                textfile.write(structure_info + "\n")
            
    if len(df_collector) > 0:
        found_structure_list = pd.concat(df_collector) 
        
    else:
        found_structure_list = []
        
    return found_structure_list


#============================================================================#
# 3)  combine the input 1) and the structures 2)
#============================================================================#
# This step consist of one function. 
#
#   combine_structure_dfs(found_structures, input_dataframe)
#
#   The aim is to take the input dataframe, and the dataframe with all 
#   the structures from the prior step and combine these to a single 
#   dataframe to continue working on.


def combine_structure_dfs(found_structures, input_dataframe):
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
        
        df = pd.DataFrame([])
        
        sub_df = found_structures[found_structures.uniprot_id == input_dataframe.uniprot[i]]
        df.insert(0, "hugo_name", [input_dataframe.hugo_name[i]]*len(sub_df), True)
        df.insert(1, "uniprot_id", [input_dataframe.uniprot[i]]*len(sub_df), True)
        df.insert(2, "uniprot_isoform", [input_dataframe['uniprot_isoform'][i]]*len(sub_df), True)
        df.insert(3, "mutations", [input_dataframe.mutations[i]]*len(sub_df), True)
        df.insert(4, "cluster_id", [input_dataframe.cluster_id[i]]*len(sub_df), True)
        df.insert(5, "structure_id", list(sub_df.index), True)
        df.insert(6, "deposition_date", list(sub_df.deposition_date), True)
        df.insert(7, "experimental_method", list(sub_df.experimental_method), True)
        df.insert(8, "resolution", list(sub_df.resolution), True)
        
        df_collector.append(df)
    
    final_df = pd.concat(df_collector)    
    
    final_df = final_df.reset_index(drop=True)
    
    return final_df

#============================================================================#
# 4)  For each sequence of structure (PDBid) an alignment to the fasta of 
#     the specified isoform. Here the differencies between the PDB and
#     fasta sequence is identified, the area of the sequence covered
#     by the PDB is annotated and the amino acids at the mutational sites
#     are found. NB! this does not account for quality of the structure. 
#============================================================================#
# This step consist of three functions. 
#
#   to_ranges(iterable)
#   align_uniprot_pdb(pdb_id, uniprot_id, isoform, mut_pos, path)
#   align_alphafold(alphafold_id, mutation_positions)
#   align(combined_structure, path)
#
#   The first functions is called by the second and the second and third 
#   function is called through the forth function.
#   The aim is to take the dataframe created in the prior step, and the 
#   relative path as input, and output an identical dataframe including 
#   information of structural coverage and amino acids in mutational sites.
#

def to_ranges(iterable):
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
        
def align_uniprot_pdb(pdb_id, uniprot_id, isoform, mut_pos, path):
    """
    ...

    Parameters
    ----------
    pdb_id : A particular string containg the foru letter code to a pdb id
    uniprot_id : A string containing the uniprot id
    isoform : an interger 
    mut_pos : an array of positions the user wish to cover with the structure.
    path : directory string.

    Returns
    -------
    output_array: An np.array containing: 
        1. chains_string e.g. 'A';'B', description of each related chain
        2. coverage_string, e.g. [(1,123);(4,74)] area of alignment per chain 
        3. mutational_column_string e.g. [E17K;T74E],[E17E;T74T]
        4. mutation_list_string, eg., [P78K];[], these are often mutations
            introduced in the experiment to keep the protein stable or to 
            investigate a particular phenomenon.

    """
    
    #Create a structure directory. 
    if os.path.exists("structure") == False:
        os.mkdir("structure")
    
    os.chdir("structure")
    
    ########################################
    
    #Get the uniprot sequence for alignment
    
    ########################################
    
    #setting the fasta path for the API (updated may 2022, novel Uniprot format)
    uniprot_Url="https://rest.uniprot.org/uniprotkb/"
    
    #Alignment to correct isoform:
    if isoform == 1:
    
        fasta_url=uniprot_Url+uniprot_id+".fasta"
        response = requests.post(fasta_url)
        sequence_data=''.join(response.text)
        Seq=StringIO(sequence_data)
        pSeq=list(SeqIO.parse(Seq,'fasta'))
        uniprot_sequence = str(pSeq[0].seq)
        uniprot_numbering = list(range(1,len(uniprot_sequence)+1,1)) 
        #gaining a list of aminoacids and a list of numberical values 
        #indicating the number of each position. 
    
    else:
    
        fasta_url=uniprot_Url+uniprot_id+"-"+str(isoform)+".fasta"
        response = requests.post(fasta_url)
        sequence_data=''.join(response.text)
        Seq=StringIO(sequence_data)
        pSeq=list(SeqIO.parse(Seq,'fasta'))
        uniprot_sequence = str(pSeq[0].seq)
        uniprot_numbering = list(range(1,len(uniprot_sequence)+1,1)) 
    
    ########################################
    
    #Prepare empty lists
    
    ########################################    
        
    # Empty lists for popultion by function
    pdb = PDBList()

    alignment_score = []
    coverage = []
    muts = []
    chains = []
    sequences = []
    alignment_score = []
    mutational_column = []
    mutation_list = []
    
    ########################################
    
    #Get the PDB sequence for the alignment
    
    ########################################
    #retrive the PDB structure, downloading the pdb files
    pdb.retrieve_pdb_file(pdb_id, file_format="mmCif")
    pdb_dir = pdb_id[1:3]
    pdb_dir = pdb_dir.lower()
    os.chdir(pdb_dir)
        
    # for loop to retrive sequence from sturcture
    for record in SeqIO.parse(pdb_id.lower()+".cif", "cif-seqres"):
        if record.seq != '' and len(set(str(record.seq))) != 1:
            #alignment. The alignemet parameters are set as: 
                
                #1) Local alignment: Aim to find the area of the uniprot
                #   sequence the PDB covers.
                #2) Match (identical amino acids) =  1 point
                #3) Non-match (not identical) = -10 points
                #4) Opening a gap: -10 points
                #5) keeping a gap open: -1
                
                #These options have been chosen to force the best fit locally
                #and highly discourage gaps, as they do not make sense in 
                #this structral space. 
            
            ########################################
    
            # Align
    
            ########################################    
            alignments = pairwise2.align.localms(uniprot_sequence, record.seq, 1, -10, -10, -1)  
            
            #sometimes a stucture where the wanted uniprot id is only included 
            #as a ligand or very small sequence is included in the structural list
            #these are removed if the overall alignment score is less than 10. 
            #10 is an arbitrary number chosen based on observations. 
            if alignments != [] and alignments[0][2] > 10:
                chains.append(record.annotations['chain'])
                alignment_score.append(alignments[0][2])
                uniprot_aligned = alignments[0][0]
                pdb_aligned = alignments[0][1]
                
                #trim first part of alignment strings if pdb contains
                #amino acids that are not in the uniprot at the begining.         
                while uniprot_aligned[0] == "-":
                    uniprot_aligned = uniprot_aligned[1:]
                    pdb_aligned = pdb_aligned[1:]
                
                #convert strings to list to get correct uniprot numbering
                list_pdb_align=[]
                list_pdb_align[:0]=pdb_aligned
                
                list_uni_align=[]
                list_uni_align[:0]=uniprot_aligned
                
                ########################################
    
                #Capture information from alignment
                
                ########################################
                            
                if len(list_pdb_align) == len(uniprot_numbering):
                    #create a dataframe with three columns and remove all where the PDB does not cover. 
                    df = pd.DataFrame({"uni": list_uni_align, "num": uniprot_numbering, "seq": list_pdb_align})
                    df = df[df.seq!="-"]
                    
                    muts = []
                    
                    #capture the amino acid in all the mutational sites of interest. 
                    for mutational_positon in list(set(mut_pos)):
                        if mutational_positon in list(df.num):
                            new_df = df[df.num == mutational_positon]
                            mutation = list(new_df.uni)[0]+str(list(new_df.num)[0])+list(new_df.seq)[0]
                            muts.append(str(mutation))
                        else:
                            muts.append("Mutation not in range")
                    
                    #capture the area the PDB covers accourding to the alignment. 
                    ranges_covered = list(to_ranges(list(df.num)))
                    coverage.append(ranges_covered)
                    mutational_column.append(muts)
                    
                    df = df.reset_index(drop=True)
                    
                    #capture all mutations in PDBfile compared to the specified 
                    #protein isoform. 
                    mutations_in_all = []
                    for i in range(len(df)):
                        if df["uni"][i] != df["seq"][i]:
                            mutations_in_all.append(f"{list(df.uni)[i]}{str(list(df.num)[i])}{list(df.seq)[i]}")
                    
                    mutation_list.append(mutations_in_all)
        
                else:
                    coverage.append("Mismatch in alignment")
                    mutational_column.append("Mutations not covered")
            
    chains_string = ';'.join([str(elem) for elem in chains])
    #alignment_score_string = ';'.join([str(elem) for elem in alignment_score])
    coverage_string = ';'.join([str(elem) for elem in coverage])
    mutational_column_string = ';'.join([str(elem) for elem in mutational_column])        
    mutation_list_string = ';'.join([str(elem) for elem in mutation_list])
                                                    
    #directory administration to return to previous directory and 
    #avoid creating a russian doll of directories
    os.chdir(path) 
    
    output_array = np.array([chains_string, coverage_string, mutational_column_string, mutation_list_string])
    
    return output_array

def align_alphafold(alphafold_id, mutation_positions):
    """
    This function takes an alphafold ID and the mutational positions. 
    The aim is to find the high quality areas of the protein, and set these in 
    relation to the mutations. 

    Parameters
    ----------
    alphafold_id : String of the ID
    mutation_positions : List of mutations

    Returns
    -------
    coverage : [(x, y)] per chain high quality areas (pDDLT > 70)
    AA_in_PDB : If the high quality portions of the AF structure covers 
    mutations. 

    """
    
    url = (f"'https://alphafold.ebi.ac.uk/files/{alphafold_id}.pdb'")
    os.system(f"wget {url}")
    
    ppdb = PandasPdb().read_pdb(f"{alphafold_id}.pdb")
    ATOMS = ppdb.df['ATOM']
    ATOMS = ATOMS[ATOMS.atom_name == "CA"]
    confidence_list = list(ATOMS['b_factor'])
    positions = np.array(range(1,len(ATOMS)+1))
    sequence = list(seq1(''.join(list(ATOMS.residue_name))))

    confidence_categories = []
    for i in range(len(confidence_list)):
        if confidence_list[i] > 70:
            confidence_categories.append("high")
        else:
            confidence_categories.append("low")
            
    #create an intermediate df containing the quality estimates
    df = pd.DataFrame({'position':positions,'sequence':sequence,'PDDLT':confidence_list,'category':confidence_categories})    
    confident_seq = np.array(df[df.category == "high"].position)
        
    #create coverage string (PDBminer output style)
    f = []
    f.append(confident_seq[0])
    for i in range(len(confident_seq)-1):
        if confident_seq[i]+1 != confident_seq[i+1]:
            f.append(confident_seq[i])
            f.append(confident_seq[i+1])
    f.append(confident_seq[-1])
    
    f = np.array(f)
    coverage = []
    for i in range(len(f[::2])):
        p = f[::2][i],f[1::2][i] 
        coverage.append(p)
        
    coverage = str(coverage)
    coverage = coverage.replace("), (", ");(" )
    
    AA_in_PDB = []
    
    for i in range(len(mutation_positions)):
        if mutation_positions[i] == "N/A":
            mutation = "N/A"
        else:
            if df.category[df.position == mutation_positions[i]].values == "high":
                mutation = f"{df.sequence[df.position == mutation_positions[i]].values[0]}{mutation_positions[i]}{df.sequence[df.position == mutation_positions[i]].values[0]}"
            else:
                mutation = 'Mutation not in range'
            
        AA_in_PDB.append(mutation)
    
    AA_in_PDB = ",".join(AA_in_PDB)
    #output coverage string
    return coverage, AA_in_PDB

def align(combined_structure, path):
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
    
    for i in range(len(combined_structure)):
        
        if type(combined_structure['mutations'][i]) != str:
            combined_structure['mutation_positions'] = "N/A"
        else:
            combined_structure['mutation_positions'] = combined_structure['mutations'].str.split(';').apply(lambda x: [int(y[1:-1]) for y in x])
        
        if combined_structure['structure_id'][i].startswith("AF"):
            alignment_info = align_alphafold(combined_structure['structure_id'][i], combined_structure['mutation_positions'][i])
            
            combined_structure.at[i, 'chains'] = "A"  
            combined_structure.at[i, 'coverage'] = alignment_info[0] 
            combined_structure.at[i, 'AA_in_PDB'] = alignment_info[1] 
            combined_structure.at[i, 'mutations_in_pdb'] = "[]" 
            
        else:    
                 
            alignment_info = align_uniprot_pdb(combined_structure.structure_id[i],
                                           combined_structure.uniprot_id[i], 
                                           int(combined_structure['uniprot_isoform'][i]),
                                           combined_structure['mutation_positions'][i],
                                           path)
                
            if alignment_info[0] != '':
                combined_structure.at[i, 'chains'] = alignment_info[0]  
                combined_structure.at[i, 'coverage'] = alignment_info[1] 
                combined_structure.at[i, 'AA_in_PDB'] = alignment_info[2] 
                combined_structure.at[i, 'mutations_in_pdb'] = alignment_info[3] 
        
            #drop missing values. 
            else:
                combined_structure = combined_structure.drop([i])
    
    combined_structure = combined_structure.reset_index(drop=True)  
                        
    return combined_structure


#============================================================================#
# 5)  The PDB files are then analyzed in terms of other present proteins, 
#     indicating a complex, ligands and other molecules present. 
#============================================================================#
# This step consist of one function. 
#
#   get_complex_information(pdb_id)
#   collect_complex_info(structural_df)
#
#   The first function is called through collect_complex_info.
#   The aim is to take the structural dataframe as input and output
#   a structural dataframe including complex and ligand columns. 
#

def get_complex_information(pdb_id):
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
                3) nucleotide_complex_list: binary 
                4) nuleotide_info: description of complex if any 
                5) ligand_complex_list: binart 
                6) ligand_info: description of complex if any. Include metal.

    """    
    #finding protein complexes and their related uniprot_ids 
    #this step also serves as a quality control of uniprot id's and chains.
    mapping_url_proteins = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{:s}"
        
    response = requests.get(mapping_url_proteins.format(pdb_id))
    response_text = json.loads(response.text)

    if response_text != {}:
    
        protein_segment_dictionary = response_text[pdb_id.lower()]
    
        if len(protein_segment_dictionary['UniProt']) > 1:
            protein_complex_list = ['protein complex']

            info = []
        
            for i in protein_segment_dictionary['UniProt']:
                info.append(f"{protein_segment_dictionary['UniProt'][i]['identifier']}, {i}, chain_{protein_segment_dictionary['UniProt'][i]['mappings'][0]['chain_id']}")
        
            info = ';'.join(info)
            protein_info = [info]

        else: 
            protein_complex_list = ["NA"]
            protein_info = ["NA"]    

    else:
        protein_complex_list = ["NA"]
        protein_info = ["NA"]
    
    #finding complexes with other ligands by identifying other molecules
    mapping_url_molecules = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{:s}"
    
    response = requests.get(mapping_url_molecules.format(pdb_id))
    response_text = json.loads(response.text)
    molecule_dictionary = response_text[pdb_id.lower()]

    for i in range(len(molecule_dictionary)):
        
        n_info = []
        
        if "nucleotide" in molecule_dictionary[i]['molecule_type']:
            nucleotide_complex_list = ['nucleotide complex']
            
            n_info.append(f"[{molecule_dictionary[i]['molecule_type']}, {molecule_dictionary[i]['in_chains']}]")

            
            nuleotide_info = n_info
        
        else: 
            nucleotide_complex_list = ["NA"]
            nuleotide_info = ["NA"]
        
        
    for i in range(len(molecule_dictionary)):
        
        l_info = []
        
        if molecule_dictionary[i]['molecule_type'] != 'polypeptide(L)':
            if molecule_dictionary[i]['molecule_type'] != 'water':
                
                ligand_complex_list = ["Other ligands"]
                
                l_info.append(f"[{molecule_dictionary[i]['molecule_name'][0]}, {molecule_dictionary[i]['in_chains']}]")
                
                ligand_info = l_info
        
        else: 
            ligand_complex_list = ["NA"]
            ligand_info = ["NA"]
            
    output_array = np.array([protein_complex_list, protein_info, nucleotide_complex_list, 
                             nuleotide_info, ligand_complex_list, ligand_info])
    
    return output_array


def collect_complex_info(structural_df):
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
    
    df = pd.DataFrame(columns=['structure_id','complex_protein',
                           'complex_protein_details', 'complex_nucleotide',
                           'complex_nucleotide_details','complex_ligand', 
                           'complex_ligand_details'])

    for pdb in range(len(structural_df)):
        
        if structural_df['structure_id'][pdb].startswith("AF"):
            list_of_values = [structural_df.structure_id[pdb],['NA'],
                              ['NA'],['NA'],
                              ['NA'],['NA'],
                              ['NA']]
            
        else:

            complex_info = get_complex_information(structural_df.structure_id[pdb])
    
            list_of_values = [structural_df.structure_id[pdb], complex_info[0], 
                              complex_info[1], complex_info[2], 
                              complex_info[3], complex_info[4], 
                              complex_info[5]]

        df.loc[pdb] = np.array(list_of_values, dtype="object") 
    
    df_combined = pd.merge(structural_df, df, how='inner', on = 'structure_id')
    
    return df_combined


#============================================================================#
# 6)  All these informations is reported in all_{uniprot_id}_structural_df.csv
#     This is done in PDBminer_run

#============================================================================#


#============================================================================#
# 7)  A cleanup of the path removing structures and if structures.
#     This is done in PDBminer_run
#============================================================================#        


#============================================================================#
# 8)  The structural_df if cleaned only keeping the PDB files that at 
#     least cover one of the specifed mutations. This is reported as
#     clean_{uniprot_id}_structural_df.csv. If there is no structures 
#     that cover the mutations, a txt file is reported indicating 
#     that an alphafold structure may be the next path. 
#============================================================================#        
# This step consist of two functions.
# 
#   cleanup_all(structural_df)
#   filter_all(structural_df)
#
#   The aim is to take the structural dataframe as input and output
#   he structural dataframe with simple ammendments. This step sould
#   become obsolete in time, when prior functions are improved. 

def cleanup_all(structural_df):

    structural_df = structural_df.drop(columns=['AA_in_PDB', 'mutation_positions'])
    structural_df.index.name = 'structure_rank'

    return structural_df

def filter_all(structural_df):
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
    
    #remove sturctures where no mutations are within range    
    for i in range(len(structural_df.AA_in_PDB)):
        chains = structural_df.AA_in_PDB[i].replace(";",",")
        chains = chains.replace("[", "")
        chains = chains.replace("]", "")
        chains = chains.replace("'", "")
        chains = chains.replace(" M", "M")
        chains = chains.split(',')
            
        if list(set(chains)) == ['Mutation not in range']:
                
            structural_df = structural_df.drop([i])
    
    structural_df = structural_df.reset_index(drop=True)        
    
    #Remove structures with only mismatch in alignment
    for i in range(len(structural_df.coverage)):
        if list(set(structural_df.coverage[i].split(";"))) == ['Mismatch in alignment']:
            structural_df = structural_df.drop([i])
    
    structural_df = structural_df.reset_index(drop=True)
    
    for i in range(len(structural_df.mutations_in_pdb)):
        if set(structural_df.mutations_in_pdb[i].split(";")) == {'[]'}:
            structural_df.iloc[i, structural_df.columns.get_loc('mutations_in_pdb')] = '[]'
            
    structural_df = structural_df.drop(columns=['AA_in_PDB', 'mutation_positions'])
    
    structural_df.index.name = 'structure_rank'

    return structural_df










