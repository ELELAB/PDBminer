#Find structure Lists

##############
# Import packs
##############

import pandas as pd

#############

from SLiM_codes import get_pdbs, get_structure_metadata

#############

def get_structure_df(uniprot_id, cluster_id):
    """ This function takes a single uniprot ID and outputs a 
    dataframe containing a sorted list of PDB ids and their metadata 
    associated to a singular Uniprot ID. 
    
    Parameters
    -------------
    uniprot_id      A sting, e.g. 'P04637'
    
    Returns
    -------------
    structure_df    A pandas dataframe containing all the names and details
                    regarding the solved structures related to the 
                    uniprot id in terms of pdb files.  
    
    """
    #find all pdbs for a uniprot id
    pdb_list = get_pdbs(uniprot_id)
    
    #empty list for collection of PDBs passing the first quality checkpoint. 
    currated_pdb_list = []
    
    if len(pdb_list) == 0:
        return uniprot_id
    
    else: 
        #create an empty df
        structure_df = pd.DataFrame()
    
        #make a for loop to populate the empty df with metadata for each pdb
        for pdb in pdb_list:
            
            #retrive all metadata        
            structure_metadata = get_structure_metadata(pdb)
    
            #Removing empty files such as PDBid 1OVC that where captured.
            if type(structure_metadata) == dict:
                #input the values into a dataframe
                structure_metadata_df = pd.DataFrame.from_dict(dict([(x, [k[x] for k in [structure_metadata]]) for x in structure_metadata]))

                #append values
                structure_df = structure_df.append(structure_metadata_df)
                currated_pdb_list.append(pdb)
        
            else: 
                continue
                
        #Add pdb ID & uniprot ID as columns
        structure_df.index = currated_pdb_list
        structure_df.insert(0, "Uniprot_ID", [uniprot_id]*len(structure_df), True)
        structure_df.insert(1, "ClusterID", [cluster_id]*len(structure_df), True)
    
        #Annotated structures as mutated or not-mutated
        structure_df["model_mutated"] = structure_df["model_mutated"].apply(lambda x: 'Not-Mutated' if x is None else 'Mutated')
        
        #Experimental methods
        #methods = structure_df['experimental_method'].unique()
    
        #Split based on experimental method
        method_dfs = [x for _, x in structure_df.groupby(structure_df['experimental_method'])]
        
        #Sort within the method
        for method in range(len(method_dfs)):
            method_dfs[method] = method_dfs[method].sort_values(["resolution", "deposition_date"], ascending = [True, False])
            
            #add numerical values in each of the dataframes to ensure that the sorting is kept in place 
            if method_dfs[method].experimental_method[0] == 'X-RAY DIFFRACTION':
                method_dfs[method].insert(0, "Rank", sorted(range(0,len(method_dfs[method]))), True)
                
            elif method_dfs[method].experimental_method[0] == 'ELECTRON MICROSCOPY':
                method_dfs[method].insert(0, "Rank", sorted(range(1000,1000+len(method_dfs[method]))), True)
                
            elif method_dfs[method].experimental_method[0] == 'SOLUTION NMR':
                method_dfs[method].insert(0, "Rank", sorted(range(2000,2000+len(method_dfs[method]))), True)
                                  
        #concatenate all the dataframes
        structure_df = pd.concat(method_dfs)
        
        #Sort based on rank & remove rank when sorted
        structure_df = structure_df.sort_values(by="Rank")
        structure_df = structure_df.drop(['Rank'], axis=1)
                                                            
        return structure_df          
                                          
                                            
def find_structure_list(input_dataframe):
    """Takes the input file and the path where it is placed and outputs
    a directory with a csv file for each Uniprot ID input and a txt file 
    including all the UNIprot IDs that does not have any solved structures.
    
    parameters
    ------------
    input_dataframe  The input df, as described in the readme file  
    
    Returns          
    --------------
    txt file         Containing a list of stings with all the uniprot IDs to
                     which there are no solved structurs.
    
    pandas df        A pd where each solved stuture
                     and a number of describtors are detailed. 
    """

    df_collector = []
    
    #take all uniprot id's from the input file
    all_uniprot_ids = list(input_dataframe.Uniprot)
    
    with open("missing_ID.txt", "w") as textfile:
        
        for row in range(len(all_uniprot_ids)):
            
            print(all_uniprot_ids[row])
            
            structure_info = get_structure_df(all_uniprot_ids[row], input_dataframe.ClusterID[row])
        
            if type(structure_info) != str: 
                df_collector.append(structure_info)
            
            else:
                textfile.write(structure_info + "\n")
            
    found_structure_list = pd.concat(df_collector)    
            
    return found_structure_list
