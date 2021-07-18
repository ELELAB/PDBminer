# -*- coding: utf-8 -*-
"""
This file contains the code required to answer
the question: 

Does a structure exist?

###################################################################

This file contains the following inport from slim, that will be imported
as a pack in future 
    
    get_pdbs(uniprot_id) 
    STRUCTURE_METADATA
    get_structure_metadata(pdb_id)
   
####################################################################
    
This file also contains: 
    
    get_structure_df(uniprot_id)        
    find_structure_list(input_dataframe, path)

"""

##############
# Import packs
##############

import pandas as pd
import requests
from tempfile import NamedTemporaryFile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os

##############

def get_pdbs(uniprot_id):
    """
    Code is taken from SLiMfast, documentation and comments there.
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

STRUCTURE_METADATA = [
    # experimental method used for structure determination
    ("_exptl.method",
     str, 
     "experimental_method"),
    # resolution
    ("_refine.ls_d_res_high",
     float,
     "resolution"),
    # deposition date
    ("_pdbx_database_status.recvd_initial_deposition_date",
     str,
     "deposition_date"),
    ]

def get_structure_metadata(pdb_id):
    """
    Code is taken from SLiMfast, documentation and comments there.
    """

    mmcif_url = "https://files.rcsb.org/view/{:s}.cif"
    metadata = {}

    with NamedTemporaryFile("w") as tmp:
        
        # Get the mmCIF file
        response = requests.get(mmcif_url.format(pdb_id))
        # Write the mmCIF file to the temporary file
        tmp.write(response.text)
         
        # Create a dictionary from the header data
        mmcif_dict = MMCIF2Dict(tmp.name)

        # For each field of interest
        for field, datatype, col in STRUCTURE_METADATA:
            
            # If the field is present in the mmCIF file
            if field in mmcif_dict.keys():

                # Retrieve data for that field
                data = mmcif_dict[field]
                
                # If the data type is not tuple
                if datatype is not tuple:
                    # If the piece of data is available
                    if data[0] != "?":
                        # Convert it to the appropriate data type
                        data = str(datatype(data[0]))
                    else:
                        data = None
                
                # If the data type is tuple   
                else:
                    # Join the pieces of data in a string,
                    # separated by a semicolon
                    data = \
                        ":".join(datatype(\
                            [item for item in data if item != "?"]))        
            else:
                # The field is not present
                data = None

            # Append data to the list
            metadata[col] = data

        # Return the dictionary of metadata
        return metadata


###########################################################################
    
def get_structure_df(uniprot_id):
    """ explain function
    """
    
    #find all pdbs for a uniprot id
    pdb_list = get_pdbs(uniprot_id)
    
    #create an empty df
    structure_df = pd.DataFrame([])
    
    #make a for loop to populate the empty df with metadata for each pdb
    for i in range(len(pdb_list)):
        
        #retrive all metadata        
        structure_metadata = get_structure_metadata(pdb_list[i])
        
        #input the values into a dataframe
        structure_metadata_df = pd.DataFrame(structure_metadata.values())
        #transpose dataframe for easy reading
        structure_metadata_df = structure_metadata_df.transpose()
        
        #name columns
        structure_metadata_df.columns = structure_metadata.keys()
        
        structure_df = structure_df.append(structure_metadata_df) 
        
    #name_rows
    structure_df.index = pdb_list
    
    #sort based on a defined hierarchy (quick and dirty version to be updated)
    structure_df.sort_values(["experimental_method", "resolution", "deposition_date"], 
                             ascending = [False, True, False], inplace = True)
    
    return structure_df

def find_structure_list(input_dataframe, path):
    """ explain function
    """
    
    dir = 'structure_lists'
    if not os.path.exists(dir):
        os.makedirs(dir)
    
    os.chdir("structure_lists")
    
    #take all uniprot id's from the input file
    all_uniprot_ids = list(input_dataframe.Uniprot)
    #print(all_uniprot_ids)
    
    #collection list for all the uniprot id's that does not
    #have any pdb structures associated to them
    missing_ID = []
    
    #loop over each uniprot ID    
    for i in range(len(all_uniprot_ids)):
        print(i)
        
        try:
            #create the dataframe
            df = get_structure_df(all_uniprot_ids[i])
            
            
            #print the dataframe as a csv
            df.to_csv(all_uniprot_ids[i]+".csv")
        
        except:
            #add to a list of uniprot ids with no pdbs
            missing_ID.append(all_uniprot_ids[i])
            
    #write into a file
    textfile = open("missing_IDs.txt", "w")
    for element in missing_ID:
        textfile.write(element + "\n")
        textfile.close()
                
    return
    