#SliM codes

##############
# Import packs
##############

import requests
from tempfile import NamedTemporaryFile
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import json
import logging as log

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
    # Mutated model
    ("_entity.pdbx_mutation",
     str,
     "model_mutated"),
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
                    #"any reason you're just not keeping the tuple?"¨
                    # reluctant to make changes in the SliM code. 
            else:
                # The field is not present
                data = None

            # Append data to the list
            metadata[col] = data

        # Return the dictionary of metadata
        return metadata
    

# Get the module logger
logger = log.getLogger(__name__)
    
def get_uniprot_segments(pdb_id, uniprot_id):
    """Get which portions (segments) of a protein, identified by its
    UniProt ID, are covered in a given PDB file, given its PDB ID. SLiMfast code.
    """

    # URL where the mapping between PDB entities in a PDB entry
    # and their corresponding UniProt data is stored on PDBe
    mapping_url = \
        "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{:s}"

    # Get the data from PDBe
    response = requests.get(mapping_url.format(pdb_id))
    
    # Convert to text and read the resulting dictionary through JSON
    response_text = json.loads(response.text)

    # Get the data about the mappings 
    data = \
        response_text[pdb_id.lower()]["UniProt"][uniprot_id]["mappings"]
    
    # Create an empty dictionary to store the mapping between
    # the PDB ID and the chains corresponding to the protein
    # identified by the UniProt ID
    mappings = {}

    # For each mapping found
    for m in data:

        # Get the chain ID
        pdb_chain = m["chain_id"]

        # Get the start and end of the chain as stored in the PDB
        pdb_start = m["start"]["author_residue_number"]
        pdb_end = m["end"]["author_residue_number"]

        # Get the start and end of the chain as they would be
        # in UniProt numbering
        uniprot_start = m["unp_start"]
        uniprot_end = m["unp_end"]

        if pdb_start == "null" or pdb_end == "null":
            warnstr = \
                f"Mapping for UniProt ID {uniprot_id} residues " \
                f"{unp_start}-{unp_end} not available in PDB " \
                f"{pdb_id} chain {pdb_chain}. Therefore, this " \
                f"segment will not be considered."
            logger.warning(warnstr)
            continue

        # Create a tuple of tuples to store those ranges
        ranges = ((pdb_start, pdb_end), (uniprot_start, uniprot_end))
    
        # If the chain ID was already encountered (i.e. the PDB
        # chain corresponds to multiple discontinuous protein
        # segments)
        if pdb_chain in mappings.keys():
            # Add the ranges to the mapping
            mappings[pdb_chain].append(ranges)

        # If it is the first time that the chain ID is encountered
        else:
            # Create a new list and add the mapping as first element
            mappings[pdb_chain] = [ranges]

    # Return the mappings
    return mappings