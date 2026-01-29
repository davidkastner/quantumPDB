"""Parse Protoss output logs and detect residue type changes."""

from Bio.PDB import PDBParser
from qp.structure.missing import get_chain_order

def parse_log(log_path, pdb_path, AA):
    """
    Parse the Protoss log file to identify residues with clashes and their chain indices.

    Parameters
    ----------
    log_path: str
        Path to the Protoss log file.
    pdb_path: str
        Path to the PDB file.
    AA: dict
        Dictionary mapping three-letter codes to one-letter codes.

    Returns
    -------
    residues_with_clashes: set of tuples
        Set of residues with clashes in the format compatible with `residues`.
        Each tuple contains ((res_id, ' '), one_letter_code, 'R') and the chain index.
    """
    chain_order = get_chain_order(pdb_path)
    residues_with_clashes = set()

    with open(log_path, "r") as f:
        for line in f:
            if "ATOM" in line:
                # Remove everything before "ATOM"
                atom_line = line[line.index("ATOM"):]

                # Extract information based on PDB format columns
                res_name = atom_line[17:20].strip()
                chain = atom_line[21].strip()
                res_id = int(atom_line[22:26].strip())
                one_letter_code = AA.get(res_name, 'X')  # Use 'X' for unknown residues
                chain_index = chain_order.get(chain, -1)
                residues_with_clashes.add((res_id, one_letter_code, chain_index, chain, res_name))

    return residues_with_clashes



def record_type(modeller_pdb_path, protoss_pdb_path):
    '''
    Checks if protoss changed HETATM or ATOM classifications.

    Parameters
    ----------
    modeller_pdb_path: string
        The path to the original modeller PDB.
    protoss_pdb_path: string
        The path to the PDB processed by protoss.

    Returns
    -------
    changed_residues: list
        List of residues that have been changed from ATOM to HETATM or vice versa.

    Note
    ----
    changed_residues = [
    ('A', (' ', 123, ' '), 'ATOM', 'HETATM'),
    ('B', (' ', 45, ' '), 'HETATM', 'ATOM'),
    ('C', ('H', 67, ' '), 'HETATM', 'ATOM')]
    '''
    
    # Initialize the parser
    parser = PDBParser(QUIET=True)
    
    # Parse the PDB files
    modeller_structure = parser.get_structure('modeller', modeller_pdb_path)
    protoss_structure = parser.get_structure('protoss', protoss_pdb_path)
    
    changed_residues = []
    
    # Iterate through all chains and residues in the modeller structure
    for modeller_chain in modeller_structure.get_chains():
        chain_id = modeller_chain.get_id()
        
        if chain_id not in protoss_structure[0]:
            print(f"> WARNING: Chain {chain_id} not found in protoss structure.")
            continue
        
        protoss_chain = protoss_structure[0][chain_id]  # Assume the chains are the same
        
        for modeller_residue in modeller_chain.get_residues():
            res_num = modeller_residue.get_id()[1]
            resname = modeller_residue.get_resname()
            
            # Find the corresponding residue in the protoss chain by residue number and name
            protoss_residue = None
            for res in protoss_chain.get_residues():
                if res.get_id()[1] == res_num and res.get_resname() == resname:
                    protoss_residue = res
                    break
            
            if protoss_residue is None:
                print(f"> WARNING: Residue {modeller_residue.get_id()} in chain {chain_id} not found in protoss structure.")
                continue
            
            # Compare the residue types (HETATM vs ATOM)
            modeller_type = 'HETATM' if modeller_residue.id[0].startswith('H_') else 'ATOM'
            protoss_type = 'HETATM' if protoss_residue.id[0].startswith('H_') else 'ATOM'
            
            if modeller_type != protoss_type:
                changed_residues.append((chain_id, resname, res_num, modeller_type, protoss_type))
    
    # Alert the user if any residues have changed
    for chain_id, resname, res_num, old_type, new_type in changed_residues:
        print(f"> WARNING: {resname} {res_num} in chain {chain_id} was changed from {old_type} to {new_type}.")

    return changed_residues