from Bio.PDB import PDBParser
def search_backbone_atoms(pdb_file):
    """
    Search for backbone atoms (C, CA, N, O) and return atom name, 
    residue number, and atom ID.
    
    Parameters
    ----------
    pdb_file : str
        Path to the PDB file
    
    Returns
    -------
    list
        List of tuples: (atom_name, residue_num, atom_id)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)
    target_atoms = {'C', 'CA', 'N', 'O'}
    backbone_data = []
    atom_index=0
    
    for atom in structure[0].get_atoms():
        atom_index += 1
        atom_name = atom.get_name().strip()
        if atom_name in target_atoms:
            
            backbone_data.append(atom_index)
    
    return backbone_data
