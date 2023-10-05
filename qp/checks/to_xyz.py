"""Converts PDB files to XYZ"""

from Bio.PDB import PDBParser

def format_element(element):
    return element[0].upper() + element[1:].lower()

def to_xyz(out, *paths):
    """
    Concatenates and writes PDB files to an XYZ file

    Parameters
    ----------
    out: str
        Path to output XYZ file
    *paths: str
        Paths to input PDB files
    """
    atoms = []
    parser = PDBParser(QUIET=True)
    for p in paths:
        structure = parser.get_structure("PDB", p)
        for atom in structure[0].get_atoms():
            atoms.append((atom.element, atom.get_coord()))

    with open(out, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for element, (x, y, z) in atoms:
            f.write(f"{format_element(element):>3}  {x:8.3f} {y:8.3f} {z:8.3f}\n")
