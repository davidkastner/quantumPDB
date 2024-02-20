"""Manipulate standard molecular files"""

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

def combine_pdbs(out_path, metals, *input_paths):
    """
    Combines multiple PDB files into a single PDB file, excluding atoms of type "H1",
    while maintaining the order of residues and atoms as they appear in the input files,
    and ensuring only the final "END" line is kept.

    Parameters
    ----------
    out_path: str
        Path to the output PDB file.
    input_paths: str
        Paths to the input PDB files.
    """
    with open(out_path, 'w') as outfile:
        for path in input_paths:
            with open(path) as infile:
                for line in infile:
                    if line.startswith("HETATM") and all(metal not in line for metal in metals):
                        continue
                    if line.startswith("END"):
                        continue
                    outfile.write(line)
        # Write the final "END" line
        outfile.write("END\n")
