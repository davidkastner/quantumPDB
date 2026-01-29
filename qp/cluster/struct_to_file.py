"""Manipulate standard molecular files"""

from Bio.PDB import PDBParser

def format_element(element):
    """Capitalize an element symbol (e.g., ``'FE'`` -> ``'Fe'``)."""
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
        if structure:
            for atom in structure[0].get_atoms():
                atoms.append((atom.element, atom.get_coord()))

    with open(out, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for element, (x, y, z) in atoms:
            f.write(f"{format_element(element):>3}  {x:8.3f} {y:8.3f} {z:8.3f}\n")

def combine_pdbs(out_path, metals, *input_paths, hetero_pdb=False):
    """Combine multiple sphere PDB files into a single cluster PDB.

    Concatenates the contents of multiple PDB files (typically the numbered
    sphere files ``0.pdb``, ``1.pdb``, etc.) while ensuring only one ``END``
    record appears at the end.

    Parameters
    ----------
    out_path : str
        Path to the output combined PDB file.
    metals : CenterResidue
        Center residue definition (currently unused but kept for API).
    *input_paths : str
        Paths to input PDB files to concatenate.
    hetero_pdb : bool, optional
        If True, include HETATM records (default False, currently unused).
    """
    with open(out_path, 'w') as outfile:
        for path in input_paths:
            with open(path) as infile:
                for line in infile:
                    if line.startswith("END"):
                        continue
                    outfile.write(line)
        # Write the final "END" line
        outfile.write("END\n")
