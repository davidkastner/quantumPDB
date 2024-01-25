"""Converts PDB files to XYZ"""

from Bio.PDB import PDBParser

def format_element(element):
    return element[0].upper() + element[1:].lower()

def to_xyz(out, paths):
    """
    Concatenates and writes PDB files to an XYZ file.

    Generates multiple xyz files from generated PDB files.
    Depending on how many spheres were computed, multiple xyz will be created.
    
    Notes
    -----
    For example, if three spheres were computed: 1.pdb, 2.pdb, and 3.pdb,
    then three different xyz's would be generated.
    They would contain 1) 0.pdb, 1.pdb, 2) 0.pdb, 1.pdb, 2.pdb,
    3) and 0.pdb, 1.pdb, 2.pdb, 3.pdb.

    Parameters
    ----------
    out: str
        Path to output XYZ file
    paths: list
        Paths to input PDB files
    """

    parser = PDBParser(QUIET=True)

    for i in range(2, len(paths) + 1):
        atoms = []

        for path in paths[:i]:
            structure = parser.get_structure("PDB", path)
            for atom in structure[0].get_atoms():
                atoms.append((atom.element, atom.get_coord()))

        # Write out the xyz file for each sphere combination
        xyz_filename = f"{out}_{i - 1}.xyz"
        with open(xyz_filename, "w") as xyz_file:
            xyz_file.write(f"{len(atoms)}\n\n")
            for element, (x, y, z) in atoms:
                formatted_element = format_element(element)
                xyz_file.write(f"{formatted_element:>3}  {x:8.3f} {y:8.3f} {z:8.3f}\n")

