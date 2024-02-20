"""Combines PDB files"""

from Bio.PDB import PDBParser, PDBIO, Select

def combine_pdbs(out, *paths):
    """
    Concatenates PDB files into a single PDB file.

    Parameters
    ----------
    out: str
        Path to output PDB file
    *paths: str
        Paths to input PDB files
    """
    class MultiPDBSelector(Select):
        def __init__(self, structures):
            self.structures = structures  # A list of Structure objects to include

        def accept_model(self, model):
            # Accept only the first model of each structure
            return model.parent.id == 0

        def accept_chain(self, chain):
            # Accept all chains, could be filtered further if needed
            return True

        def accept_residue(self, residue):
            # Accept all residues, assuming they are correctly numbered and non-overlapping
            return True

    # Load the PDB files into Structure objects based on the paths provided
    parser = PDBParser(QUIET=True)
    structures = [parser.get_structure("PDB_{}".format(i), path) for i, path in enumerate(paths)]

    # Initialize PDBIO
    io = PDBIO()

    # For each structure, set it and then selectively save using the MultiPDBSelector
    combined_structure = structures[0]  # Start with the first structure
    for structure in structures[1:]:  # Add models from subsequent structures
        for model in structure:
            if not combined_structure.has_id(model.id):
                combined_structure.add(model.copy())

    io.set_structure(combined_structure)
    io.save(out, MultiPDBSelector(structures))

