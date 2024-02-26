"""Build missing residues and atoms using MODELLER

**Usage**::

    >>> from qp.structure import missing_loops
    >>> pdb = "1lnh"

    # Parse PDB file, filling in missing residues
    >>> residues = missing_loops.get_residues("path/to/PDB.pdb")

    # Write alignment file
    >>> missing_loops.write_alignment(
    ...     residues, 
    ...     pdb, 
    ...     "path/to/PDB.pdb", 
    ...     "path/to/ALI.ali"
    ... )
    >P1;1lnh (template from original PDB file)
    structureX:1lnh.pdb:FIRST:@ END::::::
    --------GHKIKGTVVLMRKNVLDVNSVTSV-------------TLDTLTAFLGRSVSLQLISAT...
    >P1;1lnh_fill (full sequence)
    sequence:::::::::
    MLGGLLHRGHKIKGTVVLMRKNVLDVNSVTSVGGIIGQGLDLVGSTLDTLTAFLGRSVSLQLISAT...

    # Run MODELLER with the given alignment file
    >>> missing_loops.build_model(
    ...     residues, 
    ...     pdb, 
    ...     "path/to/PDB.pdb", 
    ...     "path/to/ALI.ali", 
    ...     "path/to/OUT.pdb"
    ... )

Optimization level (``optimize`` argument in ``missing_loops.build_model``): 

* 0. No optimization. Missing coordinates filled in using MODELLER's topology library.
* 1. Optimize missing residues and residues with missing atoms only. (Default)
* 2. Optimize the entire structure. Hetero atoms are included but treated as rigid bodies. 
"""

import os
from modeller import log, Environ, Selection
from modeller.automodel import AutoModel
from modeller import alignment, model


log.none()


#: Amino acid 3 to 1 lookup table
def define_residues():
    """
    Access a look up table of standard and non-standard amino acids.

    Returns
    -------
    AA: dict
        Dictionary of standard and non-standard amino acids.
    """
    AA = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        "MSE": "M",
    }

    return AA


def get_residues(path, AA):
    """
    Extracts residues from a PDB file, filling in missing residues based on
    sequence number

    Parameters
    ----------
    path: str
        Path to PDB file

    Returns
    -------
    residues: list of list
        Residues separated by chain. Stored as a tuple of
        ``((sequence number, insertion code), one letter code, flag)``,
        where flag is R for completely absent, A for missing atoms, empty otherwise
    """
    with open(path, "r") as f:
        p = f.read().splitlines()

    missing = {}
    ind = {}
    seen = {}
    residues = []
    cur = None

    for line in p:
        if line.startswith("ENDMDL"):
            break

        elif line.startswith("REMARK 465   "):  # missing residues
            res = line[15:18]
            if res in AA and (line[13] == " " or line[13] == "1"): # want only the first model
                chain = line[19]
                rid = int(line[21:26])
                ic = line[26]
                missing.setdefault(chain, set()).add(((rid, ic), AA[res], "R"))

        elif line.startswith("REMARK 470   "):  # missing atoms
            res = line[15:18]
            if res in AA and (line[13] == " " or line[13] == "1"):
                chain = line[19]
                rid = int(line[20:24])
                ic = line[24]
                missing.setdefault(chain, set()).add(((rid, ic), AA[res], "A"))

        elif line.startswith("ATOM") or line.startswith("HETATM"):
            res = line[17:20]
            chain = line[21]
            rid = int(line[22:26])
            ic = line[26]

            if chain != cur:
                residues.append([])
                if chain in missing and chain not in ind:
                    missing[chain] = sorted(missing[chain])
                    ind[chain] = 0
                if chain not in seen:
                    seen[chain] = set()
                cur = chain

            if chain in missing:
                while (
                    ind[chain] < len(missing[chain])
                    and (rid, ic) >= missing[chain][ind[chain]][0]
                ):
                    residues[-1].append(missing[chain][ind[chain]])
                    ind[chain] += 1
                    seen[chain].add(residues[-1][-1][:2])

            if line.startswith("HETATM") and res != "MSE":
                resname = "w" if res == "HOH" else "."
            else:
                resname = AA.get(res, ".")
            if ((rid, ic), resname) not in seen[chain]:
                residues[-1].append(((rid, ic), resname, ""))
                seen[chain].add(residues[-1][-1][:2])

        elif line.startswith("TER"):
            chain = line[21]
            if chain in missing:
                while ind[chain] < len(missing[chain]):
                    residues[-1].append(missing[chain][ind[chain]])
                    ind[chain] += 1

    return residues

def clean_termini(residues):
    """
    Remove unresolve amino acids from the N- and C-termini.

    Most proteins will have long sequences of missing residues at the ends.
    Adding these with Modeller is unreliable and are removed with this function.

    Parameters
    ----------
    residues: list of lists of tuples
        Residues by chain [[((),)],[((),)]]
    
    Returns
    -------
    residues: list of lists of tuples
        Residues by chain. Stored as a tuple without unresolved termini

    """
    stripped_residues = []
    for chain in residues:
        # Strip 'R' tuples from the beginning of the chain
        while chain and chain[0][-1] == 'R':
            chain.pop(0)

        # Strip 'R' tuples from the end of the chain
        while chain and chain[-1][-1] == 'R':
            chain.pop()

        # Add the processed chain to the stripped protein
        if chain:  # Only add non-empty chains
            stripped_residues.append(chain)

    return stripped_residues
    

def write_alignment(residues, pdb, path, out):
    """
    Writes MODELLER alignment file for missing residues, according to
    https://salilab.org/modeller/10.4/manual/node501.html and
    https://salilab.org/modeller/wiki/Missing_residues

    Parameters
    ----------
    residues: list of list
        Residues separated by chain. Stored as a tuple of
        ``((sequence number, insertion code), one letter code, flag)``
    pdb: str
        PDB code
    path: str
        Path to PDB file
    out: str
        Path to output file
    """
    seq = "/".join(
        "".join(res[1] if res[2] != "R" else "-" for res in chain) for chain in residues
    )
    seq_fill = "/".join("".join(res[1] for res in chain) for chain in residues)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as f:
        f.write(f">P1;{pdb}\nstructureX:{path}:FIRST:@ END::::::\n{seq}*\n")
        f.write(f">P1;{pdb}_fill\nsequence:::::::::\n{seq_fill}*\n")


def transfer_numbering(e, ali, path, out):
    """
    Transfer residue numbers and chain ids from reference model to the built model.

    By default, Modeller will restart the number at 1 and rename the chains.
    This function, we transfer the numbering over from the template PDB.
    However, it can cause some issues with the numbering in some cases.
    We use fix_numbering() to fix those issues.

    See Also
    --------
    qp.structure.missing_loops.fix_numbering()

    """
    # Read the alignment for the transfer
    aln = alignment(e, file=ali)

    # Read the template and target models:
    mdl_reference = model(e, file=path)
    mdl_built = model(e, file=os.path.basename(out))

    # Transfer the residue and chain ids and write out the modified MODEL:
    mdl_built.res_num_from(mdl_reference, aln)
    file=os.path.basename(out)
    mdl_built.write(file)

    with open(file, 'r') as f:
        content = f.read()

    corrected_content = fix_numbering(content)
    with open(file, 'w') as f:
        f.write(corrected_content)


def fix_numbering(pdb_content):
    """
    Handles insertion codes from Modeller.

    Modeller resets the residue indices, which we fix with transfer_numbering().
    However, transfer_numbering() names residues with insertion codes.
    For example, 5,5A,5B,5C instead of 5,6,7,8.
    This function fixes is residue numbering behavior.
    Takes the generated PDB, fixes the residue numbers and overwrites the PDB.

    Parameters
    ----------
    pdb_content : str
        Contents of a PDB as a str with line breaks.

    Returns
    -------
    pdb_content :str
        Contents of a PDB as a str with line breaks.
    """

    lines = pdb_content.split('\n')
    new_lines = []

    prev_chain = None
    prev_residue_number = None
    seen_residues = set()

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            current_chain = line[21]
            residue_number_with_letter = line[22:27].strip()  # Residue number along with potential insertion code
            has_char = line[26] != ' '

            # If we encounter a new chain or a residue we haven't seen, reset/update values
            if current_chain != prev_chain:
                prev_chain = current_chain
                prev_residue_number = None
                seen_residues = set()
            if residue_number_with_letter not in seen_residues:
                if has_char:
                    prev_residue_number = prev_residue_number + 1 if prev_residue_number is not None else int(residue_number_with_letter[:-1])
                else:
                    prev_residue_number = int(residue_number_with_letter)
                seen_residues.add(residue_number_with_letter)

            # If there's an insertion code, replace the residue number and remove the insertion code
            if has_char:
                line = line[:22] + f"{prev_residue_number:4} " + line[27:]

            new_lines.append(line)
        else:
            new_lines.append(line)
    
    pdb_content = '\n'.join(new_lines)
    return pdb_content


def build_model(residues, pdb, path, ali, out, optimize=1):
    """
    Runs MODELLER for the given alignment file

    Parameters
    ----------
    residues: list of list
        Residues separated by chain. Stored as a tuple of
        ``((sequence number, insertion code), one letter code, flag)``
    pdb: str
        PDB code
    path: str
        Path to PDB file
    ali: str
        Path to alignment file
    out: str
        Path to output PDB file
    optimize: int
        Flag for level of optimization to use
        (0 - no optimization,
        1 - only missing residues or residues with missing atoms,
        2 - everything)
    """
    ali = os.path.abspath(ali)
    cwd = os.getcwd()
    dir = os.path.dirname(os.path.abspath(out))
    os.makedirs(dir, exist_ok=True)
    os.chdir(dir) # MODELLER only supports writing files to the current working directory

    class CustomModel(AutoModel):
        def get_model_filename(self, root_name, id1, id2, file_ext):
            return os.path.basename(out)

        def special_patches(self, aln): # renumber residues to avoid hybrid-36 notation with large models
            chain_ids = []
            for i in range(26):
                chain_ids.append(chr(ord("A") + i))
            for i in range(10):
                chain_ids.append(str(i))
            for i in range(26):
                chain_ids.append(chr(ord("a") + i))
            n = len(residues)
            self.rename_segments(chain_ids[:n], [1] * n)

    missing = []
    for i, c in enumerate(residues):
        chain = chr(ord("A") + i)
        ind = None
        for j, res in enumerate(c):
            if res[2] and ind is None:
                ind = j + 1
            elif not res[2] and ind is not None:
                missing.append((f"{ind}:{chain}", f"{j}:{chain}"))
                ind = None
        if ind is not None:
            missing.append((f"{ind}:{chain}", f"{len(c)}:{chain}"))

    if optimize == 1 and missing:
        CustomModel.select_atoms = lambda self: Selection(
            *[self.residue_range(x, y) for x, y in missing]
        )

    e = Environ()
    e.io.hetatm = True
    e.io.water = True
    a = CustomModel(e, alnfile=ali, knowns=pdb, sequence=f"{pdb}_fill")
    a.starting_model = 1
    a.ending_model = 1

    if not optimize or (optimize == 1 and not missing):
        a.make(exit_stage=2)
        os.rename(f"{pdb}_fill.ini", os.path.basename(out))
        transfer_numbering(e, ali, pdb, out)
    else:
        a.make()
        for ext in ["ini", "rsr", "sch", "D00000001", "V99990001"]:
            os.remove(f"{pdb}_fill.{ext}")
        transfer_numbering(e, ali, path, out)

    os.chdir(cwd)

