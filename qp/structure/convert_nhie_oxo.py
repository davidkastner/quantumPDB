"""Convert NHIE active sites to reactive OXO/succinate state.

This module provides specialized functions for converting non-heme iron
enzyme (NHIE) active sites from the resting state (with alpha-ketoglutarate
and water/NO ligands) to the reactive ferryl-oxo state (with succinate
and OXO ligands). This is useful for modeling the reactive intermediate
in NHIE catalysis.

Note: This is a specialized module for NHIE enzyme studies and is not
part of the standard QuantumPDB workflow.
"""

import numpy as np
import shutil
from Bio.PDB import PDBParser, PDBIO


def read_pdb(pdb_path):
    """Parse a PDB file and extract all ATOM/HETATM records.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.

    Returns
    -------
    list of dict
        Each dict contains PDB atom fields: record, serial, name, altLoc,
        resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor,
        element, charge.
    """
    atoms = []
    with open(pdb_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = {
                    'record': line[0:6].strip(),
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'altLoc': line[16:17].strip(),
                    'resName': line[17:20].strip(),
                    'chainID': line[21:22].strip(),
                    'resSeq': int(line[22:26]),
                    'iCode': line[26:27].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': float(line[54:60]),
                    'tempFactor': float(line[60:66]),
                    'element': line[76:78].strip(),
                    'charge': line[78:80].strip()
                }
                atoms.append(atom)
    return atoms

def write_pdb(atoms, pdb_path):
    """Write a list of atom dictionaries to a PDB file.

    Atoms with ``record == 'DELETE'`` are skipped.

    Parameters
    ----------
    atoms : list of dict
        Atom records as returned by :func:`read_pdb`.
    pdb_path : str
        Path to the output PDB file.
    """
    with open(pdb_path, 'w') as file:
        for atom in atoms:
            if atom['record'] != 'DELETE':
                file.write(f"{atom['record']:<6}{atom['serial']:>5} {atom['name']:<4}{atom['altLoc']:<1}{atom['resName']:<3} {atom['chainID']:<1}{atom['resSeq']:>4}{atom['iCode']:<1}   {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}{atom['occupancy']:>6.2f}{atom['tempFactor']:>6.2f}          {atom['element']:>2}{atom['charge']:>2}\n")

def find_and_convert_ligands_to_oxo(atoms, iron_names, distance_cutoff=2.8):
    """Convert water or NO ligands near iron to OXO residues.

    Searches for water (HOH, WAT) or nitric oxide (HOA, NO) ligands within
    the distance cutoff of any iron atom and converts them to OXO (oxo)
    ligands with the Fe-O bond set to 1.65 A.

    Parameters
    ----------
    atoms : list of dict
        Atom records from :func:`read_pdb`.
    iron_names : list of str
        Iron atom names to search for (e.g., ``['FE', 'FE2']``).
    distance_cutoff : float, optional
        Maximum Fe-ligand distance in angstroms (default 2.8).

    Returns
    -------
    set
        Serial numbers of iron atoms that received an OXO ligand.
    """
    ligands_to_convert = {'HOH': 'O', 'WAT': 'O', 'HOA': 'N', 'NO': 'N'}
    oxo_placed = set()

    for iron in atoms:
        if iron['name'] in iron_names:
            iron_coord = np.array([iron['x'], iron['y'], iron['z']])
            for ligand in atoms:
                if ligand['resName'] in ligands_to_convert and ligand['name'] == ligands_to_convert[ligand['resName']]:
                    ligand_coord = np.array([ligand['x'], ligand['y'], ligand['z']])
                    distance = np.linalg.norm(iron_coord - ligand_coord)
                    if distance <= distance_cutoff:
                        ligand['resName'] = 'OXO'
                        ligand['name'] = 'O'
                        ligand['element'] = 'O'
                        ligand['x'], ligand['y'], ligand['z'] = iron_coord + 1.65 * (ligand_coord - iron_coord) / distance
                        for other_atom in atoms:
                            if other_atom['resSeq'] == ligand['resSeq'] and other_atom['chainID'] == ligand['chainID'] and other_atom['serial'] != ligand['serial']:
                                other_atom['record'] = 'DELETE'
                        oxo_placed.add(iron['serial'])
                        print(f'   > Converted {ligand["serial"]} {ligand["resName"]}/{ligand["resName"]} to OXO at distance {round(distance, 2)}')
                        break
    return oxo_placed

def find_histidines(atoms, iron):
    """Find histidine NE2 atoms coordinating an iron center.

    Parameters
    ----------
    atoms : list of dict
        Atom records from :func:`read_pdb`.
    iron : dict
        Iron atom record.

    Returns
    -------
    list of dict
        Histidine NE2 atom records within 3.0 A of the iron.
    """
    iron_coord = np.array([iron['x'], iron['y'], iron['z']])
    histidines = []
    for residue in atoms:
        if residue['resName'] == 'HIS' and residue['name'] == 'NE2':
            ne2_coord = np.array([residue['x'], residue['y'], residue['z']])
            distance = np.linalg.norm(iron_coord - ne2_coord)
            coordination_cutoff_for_his = 3.0
            if distance <= coordination_cutoff_for_his:
                histidines.append(residue)
    print(f"   > Found {len(histidines)} His coordinating iron {iron['serial']}")
    return histidines

def calculate_dihedral(p1, p2, p3, p4):
    """Calculate the dihedral angle defined by four points.

    Parameters
    ----------
    p1, p2, p3, p4 : numpy.ndarray
        3D coordinates of the four points.

    Returns
    -------
    float
        Dihedral angle in degrees.
    """
    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.degrees(np.arctan2(y, x))

def place_oxo_manually(atoms, iron):
    """Place an OXO ligand opposite the coordinating histidines.

    For NHIE enzymes with two facial histidines, this function places
    the OXO ligand at the position opposite to one of the histidines,
    selecting the position with a dihedral angle closest to 90 degrees
    relative to the AKG C2-O5-Fe axis and avoiding steric clashes.

    Parameters
    ----------
    atoms : list of dict
        Atom records from :func:`read_pdb` (modified in place).
    iron : dict
        Iron atom record.

    Returns
    -------
    bool
        True if OXO was successfully placed, False otherwise.
    """
    histidines = find_histidines(atoms, iron)
    if len(histidines) != 2:
        print(f"> Error: Expected two His coordinating iron {iron['serial']}")
        return False

    iron_coord = np.array([iron['x'], iron['y'], iron['z']])
    potential_oxos = []

    for ne2 in histidines:
        ne2_coord = np.array([ne2['x'], ne2['y'], ne2['z']])
        oxo_coord = iron_coord + 1.65 * (iron_coord - ne2_coord) / np.linalg.norm(iron_coord - ne2_coord)
        print(f"> Checking potential oxo site opposite His {ne2['serial']}")

        clash_threshold = 1.2
        clash = False
        for atom in atoms:
            if atom['serial'] != iron['serial'] and not (atom['resName'] == 'AKG' and atom['name'] in ['C1', 'O1', 'O2']):  # Exclude specific AKG atoms from clash detection
                distance = np.linalg.norm(np.array([atom['x'], atom['y'], atom['z']]) - oxo_coord)
                if distance < clash_threshold:
                    clash = True
                    print(f"> Clash found with atom {atom['serial']} ({atom['name']} {atom['resName']}) at distance {round(distance, 2)}")
                    break
        if not clash:
            potential_oxos.append((oxo_coord, ne2['serial']))
            print(f"> Potential oxo site opposite His {ne2['serial']}.")

    # Find AKG coordinates for dihedral calculation
    try:
        akg_c2 = next(atom for atom in atoms if atom['resName'] == 'AKG' and atom['name'] == 'C2')
        akg_o5 = next(atom for atom in atoms if atom['resName'] == 'AKG' and atom['name'] == 'O5')
    except StopIteration:
        print("> Error: Could not find AKG atoms (C2 and O5) for dihedral calculation.")
        return False
    
    c2_coord = np.array([akg_c2['x'], akg_c2['y'], akg_c2['z']])
    o5_coord = np.array([akg_o5['x'], akg_o5['y'], akg_o5['z']])

    best_site = None
    best_dihedral = None

    # Calculate dihedral for all potential sites, including clashing ones
    for ne2 in histidines:
        ne2_coord = np.array([ne2['x'], ne2['y'], ne2['z']])
        oxo_coord = iron_coord + 1.65 * (iron_coord - ne2_coord) / np.linalg.norm(iron_coord - ne2_coord)
        dihedral = abs(calculate_dihedral(c2_coord, o5_coord, iron_coord, oxo_coord))
        print(f"> Dihedral angle for potential OXO site opposite His {ne2['serial']}: {round(dihedral, 2)}")

    # Select the best site among the non-clashing ones
    for oxo_coord, _ in potential_oxos:
        dihedral = abs(calculate_dihedral(c2_coord, o5_coord, iron_coord, oxo_coord))
        print(f"> Dihedral angle for non-clashing site: {round(dihedral, 2)}")
        if best_dihedral is None or abs(dihedral - 90) < abs(best_dihedral - 90):
            best_site = oxo_coord
            best_dihedral = dihedral

    if best_site is None:
        print("> Error: Could not determine the best site for OXO")
        return False

    oxo = {
        'record': 'HETATM',
        'serial': max(atom['serial'] for atom in atoms) + 1,
        'name': 'O',
        'altLoc': '',
        'resName': 'OXO',
        'chainID': iron['chainID'],
        'resSeq': iron['resSeq'],
        'iCode': '',
        'x': best_site[0],
        'y': best_site[1],
        'z': best_site[2],
        'occupancy': 1.00,
        'tempFactor': 0.00,
        'element': 'O',
        'charge': ''
    }
    atoms.append(oxo)
    rounded_best_site = [round(num, 2) for num in best_site]
    print(f'> OXO placed at {rounded_best_site} with dihedral {round(best_dihedral, 2)}\n')
    return True

def convert_akg_to_sin(atoms):
    """Convert alpha-ketoglutarate (AKG) to succinate (SIN).

    Removes the C1 carbonyl oxygens and converts C1 to an oxo group,
    simulating the decarboxylation that occurs during NHIE catalysis.
    Atoms are modified in place; deleted atoms have ``record = 'DELETE'``.

    Parameters
    ----------
    atoms : list of dict
        Atom records from :func:`read_pdb` (modified in place).
    """
    for atom in atoms:
        if atom['resName'] == 'AKG':
            if atom['name'] == 'O1' or atom['name'] == 'O2':
                atom['record'] = 'DELETE'
            elif atom['name'] == 'C1':
                atom['name'] = 'O1'
                atom['element'] = 'O'
                c2_atom = next(a for a in atoms if a['name'] == 'C2' and a['resSeq'] == atom['resSeq'] and a['chainID'] == atom['chainID'])
                c2_coord = np.array([c2_atom['x'], c2_atom['y'], c2_atom['z']])
                atom_coord = np.array([atom['x'], atom['y'], atom['z']])
                bond_vector = atom_coord - c2_coord
                bond_length = np.linalg.norm(bond_vector)
                new_bond_vector = bond_vector * (1.2 / bond_length)
                atom['x'], atom['y'], atom['z'] = c2_coord + new_bond_vector
                atom['resName'] = 'SIN'
            else:
                atom['resName'] = 'SIN'

def add_oxo_sin(pdb_path):
    """Convert an NHIE structure from resting to reactive state.

    Main orchestration function that:
    1. Converts water/NO ligands near iron to OXO
    2. Places OXO manually if no suitable ligand found
    3. Converts AKG to succinate (SIN)

    The original file is backed up with ``_akg`` suffix.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file (modified in place).
    """
    # Start by making a copy of the file
    shutil.copy(pdb_path, f'{pdb_path[:-4]}_akg.pdb')
    atoms = read_pdb(pdb_path)
    
    iron_names = ['FE', 'FE2', 'FE3']
    
    oxo_placed = find_and_convert_ligands_to_oxo(atoms, iron_names)
    
    for iron in atoms:
        if iron['name'] in iron_names and iron['serial'] not in oxo_placed:
            oxo_placed_manually = place_oxo_manually(atoms, iron)
            if oxo_placed_manually:
                oxo_placed.add(iron['serial'])

    if not oxo_placed:
        print("> Error: Could not place OXO")

    convert_akg_to_sin(atoms)
    atoms = [atom for atom in atoms if atom['record'] != 'DELETE']
    write_pdb(atoms, pdb_path)

def remove_oxo_hydrogens(protoss_pdb):
    """Remove hydrogen atoms from OXO residues.

    Protoss may add hydrogens to the OXO oxygen; this function removes
    them to maintain the correct oxo (O2-) charge state.

    Parameters
    ----------
    protoss_pdb : str
        Path to the PDB file (modified in place).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', protoss_pdb)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == 'OXO':
                    to_remove = []
                    for atom in residue:
                        if atom.element == 'H':
                            to_remove.append(atom)
                    for atom in to_remove:
                        residue.detach_child(atom.id)

    io = PDBIO()
    io.set_structure(structure)
    io.save(protoss_pdb)

def update_oxo_sdf(pdb_path, sdf_path):
    """Append OXO entries to the ligand SDF file.

    Creates SDF entries for OXO residues with -2 formal charge so that
    the charge embedding correctly accounts for the oxo ligand.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file containing OXO residues.
    sdf_path : str
        Path to the SDF file to append to.
    """
    atoms = read_pdb(pdb_path)
    oxo_atoms = [atom for atom in atoms if atom['resName'] == 'OXO']

    with open(sdf_path, 'a') as sdf_file:
        for oxo in oxo_atoms:
            sdf_file.write(f"{oxo['resName']}_{oxo['chainID']}_{oxo['resSeq']}\n")
            sdf_file.write("Generated by script\n")
            sdf_file.write(" QuantumPDB\n")
            sdf_file.write(" 1 0  0  0  0  0            999 V2000\n")
            sdf_file.write(f"   {oxo['x']:>10.4f}{oxo['y']:>10.4f}{oxo['z']:>10.4f} O   0  0  0  0  0  0  0  0  0  0  0  0\n")
            sdf_file.write(f"M  CHG  1   1  -2\n")
            sdf_file.write("M  END\n$$$$\n")

def update_sin_sdf(akg_sdf_path, sin_sdf_path):
    """Replace AKG entries with SIN entries in the ligand SDF file.

    Removes all AKG (alpha-ketoglutarate) entries from the original SDF
    and appends the SIN (succinate) entries from the Protoss-processed
    temporary file.

    Parameters
    ----------
    akg_sdf_path : str
        Path to the original SDF file (modified in place).
    sin_sdf_path : str
        Path to the temporary SDF file with SIN entries from Protoss.
    """
    # Read the old SDF file and remove AKG entries
    with open(akg_sdf_path, 'r') as old_sdf_file:
        old_sdf_lines = old_sdf_file.readlines()

    filtered_old_sdf_lines = []
    i = 0
    while i < len(old_sdf_lines):
        if old_sdf_lines[i].startswith("AKG"):
            while i < len(old_sdf_lines) and not old_sdf_lines[i].strip() == "$$$$":
                i += 1
            i += 1  # Skip the "$$$$"
        else:
            filtered_old_sdf_lines.append(old_sdf_lines[i])
            i += 1

    # Read the temporary SIN SDF file
    with open(sin_sdf_path, 'r') as temp_sin_sdf_file:
        sin_sdf_lines = temp_sin_sdf_file.readlines()

    # Filter SIN entries
    filtered_sin_sdf_lines = []
    i = 0
    while i < len(sin_sdf_lines):
        if sin_sdf_lines[i].startswith("SIN"):
            while i < len(sin_sdf_lines) and not sin_sdf_lines[i].strip() == "$$$$":
                filtered_sin_sdf_lines.append(sin_sdf_lines[i])
                i += 1
            filtered_sin_sdf_lines.append(sin_sdf_lines[i])  # Append the "$$$$"
        i += 1

    # Write the updated SDF file with SIN entries
    with open(akg_sdf_path, 'w') as sdf_file:
        sdf_file.writelines(filtered_old_sdf_lines)
        sdf_file.writelines(filtered_sin_sdf_lines)