import numpy as np
import shutil
from Bio.PDB import PDBParser, PDBIO

def read_pdb(pdb_path):
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
    with open(pdb_path, 'w') as file:
        for atom in atoms:
            if atom['record'] != 'DELETE':
                file.write(f"{atom['record']:<6}{atom['serial']:>5} {atom['name']:<4}{atom['altLoc']:<1}{atom['resName']:<3} {atom['chainID']:<1}{atom['resSeq']:>4}{atom['iCode']:<1}   {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}{atom['occupancy']:>6.2f}{atom['tempFactor']:>6.2f}          {atom['element']:>2}{atom['charge']:>2}\n")

def find_and_convert_ligands_to_oxo(atoms, iron_names, distance_cutoff=2.5):
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
                        ligand['x'], ligand['y'], ligand['z'] = iron_coord + 1.6 * (ligand_coord - iron_coord) / distance
                        for other_atom in atoms:
                            if other_atom['resSeq'] == ligand['resSeq'] and other_atom['chainID'] == ligand['chainID'] and other_atom['serial'] != ligand['serial']:
                                other_atom['record'] = 'DELETE'
                        oxo_placed.add(iron['serial'])
                        print(f'   > Converted {ligand["serial"]} {ligand["resName"]}/{ligand["resName"]} to OXO at distance {round(distance, 2)}')
                        break
    return oxo_placed

def find_histidines(atoms, iron):
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
    histidines = find_histidines(atoms, iron)
    if len(histidines) != 2:
        print(f"   > Error: Expected two His coordinating iron {iron['serial']}")
        return False

    iron_coord = np.array([iron['x'], iron['y'], iron['z']])
    potential_oxos = []

    for ne2 in histidines:
        ne2_coord = np.array([ne2['x'], ne2['y'], ne2['z']])
        oxo_coord = iron_coord + 1.6 * (iron_coord - ne2_coord) / np.linalg.norm(iron_coord - ne2_coord)
        print(f"   > Checking potential oxo site opposite His {ne2['serial']}")

        clash_threshold = 1.8
        clash = False
        for atom in atoms:
            if atom['serial'] != iron['serial'] and not (atom['resName'] == 'AKG' and atom['name'] in ['C1', 'O1', 'O2']):  # Exclude specific AKG atoms from clash detection
                distance = np.linalg.norm(np.array([atom['x'], atom['y'], atom['z']]) - oxo_coord)
                if distance < clash_threshold:
                    clash = True
                    print(f"   > Clash found with atom {atom['serial']} ({atom['name']} {atom['resName']}) at distance {round(distance, 2)}")
                    break
        if not clash:
            potential_oxos.append((oxo_coord, ne2['serial']))
            print(f"   > Potential oxo site opposite His {ne2['serial']}.")

    # Find AKG coordinates for dihedral calculation
    try:
        akg_c2 = next(atom for atom in atoms if atom['resName'] == 'AKG' and atom['name'] == 'C2')
        akg_o5 = next(atom for atom in atoms if atom['resName'] == 'AKG' and atom['name'] == 'O5')
    except StopIteration:
        print("   > Error: Could not find AKG atoms (C2 and O5) for dihedral calculation.")
        return False
    
    c2_coord = np.array([akg_c2['x'], akg_c2['y'], akg_c2['z']])
    o5_coord = np.array([akg_o5['x'], akg_o5['y'], akg_o5['z']])

    best_site = None
    best_dihedral = None

    # Calculate dihedral for all potential sites, including clashing ones
    for ne2 in histidines:
        ne2_coord = np.array([ne2['x'], ne2['y'], ne2['z']])
        oxo_coord = iron_coord + 1.6 * (iron_coord - ne2_coord) / np.linalg.norm(iron_coord - ne2_coord)
        dihedral = abs(calculate_dihedral(c2_coord, o5_coord, iron_coord, oxo_coord))
        print(f"   > Dihedral angle for potential OXO site opposite His {ne2['serial']}: {round(dihedral, 2)}")

    # Select the best site among the non-clashing ones
    for oxo_coord, _ in potential_oxos:
        dihedral = abs(calculate_dihedral(c2_coord, o5_coord, iron_coord, oxo_coord))
        print(f"   > Dihedral angle for non-clashing site: {round(dihedral, 2)}")
        if best_dihedral is None or abs(dihedral - 90) < abs(best_dihedral - 90):
            best_site = oxo_coord
            best_dihedral = dihedral

    if best_site is None:
        print("   > Error: Could not determine the best site for OXO")
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
    print(f'   > OXO placed at {rounded_best_site} with dihedral {round(best_dihedral, 2)}\n')
    return True

def convert_akg_to_sin(atoms):
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

def add_oxo_and_sin(pdb_path):
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
        print("   > Error: Could not place OXO")

    convert_akg_to_sin(atoms)
    atoms = [atom for atom in atoms if atom['record'] != 'DELETE']
    write_pdb(atoms, pdb_path)

def remove_oxo_hydrogens(protoss_pdb):
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

def update_ligands_sdf(pdb_path, sdf_path):
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
