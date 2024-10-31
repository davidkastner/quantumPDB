"""Calculate QM properties of the ligands

**Usage**

>>> protonate.compute_charge("path/to/OUT.sdf")
{"AAG_A331": 0, "SO4_A901": -2, "SO4_A325": -2, 
"SO4_A903": -2, "AKG_A330": -2,  "GOL_A328": 0, 
"GOL_A329": 0, "GOL_A900": 0, "GOL_A904": 0}

"""

import numpy as np
from collections import defaultdict
from typing import List, Tuple, Dict, Any


def read_ligands(path_ligand: str) -> List[List[str]]:
    with open(path_ligand, "r") as f:
        sdf = f.read()
    return [
        [t for t in s.splitlines() if t != ""]
        for s in sdf.split("$$$$")
    ]


def parse_ligand_name(ligand_block: List[str]) -> Tuple[List[Tuple[str, str, str]], str]:
    title = ligand_block[0].split("_")
    name_id = [
        (res_name, chain_id, res_id) 
        for res_name, chain_id, res_id 
        in zip(title[::3], title[1::3], title[2::3])
    ]
    name = " ".join([f"{res_name}_{chain_id}{res_id}" for res_name, chain_id, res_id in name_id])
    return name_id, name


def compute_charge(path_ligand: str, path_pdb: str) -> Dict[str, int]:
    """
    Computes the total charge of each ligand

    Parameters
    ----------
    path: str
        Path to ligand SDF file

    Returns
    -------
    charge: dict
        Keyed by ligand ID
    """
    with open(path_pdb, "r") as f:
        pdb_lines = f.readlines()
    ligands = read_ligands(path_ligand)

    charge = {}
    for l in ligands:
        if not l:
            break
        name_id, name = parse_ligand_name(l)

        c = 0
        n_atom = 0
        for line in l:
            if "V2000" in line:
                n_atom = int(line[:3]) # https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
            if line.startswith("M  RGP"):
                # R# atom is found in the addition between CSO and TAN
                # It's not a real atom and recorded in the RGP entry
                n_atom -= int(line.split()[2])
            if line.startswith("M  CHG"):
                c += sum([int(x) for x in line.split()[4::2]])
        charge[name] = c
        cnt = 0
        for res_name, chain_id, res_id in name_id:
            if res_name == "NO":
                cnt = n_atom
                break
            # to detect removed atoms
            # NO is corrected from NH2OH, thus the deleted 3 hydrogen shouldn't affect the charge
            for line in pdb_lines:
                if line[17:20].strip() == res_name and line[21] == chain_id and line[22:26].strip() == res_id and line[0:3] != "TER":
                    cnt += 1
        charge[name] -= (n_atom - cnt)

    # check Na+, K+, Mg2+, Ca2+
    for line in pdb_lines:
        if len(line) > 26 and line[0:3] != "TER":
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_id = line[22:26].strip()
            name = f"{res_name}_{chain_id}{res_id}"
            if res_name in {"NA", "K"} and name not in charge:
                print(f"Find extra {res_name} in protein, charge + 1")
                charge[name] = 1
            elif res_name in {"MG", "CA"} and name not in charge:
                print(f"Find extra {res_name} in protein, charge + 2")
                charge[name] = 2
    return charge


def compute_spin(path_ligand):
    ligands = read_ligands(path_ligand)
    spin = {}
    for l in ligands:
        if not l:
            break
        name_id, name = parse_ligand_name(l)
        for res_name, _, _ in name_id:
            if res_name == "NO":
                spin[name] = 1
            elif res_name == "OXY":
                spin[name] = 2
    return spin


def get_coord(line_segments: List[str]) -> np.ndarray:
    return np.array([float(x) for x in line_segments[:3]])


def collect_RGP_atoms(path_ligand: str) -> Dict[str, Dict[int, Dict[str, Any]]]:
    ligands = read_ligands(path_ligand)
    RGP_atoms = defaultdict(dict)
    for l in ligands:
        if not l:
            break
        _, name = parse_ligand_name(l)
        n_atom = 0 # initialize n_atom
        RGP_indices = set()
        for i, line in enumerate(l):
            line_segments = line.split()
            
            if i == 3:
                n_atom = int(line_segments[0])
                n_bond = int(line_segments[1])
            if len(line_segments) > 3 and line_segments[3] == "R#":
                RGP_atoms[name][i - 3] = {
                    "coord": get_coord(line_segments),
                }
                RGP_indices.add(i - 3)
            if i >= n_atom + 4 and i < n_atom + 4 + n_bond:
                begin_atom_idx = int(line_segments[0])
                end_atom_idx = int(line_segments[1])
                if begin_atom_idx in RGP_indices:
                    RGP_atoms[name][begin_atom_idx]["linking_atom_coord"] = get_coord(l[end_atom_idx + 3].split())
                elif end_atom_idx in RGP_indices:
                    RGP_atoms[name][end_atom_idx]["linking_atom_coord"] = get_coord(l[begin_atom_idx + 3].split())
    return RGP_atoms
