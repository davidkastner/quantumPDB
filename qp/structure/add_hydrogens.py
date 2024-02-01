"""Add hydrogens using Protoss

**Usage**

#. Submitting existing or custom PDB file::

    >>> from qp.structure import add_hydrogens
    >>> pid = add_hydrogens.upload("path/to/PDB.pdb")
    >>> job = add_hydrogens.submit(pid)
    >>> add_hydrogens.download(job, "path/to/OUT.pdb")

#. Submitting PDB code::

    >>> from qp.structure import add_hydrogens
    >>> pdb = "1dry"
    >>> job = add_hydrogens.submit(pdb)
    >>> add_hydrogens.download(job, "path/to/OUT.pdb")
    >>> add_hydrogens.download(job, "path/to/OUT.sdf", "ligands")
    >>> add_hydrogens.compute_charge("path/to/OUT.sdf")
    {"AAG_A331": 0, "SO4_A901": -2, "SO4_A325": -2, 
     "SO4_A903": -2, "AKG_A330": -2,  "GOL_A328": 0, 
     "GOL_A329": 0, "GOL_A900": 0, "GOL_A904": 0}

Protoss automatically removes alternative conformations and overlapping entries. 
Download the log file (``key="log"`` in ``add_hydrogens.download``) to see affected atoms. 

Some metal-coordinating residues may be incorrectly protonated. Use 
``add_hydrogens.adjust_activesites(path, metals)`` with the metal IDs to deprotonate
these residues. 
"""

import os
import requests
import json
import time
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom


# NON_STANDARD_AA_PROTOSS = {
#     "CSO": "CSD"
# }


# def replace(entry):
#     for AA1, AA2 in NON_STANDARD_AA_PROTOSS.items():
#         entry = entry.replace(AA1, AA2)
#     return entry
    

# def fix(path):
#     """
#     Replace the non-standard AA that can't be identified by Protoss API with valid ones
#     """
#     with open(path, "r") as f:
#         lines = f.readlines()
    
#     with open(path, "w") as f:
#         f.writelines(map(replace, lines))


def upload(path):
    """
    Uploads a PDB file to the ProteinsPlus web server

    Parameters
    ----------
    path: str
        Path to PDB file

    Returns
    -------
    pid: str
        ProteinsPlus ID
    """
    pp = requests.post(
        "https://proteins.plus/api/pdb_files_rest",
        files={"pdb_file[pathvar]": open(path, "rb")},
    )
    if pp.status_code == 400:
        raise ValueError("Bad request")
    loc = json.loads(pp.text)["location"]

    r = requests.get(loc)
    while r.status_code == 202:
        time.sleep(1)
        r = requests.get(loc)
    return json.loads(r.text)["id"]


def submit(pid):
    """
    Submits a PDB code to the Protoss web API

    Parameters
    ----------
    pid: str
        PDB code or ProteinsPlus ID

    Returns
    -------
    job: str
        URL of the Protoss job location
    """
    protoss = requests.post(
        "https://proteins.plus/api/protoss_rest",
        json={"protoss": {"pdbCode": pid}},
        headers={"Accept": "application/json"},
    )
    if protoss.status_code == 400:
        raise ValueError("Invalid PDB code")
    return json.loads(protoss.text)["location"]


def download(job, out, key="protein"):
    """
    Downloads a Protoss output file

    Parameters
    ----------
    job: str
        URL of the Protoss job location
    out: str
        Path to output file
    key: str
        Determines which file to download ("protein", "ligand", or "log")
    """
    r = requests.get(job)
    while r.status_code == 202:
        time.sleep(1)
        r = requests.get(job)

    protoss = requests.get(json.loads(r.text)[key])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as f:
        f.write(protoss.text)


def flip_coordiating_HIS(points, res):
    flip_flag = False
    try:
        CE1 = res["CE1"]
        NE2 = res["NE2"]
    except IndexError:
        print("Non standard atom names of HIS")
        return
    for p in points:
        dist_CE1_metal = p - CE1
        dist_NE2_metal = p - NE2
        if dist_CE1_metal < dist_NE2_metal and dist_CE1_metal < 3.5:
            flip_flag = True
            break
    if flip_flag:
        old_coords = dict()
        try:
            for atom_name in ["HD2", "HE1", "CD2", "ND1", "CE1", "NE2", "CG", "CB"]:
                old_coords[atom_name] = res[atom_name].get_coord()
        except IndexError:
            print("Non standard atom names of HIS")
            return
        
        coord_CB, coord_CG = old_coords["CB"], old_coords["CG"]
        axis = coord_CG - coord_CB
        for atom_name in ["HD2", "HE1", "CD2", "ND1", "CE1", "NE2"]:
            coord = old_coords[atom_name]
            loc = coord - coord_CB
            proj = axis * np.dot(axis, loc) / np.dot(axis, axis)
            flipped_loc = 2 * proj - loc
            flipped_coord = flipped_loc + coord_CB
            res[atom_name].set_coord(flipped_coord)

        if "HE2" in res and "HD1" not in res:
            coord_CE1 = res["CE1"].get_coord()
            coord_ND1 = res["ND1"].get_coord()
            vec_ND1_CE1 = coord_CE1 - coord_ND1
            vec_ND1_CG = coord_CG - coord_ND1
            bisector = vec_ND1_CE1 + vec_ND1_CG
            vec_ND1_HD1 = - bisector / np.linalg.norm(bisector) # * 1.00
            res["HE2"].set_coord(vec_ND1_HD1 + coord_ND1)
            res["HE2"].name = "HD1"


def add_hydrogen_CSO(res):
    coords = dict()
    try:
        for atom_name in ["C", "N", "CA", "CB", "SG"]:
            coords[atom_name] = res[atom_name].get_coord()
    except IndexError:
        print("Non standard atom names of CSO")
        return
    if "HA" not in res:
        vec_CA_C = coords["C"] - coords["CA"]
        vec_CA_N = coords["N"] - coords["CA"]
        vec_CA_CB = coords["CB"] - coords["CA"]
        vec_sum = vec_CA_C + vec_CA_N + vec_CA_CB
        vec_CA_HA = - vec_sum / np.linalg.norm(vec_sum) # * 1.00
        coord_HA = vec_CA_HA + coords["CA"]
        res.add(Atom("HA", coord_HA, 0, 1, " ", "HA", None, "H"))
    if "HB1" not in res and "HB2" not in res:
        ANGLE_HCH = 109.51 / 180 * np.pi
        vec_CB_SG = coords["SG"] - coords["CB"]
        vec_CB_CA = coords["CA"] - coords["CB"]
        vec_xy = vec_CB_SG + vec_CB_CA
        vec_xy = - vec_xy / np.linalg.norm(vec_xy) * np.cos(ANGLE_HCH / 2) # * 1.00
        vec_z = np.cross(vec_CB_SG, vec_CB_CA)
        vec_z = vec_z / np.linalg.norm(vec_z) * np.sin(ANGLE_HCH / 2) # * 1.00
        coord_HB1 = coords["CB"] + vec_xy + vec_z
        coord_HB2 = coords["CB"] + vec_xy - vec_z
        res.add(Atom("HB1", coord_HB1, 0, 1, " ", "HB1", None, "H"))
        res.add(Atom("HB2", coord_HB2, 0, 1, " ", "HB2", None, "H"))
    return res


def adjust_activesites(path, metals):
    """
    Deprotonates metal-coordinating residues that are (most likely) incorrectly
    protonated by Protoss. Removes hydrogens from coordinating tyrosines and
    cysteines, N-ligands and backbones using a distance cutoff of 3 A. Removes hydrogens
    from NO, which is protonated as an hydroxylamine.

    Parameters
    ----------
    path: str
        Path to existing Protoss output file, will be overwritten
    metals: list
        List of active site metal IDs
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", path)

    points = []
    for res in structure[0].get_residues():
        if res.get_resname() in metals:
            points.append(res.get_unpacked_list()[0])

    for res in structure[0].get_residues():
        if res.get_resname() == "HIS":
            flip_coordiating_HIS(points, res)
        if res.get_resname() == "CSO":
            add_hydrogen_CSO(res)

    class AtomSelect(Select):
        def accept_atom(self, atom):
            res = atom.get_parent()
            if res.get_resname() == "NO":
                if "H" in atom.get_name():
                    return False
                
            coord = None
            if atom.get_name() == "HH" and res.get_resname() == "TYR":
                coord = res["OH"]
            elif atom.get_name() == "HG" and res.get_resname() == "CYS":
                coord = res["SG"]
            elif atom.get_name() in ["H", "H2"] and "N" in res:
                coord = res["N"]

            if coord:
                for p in points:
                    if p - coord < 3:
                        return False
            return True

    io = PDBIO()
    io.set_structure(structure)
    io.save(path, AtomSelect())


def compute_charge(path):
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
    with open(path, "r") as f:
        sdf = f.read()
    ligands = [
        [t for t in s.splitlines() if t != ""]
        for s in sdf.split("$$$$")
        if s != "\n" and s != ""
    ]

    charge = {}
    for l in ligands:
        n = l[0].split("_")
        name = " ".join([f"{a}_{b}{c}" for a, b, c in zip(n[::3], n[1::3], n[2::3])])
        c = 0
        for line in l:
            if line.startswith("M  CHG"):
                c += sum([int(x) for x in line.split()[4::2]])
                break
        charge[name] = c
    return charge
