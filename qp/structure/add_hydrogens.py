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
from warnings import warn
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, Polypeptide
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.NeighborSearch import NeighborSearch


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
    retries = 5
    delay = 60  # seconds

    for _ in range(retries):
        try:
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

            return json.loads(r.text)["id"]  # Exit if successful
        except KeyError:
            print(f"KeyError encountered. Retrying in {delay} seconds...")
            time.sleep(delay)

    raise KeyError(f"Failed to upload the file and retrieve 'id' after {retries} attempts.")


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
    retries = 5
    delay = 60  # seconds

    for _ in range(retries):
        try:
            protoss = requests.post(
                "https://proteins.plus/api/protoss_rest",
                json={"protoss": {"pdbCode": pid}},
                headers={"Accept": "application/json"},
            )
            if protoss.status_code == 400:
                raise ValueError("Invalid PDB code")
            return json.loads(protoss.text)["location"]  # Exit if successful
        except KeyError:
            print(f"KeyError encountered. Retrying in {delay} seconds...")
            time.sleep(delay)

    raise KeyError(f"Failed to submit the PDB code and retrieve 'location' after {retries} attempts.")



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
    # Sometimes the Protoss server doesn't respond correctly with the first query
    retries = 5
    delay = 60  # seconds

    for _ in range(retries):
        try:
            r = requests.get(job)
            while r.status_code == 202:
                time.sleep(1)
                r = requests.get(job)

            protoss = requests.get(json.loads(r.text)[key])
            os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
            with open(out, "w") as f:
                f.write(protoss.text)
            return  # Exit if successful
        except KeyError:
            time.sleep(delay)
            continue  # Retry on failure

    raise KeyError(f"Failed to download the file with key '{key}' after {retries} attempts.")


def repair_ligands(path, orig):
    parser = PDBParser(QUIET=True)
    prot_structure = parser.get_structure("Prot", path)
    orig_structure = parser.get_structure("Orig", orig)

    for res in prot_structure[0].get_residues():
        if res.get_resname() == "MOL":
            resid = res.get_id()
            chain = res.get_parent()
            chain.detach_child(resid)

            missing = []
            found = False
            for r in orig_structure[0][chain.get_id()].get_residues():
                if r.get_id()[1] == resid[1]:
                    found = True
                if found:
                    if r.get_id() not in chain:
                        chain.add(r)
                        missing.append(r)
                    else:
                        break

            for r in missing:
                for a in r.get_unpacked_list():
                    if a.element == "H":
                        r.detach_child(atom.get_id())
            
            for atom in res.get_unpacked_list():
                if atom.element != "H":
                    continue
                closest = None
                for r in missing:
                    for a in r.get_unpacked_list():
                        if closest is None or atom - a < atom - closest:
                            closest = a
                closest.get_parent().add(atom)

    io = PDBIO()
    io.set_structure(prot_structure)
    io.save(path)


def flip_coordinated_HIS(points, res):
    """
    Flip the imidazole ring of HIS if it can be coordinated with the metal

    The conformation and tautomerism states of HIS is adjusted by Protoss. 
    However, sometimes the coordinated nitrogen will be flipped away from 
    where it's supposed to be. This procedure detect if the nitrogen (NE2/ND1) is
    more favorable to be nearer to the metal than its neighbor carbon (CE1/CD2) and
    re-flip the ring and re-arrange the hydrogen.

    Parameters
    ----------
    points: list of Bio.PDB.Atom.Atom
        Metal atoms
    res: Bio.PDB.Residue.Residue
        A HIS residue
    """
    flip_flag = ""
    try:
        CE1 = res["CE1"]
        NE2 = res["NE2"]
        CD2 = res["CD2"]
        ND1 = res["ND1"]
        CB = res["CB"]
        CG = res["CG"]
    except IndexError:
        print("Non standard atom names of HIS")
        return

    closest = sorted(points, key=lambda p: min(p - x for x in [CE1, NE2, CD2, ND1]))[0]
    dist_CE1_metal = closest - CE1
    dist_NE2_metal = closest - NE2
    dist_CD2_metal = closest - CD2
    dist_ND1_metal = closest - ND1
    if dist_CE1_metal < dist_NE2_metal and dist_CE1_metal < 3.5 and dist_CE1_metal < dist_ND1_metal:
        flip_flag = "E"
    elif dist_CD2_metal < dist_ND1_metal and dist_CD2_metal < 3.5 and dist_CD2_metal < dist_NE2_metal:
        flip_flag = "D"

    if flip_flag:
        old_coords = dict()
        for atom_name in ["HD2", "HE1", "HD1", "HE2", "CD2", "ND1", "CE1", "NE2", "CG", "CB"]:
            if atom_name in res:
                old_coords[atom_name] = res[atom_name].get_coord()
        
        coord_CB, coord_CG = old_coords["CB"], old_coords["CG"]
        axis = coord_CG - coord_CB
        for atom_name in ["HD2", "HE1", "HD1", "HE2", "CD2", "ND1", "CE1", "NE2"]:
            if atom_name in res:
                coord = old_coords[atom_name]
                loc = coord - coord_CB
                proj = axis * np.dot(axis, loc) / np.dot(axis, axis)
                flipped_loc = 2 * proj - loc
                flipped_coord = flipped_loc + coord_CB
                res[atom_name].set_coord(flipped_coord)

        if flip_flag == "E" and "HE2" in res and "HD1" not in res:
            coord_CE1 = res["CE1"].get_coord()
            coord_ND1 = res["ND1"].get_coord()
            vec_ND1_CE1 = coord_CE1 - coord_ND1
            vec_ND1_CG = coord_CG - coord_ND1
            bisector = vec_ND1_CE1 + vec_ND1_CG
            vec_ND1_HD1 = - bisector / np.linalg.norm(bisector) # * 1.00
            res["HE2"].set_coord(vec_ND1_HD1 + coord_ND1)
            res["HE2"].name = "HD1"

        if flip_flag == "D" and "HD1" in res and "HE2" not in res:
            coord_CE1 = res["CE1"].get_coord()
            coord_NE2 = res["NE2"].get_coord()
            coord_CD2 = res["CD2"].get_coord()
            vec_NE2_CE1 = coord_CE1 - coord_NE2
            vec_NE2_CD2 = coord_CD2 - coord_NE2
            bisector = vec_NE2_CE1 + vec_NE2_CD2
            vec_NE2_HE2 = - bisector / np.linalg.norm(bisector) # * 1.00
            res["HD1"].set_coord(vec_NE2_HE2 + coord_NE2)
            res["HD1"].name = "HE2"


def add_hydrogen_CSO(res, structure):
    """
    Add hydrogens to CSO

    CSO is a non-standard amino acid which cannot be identified by Protoss REST API.
    This procedure add hydrogens by geometric rules to its carbons 
    and the hydroxyl oxygen (not added in the special case where a TAN ligand
    is reacting with CSO's hydroxyl group).

    Parameters
    ----------
    res: Bio.PDB.Residue.Residue
        A CSO residue
    structure: Bio.PDB.Structure.Structure
        The whole structure of the protein
    """
    coords = dict()
    try:
        for atom_name in ["C", "N", "CA", "CB", "SG", "OD"]:
            coords[atom_name] = res[atom_name].get_coord()
    except IndexError:
        print("Non standard atom names of CSO")
        return
    if "H_1" not in res and "HA" not in res:
        vec_CA_C = coords["C"] - coords["CA"]
        vec_CA_N = coords["N"] - coords["CA"]
        vec_CA_CB = coords["CB"] - coords["CA"]
        vec_sum = vec_CA_C + vec_CA_N + vec_CA_CB
        vec_CA_HA = - vec_sum / np.linalg.norm(vec_sum) # * 1.00
        coord_HA = vec_CA_HA + coords["CA"]
        res.add(Atom("HA", coord_HA, 0, 1, " ", "HA", None, "H"))
    if "H_2" not in res and "H_3" not in res and "HB1" not in res and "HB2" not in res:
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
    if "H_4" not in res and "HD" not in res:
        for res2 in structure[0].get_residues():
            if res2.get_resname() == "TAN":
                try:
                    C = res2["C"]
                except:
                    print("Non standard atom names of TAN")
                if C - res["OD"] < 1.6 and "H1" in res2:
                    # TAN reacting with CSO, ref. 3X20, 3X24, 3X25
                    print("TAN reacting with CSO!")
                    return
        vec_CB_SG = coords["SG"] - coords["CB"]
        vec_OD_HD = vec_CB_SG / np.linalg.norm(vec_CB_SG) # * 1.00
        coord_HD = coords["OD"] + vec_OD_HD
        res.add(Atom("HD", coord_HD, 0, 1, " ", "HD", None, "H"))


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
            atoms = res.get_unpacked_list()
            if len(atoms) > 1:
                continue # skip if not metal-like
            points.append(atoms[0])

    for res in structure[0].get_residues():
        if res.get_resname() == "HIS" and points:
            flip_coordinated_HIS(points, res)
        if res.get_resname() == "CSO":
            add_hydrogen_CSO(res, structure)

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
            elif is_aa(res) and atom.get_name() in ["H", "H2"] and "N" in res:
                coord = res["N"]

            if coord:
                for p in points:
                    if p - coord < 3:
                        return False
            if atom.element == "H" and res.get_resname() != "HOH" and not is_aa(res):
                for p in points:
                    if p - atom < 1.5:
                        return False
            return True

    io = PDBIO()
    io.set_structure(structure)
    io.save(path, AtomSelect())


def res_priority(res, res_info, center_residues):
    if Polypeptide.is_aa(res, standard=True) or res.get_resname() in center_residues:
        return 4e7
    else:
        return res_info[res]["num_atom"] + \
            int(res_info[res]["avg_occupancy"] * 100) * 1e5 + \
            int(Polypeptide.is_aa(res)) * 2e7


def partial_res_info(res, res_info):
    _, _, chain, code = res.get_full_id()
    _, resid, _ = code
    resname = res.get_resname()
    return f'resname {resname}_{chain}{resid} with averaged occupancy {res_info[res]["avg_occupancy"]:.2f} and {res_info[res]["num_atom"]} atoms'


def clean_partial_occupancy(path, path_new, center_residues):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure("PDB", path)
    io.set_structure(structure)
    kept_partial_res = dict()
    partial_atom_list = []
    for res in structure[0].get_residues():
        for atom in res:
            atom_occupancy = atom.get_occupancy()
            if atom_occupancy < 1.0:
                partial_atom_list.append(atom)
                if res not in kept_partial_res:
                    kept_partial_res[res] = {
                        "avg_occupancy": atom_occupancy,
                        "num_atom": len(res),
                        "kept": True
                    }
                else:
                    kept_partial_res[res]["avg_occupancy"] += atom_occupancy
    if not partial_atom_list:
        return False
    for res in kept_partial_res:
        kept_partial_res[res]["avg_occupancy"] /= kept_partial_res[res]["num_atom"]
    for res in kept_partial_res:
        if not kept_partial_res[res]["kept"]:
            continue
        search = NeighborSearch(partial_atom_list)
        for atom in res:
            if not kept_partial_res[res]["kept"]:
                break
            for clash_res in search.search(center=atom.get_coord(), radius=3, level="R"):
                if clash_res == res or not kept_partial_res[clash_res]["kept"]:
                    continue
                p1 = res_priority(res, kept_partial_res, center_residues)
                p2 = res_priority(clash_res, kept_partial_res, center_residues)
                if p1 > p2:
                    kept_partial_res[clash_res]["kept"] = False
                    warn(f"{partial_res_info(clash_res, kept_partial_res)} is deleted in comparison with {partial_res_info(res, kept_partial_res)}")
                elif p1 < p2:
                    kept_partial_res[res]["kept"] = False
                    warn(f"{partial_res_info(res, kept_partial_res)} is deleted in comparison with {partial_res_info(clash_res, kept_partial_res)}")
                    break

    class ResSelect(Select):
        def accept_residue(self, residue):
            return (residue not in kept_partial_res) or kept_partial_res[residue]["kept"]

    io.save(path_new, ResSelect())
    return True


def compute_charge(path_ligand, path_pdb):
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
    with open(path_ligand, "r") as f:
        sdf = f.read()
    with open(path_pdb, "r") as f:
        pdb_lines = f.readlines()
    ligands = [
        [t for t in s.splitlines() if t != ""]
        for s in sdf.split("$$$$")
        if s != "\n" and s != ""
    ]

    charge = {}
    for l in ligands:
        title = l[0].split("_")
        name_id = [
            (res_name, chain_id, res_id) 
            for res_name, chain_id, res_id 
            in zip(title[::3], title[1::3], title[2::3])
        ]
        name = " ".join([f"{res_name}_{chain_id}{res_id}" for res_name, chain_id, res_id in name_id])

        c = 0
        n_atom = 0
        for line in l:
            if "V2000" in line:
                n_atom = int(line.split()[0][:3])
            if line.startswith("M  RGP"):
                # R# atom is found in the addition between CSO and TAN
                # It's not a real atom and recorded in the RGP entry
                n_atom -= sum([int(x) for x in line.split()[4::2]])
            if line.startswith("M  CHG"):
                c += sum([int(x) for x in line.split()[4::2]])
                break
        charge[name] = c
        cnt = 0
        for res_name, chain_id, res_id in name_id:
            if res_name == "NO":
                cnt = n_atom
                break
            # to detect removed atoms
            # NO is corrected from NH2OH, thus the deleted 3 hydrogen shouldn't affect the charge
            for line in pdb_lines:
                if line[17:20].strip() == res_name and line[21] == chain_id and line[22:26].strip() == res_id:
                    cnt += 1
        charge[name] -= (n_atom - cnt)
    return charge


def compute_spin(path_ligand):
    with open(path_ligand, "r") as f:
        sdf = f.read()

    ligands = [
        [t for t in s.splitlines() if t != ""]
        for s in sdf.split("$$$$")
        if s != "\n" and s != ""
    ]
    spin = {}
    for l in ligands:
        title = l[0].split("_")
        name_id = [
            (res_name, chain_id, res_id) 
            for res_name, chain_id, res_id 
            in zip(title[::3], title[1::3], title[2::3])
        ]
        name = " ".join([f"{res_name}_{chain_id}{res_id}" for res_name, chain_id, res_id in name_id])
        for res_name, chain_id, res_id in name_id:
            if res_name == "NO":
                spin[name] = 1
            elif res_name == "OXY":
                spin[name] = 2
    return spin