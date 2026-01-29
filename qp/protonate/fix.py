"""Fix protonation errors in metal coordination sites."""

import numpy as np
from typing import List
from queue import Queue
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser, PDBIO, Select, Polypeptide
from qp.cluster.spheres import CenterResidue


def res_priority(res, res_info, center_residues):
    """Compute a priority score for a residue during partial occupancy resolution.

    Higher scores mean higher priority (kept over lower-priority residues).
    Center residues get the highest priority, followed by standard amino acids,
    then non-standard residues scored by atom count and average occupancy.

    Parameters
    ----------
    res : Bio.PDB.Residue.Residue
        The residue to score.
    res_info : dict
        Partial occupancy info dict with ``'avg_occupancy'`` and ``'num_atom'`` keys.
    center_residues : CenterResidue
        The center residue definition for the current structure.

    Returns
    -------
    float
        Priority score.
    """
    if res in center_residues:
        return 8e7
    elif Polypeptide.is_aa(res, standard=True):
        return 4e7
    else:
        return res_info[res]["num_atom"] + \
            int(res_info[res]["avg_occupancy"] * 100) * 1e5 + \
            int(Polypeptide.is_aa(res)) * 2e7


def partial_res_info(res, res_info):
    """Format a residue's partial occupancy information as a human-readable string.

    Parameters
    ----------
    res : Bio.PDB.Residue.Residue
        The residue to describe.
    res_info : dict
        Partial occupancy info dict with ``'avg_occupancy'`` and ``'num_atom'`` keys.

    Returns
    -------
    str
        Formatted string, e.g. ``'residue FE_A199 with avg. occupancy 0.75 and 1 atoms'``.
    """
    _, _, chain, code = res.get_full_id()
    _, resid, _ = code
    resname = res.get_resname()
    return f'residue {resname}_{chain}{resid} with avg. occupancy {res_info[res]["avg_occupancy"]:.2f} and {res_info[res]["num_atom"]} atoms'


def clean_occupancy(path: str, center_residues: CenterResidue) -> bool:
    """Resolve partial occupancy conflicts by selecting a self-consistent coordinate set.

    Iterates over all residues with partial occupancy and removes lower-priority
    residues that clash spatially. Priority order: center residues > standard
    amino acids > non-standard residues > higher average occupancy > more atoms.
    Distance criteria are 1.0 A for center residues (allowing coordination bonds)
    and 2.5 A for all others.

    The PDB file at ``path`` is overwritten in place with the cleaned structure.

    Parameters
    ----------
    path : str
        Path to the PDB file (modified in place).
    center_residues : CenterResidue
        The center residue definition for priority assignment.

    Returns
    -------
    bool
        True if partial occupancy was detected and cleaned, False otherwise.
    """
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure("PDB", path)
    io.set_structure(structure)
    kept_partial_res = dict()
    partial_atom_list = []

    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                partial_in_res = False
                total_occupancy = 0

                # Check if there are partial occupied atoms in res
                for atom in res:
                    atom_occupancy = atom.get_occupancy()
                    total_occupancy += atom_occupancy
                    if atom_occupancy < 1.0:
                        partial_in_res = True
                        partial_atom_list.append(atom)
                
                # If so, add res to the partial res dict
                if partial_in_res:
                    num_atom = len(res)
                    kept_partial_res[res] = {
                        "avg_occupancy": total_occupancy / num_atom,
                        "num_atom": num_atom,
                        "kept": True,
                        "search_radius": 1.0 if res in center_residues else 2.5 # allow for bonding to center, otherwise vdw
                    }
    
    # no partial occupancy detected
    if not partial_atom_list:
        return False
    
    check_queue = Queue()
    for res in kept_partial_res:
        check_queue.put(res)

    search = NeighborSearch(partial_atom_list)
    while not check_queue.empty():
        res = check_queue.get()

        # don't check deleted res
        if not kept_partial_res[res]["kept"]:
            continue

        recheck_res = set()
        for atom in res:
            for clash_res in search.search(center=atom.get_coord(), radius=kept_partial_res[res]["search_radius"], level="R"):
                # don't delete self and deleted res
                if (
                    clash_res == res or # self
                    not kept_partial_res[clash_res]["kept"] # deleted res
                ):
                    continue

                # if center is in clash with current atom
                # reconsider this clash from the center's point of view with tighter search radius
                if clash_res in center_residues:
                    # if center doesn't contain partial occupancy, put it into the queue for rechecking,
                    if clash_res not in kept_partial_res:
                        kept_partial_res[clash_res] = {
                            "avg_occupancy": 1.0,
                            "num_atom": len(clash_res),
                            "kept": True,
                            "search_radius": 1.0
                        }
                        recheck_res.add(clash_res)

                    # otherwise the center would have been checked
                    continue

                # compare priority and delete the lower one
                p1 = res_priority(res, kept_partial_res, center_residues)
                p2 = res_priority(clash_res, kept_partial_res, center_residues)
                if p1 > p2:
                    kept_partial_res[clash_res]["kept"] = False
                    print(f"> Deleting {partial_res_info(clash_res, kept_partial_res)}, keeping {partial_res_info(res, kept_partial_res)} in search radius {kept_partial_res[res]['search_radius']}Å")
                elif p1 < p2:
                    kept_partial_res[res]["kept"] = False
                    print(f"> Deleting {partial_res_info(res, kept_partial_res)}, keeping {partial_res_info(clash_res, kept_partial_res)} in search radius {kept_partial_res[res]['search_radius']}Å")
                    break
        
        for res in recheck_res:
            print("Recheck", res)
            check_queue.put(res)

    class ResSelect(Select):
        def accept_residue(self, residue):
            return (residue not in kept_partial_res) or kept_partial_res[residue]["kept"]
    
    # Creating a backup copy is unnecessary although it could be useful so I'm commenting it
    # shutil.copy(path, f"{path[:-4]}_old.pdb")
    io.save(path, ResSelect())
    print("> Partial occupancy cleaning finished")
    return True


def flip_coordinated_HIS(points, res):
    """Flip a histidine imidazole ring to orient nitrogen toward a metal.

    Protoss assigns histidine tautomers and conformations, but sometimes
    the coordinating nitrogen (NE2 or ND1) is oriented away from the metal.
    This function detects when a carbon (CE1 or CD2) is closer to the metal
    than its adjacent nitrogen and flips the ring around the CB--CG axis
    to restore proper coordination geometry.

    Parameters
    ----------
    points : list of Bio.PDB.Atom.Atom
        Metal atoms to check coordination against.
    res : Bio.PDB.Residue.Residue
        A histidine residue (modified in place).

    Notes
    -----
    The residue is modified in place. Hydrogen positions are adjusted
    to match the new ring orientation.
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
        print("> Non standard atom names of HIS")
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
    """Add hydrogens to cysteine sulfenic acid (CSO) residues.

    CSO is a post-translational modification not recognized by Protoss.
    This function adds hydrogens to CA (HA), CB (HB1, HB2), and the
    hydroxyl oxygen (HD) using ideal geometry. The hydroxyl hydrogen
    is omitted if a TAN (2,3,3-trichloroallyl alcohol) ligand is
    covalently bonded to the OD atom.

    Parameters
    ----------
    res : Bio.PDB.Residue.Residue
        A CSO residue (modified in place).
    structure : Bio.PDB.Structure.Structure
        The full protein structure (used to detect TAN--CSO bonds).

    Notes
    -----
    The residue is modified in place. Existing hydrogens are not replaced.
    """
    coords = dict()
    try:
        for atom_name in ["C", "N", "CA", "CB", "SG", "OD"]:
            coords[atom_name] = res[atom_name].get_coord()
    except IndexError:
        print("> Non standard atom names of CSO")
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
                    print("> Non standard atom names of TAN")
                if C - res["OD"] < 1.6 and "H1" in res2:
                    # TAN reacting with CSO, ref. 3X20, 3X24, 3X25
                    print("> TAN reacting with CSO!")
                    return
        vec_CB_SG = coords["SG"] - coords["CB"]
        vec_OD_HD = vec_CB_SG / np.linalg.norm(vec_CB_SG) # * 1.00
        coord_HD = coords["OD"] + vec_OD_HD
        res.add(Atom("HD", coord_HD, 0, 1, " ", "HD", None, "H"))


def fix_OXT(path: str) -> None:
    """Fix corrupted C-terminal OXT atoms by reflecting their coordinates.

    When the OXT and O atoms of a C-terminal residue are too close (< 1.8 A),
    the O atom is reflected across the CA--C bond axis to restore correct
    geometry. The PDB file is overwritten in place.

    Parameters
    ----------
    path : str
        Path to the PDB file (modified in place).
    """
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure("PDB", path)
    io.set_structure(structure)
    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                res: Residue
                OXT_flag = True
                for name in {"OXT", "O", "CA", "C"}:
                    if name not in res:
                        OXT_flag = False
                        break
                if not OXT_flag:
                    continue
                OXT: Atom = res["OXT"]
                O: Atom = res["O"]
                CA: Atom = res["CA"]
                C: Atom = res["C"]
                if OXT - O > 1.8:
                    break
                print("> Fixing corrupted C terminus")
                vec_C_CA = CA.get_coord() - C.get_coord()
                vec_C_OXT = OXT.get_coord() - C.get_coord()
                proj_vec_C_OXT = np.dot(vec_C_CA, vec_C_OXT) / np.dot(vec_C_CA, vec_C_CA) * vec_C_CA
                dist_vec_C_OXT = proj_vec_C_OXT - vec_C_OXT
                reflect_vec_C_OXT = vec_C_OXT + 2 * dist_vec_C_OXT
                O.set_coord(C.get_coord() + reflect_vec_C_OXT)
    
    io.save(path, Select())


def adjust_activesites(path, center_residue: CenterResidue):
    """Deprotonate metal-coordinating residues incorrectly protonated by Protoss.

    Removes hydrogens from residues coordinating to metal centers:

    - Tyrosine OH groups within 3 A of a metal
    - Cysteine SG groups within 3 A of a metal
    - Backbone N atoms within 3 A of a metal
    - All hydrogens on NO ligands (protonated as hydroxylamine by Protoss)
    - Any ligand H atoms within 1.5 A of a metal

    Also flips incorrectly oriented histidine rings and adds hydrogens to
    CSO (cysteine sulfenic acid) residues.

    Parameters
    ----------
    path : str
        Path to the Protoss output PDB file (modified in place).
    center_residue : CenterResidue
        Center residue definition identifying metal atoms.

    Notes
    -----
    The file at ``path`` is overwritten with the corrected structure.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", path)

    points = []
    for res in structure[0].get_residues():
        if res in center_residue:
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