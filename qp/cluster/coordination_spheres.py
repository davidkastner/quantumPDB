"""Extract active site coordination sphere clusters

**Usage**::

    >>> from qp.cluster import coordination_spheres

    >>> coordination_spheres.extract_clusters(
    ...     "path/to/PDB.pdb", 
    ...     "path/to/out/dir/", 
    ...     center_residues=["FE", "FE2"], # List of resnames of the residues to use as the cluster center
    ...     sphere_count=2,              # Number of spheres to extract
    ...     ligands=["AKG"]       # PDB IDs of additional ligands
    ... )

Extracting clusters leaves open valences in the outermost sphere. Capping may be
performed by specifying ``capping`` in ``coordination_spheres.extract_clusters``:

* 0. No capping. (Default)
* 1. Cap with hydrogens.
* 2. Cap with ACE/NME groups. 
"""

import os
from typing import Set, Literal, Optional
import numpy as np
from Bio.PDB import PDBParser, Polypeptide, PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.NeighborSearch import NeighborSearch
from scipy.spatial import Voronoi
from qp.structure import struct_to_file
from sklearn.cluster import DBSCAN


RANDOM_SEED = 66265


def get_grid_coord_idx(coord, coord_min, mean_distance):
    """
    Compute a point's position in a 1D grid
    
    Parameters
    ----------
    coord: float
        One coordinate of a point
    coord_min: float
        The minimum coordinate of the grid
    mean_distance: float
        The distance between neighbors in the grid

    Returns
    -------
    idx: int
        The idx of the point in the grid
    """
    return int((coord - coord_min + mean_distance * 0.5) // mean_distance)


def get_grid(coords, mean_distance):
    """
    Compute the grid's parameter for a given 1D point list

    Parameters
    ----------
    coords: numpy.array
        The coordinates of the 1D point list
    mean_distance: float
        The distance between neighbors in the grid
    
    Returns
    -------
    coord_min: float
        The minimum coordinate of the grid
    coord_max: float
        The maximum coordinate of the grid
    grid: numpy.array
        The 1D grid
    """
    coord_min, coord_max = coords.min() - mean_distance, coords.max() + mean_distance
    npoints = get_grid_coord_idx(coord_max, coord_min, mean_distance)
    return coord_min, coord_max, coord_min + np.linspace(0, npoints - 1, npoints) * mean_distance


def _visualize_dummy(dummy, path):
    # debug function
    with open(path + "/dummy.xyz", "w") as f:
        f.write(f"{len(dummy)}\n\n")
        for coord in dummy:
            f.write(f"He {coord[0]} {coord[1]} {coord[2]}\n")


def fill_dummy(points, mean_distance=3, noise_amp=0.2):
    """
    Fill dummy atoms in a point cloud

    Parameters
    ----------
    points: numpy.array
        The 3D coordinates of the point cloud
    mean_distance: float
        The distance between neighbors in the grid
    noise_amp: float
        The amplitude of the noise of dummy atoms' position
    
    Returns
    -------
    points:
        The 3D coordinates of the point cloud filled with dummy atoms
    """
    conf = np.stack(points, axis=0)
    x_min, _, x_grids = get_grid(conf[:,0], mean_distance)
    y_min, _, y_grids = get_grid(conf[:,1], mean_distance)
    z_min, _, z_grids = get_grid(conf[:,2], mean_distance)
    flags = np.ones((len(x_grids), len(y_grids), len(z_grids)), dtype=bool)
    for point in points:
        x, y, z = point
        flags[
            get_grid_coord_idx(x, x_min, mean_distance),
            get_grid_coord_idx(y, y_min, mean_distance), 
            get_grid_coord_idx(z, z_min, mean_distance)
        ] = False
    flags = flags.flatten()
    np.random.seed(RANDOM_SEED)
    noise = mean_distance * noise_amp * (np.random.rand(len(x_grids), len(y_grids), len(z_grids), 3) - 0.5)
    dummy = np.stack(np.meshgrid(x_grids, y_grids, z_grids, indexing="ij"), axis=-1)
    dummy = (dummy + noise).reshape(-1, 3)
    dummy = dummy[flags, :]
    return np.concatenate([points, dummy], axis=0)


def calc_dist(point_a, point_b):
    """
    Calculate the Euclidean distance between two points

    Parameters
    ----------
    point_a: numpy.array
        Point A
    point_b: numpy.array
        Point B

    Returns
    -------
    dist: float
        The Euclidean distance between Point A and Point B    

    """
    return np.linalg.norm(point_a - point_b)


def voronoi(model, center_residues, ligands, smooth_method, **smooth_params):
    """
    Computes the Voronoi tessellation of a protein structure.

    Parameters
    ----------
    model: Bio.PDB.Model
        Protein structure model

    Returns
    -------
    neighbors: dict
        Adjacency list corresponding to neighboring Voronoi cells
    """
    atoms = []
    points = []
    for res in model.get_residues():
        # if Polypeptide.is_aa(res) or \
        #     res.get_resname() in center_residues or \
        #     (res.get_resname() in ligands):
        for atom in res.get_unpacked_list(): # includes atoms from multiple conformations
            atoms.append(atom)
            points.append(atom.get_coord())

    points_count = len(points)
    if smooth_method == "dummy_atom":
        new_points = fill_dummy(points, **smooth_params)
        vor = Voronoi(new_points)
    else:
        vor = Voronoi(points)

    neighbors = {}
    for a, b in vor.ridge_points:
        if a < points_count and b < points_count:
            dist = calc_dist(points[a], points[b])
            neighbors.setdefault(atoms[a], []).append((atoms[b], dist))
            neighbors.setdefault(atoms[b], []).append((atoms[a], dist))
    return neighbors


def merge_centers(cur, search, seen, radius=0.0):
    if radius == 0.0:
        return {cur}

    center = {cur}
    seen.add(cur)
    nxt = set()
    for atom in cur.get_unpacked_list():
        if atom.element == "H":
            continue
        for res in search.search(atom.get_coord(), radius, "R"):
            nxt.add(res)
            
    for res in nxt:
        if res not in seen:
            center |= merge_centers(res, search, seen, radius)
    return center


def get_center_residues(model, center_residues, merge_cutoff=0.0):
    found = set()
    center_atoms = []
    for res in model.get_residues():
        # We will assume the ligand is labeled as a heteroatom with res.id[0] != ' '
        if res.get_resname() in center_residues and res.id[0] != ' ':
            found.add(res)
            center_atoms.extend([atom for atom in res.get_unpacked_list() if atom.element != "H"])
    if not len(center_atoms):
        print("> WARNING: No matching cluster center found. Skipping cluster.")
        return []
    
    search = NeighborSearch(center_atoms)
    seen = set()
    centers = []
    for res in found:
        if res not in seen:
            centers.append(merge_centers(res, search, seen, merge_cutoff))
    return centers


def box_outlier_thres(data, coeff=1.5):
    """
    Compute the threshold for the boxplot outlier detection method

    Parameters
    ----------
    data: list
        The data for the boxplot statistics
    coeff: float
        The coefficient for the outlier criterion from quartiles of the data

    Returns
    -------
    lb: float
        the lower bound of non-outliers
    ub: float
        the upper bound of non-outliers

    """
    Q3 = np.quantile(data, 0.75)
    Q1 = np.quantile(data, 0.25)
    IQR = Q3 - Q1
    return Q1 - coeff * IQR, Q3 + coeff * IQR


def check_NC(atom, metal):
    """
    Check if a nitrogen / carbon atom in the first sphere is coordinated.

    If the nitrogen / carbon atom is in the backbone, or it's the nearest atom to the metal
    among all atoms in its residue, it's considered coordinated.

    Parameters
    ----------
    atom: Bio.PDB.Atom
        The nitrogen / carbon atom to be checked
    metal:
        The metal atom (coordination center)

    Returns
    -------
    flag: bool
        Whether the atom is coordinated or not
    """
    if atom.get_name() == "N":
        # TODO: check special coordinated nitrogens in the backbone
        return True
    else:
        ref_dist = calc_dist(atom.get_coord(), metal.get_coord())
        res = atom.get_parent()
        for atom in res.get_unpacked_list():
            if calc_dist(atom.get_coord(), metal.get_coord()) < ref_dist:
                return False
        return True


def get_next_neighbors(
    start, neighbors, sphere_count, ligands,
    first_sphere_radius=3,
    smooth_method="boxplot", 
    include_ligands=2,
    **smooth_params):
    """
    Iteratively determines spheres around a given starting atom

    Parameters
    ----------
    start: Bio.PDB.Residue
        Starting metal atom
    neighbors: dict
        Adjacency list of neighboring atoms
    sphere_count: int
        Number of spheres to extract
    ligands: list
        A list of ligands to include
    smooth_method: ("boxplot" | "dbscan" | "dummy_atom")
        The method used to smoothen the spheres
    include_ligands: int
        the mode of including ligands in the sphere
    smooth_params:
        params of the specific smooth method

    Returns
    -------
    metal_id: str
        Active site identifier
    seen: set
        Set of residues from all spheres
    spheres: list of sets
        Sets of residues separated by spheres
    """
    seen = start.copy()
    spheres = [start]
    lig_adds = [set()]
    for i in range(0, sphere_count):
        # get candidate atoms in the new sphere
        nxt = set()
        lig_add = set()
        if i == 0 and first_sphere_radius > 0:
            start_atoms = []
            for res in start:
                start_atoms.extend(res.get_unpacked_list())

            is_metal_like = (len(start_atoms) == len(start))
            search = NeighborSearch([atom for atom in start_atoms[0].get_parent().get_parent().get_parent().get_atoms() if atom.element != "H" and atom not in start_atoms])
            first_sphere = []
            for center in start_atoms:
                first_sphere += search.search(center=center.get_coord(), radius=first_sphere_radius, level="A")
            for atom in first_sphere:
                if atom.get_parent() not in seen:
                    element = atom.element
                    if (
                        not is_metal_like 
                        or element in "OS" 
                        or (element in "NC" and any(check_NC(atom, center) for center in start_atoms))
                    ): # only consider coordinated atoms
                        res = atom.get_parent()
                        seen.add(res)
                        if Polypeptide.is_aa(res):
                            nxt.add(res)
                        else:
                            if (
                                include_ligands != 1 or
                                res.get_resname() != "HOH" # mode 1: exclude all waters
                            ):
                                lig_add.add(res)
        else:
            candidates = []
            frontiers = spheres[-1] if spheres[-1] else spheres[0]
            for res in frontiers:
                for atom in res.get_unpacked_list():
                    for n, dist in neighbors[atom]:
                        par = n.get_parent()
                        candidates.append((n, dist, par))

            # screen candidates
            if smooth_method == "box_plot":
                dist_data = [dist for n, dist, par in candidates]
                lb, ub = box_outlier_thres(dist_data, **smooth_params)
                screened_candidates = [(dist, par) for n, dist, par in candidates if dist < ub]
            elif smooth_method == "dbscan":
                optics = DBSCAN(**smooth_params)
                X = [n.get_coord() for n, dist, par in candidates]
                for res in spheres[-1]:
                    for atom in res.get_unpacked_list():
                        X.append(atom.get_coord())
                cluster_idx = optics.fit_predict(X) + 1
                largest_idx = np.bincount(cluster_idx).argmax()
                screened_candidates = [(dist, par) for i, (n, dist, par) in enumerate(candidates) if cluster_idx[i] == largest_idx]
            else:
                screened_candidates = [(dist, par) for n, dist, par in candidates]

            screened_candidates = [par for dist, par in screened_candidates if i != 0 or dist < 2.7]

            # build the new sphere
            for par in screened_candidates:
                if par not in seen:
                    seen.add(par)
                ## Added to include ligands in the first coordination sphere
                ## Produced unweildy clusters as some ligands are large
                ## Potentially useful in the future
                    if Polypeptide.is_aa(par):
                        nxt.add(par)
                    else:
                        if (
                            (include_ligands == 0 and par.get_resname() in ligands) or 
                            # mode 0: only include ligands in the first sphere unless specified
                            (include_ligands == 1 and par.get_resname() != "HOH") or # mode 1: exclude all waters
                            include_ligands == 2 # mode 2: include everything
                        ):
                            lig_add.add(par)

        spheres.append(nxt)
        lig_adds.append(lig_add)

    metal_id = []
    for res in start:
        res_id = res.get_full_id()
        chain_name = res_id[2]
        metal_index = str(res_id[3][1])
        metal_id.append(chain_name + metal_index)
    for i in range(len(spheres)):
        spheres[i] = spheres[i] | lig_adds[i]
    
    return "_".join(sorted(metal_id)), seen, spheres


def prune_atoms(center, residues, spheres, max_atom_count, ligands):
    """
    Prune residues from the cluster model to meet the max atom count constraint,
    while keeping specified ligands and co-factors.

    Parameters
    ----------
    center: set
        Set of central residues
    residues: set
        Set of residues in the cluster
    spheres: list of sets
        List of residue sets, each corresponding to a coordination sphere
    max_atom_count: int
        Maximum allowed atom count in the cluster
    ligands_to_keep: list
        List of ligand names to keep in the cluster
    """

    atom_cnt = 0
    for res in residues:
        atom_cnt += len(res)
    if atom_cnt <= max_atom_count:
        return

    center_atoms = []
    for c in center:
        center_atoms.extend(c.get_unpacked_list())
    def dist(res):
        return min(atom - x for x in center_atoms for atom in res.get_unpacked_list())
                   
    prune = set()
    for res in sorted(residues, key=dist, reverse=True):
        # Check if the residue is in the ligands_to_keep list
        if res.get_resname() not in ligands:
            prune.add(res)
            atom_cnt -= len(res)
            if atom_cnt <= max_atom_count:
                break

    residues -= prune
    for s in spheres:
        s -= prune
    while not spheres[-1]:
        spheres.pop()


def scale_hydrogen(a, b, scale):
    """
    Replaces an atom with hydrogen, rescaling the original bond length

    Parameters
    ----------
    a: Bio.PDB.Atom
        Bonded atom to keep
    b: Bio.PDB.Atom
        Bonded atom to replace
    scale: float
        Bond length scale

    Returns
    -------
    pos: array of float
        Coordinates of new hydrogen atom
    """
    p = a.get_coord()
    q = b.get_coord()
    return scale * (q - p) + p


def get_normalized_vector(atom1: Atom, atom2: Atom) -> np.array:
    v = atom2.get_coord() - atom1.get_coord()
    return v / np.linalg.norm(v)


def build_hydrogen(parent: Residue, template: Optional[Residue], atom: Literal["N", "C"]):
    """
    Cap with hydrogen, building based on the upstream or downstream residue

    Parameters
    ----------
    chain: Bio.PDB.Chain
        Chain with desired residue
    parent: Bio.PDB.Residue
        Residue to cap
    template: Bio.PDB.Residue
        Upstream or downstream residue
    atom: str
        Flag for adding to the 'N' or 'C' side of the residue

    Returns
    -------
    res: Bio.PDB.Residue
        Residue containing added hydrogen
    """
    if template is not None:
        if atom == "N":
            pos = scale_hydrogen(parent["N"], template["C"], 1 / 1.32)
        else:
            pos = scale_hydrogen(parent["C"], template["N"], 1.09 / 1.32)
    else:
        CA = parent["CA"]
        if atom == "N":
            N = parent["N"]
            H = parent["H"]
            bis = get_normalized_vector(N, CA) + get_normalized_vector(N, H)
            bis /= np.linalg.norm(bis)
            pos = N.get_coord() - bis
        else:
            C = parent["C"]
            O = parent["O"]
            bis = get_normalized_vector(C, CA) + get_normalized_vector(C, O)
            bis /= np.linalg.norm(bis)
            pos = C.get_coord() - bis * 1.09

    for name in ["H1", "H2", "H3"]:
        if name not in parent:
            break
    atom = Atom(name, pos, 0, 1, " ", name, None, "H")
    parent.add(atom)
    return atom


def build_heavy(chain, parent, template, atom):
    """
    Cap with ACE/NME, building based on the upstream or downstream residue

    Parameters
    ----------
    chain: Bio.PDB.Chain
        Chain with desired residue
    parent: Bio.PDB.Residue
        Residue to cap
    template: Bio.PDB.Residue
        Upstream or downstream residue
    atom: str
        Flag for adding to the "N" (ACE) or "C" (NME) side of the residue

    Returns
    -------
    res: Bio.PDB.Residue
        Residue containing added group
    """

    pos = {"CH3": template["CA"].get_coord()}
    if template.get_resname() == "GLY":
        pos["HH31"] = template["HA2"].get_coord()
        pos["HH32"] = template["HA3"].get_coord()
    else:
        pos["HH31"] = template["HA"].get_coord()
        pos["HH32"] = scale_hydrogen(template["CA"], template["CB"], 1.09 / 1.54)
    
    if atom == "N":
        pos["C"] = template["C"].get_coord()
        pos["O"] = template["O"].get_coord()
        pos["HH33"] = scale_hydrogen(template["CA"], template["N"], 1.09 / 1.46)
    else:
        pos["N"] = template["N"].get_coord()
        if template.get_resname() == "PRO":
            pos["H"] = scale_hydrogen(template["N"], template["CD"], 1 / 1.46)
        else:
            pos["H"] = template["H"].get_coord()
        pos["HH33"] = scale_hydrogen(template["CA"], template["C"], 1.09 / 1.51)

    adj_id = ("H_" + ("NME" if atom == "N" else "ACE"), template.get_id()[1], " ")
    skip = chain.has_id(adj_id) # skip building methyl if already present in adjacent cap
    if skip:
        chain[adj_id].detach_child("HH33")

    name = "ACE" if atom == "N" else "NME"
    res_id = ("H_" + name, template.get_id()[1], " ")
    res = Residue(res_id, name, " ")
    for k, v in pos.items():
        if skip and k in ["CH3", "HH31", "HH32", "HH33"]:
            continue
        res.add(Atom(k, v, 0, 1, " ", k, None, k[0]))
    chain.add(res)
    return res


def check_atom_valence(res: Residue, tree: NeighborSearch, atom: Literal["N", "C"], cn: int) -> bool:
    neighbors = tree.search(res[atom].get_coord(), radius=1.8)
    if len(neighbors) > cn:
        return True
    else:
        for neighbor in neighbors:
            if neighbor.get_name() == "C" and atom == "N":
                return True
            elif neighbor.get_name() == "N" and atom == "C":
                return True
    return False 


def cap_chains(model: Model, residues: Set[Residue], capping: int) -> Set[Residue]:
    """
    Cap chain breaks for a set of extracted residues

    Parameters
    ----------
    model: Bio.PDB.Model
        Protein structure model
    residues: set
        Set of residues
    capping: int
        Flag for capping group, H (1) or ACE/NME (2)

    Returns
    -------
    cap_residues: set
        Set of residues containing added groups
    """
    orig_chains = {}
    for chain in model:
        orig_chains[chain.get_id()] = chain.get_unpacked_list()

    cap_residues = set()

    for res in list(sorted(residues)):
        if not Polypeptide.is_aa(res) or res.get_id()[0] != " ":
            continue

        res_id = res.get_full_id()
        chain = model[res_id[2]]
        chain_tree = NeighborSearch(list(chain.get_atoms()))
        chain_list = orig_chains[chain.get_id()]
        ind = chain_list.index(res)

        if ind > 0:
            pre = chain_list[ind - 1]
            if (
                pre.get_id()[1] == res_id[3][1] - 1
                and pre.get_id()[0] == " "
                and pre not in residues
                and Polypeptide.is_aa(pre)
            ):  # ignores hetero residues
                if capping == 1:
                    cap_residues.add(build_hydrogen(res, pre, "N"))
                else:
                    cap_residues.add(build_heavy(res, pre, "N"))
            elif not check_atom_valence(res, chain_tree, "N", 3):
                cap_residues.add(build_hydrogen(res, None, "N"))

        if ind < len(chain_list) - 1:
            nxt = chain_list[ind + 1]
            if (
                nxt.get_id()[1] == res_id[3][1] + 1
                and nxt.get_id()[0] == " "
                and nxt not in residues
                and Polypeptide.is_aa(nxt)
            ):
                if capping == 1:
                    cap_residues.add(build_hydrogen(res, nxt, "C"))
                else:
                    cap_residues.add(build_heavy(res, nxt, "C"))
            elif not check_atom_valence(res, chain_tree, "C", 3):
                cap_residues.add(build_hydrogen(res, None, "C"))

    return cap_residues


def write_pdb(io, sphere, out):
    """
    Write coordination sphere to PDB file.

    Parameters
    ----------
    io: Bio.PDB.PDBIO
        Bio.PDB writer
    sphere: list
        List of coordination sphere residues
    out: str
        Path to output PDB file
    """

    class ResSelect(Select):
        def accept_residue(self, residue):
            return residue in sphere

    io.save(out, ResSelect())


def residue_in_ligands(resname, resid, res_is_aa, ligand_keys):
    res_key = f"{resname}_{resid[2]}{resid[3][1]}"
    if res_is_aa:
        return res_key in ligand_keys
    else:
        for ligand_key in ligand_keys:
            ligand_res_keys = ligand_key.split()
            if res_key in ligand_res_keys:
                return True
        return False


def compute_charge(spheres, structure, ligand_charge):
    """
    Computes the total charge of coordinating AAs

    Parameters
    ----------
    spheres: list of sets
        Sets of residues separated by spheres
    structure: Bio.PDB.Structure
        The protein structure
    ligand_charge: dict
        Key, value pairs of ligand names and charges

    Returns
    -------
    charge: list
        Total charge of AAs in each sphere
    """
    # Identifying N-terminal and C-terminal residues for each chain
    n_terminals = set()
    c_terminals = set()
    # Loop over the residues to get first and last as indices may be different
    for chain in structure.get_chains():
        chain_residues = list(chain.get_residues())
        if chain_residues:
            n_terminals.add(chain_residues[0].get_full_id())
            c_terminals.add(chain_residues[-1].get_full_id())

    pos = {
        "ARG": ["HE", "HH11", "HH12", "HH21", "HH22"],
        "LYS": ["HZ1", "HZ2", "HZ3"],
        "HIS": ["HD1", "HD2", "HE1", "HE2"],
        "MLZ": [],
        "M3L": []
    }
    neg = {
        "ASP": ["HD2", "HOD1", "HOD2"],
        "GLU": ["HE2", "HOE1", "HOE2"],
        "CYS": ["HG"],
        "TYR": ["HH"],
        "OCS": [],
        "CSD": ["HD1", "HD2"],
        "KCX": ["HQ1", "HQ2", "HOQ1", "HOQ2"],
        "HIS": ["HD1", "HE2"]
    }

    charge = []
    for s in spheres[1:]:
        c = 0
        for res in s:
            res_id = res.get_full_id()
            resname = res.get_resname()
            res_is_aa = Polypeptide.is_aa(res)
            if not residue_in_ligands(resname, res_id, res_is_aa, ligand_charge.keys()):
                if resname in pos and all(res.has_id(h) for h in pos[resname]):
                    c += 1
                elif resname in neg and all(not res.has_id(h) for h in neg[resname]):
                    c -= 1
                if res_is_aa and resname != "PRO" and all(not res.has_id(h) for h in ["H", "H2"]):
                    # TODO: termini
                    c -= 1

                # Check for charged N-terminus
                if res_id in n_terminals \
                    and res.has_id("N"): # exclude sugar chain terminus
                    c += 1

                # Check for charged C-terminus
                if res.has_id("OXT"):
                    c -= 1

        charge.append(c)
    return charge


def count_residues(spheres):
    """
    Counts the frequency of coordinating residues

    Parameters
    ----------
    spheres: list of sets
        Sets of residues separated by spheres

    Returns
    -------
    count: list of dicts
        Frequency table by sphere
    """
    count = []
    for s in spheres[1:]:
        c = {}
        for res in s:
            c[res.get_resname()] = c.get(res.get_resname(), 0) + 1
        count.append(c)
    return count


def make_res_key(res):
    resname = res.get_resname()
    resid = res.get_id()[1]
    chainid = res.get_parent().get_id()
    return f"{resname}_{chainid}{resid}"    


def complete_oligomer(ligand_keys, model, residues, spheres, include_ligands):
    ligand_res_found = dict()
    oligomer_found = dict()
    for ligand_key in ligand_keys:
        ligand_res_keys = ligand_key.split()
        if len(ligand_res_keys) == 1:
            continue
        oligomer_found[ligand_key] = dict()
        for ligand_res_key in ligand_res_keys:
            ligand_res_found[ligand_res_key] = {
                "sphere": -1,
                "oligomer": ligand_key
            }
            oligomer_found[ligand_key][ligand_res_key] = False
    if not oligomer_found:
        return
    for i, sphere in enumerate(spheres):
        if include_ligands == 0 and i > 0:
            break
        for res in sphere:
            res_key = make_res_key(res)
            if res_key in ligand_res_found and not Polypeptide.is_aa(res):
                ligand_res_found[res_key]["sphere"] = i
                oligomer = ligand_res_found[res_key]["oligomer"]
                oligomer_found[oligomer][res_key] = True
    for chain in model:
        for res in chain.get_unpacked_list():
            res_key = make_res_key(res)
            if (
                res_key in ligand_res_found and 
                not Polypeptide.is_aa(res)
            ):
                oligomer = ligand_res_found[res_key]["oligomer"]
                found_sphere = ligand_res_found[res_key]["sphere"]
                if found_sphere < 0 and any(oligomer_found[oligomer].values()):
                    if include_ligands == 0:
                        spheres[0].add(res)
                    else:
                        spheres[-1].add(res)
                    residues.add(res)
                    print(f"To avoid unpredictable charge error, {res_key} in {oligomer} is added to spheres")


def extract_clusters(
    path,
    out,
    center_residues,
    sphere_count=2,
    first_sphere_radius=4.0,
    max_atom_count=None,
    merge_cutoff=0.0,
    smooth_method="box_plot",
    ligands=[],
    capping=1,
    charge=True,
    ligand_charge=dict(),
    count=True,
    xyz=True,
    hetero_pdb=False,
    include_ligands=2,
    **smooth_params
):
    """
    Extract active site coordination spheres. Neighboring residues determined by
    Voronoi tessellation.

    Parameters
    ----------
    path: str
        Path to PDB file
    out: str
        Path to output directory
    center_residues: list
        List of resnames of the residues to use as the cluster center
    sphere_count: int
        Number of coordinations spheres to extract
    first_sphere_radius: float
        the radius cutoff of the first sphere
    max_atom_count: int
        the maximum number of atoms in the whole cluster
    merge_cutoff: int
        the distance cutoff when merging two centers of spheres
    smooth_method: ("boxplot" | "dbscan" | "dummy_atom")
        The method used to smoothen the spheres
    ligands: list
        Other ligand IDs to include, in addition to AAs and waters
    capping: int
        Whether to cap chains with nothing (0), H (1), or ACE/NME (2)
    charge: bool
        If true, total charge of coordinating AAs will be written to out/charge.csv
    ligand_charge: dict
        The dict containing each ligand's name and charge
    count: bool
        If true, residue counts will be written to out/count.csv
    xyz: bool
        If true, XYZ files will be written according to the output PDBs
    hetero_pdb: bool
        If true, keep all the heteroatoms in the cluster PDB output
    include_ligands: int
        the mode of including ligands in the sphere
    smooth_params:
        params of the specific smooth method 
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", path)
    io = PDBIO()
    io.set_structure(structure)

    model = structure[0]
    neighbors = voronoi(model, center_residues, ligands, smooth_method, **smooth_params)

    centers = get_center_residues(model, center_residues, merge_cutoff)

    aa_charge = {}
    res_count = {}
    cluster_paths = []
    for c in centers:
        metal_id, residues, spheres = get_next_neighbors(
            c, neighbors, sphere_count, ligands, first_sphere_radius, smooth_method, include_ligands, **smooth_params
        )
        complete_oligomer(ligand_charge, model, residues, spheres, include_ligands)
        cluster_path = f"{out}/{metal_id}"
        cluster_paths.append(cluster_path)
        os.makedirs(cluster_path, exist_ok=True)

        if max_atom_count is not None:
            prune_atoms(c, residues, spheres, max_atom_count, ligands)
        if charge:
            aa_charge[metal_id] = compute_charge(spheres, structure, ligand_charge)
        if count:
            res_count[metal_id] = count_residues(spheres)
        if capping:
            cap_residues = cap_chains(model, residues, capping)
            if capping == 2:
                spheres[-1] |= cap_residues

        sphere_paths = []
        for i, s in enumerate(spheres):
            sphere_path = f"{cluster_path}/{i}.pdb"
            sphere_paths.append(sphere_path)
            write_pdb(io, s, sphere_path)
        if capping:
            for cap in cap_residues:
                cap.get_parent().detach_child(cap.get_id())
        if xyz:
            struct_to_file.to_xyz(f"{cluster_path}/{metal_id}.xyz", *sphere_paths)
            struct_to_file.combine_pdbs(f"{cluster_path}/{metal_id}.pdb", center_residues, *sphere_paths, hetero_pdb=hetero_pdb)

    if charge:
        with open(f"{out}/charge.csv", "w") as f:
            f.write(f"Name,{','.join(str(i + 1) for i in range(sphere_count))}\n")
            for k, v in sorted(aa_charge.items()):
                f.write(k)
                for s in v:
                    f.write(f",{s}")
                f.write(f"\n")

    if count:
        with open(f"{out}/count.csv", "w") as f:
            f.write(f"Name,{','.join(str(i + 1) for i in range(sphere_count))}\n")
            for k, v in sorted(res_count.items()):
                f.write(k)
                for sphere in v:
                    s = ", ".join(f"{r} {c}" for r, c in sorted(sphere.items()))
                    f.write(f',"{s}"')
                f.write("\n")

    return cluster_paths