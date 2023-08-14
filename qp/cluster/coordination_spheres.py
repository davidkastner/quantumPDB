import os
from Bio.PDB import PDBParser, Polypeptide, PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from scipy.spatial import Voronoi


def voronoi(model):
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
        for atom in res.get_unpacked_list(): # includes atoms from multiple conformations
            atoms.append(atom)
            points.append(atom.get_coord())

    vor = Voronoi(points)

    neighbors = {}
    for a, b in vor.ridge_points:
        neighbors.setdefault(atoms[a], []).append(atoms[b])
        neighbors.setdefault(atoms[b], []).append(atoms[a])
    return neighbors


def get_next_neighbors(start, neighbors, limit, ligands):
    """
    Iteratively determines spheres around a given starting atom

    Parameters
    ----------
    start: Bio.PDB.Residue
        Starting metal atom
    neighbors: dict
        Adjacency list of neighboring atoms
    limit: int
        Number of spheres to extract
    
    Returns
    -------
    metal_id: str
        Active site identifier
    seen: set
        Set of residues from all spheres
    spheres: list of sets
        Sets of residues separated by spheres
    """
    seen = {start}
    spheres = [{start}]
    for i in range(limit):
        nxt = set()
        for res in spheres[-1]:
            for atom in res.get_unpacked_list():
                for n in neighbors[atom]:
                    par = n.get_parent()
                    if (par not in seen and (Polypeptide.is_aa(par) or
                        par.get_resname() == "HOH" or 
                        par.get_resname() in ligands)):
                        nxt.add(par)
                        seen.add(par)
        spheres.append(nxt)
    
    res_id = start.get_full_id()
    metal_id = res_id[2] + str(res_id[3][1])
    return metal_id, seen, spheres


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


def construct_hydrogen(chain, parent, template, atom):
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
    if atom == "N":
        pos = scale_hydrogen(parent["N"], template["C"], 1 / 1.32) # TODO verify scaling
    else:
        pos = scale_hydrogen(parent["C"], template["N"], 1.09 / 1.32)

    res_id = ("H_" + parent.get_resname(), template.get_id()[1], " ")
    if not chain.has_id(res_id):
        res = Residue(res_id, parent.get_resname(), " ")
        res.add(Atom("H1", pos, 0, 1, " ", "H1", None, "H"))
        chain.add(res)
    else:
        chain[res_id].add(Atom("H2", pos, 0, 1, " ", "H2", None, "H"))
    return chain[res_id]


def construct_heavy(chain, parent, template, atom):
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
        pos["HH32"] = scale_hydrogen(template["CA"], template["CB"], 1.09 / 1.54) # TODO verify scaling
    
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


def cap_chains(model, residues, capping):
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
    cap_residues = set()

    for res in list(sorted(residues)):
        if not Polypeptide.is_aa(res):
            continue

        res_id = res.get_full_id()
        chain = model[res_id[2]]
        chain_list = chain.get_unpacked_list()
        ind = chain_list.index(res)
        
        pre = chain_list[ind - 1]
        if (pre.get_id()[1] == res_id[3][1] - 1 and pre.get_id()[0] == " " and 
            pre not in residues and Polypeptide.is_aa(pre)):
            if capping == 1:
                cap_residues.add(construct_hydrogen(chain, res, pre, "N"))
            else:
                cap_residues.add(construct_heavy(chain, res, pre, "N"))
        
        nxt = chain_list[ind + 1]
        if (nxt.get_id()[1] == res_id[3][1] + 1 and nxt.get_id()[0] == " " and 
            nxt not in residues and Polypeptide.is_aa(nxt)):
            if capping == 1:
                cap_residues.add(construct_hydrogen(chain, res, nxt, "C"))
            else:
                cap_residues.add(construct_heavy(chain, res, nxt, "C"))

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


def extract_clusters(path, out, metals, limit=2, ligands=[], capping=0):
    """
    Extract active site coordination spheres. Neighboring residues determined by
    Voronoi tessellation.

    Parameters
    ----------
    path: str
        Path to PDB file
    out: str
        Path to output directory
    metals: list
        List of active site metal IDs
    limit: int
        Number of coordinations spheres to extract
    ligands: list
        Other ligand IDs to include, in addition to AAs and waters
    capping: int
        Whether to cap chains with nothing (0), H (1), or ACE/NME (2)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", path)
    io = PDBIO()
    io.set_structure(structure)

    model = structure[0]
    neighbors = voronoi(model)

    for res in model.get_residues():
        if res.get_resname() in metals:
            metal_id, residues, spheres = get_next_neighbors(res, neighbors, limit, ligands)
            os.makedirs(f"{out}/{metal_id}", exist_ok=True)
            if capping:
                cap_residues = cap_chains(model, residues, capping)
                spheres[-1] |= cap_residues
            for i in range(limit + 1):
                write_pdb(io, spheres[i], f"{out}/{metal_id}/{i}.pdb")
            if capping:
                for cap in cap_residues:
                    cap.get_parent().detach_child(cap.get_id())
