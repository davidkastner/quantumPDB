"""Temporary mmCIF → classic PDB converter.

This is a stopgap so QuantumPDB's PDB-centric pipeline can accept local
``.cif``/``.mmcif`` inputs and RCSB entries that ship mmCIF only
(``pdb_format_compatible = N``).

Known limits (classic PDB format):
- At most 99,999 atoms and 62 chains
- Residue names truncated/remapped to 3 characters
- Chain IDs remapped to a single character when needed
- B-factors capped at 999.99

If center residues use a remapped ligand/chain name, consult the sidecar
``{stem}_mmcif_remap.json`` written next to the output PDB.
"""

from __future__ import annotations

import json
import os
import string
import warnings
from typing import Dict, List

from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.PDBExceptions import PDBIOException

MAX_PDB_ATOMS = 99999
MAX_PDB_CHAINS = 62
MAX_BFACTOR = 999.99
CHAIN_ALPHABET = string.ascii_uppercase + string.ascii_lowercase + string.digits
RESNAME_ALPHABET = string.ascii_uppercase + string.digits


class OversizedStructureError(ValueError):
    """Raised when a structure cannot fit classic PDB size limits.

    Batch runs should catch this, warn, and skip the entry rather than abort.
    """


def convert_mmcif_to_pdb(cif_path: str, pdb_path: str) -> dict:
    """Convert an mmCIF file to a classic PDB file with format sanitization.

    Parameters
    ----------
    cif_path : str
        Path to the input ``.cif`` / ``.mmcif`` file.
    pdb_path : str
        Path to the output ``.pdb`` file.

    Returns
    -------
    dict
        Remap metadata written to the sidecar JSON (chain/resname maps,
        warnings, atom count).

    Raises
    ------
    OversizedStructureError
        If the structure exceeds classic PDB atom/chain limits after sanitization.
    ValueError
        If the written PDB fails smoke re-parse.
    PDBIOException
        If Bio.PDB cannot write the structure.
    """
    cif_path = os.path.abspath(cif_path)
    pdb_path = os.path.abspath(pdb_path)
    if not os.path.isfile(cif_path):
        raise FileNotFoundError(f"mmCIF file not found: {cif_path}")

    parser = MMCIFParser(QUIET=True)
    structure_id = os.path.splitext(os.path.basename(pdb_path))[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(structure_id, cif_path)

    remap_info = _sanitize_structure(structure)
    n_atoms = sum(1 for _ in structure.get_atoms())
    if n_atoms > MAX_PDB_ATOMS:
        raise OversizedStructureError(
            f"Structure has {n_atoms} atoms; classic PDB supports at most "
            f"{MAX_PDB_ATOMS}. Native mmCIF support is required."
        )

    n_chains = _count_chains(structure)
    if n_chains > MAX_PDB_CHAINS:
        raise OversizedStructureError(
            f"Structure has {n_chains} chains; classic PDB supports at most "
            f"{MAX_PDB_CHAINS}. Native mmCIF support is required."
        )

    os.makedirs(os.path.dirname(pdb_path) or ".", exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    try:
        io.save(pdb_path)
    except PDBIOException as exc:
        raise PDBIOException(
            f"Failed to write classic PDB for {cif_path}: {exc}"
        ) from exc

    _smoke_reparse(pdb_path)

    remap_info["cif_path"] = cif_path
    remap_info["pdb_path"] = pdb_path
    remap_info["n_atoms"] = n_atoms
    remap_info["n_chains"] = n_chains
    _write_remap_sidecar(pdb_path, remap_info)
    return remap_info


def _count_chains(structure) -> int:
    count = 0
    for model in structure:
        count = max(count, sum(1 for _ in model.get_chains()))
    return count


def _sanitize_structure(structure) -> dict:
    """Remap multi-char chains / long resnames and clamp B-factors in place."""
    warnings_list: List[str] = []
    chain_map: Dict[str, str] = {}
    resname_map: Dict[str, str] = {}

    used_chains = set()
    for model in structure:
        for chain in model:
            if len(chain.id) == 1:
                used_chains.add(chain.id)

    for model in structure:
        for chain in list(model.get_chains()):
            if len(chain.id) <= 1:
                continue
            old_id = chain.id
            if old_id in chain_map:
                new_id = chain_map[old_id]
            else:
                new_id = _allocate_chain_id(used_chains)
                chain_map[old_id] = new_id
                used_chains.add(new_id)
                warnings_list.append(f"Remapped chain {old_id!r} → {new_id!r}")
            _rename_chain(chain, new_id)

    used_resnames = {
        residue.resname for residue in structure.get_residues() if len(residue.resname) <= 3
    }

    for residue in list(structure.get_residues()):
        old_name = residue.resname
        if len(old_name) <= 3:
            continue
        if old_name in resname_map:
            new_name = resname_map[old_name]
        else:
            new_name = _allocate_resname(old_name, used_resnames)
            resname_map[old_name] = new_name
            used_resnames.add(new_name)
            warnings_list.append(f"Remapped residue name {old_name!r} → {new_name!r}")
        _rename_residue(residue, new_name)

    for atom in structure.get_atoms():
        bfactor = atom.get_bfactor()
        if bfactor is not None and bfactor > MAX_BFACTOR:
            warnings_list.append(
                f"Clamped B-factor {bfactor:.2f} → {MAX_BFACTOR} for atom {atom.full_id}"
            )
            atom.set_bfactor(MAX_BFACTOR)

    return {
        "chain_map": chain_map,
        "resname_map": resname_map,
        "warnings": warnings_list,
    }


def _allocate_chain_id(used: set) -> str:
    for candidate in CHAIN_ALPHABET:
        if candidate not in used:
            return candidate
    raise OversizedStructureError(
        f"Cannot allocate a single-character chain ID; classic PDB supports "
        f"at most {MAX_PDB_CHAINS} chains."
    )


def _allocate_resname(original: str, used: set) -> str:
    candidate = original[:3]
    if candidate not in used:
        return candidate
    for c1 in RESNAME_ALPHABET:
        for c2 in RESNAME_ALPHABET:
            for c3 in RESNAME_ALPHABET:
                name = f"{c1}{c2}{c3}"
                if name not in used:
                    return name
    raise ValueError("Exhausted 3-character residue name space for PDB remapping.")


def _rename_chain(chain, new_id: str) -> None:
    model = chain.get_parent()
    old_id = chain.id
    if old_id == new_id:
        return
    model.detach_child(old_id)
    chain.id = new_id
    model.add(chain)


def _rename_residue(residue, new_name: str) -> None:
    chain = residue.get_parent()
    old_id = residue.id
    hetflag, seq, icode = old_id
    if hetflag.startswith("H_"):
        new_hetflag = f"H_{new_name}"
    else:
        new_hetflag = hetflag
    new_id = (new_hetflag, seq, icode)
    chain.detach_child(old_id)
    residue.id = new_id
    residue.resname = new_name
    chain.add(residue)


def _smoke_reparse(pdb_path: str) -> None:
    parser = PDBParser(QUIET=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure("smoke", pdb_path)
        n_atoms = sum(1 for _ in structure.get_atoms())
        if n_atoms == 0:
            raise ValueError("Re-parsed PDB contains no atoms")
    except Exception as exc:
        raise ValueError(f"Written PDB failed smoke re-parse ({pdb_path}): {exc}") from exc

    # Column sanity: occupancy and B-factor must be separable on ATOM/HETATM lines
    with open(pdb_path, "r") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 66:
                raise ValueError(f"PDB line too short after conversion: {line.rstrip()!r}")
            resname = line[17:20]
            if " " in resname.strip() and len(resname.strip()) > 3:
                raise ValueError(f"Invalid residue name field: {line.rstrip()!r}")
            # Detect glued occupancy/B from oversized resnames (e.g. '1.00113.59')
            occ_b = line[54:66]
            if occ_b.count(".") > 2:
                raise ValueError(
                    f"Occupancy/B-factor columns look corrupted: {line.rstrip()!r}"
                )


def _write_remap_sidecar(pdb_path: str, remap_info: dict) -> str:
    stem, _ = os.path.splitext(pdb_path)
    sidecar = f"{stem}_mmcif_remap.json"
    with open(sidecar, "w") as handle:
        json.dump(remap_info, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return sidecar


def remap_sidecar_path(pdb_path: str) -> str:
    """Return the expected remap sidecar path for a converted PDB."""
    stem, _ = os.path.splitext(os.path.abspath(pdb_path))
    return f"{stem}_mmcif_remap.json"


def load_remap_sidecar(pdb_path: str) -> dict:
    """Load mmCIF remap metadata next to a PDB, if present.

    Parameters
    ----------
    pdb_path : str
        Path to the classic PDB (sidecar is ``{stem}_mmcif_remap.json``).

    Returns
    -------
    dict
        Remap dict with at least ``resname_map`` and ``chain_map`` keys.
        Missing sidecars yield empty maps.
    """
    sidecar = remap_sidecar_path(pdb_path)
    if not os.path.isfile(sidecar):
        return {"resname_map": {}, "chain_map": {}, "warnings": []}
    with open(sidecar, "r") as handle:
        data = json.load(handle)
    return {
        "resname_map": dict(data.get("resname_map") or {}),
        "chain_map": dict(data.get("chain_map") or {}),
        "warnings": list(data.get("warnings") or []),
    }


def expand_resnames_for_matching(names, resname_map=None) -> set:
    """Expand center resname tokens with mmCIF→PDB remap aliases.

    Parameters
    ----------
    names : iterable of str
        Center residue name tokens (fuzzy mode).
    resname_map : dict, optional
        ``{original_name: mapped_name}`` from a remap sidecar.

    Returns
    -------
    set of str
        Names that should be compared against ``Residue.get_resname()``.
    """
    resname_map = resname_map or {}
    expanded = set()
    for name in names:
        expanded.add(name)
        mapped = resname_map.get(name)
        if mapped:
            expanded.add(mapped)
    return expanded


def normalize_center_key(key: str, resname_map=None) -> str:
    """Rewrite a strict ``RESNAME_CHAINSEQ`` center key through ``resname_map``.

    Parameters
    ----------
    key : str
        Strict center token such as ``A1E3R_A302``.
    resname_map : dict, optional
        ``{original_name: mapped_name}`` from a remap sidecar.

    Returns
    -------
    str
        Key using the mapped residue name when available.
    """
    resname_map = resname_map or {}
    if "_" not in key:
        return key
    resname, rest = key.split("_", 1)
    mapped = resname_map.get(resname)
    if mapped:
        resname = mapped
    return f"{resname}_{rest}"
