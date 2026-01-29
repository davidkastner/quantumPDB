"""Add hydrogens using Protoss

**Usage**

#. Submitting existing or custom PDB file::

    >>> from qp.protonate import get_protoss
    >>> pid = get_protoss.upload("path/to/PDB.pdb")
    >>> job = get_protoss.submit(pid)
    >>> get_protoss.download(job, "path/to/OUT.pdb")

#. Submitting PDB code::

    >>> from qp.protonate import get_protoss
    >>> pdb = "1dry"
    >>> job = get_protoss.submit(pdb)
    >>> get_protoss.download(job, "path/to/OUT.pdb")
    >>> get_protoss.download(job, "path/to/OUT.sdf", "ligands")

Protoss automatically removes alternative conformations and overlapping entries. 
Download the log file (``key="log"`` in ``get_protoss.download``) to see affected atoms. 

Some metal-coordinating residues may be incorrectly protonated. Use 
``get_protoss.adjust_activesites(path, metals)`` with the metal IDs to deprotonate
these residues. 
"""

import os
import json
import time
import requests
from Bio.PDB import PDBParser, PDBIO


def upload(path):
    """Upload a PDB file to the ProteinsPlus web server.

    Parameters
    ----------
    path : str
        Path to the PDB file to upload.

    Returns
    -------
    str
        ProteinsPlus ID for the uploaded structure.

    Raises
    ------
    ValueError
        If the server returns a 400 Bad Request error.
    KeyError
        If the upload fails after 5 retry attempts.
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
            print(f"> KeyError encountered. Retrying in {delay} seconds...")
            time.sleep(delay)

    raise KeyError(f"> Failed to upload the file and retrieve 'id' after {retries} attempts.")


def submit(pid):
    """Submit a PDB code or ProteinsPlus ID to the Protoss web API.

    Parameters
    ----------
    pid : str
        Four-character PDB code or ProteinsPlus ID from :func:`upload`.

    Returns
    -------
    str
        URL of the Protoss job location for status polling.

    Raises
    ------
    ValueError
        If the PDB code is invalid (server returns 400).
    KeyError
        If the submission fails after 5 retry attempts.
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
            print(f"> KeyError encountered. Retrying in {delay} seconds...")
            time.sleep(delay)

    raise KeyError(f"> Failed to submit the PDB code and retrieve 'location' after {retries} attempts.")



def download(job, out, key="protein"):
    """Download a Protoss output file.

    Polls the Protoss job URL until completion, then downloads the
    requested output file.

    Parameters
    ----------
    job : str
        URL of the Protoss job location from :func:`submit`.
    out : str
        Path to the output file (directory created if needed).
    key : str, optional
        File type to download: ``'protein'`` (protonated PDB),
        ``'ligand'`` (ligand SDF), or ``'log'`` (processing log).
        Default is ``'protein'``.

    Raises
    ------
    KeyError
        If the download fails after 5 retry attempts.
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

    raise KeyError(f"> Failed to download the file with key '{key}' after {retries} attempts.")


def repair_ligands(path, orig):
    """Repair ligands that Protoss renamed to ``MOL``.

    Protoss sometimes replaces unrecognized ligand residues with a generic
    ``MOL`` label. This function restores the original residue names and
    structures by matching them back from the pre-Protoss PDB, then
    reassigns hydrogen atoms to the closest heavy atoms.

    Parameters
    ----------
    path : str
        Path to the Protoss output PDB file (modified in place).
    orig : str
        Path to the original (pre-Protoss) PDB file.
    """
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
