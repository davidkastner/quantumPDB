"""Fetch PDB and perform quality checks"""

import requests
import os


def fetch_pdb(pdb, out):
    """
    Fetches the PDB file for a given PDB code

    Parameters
    ----------
    pdb: str
        PDB code
    out: str
        Path to output PDB file
    """
    url = f"https://files.rcsb.org/view/{pdb}.pdb"
    r = requests.get(url)
    if r.status_code != 200:
        raise ValueError("Error fetching PDB file")

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as f:
        f.write(r.text)


def parse_input(input_path, output_path):
    """
    Parses the input PDBs and returns a list of tuples

    Parameters
    ----------
    input_path: list
        List of PDB codes or paths to PDB files
    output_path: str
        Path to output directory

    Returns
    -------
    pdb_all
        List of tuples containing parsed PDB ID and path to PDB file
    """

    # Store input PDBs as a tuple of
    #   parsed ID (PDB code or filename)
    #   path to PDB file (existing or to download)
    pdb_all = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            if ext == ".pdb":
                pdb_all.append((pdb, pdb_id))
            else:
                with open(pdb_id, "r") as f:
                    pdb_all.extend(
                        [(pdb, f"{output_path}/{pdb}/{pdb}.pdb") for pdb in f.read().splitlines()]
                    )
        else:
            pdb_all.append((pdb_id, f"{output_path}/{pdb_id}/{pdb_id}.pdb"))
    
    return pdb_all
