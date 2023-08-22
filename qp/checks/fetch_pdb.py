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
