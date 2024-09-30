"""Fetch PDB and perform quality checks"""

import os
import csv
import sys
import yaml
import requests


def read_config(config_file):
    """Read in the configuration yaml file"""
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)


def parse_input(input, output, center_yaml_residues):
    """Returns system information from the provided input"""

    # If the input was a path or a single PDB, convert it to a list
    if isinstance(input, str):
        input = [input]
    output = os.path.abspath(output)
    center_csv_residues = get_centers(input)
    pdb_all = get_pdbs(input, output)

    # Determine how the user has decided to provide the center residues
    if center_csv_residues:
        use_csv_centers = True
        print("> Using residues from the input csv as centers\n")
    elif center_yaml_residues:
        use_csv_centers = False
        print("> Using residues from the input yaml as centers.\n")
    else:
        print("> No center residues were provided in the input csv or the yaml.\n")
        exit()

    # Determine the center residues for the current PDB
    if use_csv_centers:
        center_residues = center_csv_residues
    else:
        if not center_yaml_residues:
            sys.exit("> No more center residues available.")
        center_residues = center_yaml_residues

    return pdb_all, center_residues


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
        raise ValueError("> ERROR: Could not fetch PDB file")

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as f:
        f.write(r.text)


def get_pdbs(input_path, output_path):
    """
    Parses the input PDBs and returns a list of tuples

    Parameters
    ----------
    input_path: list
        List of PDB codes, paths to PDB files, or path to CSV file
    output_path: str
        Path to output directory

    Returns
    -------
    pdb_all
        List of tuples containing parsed PDB ID and path to PDB file


    Notes
    -----
    Store input PDBs as a tuple of
    Parsed ID (PDB code or filename)
    Path to PDB file (existing or to download)
    """
    pdb_all = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            if ext == ".pdb":
                pdb_all.append((pdb, pdb_id))
            elif ext == ".csv":
                with open(pdb_id, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        pdb = row['pdb_id']
                        pdb_all.append((pdb, os.path.join(output_path, pdb, f"{pdb}.pdb")))
            else:
                with open(pdb_id, "r") as f:
                    pdb_all.extend(
                        [(pdb, os.path.join(output_path, pdb, f"{pdb}.pdb")) for pdb in f.read().splitlines()]
                    )
        else:
            pdb_all.append((pdb_id, os.path.join(output_path, pdb_id, f"{pdb_id}.pdb")))

    return pdb_all


def get_centers(input_path):
    """
    Parses the input centers and returns a list of the center residue PDB IDs.

    Parameters
    ----------
    input_path: list
        List of pdbs or the input csv file.

    Returns
    -------
    centers: list
        List PDB IDs for user-defined center residues.
    """
    
    centers = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            if ext == ".pdb":
                centers.append((pdb, pdb_id))
            elif ext == ".csv":
                input_csv = pdb_id
                with open(input_csv, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    # Check if 'center' column exists
                    if 'center' not in reader.fieldnames:
                        print(f"> WARNING: The 'center' option is not being used in {input_csv}. Returning empty list.")
                        return []
                    # If the 'center' column exists, proceed to collect centers
                    for row in reader:
                        center = row.get('center', None)
                        if center:
                            centers.append(center)
    return centers
