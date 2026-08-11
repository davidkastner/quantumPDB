"""Fetch PDB and perform quality checks"""

import os
import csv
import sys
import yaml
import requests

from qp.structure.mmcif_to_pdb import OversizedStructureError, convert_mmcif_to_pdb


def read_config(config_file):
    """Read and parse a YAML configuration file.

    Parameters
    ----------
    config_file : str
        Path to the configuration YAML file.

    Returns
    -------
    dict
        Parsed configuration data with all pipeline parameters.
    """
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)


def parse_input(input, output, center_yaml_residues):
    """Parse input sources and determine center residues for each structure.

    Processes the user-provided input (PDB codes, file paths, or CSV) and
    determines the center residue definitions from either the CSV ``center``
    column or the YAML ``center_residues`` parameter.

    Parameters
    ----------
    input : str or list
        Input specification: PDB code(s), path to PDB/mmCIF file(s), or path to CSV.
    output : str
        Path to the output directory.
    center_yaml_residues : list
        Center residue definitions from the YAML config file.

    Returns
    -------
    tuple of (list, list)
        ``(pdb_all, center_residues)`` where ``pdb_all`` is a list of
        ``(pdb_id, pdb_path, source_cif)`` tuples and ``center_residues`` is a
        list of center definition strings. ``source_cif`` is a local mmCIF
        path when the user supplied ``.cif``/``.mmcif``, otherwise ``None``.

    Raises
    ------
    SystemExit
        If no center residues are provided in either the CSV or YAML.
    """

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
    Fetch a structure for a PDB code, falling back to mmCIF when needed.

    Tries the classic PDB download first. If RCSB returns 404 (common for
    entries with ``pdb_format_compatible = N``), downloads the mmCIF and
    converts it to ``out``.

    Parameters
    ----------
    pdb: str
        PDB code
    out: str
        Path to output PDB file

    Returns
    -------
    str
        ``"pdb"`` if a classic PDB was downloaded, or ``"mmcif"`` if the
        file was obtained via mmCIF conversion.

    Raises
    ------
    ValueError
        If neither PDB nor mmCIF is available for the ID.
    IOError
        If there is a network or server-side error.
    """
    pdb = pdb.strip()
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    pdb_url = f"https://files.rcsb.org/download/{pdb}.pdb"
    try:
        r = requests.get(pdb_url, timeout=30)
    except requests.exceptions.RequestException as e:
        raise IOError(f"Could not connect to the server. Details: {e}")

    if r.status_code == 200:
        with open(out, "w") as f:
            f.write(r.text)
        return "pdb"
    if r.status_code != 404:
        raise IOError(f"Server returned an error with status code {r.status_code}.")

    # Classic PDB unavailable — try mmCIF (entries marked pdb_format_compatible=N)
    print(f"> PDB unavailable for {pdb}; fetching mmCIF")
    cif_url = f"https://files.rcsb.org/download/{pdb}.cif"
    try:
        r_cif = requests.get(cif_url, timeout=60)
    except requests.exceptions.RequestException as e:
        raise IOError(f"Could not connect to the server. Details: {e}")

    if r_cif.status_code == 404:
        raise ValueError(f"PDB ID '{pdb}' is not valid or does not exist.")
    if r_cif.status_code != 200:
        raise IOError(f"Server returned an error with status code {r_cif.status_code}.")

    cif_path = os.path.join(os.path.dirname(out), f"{pdb}.cif")
    with open(cif_path, "w") as f:
        f.write(r_cif.text)
    print(f"> Converting mmCIF → PDB ({cif_path})")
    convert_mmcif_to_pdb(cif_path, out)
    return "mmcif"


def ensure_structure_pdb(pdb_id, pdb_path, source_cif=None):
    """Ensure ``pdb_path`` exists as a classic PDB, converting mmCIF if needed.

    Parameters
    ----------
    pdb_id : str
        Structure identifier (PDB code or local basename).
    pdb_path : str
        Target classic PDB path used by the rest of the pipeline.
    source_cif : str or None
        Optional local mmCIF path to convert when ``pdb_path`` is missing.

    Returns
    -------
    str
        One of ``"exists"``, ``"converted"``, ``"pdb"``, or ``"mmcif"``.
    """
    pdb_path = os.path.abspath(pdb_path)
    if os.path.isfile(pdb_path):
        return "exists"

    if source_cif:
        print(f"> Converting mmCIF → PDB ({source_cif})")
        convert_mmcif_to_pdb(source_cif, pdb_path)
        return "converted"

    return fetch_pdb(pdb_id, pdb_path)


def get_pdbs(input_path, output_path):
    """
    Parses the input PDBs and returns a list of tuples

    Parameters
    ----------
    input_path: list
        List of PDB codes, paths to PDB/mmCIF files, or path to CSV file
    output_path: str
        Path to output directory

    Returns
    -------
    pdb_all
        List of ``(pdb_id, pdb_path, source_cif)`` tuples. ``source_cif`` is
        set for local ``.cif``/``.mmcif`` inputs; otherwise ``None``.


    Notes
    -----
    Store input PDBs as a tuple of
    Parsed ID (PDB code or filename)
    Path to PDB file (existing or to download/convert)
    Optional local mmCIF source path
    """
    pdb_all = []
    for pdb_id in input_path:
        if os.path.isfile(pdb_id):
            pdb, ext = os.path.splitext(os.path.basename(pdb_id))
            pdb = pdb.replace(".", "_")
            ext_lower = ext.lower()
            if ext_lower == ".pdb":
                pdb_all.append((pdb, pdb_id, None))
            elif ext_lower in (".cif", ".mmcif"):
                target = os.path.join(output_path, pdb, f"{pdb}.pdb")
                pdb_all.append((pdb, target, os.path.abspath(pdb_id)))
            elif ext_lower == ".csv":
                with open(pdb_id, "r") as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        pdb = row['pdb_id']
                        pdb_all.append((pdb, os.path.join(output_path, pdb, f"{pdb}.pdb"), None))
            else:
                with open(pdb_id, "r") as f:
                    pdb_all.extend(
                        [
                            (pdb, os.path.join(output_path, pdb, f"{pdb}.pdb"), None)
                            for pdb in f.read().splitlines()
                            if pdb.strip()
                        ]
                    )
        else:
            pdb_all.append((pdb_id, os.path.join(output_path, pdb_id, f"{pdb_id}.pdb"), None))

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
            if ext == ".csv":
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
