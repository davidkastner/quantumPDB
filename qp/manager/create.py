"""Manager for submitting DFT single point calcultions."""

import os
import sys
import glob
import shutil
import pandas as pd
from itertools import groupby
from operator import itemgetter
from qp.manager import job_scripts
from qp.manager import charge_embedding

def compress_sequence(seq):
    """Condense a list of atom indices into range notation for TeraChem constraints.

    Converts a sorted list of integers into a compact string representation
    using ranges (e.g., ``1-10``) and comma-separated values for non-consecutive
    indices.

    Parameters
    ----------
    seq : list of int
        Sorted list of atom indices to freeze.

    Returns
    -------
    str
        TeraChem constraint block with ``xyz`` prefix for each range or group.

    Examples
    --------
    >>> compress_sequence([1, 2, 3, 5, 7, 8, 9])
    'xyz 1-3\\nxyz 5\\nxyz 7-9'
    """
    ranges = []
    buffer = []  # Stores standalone numbers temporarily
    for k, g in groupby(enumerate(seq), lambda ix: ix[0] - ix[1]):
        group = list(map(itemgetter(1), g))
        if len(group) > 1:
            if buffer:
                ranges.append(f"xyz {','.join(map(str, buffer))}")
                buffer.clear()
            ranges.append(f"xyz {group[0]}-{group[-1]}")
        else:
            buffer.append(group[0])
    
    if buffer:  # Append any remaining standalone numbers
        ranges.append(f"xyz {','.join(map(str, buffer))}")

    return '\n'.join(ranges)


def find_heavy():
    """Find heavy (non-hydrogen) atom indices for geometry optimization constraints.

    Reads the XYZ file in the current directory and identifies all atoms
    that are not hydrogen. These atoms are typically frozen during partial
    geometry optimizations to relax only the hydrogen positions.

    Returns
    -------
    str
        TeraChem constraint block specifying which atoms to freeze,
        formatted using :func:`compress_sequence`.

    Raises
    ------
    StopIteration
        If no XYZ file is found in the current directory.
    """
    xyz_file = next(f for f in os.listdir() if f.endswith('.xyz'))

    with open(xyz_file, 'r') as f:
        non_hydrogens = [i for i, line in enumerate(f.readlines()[2:], 1) if line.split()[0] != 'H']

    return compress_sequence(non_hydrogens)


def residue_exists(sphere_path, res_name, chain, res_id):
    """Check whether a specific residue exists in a PDB file.

    Parameters
    ----------
    sphere_path : str
        Path to the sphere PDB file (e.g., ``0.pdb``).
    res_name : str
        Three-letter residue name to search for.
    chain : str
        Single-character chain ID.
    res_id : int
        Residue sequence number.

    Returns
    -------
    bool
        True if a matching residue is found, False otherwise.
    """
    with open(sphere_path, 'r') as f:
        for line in f:
            if res_name == line[17:20].strip() \
                and chain == line[21] and res_id == int(line[22:26]):
                return True
    return False


def ligand_in_spheres(ligand, structure_dir, num_sphere):
    """Check whether a ligand is present in any of the sphere PDB files.

    For oligomeric ligands (multi-residue), all component residues must be
    found. A warning is printed if only some residues of an oligomer are
    present, as this can cause charge errors.

    Parameters
    ----------
    ligand : str
        Ligand key string (space-separated residue keys for oligomers).
    structure_dir : str
        Path to the cluster directory containing numbered sphere PDB files.
    num_sphere : int
        Number of spheres to search (0 through ``num_sphere``).

    Returns
    -------
    bool
        True if the ligand (all residues for oligomers) is found.
    """
    ligand_residues = ligand.split()
    partial_found = False
    for ligand_residue in ligand_residues:
        res_name, res_id_full = ligand_residue.split('_')
        chain = res_id_full[0]
        res_id = int(res_id_full[1:])

        for i in range(num_sphere + 1):
            sphere_path = os.path.join(structure_dir, f"{i}.pdb")
            if residue_exists(sphere_path, res_name, chain, res_id):
                partial_found = True
                break
        else: # if ligand residue is not found in any sphere
            if len(ligand_residues) > 1 and partial_found:
                print(f"For oligomer ligand {ligand}, {ligand_residue} is partially found in spheres. This will cause unpredictable charge error!")
                pass
            return False
    return True


def get_electronic(pdb_id, pdb_list_path):
    """Retrieve oxidation state and spin multiplicity from the input CSV.

    Looks up the specified PDB ID in the user-provided CSV file and
    extracts the ``oxidation`` and ``multiplicity`` columns for
    calculating total charge and spin.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB accession code (lowercase).
    pdb_list_path : str
        Path to the CSV file containing ``pdb_id``, ``oxidation``,
        and ``multiplicity`` columns.

    Returns
    -------
    tuple of (int, int)
        ``(oxidation_state, spin_multiplicity)``. Missing values default
        to 0 for oxidation and 1 for multiplicity.

    Raises
    ------
    SystemExit
        If the CSV cannot be read or parsed.
    """

    try:
        df = pd.read_csv(pdb_list_path)
        # Check if the row with the matching pdb_id exists
        row = df[df['pdb_id'] == pdb_id]
        if row.empty:
            print(f"> No data found for pdb_id: {pdb_id}")
            return None
        # Extract multiplicity and oxidation state directly
        multiplicity = int(row['multiplicity'].fillna(1).iloc[0])
        oxidation = int(row['oxidation'].fillna(0).iloc[0])
        return oxidation, multiplicity

    except Exception as e:
        print(f"> ERROR: {e}")
        sys.exit()


def get_charge(structure_dir=None):
    """Extract total charge and extra spin from charge.csv and spin.csv.

    Parses the charge CSV generated by the cluster extraction stage to
    sum up amino acid charges and ligand charges for residues present
    in the extracted spheres. Also reads spin.csv if present for
    radical ligands (e.g., NO, O2).

    Parameters
    ----------
    structure_dir : str, optional
        Path to the chain directory (e.g., ``output/1os7/A200``).
        If None, inferred from the current working directory.

    Returns
    -------
    tuple of (int, int)
        ``(total_charge, extra_spin)`` where ``total_charge`` is the sum
        of amino acid and ligand charges, and ``extra_spin`` is the
        additional spin contribution from radical species.
    """
    if structure_dir is None:
        current_dir = os.getcwd()
        charge_dir = os.path.abspath(os.path.join(current_dir, "../../"))
        structure_dir = os.path.abspath(os.path.join(current_dir, "../"))
    else:
        charge_dir = os.path.abspath(os.path.join(structure_dir, "../"))
    charge_csv_path = os.path.join(charge_dir, "charge.csv")
    spin_csv_path = os.path.join(charge_dir, "spin.csv")
    
    charge = 0
    section = 1
    num_sphere = 0
    chain_identifier = os.path.basename(structure_dir)
    chain = os.path.basename(structure_dir)[0]

    with open(charge_csv_path, 'r') as charge_csv_content:
        for line in charge_csv_content:
            line = line.strip()
            
            # Check if entered the next section
            if not line:
                section += 1
                continue

            # Parse charge.csv section 1
            if section == 1 and line.startswith(chain_identifier):
                parts = line.split(',')
                num_sphere = len(parts) - 2
                current_chain_identifier = parts[0]
                if current_chain_identifier == chain_identifier:
                    charge += sum([int(x) for x in parts[1:]])

            # Parse charge.csv section 2
            elif section == 2:
                ligand, value = line.split(',')
                if ligand_in_spheres(ligand, structure_dir, num_sphere):
                    charge += int(value)                  

    spin = 0
    if os.path.exists(spin_csv_path):
        with open(spin_csv_path, 'r') as spin_csv_content:
            for line in spin_csv_content:
                line = line.strip()
                if not line:
                    continue
                ligand, value = line.split(',')
                if ligand_in_spheres(ligand, structure_dir, num_sphere):
                    spin += int(value)   
    
    return charge, spin


def create_jobs(pdb_list_path, output_dir, optimization, basis, method, guess, use_charge_embedding, charge_embedding_cutoff, charge_embedding_charges, gpus, memory, scheduler, pcm_radii_file, dielectric, use_implicit_solvent=True):
    """Generate QM job input files for all extracted clusters.

    Creates TeraChem input files (``qmscript.in``) and scheduler submission
    scripts (``jobscript.sh``) for each cluster directory. Automatically
    determines charge and spin multiplicity from the input CSV and
    charge.csv files.

    Parameters
    ----------
    pdb_list_path : str
        Path to the input CSV with ``pdb_id``, ``oxidation``, ``multiplicity``.
    output_dir : str
        Base output directory containing PDB subdirectories.
    optimization : bool
        If True, run geometry optimization; otherwise single-point energy.
    basis : str
        Basis set name (e.g., ``'lacvps_ecp'``).
    method : str
        DFT method/functional (e.g., ``'wpbeh'``).
    guess : str
        Initial wavefunction guess method (e.g., ``'generate'``).
    use_charge_embedding : bool
        If True, generate point charge embedding file.
    charge_embedding_cutoff : float
        Distance cutoff in angstroms for MM point charges.
    charge_embedding_charges : str or None
        Path to a JSON file with custom partial charges. When ``None``,
        the built-in AMBER ff14SB charges are used.
    gpus : int
        Number of GPUs to request.
    memory : str
        Memory allocation string (e.g., ``'8G'``).
    scheduler : str
        Job scheduler type (``'slurm'`` or ``'sge'``).
    pcm_radii_file : str
        Path to PCM radii file for implicit solvent.
    dielectric : float
        Dielectric constant for PCM solvent model.
    use_implicit_solvent : bool, optional
        If True, include PCM implicit solvent in TeraChem input. Can be
        enabled alongside ``use_charge_embedding``. Default is True.
    """

    orig_dir = os.getcwd()
    os.chdir(output_dir)
    base_dir = os.getcwd()
    
    pdb_dirs = sorted([d for d in os.listdir() if os.path.isdir(d) and not d == 'Protoss'])
    for pdb in pdb_dirs:
        pdb_dir_path = os.path.join(base_dir, pdb)
        structure_dirs = sorted(glob.glob(os.path.join(pdb_dir_path, '[A-Z][0-9]*')))
        
        if len(structure_dirs) == 0:
            continue

        for structure in structure_dirs:
            qm_path = os.path.join(structure, method)
            os.makedirs(qm_path, exist_ok=True)
            
            xyz_files = glob.glob(os.path.join(structure, '*.xyz'))
            if len(xyz_files) != 1:
                print(f"ERROR: Expected 1 xyz file in {structure}, found {len(xyz_files)}")
                continue
            
            coord_file = os.path.join(qm_path, os.path.basename(xyz_files[0]))
            shutil.copyfile(xyz_files[0], coord_file)
            os.chdir(qm_path)
            
            oxidation, multiplicity = get_electronic(pdb.lower(), pdb_list_path)
            charge, extra_spin = get_charge()
            multiplicity += extra_spin
            total_charge = charge + oxidation
            
            # Get heavy atoms to fix if geometry optimization was requested
            heavy_list = find_heavy() if optimization else ""
            constraint_freeze = f"$constraint_freeze\n{heavy_list}\n$end" if optimization else ""
            
            coord_file = os.path.basename(coord_file)
            pdb_name = os.path.basename(pdb)
            structure_name = os.path.basename(structure)
            job_name = f"{pdb_name}{structure_name}"
            
            if scheduler == "slurm":
                qmscript = job_scripts.write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding, use_implicit_solvent)
                jobscript = job_scripts.write_slurm_job(job_name, gpus, memory)
            if scheduler == "sge":
                qmscript = job_scripts.write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding, use_implicit_solvent)
                jobscript = job_scripts.write_sge_job(job_name, gpus, memory)
            
            with open('qmscript.in', 'w') as f:
                f.write(qmscript)
            with open('jobscript.sh', 'w') as f:
                f.write(jobscript)

            if use_charge_embedding:
                charge_embedding.get_charges(charge_embedding_cutoff, charge_embedding_charges)

            print(f"> Created QM job files for {pdb}/{structure_name}/{method}/")
            
            os.chdir(base_dir)
    
    os.chdir(orig_dir)
    return
            
