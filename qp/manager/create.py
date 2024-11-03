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
    """Condenses frozen atoms to make it more readable."""
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
    """Freeze heavy atoms for the geometry optimization."""
    xyz_file = next(f for f in os.listdir() if f.endswith('.xyz'))

    with open(xyz_file, 'r') as f:
        non_hydrogens = [i for i, line in enumerate(f.readlines()[2:], 1) if line.split()[0] != 'H']

    return compress_sequence(non_hydrogens)


def residue_exists(sphere_path, res_name, chain, res_id):
    """Returns True if ligand_name exists in the pdb file."""
    with open(sphere_path, 'r') as f:
        for line in f:
            if res_name == line[17:20].strip() \
                and chain == line[21] and res_id == int(line[22:26]):
                return True
    return False


def ligand_in_spheres(ligand, structure_dir, num_sphere):
    ligand_residues = ligand.split()
    not_found = []
    for ligand_residue in ligand_residues:
        res_name, res_id_full = ligand_residue.split('_')
        chain = res_id_full[0]
        res_id = int(res_id_full[1:])

        for i in range(num_sphere + 1):
            sphere_path = os.path.join(structure_dir, f"{i}.pdb")
            if residue_exists(sphere_path, res_name, chain, res_id):
                break
        else: # if any ligand residue is not found in any sphere
            not_found.append(ligand_residue)
    if not not_found:
        return True
    else:
        if len(ligand_residues) > 1 and len(not_found) < len(ligand_residues):
            print(f"For oligomer ligand {ligand}, {', '.join(not_found)} not found in spheres {structure_dir}. This will cause unpredictable charge error!")
        return False


def get_electronic(pdb_id, pdb_list_path):
    """Query the local pdb master list CSV file."""

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
    """Extract the charge values from charge.csv"""
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
                num_sphere = len(parts) - 1
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


def create_jobs(pdb_list_path, output_dir, optimization, basis, method, guess, use_charge_embedding, charge_embedding_cutoff, gpus, memory, scheduler, pcm_radii_file, dielectric):
    """Generate and submit jobs to queueing system, returning True if any jobs were submitted."""
    
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
            shutil.copy(xyz_files[0], coord_file)
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
                qmscript = job_scripts.write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding)
                jobscript = job_scripts.write_slurm_job(job_name, gpus, memory)
            if scheduler == "sge":
                qmscript = job_scripts.write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding)
                jobscript = job_scripts.write_sge_job(job_name, gpus, memory)
            
            with open('qmscript.in', 'w') as f:
                f.write(qmscript)
            with open('jobscript.sh', 'w') as f:
                f.write(jobscript)

            if use_charge_embedding:
                charge_embedding.get_mm_charges(charge_embedding_cutoff)

            print(f"> Created QM job files for {pdb}/{structure_name}/{method}/")
            
            os.chdir(base_dir)
    
    os.chdir(orig_dir)
    return
            
