"""Manager for submitting DFT single point calcultions."""

import os
import glob
import shutil
import requests
import pandas as pd
from io import StringIO
from itertools import groupby
from operator import itemgetter

def compress_sequence(seq):
    """Condenses frozen atoms to make it more readable."""
    ranges = []
    buffer = []  # Stores standalone numbers temporarily
    for k, g in groupby(enumerate(seq), lambda ix: ix[0] - ix[1]):
        group = list(map(itemgetter(1), g))
        if len(group) > 1:  # It's a range
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

def residue_exists(first_sphere_path, ligand_name):
    """Returns True if ligand_name exists in the pdb file."""
    with open(first_sphere_path, 'r') as f:
        for line in f:
            if line.startswith("HETATM") and ligand_name in line:
                return True
    return False


def get_oxidation(pdb_id):
    """Query the local pdb master list CSV file."""

    path_to_master_list = "/home/kastner/projects/quantumPDB/protein_master_list.csv"
    try:
        # Read the data into a DataFrame
        df = pd.read_csv(path_to_master_list)

        # Check if the row with the matching pdb_id exists
        row = df[df['pdb_id'] == pdb_id]
        if row.empty:
            print(f"No data found for pdb_id: {pdb_id}")
            return None

        # Check for a value in the 'literature_oxidation' column, if not found, use 'pdb_oxidation'
        if not pd.isna(row['literature_oxidation'].iloc[0]) and row['literature_oxidation'].iloc[0] != '':
            return int(row['literature_oxidation'].iloc[0])
        else:
            return int(row['pdb_oxidation'].iloc[0])
    except Exception as e:
        print(f"Error: {e}")
        return None


def get_charge(oxidation):
    """Extract the charge values from charge.csv"""
    current_dir = os.getcwd()
    
    # Calculate relative paths
    charge_dir = os.path.abspath(os.path.join(current_dir, "../../"))
    structure_dir = os.path.abspath(os.path.join(current_dir, "../"))
    charge_csv_path = os.path.join(charge_dir, "charge.csv")
    first_sphere_path = os.path.join(structure_dir, "1.pdb")
    
    charge = 0
    section = 1
    identifier = os.path.basename(structure_dir)
    chain = os.path.basename(structure_dir)[0]

    with open(charge_csv_path, 'r') as charge_csv_content:
        for line in charge_csv_content:
            line = line.strip()
            
            # Check if entered the next section
            if not line:
                section += 1
                continue

            # Parse charge.csv section 1
            if section == 1 and line.startswith(identifier):
                parts = line.split(',')
                current_identifier = parts[0]
                if current_identifier == identifier:
                    charge += sum([int(x) for x in parts[1:]])

            # Parse charge.csv section 2
            elif section == 2 and '_' + chain in line:
                ligand, value = line.split(',')
                if residue_exists(first_sphere_path, ligand.split('_')[0]):
                    charge += int(value)
    
    # Add the in the metal oxidation state
    charge += oxidation
    
    return charge


def submit_jobs(job_count, minimization, basis, method, guess, gpus, memory):
    """Generate and submit jobs to queueing system."""

    pdb_dirs = sorted([d for d in os.listdir() if os.path.isdir(d) and not d == 'Protoss'])

    for pdb in pdb_dirs:
        structure_dirs = sorted(glob.glob(os.path.join(pdb, '[A-Z][0-9]*')))

        for structure in structure_dirs:
            qm_path = os.path.join(structure, method)
            os.makedirs(qm_path, exist_ok=True)

            # Check to see if a job has already been run here
            if os.path.exists(os.path.join(qm_path, 'qmscript.out')):
                continue
            
            # Remove previous slurm log files
            for file in glob.glob(os.path.join(qm_path, 'Z*')):
                os.remove(file)

            # Move the xyz file into the QM job directory
            xyz_files = glob.glob(os.path.join(structure, '*.xyz'))
            if len(xyz_files) != 1:
                print(f"Error: Expected 1 xyz file in {structure}, found {len(xyz_files)}")
                continue
            coord_file = os.path.join(qm_path, os.path.basename(xyz_files[0]))
            shutil.copy(xyz_files[0], coord_file)

            # Get total charge
            current_directory = os.getcwd()
            os.chdir(qm_path)

            oxidation = get_oxidation(pdb.lower())
            charge = get_charge(oxidation)

            # Get atoms to restrain if requested
            if minimization:
                heavy_list = find_heavy()
                constraint_freeze = f"$constraint_freeze\n{heavy_list}\n$end"
            else:
                constraint_freeze = ""

            # Determine multiplicity
            os.chdir(current_directory)
            if method[0] == "u":
                multiplicity = 5 if oxidation == 2 else 6
            else:
                multiplicity = 1

            coord_name = os.path.basename(coord_file)
            pdb_name = os.path.basename(pdb)
            structure_name = os.path.basename(structure)
            job_name = f"{pdb_name}{structure_name}"
            
            # Generate TeraChem job script contents
            qmscript = write_qmscript(minimization,
                                      coord_name,
                                      basis, 
                                      method, 
                                      charge, 
                                      multiplicity,
                                      guess, 
                                      constraint_freeze,)

            jobscript = write_jobscript(job_name,
                                        gpus,
                                        memory,)

            with open(os.path.join(qm_path, 'qmscript.in'), 'w') as f:
                f.write(qmscript)
            with open(os.path.join(qm_path, 'jobscript.sh'), 'w') as f:
                f.write(jobscript)
            os.system(f'cd {qm_path} && qsub jobscript.sh')

            job_count -= 1
            if job_count <= 0:
                return

def write_qmscript(minimzation,
                  coordinate_file,
                  basis, 
                  method, 
                  charge, 
                  multiplicity,
                  guess, 
                  constraint_freeze):
    """Generate TeraChem job submission scripts."""

    minimization_keywords = """new_minimizer yes\nrun minimize\n""" if minimzation == True else ""
    qmscript_content = f"""levelshift yes
levelshiftvala 0.25
levelshiftvalb 0.25
{minimization_keywords}coordinates {coordinate_file}
basis {basis}
method {method}
charge {charge}
spinmult {multiplicity}
guess {guess}
maxit 500 
dftd d3
scrdir ./scr
pcm cosmo
epsilon 10
pcm_matrix no
scf diis+a
pcm_radii read
pcm_radii_file /home/kastner/reference/pcm_radii
ml_prop yes
end

{constraint_freeze}
"""

    return qmscript_content

def write_jobscript(job_name,
                  gpus,
                  memory):
    """Generate bash submission scripts."""

    jobscript_content = f"""#!/bin/bash
#$ -N Z{job_name}
#$ -cwd
#$ -l h_rt=300:00:00
#$ -l h_rss={memory}
#$ -q (gpusnew|gpus|gpusbig)
#$ -l gpus=1
#$ -pe smp {gpus}
# -fin qmscript.in
# -fin *.xyz
# -fout scr/

module load cuda/10.0
module load terachem/071920-cuda10.0-intel16
module load intel/16.0.109
export OMP_NUM_THREADS={gpus}
terachem qmscript.in > $SGE_O_WORKDIR/qmscript.out
"""

    return jobscript_content

if __name__ == '__main__':

    job_count = int(input("Enter the number of jobs to be submitted simultaneously: "))
    method = input("What functional would you like to use (e.g. uwpbeh)? ").lower()
    minimization = False

    if method == "uwpbeh":
        basis = "lacvps_ecp"
        guess = "generate"
        gpus = 1
        memory = "16G"
    elif method == "ugfn2xtb":
        basis = "gfn2xtb"
        guess = "hcore"
        gpus = 1
        memory = "8G"
    elif method == "gfn2xtb":
        basis = "gfn2xtb"
        guess = "hcore"
        gpus = 1
        memory = "8G"

    submit_jobs(job_count, minimization, basis, method, guess, gpus, memory)
