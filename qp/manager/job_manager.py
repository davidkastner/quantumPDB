"""Manager for submitting DFT single point calcultions."""

import os
import sys
import glob
import shutil
import requests
import pandas as pd
from itertools import groupby
from operator import itemgetter
from qp.manager import failure_checkup

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


def get_electronic(pdb_id, master_list):
    """Query the local pdb master list CSV file."""

    try:
        df = pd.read_csv(master_list)
        # Check if the row with the matching pdb_id exists
        row = df[df['pdb_id'] == pdb_id]
        if row.empty:
            print(f"No data found for pdb_id: {pdb_id}")
            return None

        # Check for a value in the 'literature_oxidation' column, if not found, use 'pdb_oxidation'
        if not pd.isna(row['literature_oxidation'].iloc[0]) and row['literature_oxidation'].iloc[0] != '':
            oxidation = int(row['literature_oxidation'].iloc[0])
        else:
            oxidation = int(row['pdb_oxidation'].iloc[0])

        # Determine the multiplicity
        spin = row['spin'].iloc[0]
        if oxidation == 2 and spin == "high":
            multiplicity = 5
        elif oxidation == 2 and spin == "low":
            multiplicity = 1
        elif oxidation == 3 and spin == "high":
            multiplicity = 6
        elif oxidation == 3 and spin == "low":
            multiplicity = 2
        else:
            print(f"\033[93mWARNING: Spin {spin} and oxidation {oxidation} are not supported for {pdb_id}.\033[0m")

        return oxidation, multiplicity

    except Exception as e:
        print(f"Error: {e}")
        sys.exit()


def get_charge(spheres):
    """Extract the charge values from charge.csv"""
    current_dir = os.getcwd()
    
    # Construct relative paths
    charge_dir = os.path.abspath(os.path.join(current_dir, "../../"))
    structure_dir = os.path.abspath(os.path.join(current_dir, "../"))
    charge_csv_path = os.path.join(charge_dir, "charge.csv")
    first_sphere_path = os.path.join(structure_dir, "1.pdb")
    
    charge = 0
    section = 1
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
                current_chain_identifier = parts[0]
                if current_chain_identifier == chain_identifier:
                    charge += sum([int(x) for x in parts[1:spheres + 1]])

            # Parse charge.csv section 2
            elif section == 2 and '_' + chain in line:
                ligand, value = line.split(',')
                if residue_exists(first_sphere_path, ligand.split('_')[0]):
                    charge += int(value)
    
    return charge


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
sleep 180
"""

    return jobscript_content


def get_master_list(url):
    """Retrieve a copy of the master list csv from shared Google Spreadsheet."""

    master_list_file = "master_list.csv"

    response = requests.get(url)
    if response.status_code == 200:
        with open(master_list_file, 'wb') as file:
            file.write(response.content)
        print("> Protein master list downloaded successfully")
        return os.path.abspath(master_list_file)
    
    else:
        print("Failed to retrieve the file.")
        sys.exit(1)


def submit_jobs(job_count, spheres, master_list_path, minimization, basis, method, guess, gpus, memory):
    """Generate and submit jobs to queueing system."""

    base_dir = os.getcwd()  # Store the base directory
    pdb_dirs = sorted([d for d in os.listdir() if os.path.isdir(d) and not d == 'Protoss'])

    for pdb in pdb_dirs:
        pdb_dir_path = os.path.join(base_dir, pdb)
        structure_dirs = sorted(glob.glob(os.path.join(pdb_dir_path, '[A-Z][0-9]*')))

        for structure in structure_dirs:
            qm_path = os.path.join(structure, method)
            os.makedirs(qm_path, exist_ok=True)

            # Check if a prevous job was run and finished correctly
            qm_output_file = os.path.join(qm_path, 'qmscript.out')
            if os.path.exists(qm_output_file):
                job_status = failure_checkup.check_failure_mode(qm_output_file)
                if job_status == "done" or job_status == "running":
                    continue
            
            # Remove old log files
            for file in glob.glob(os.path.join(qm_path, 'Z*')):
                os.remove(file)

            # Copy over the xyz corresponding to the requested cluster model
            xyz_files = glob.glob(os.path.join(structure, '*.xyz'))
            coord_file = os.path.join(qm_path, xyz_files[spheres])
            shutil.copy(xyz_files[spheres], qm_path)

            os.chdir(qm_path)

            oxidation, multiplicity = get_electronic(pdb.lower(), master_list_path)
            charge = get_charge(spheres)
            total_charge = charge + oxidation

            if minimization:
                heavy_list = find_heavy()
                constraint_freeze = f"$constraint_freeze\n{heavy_list}\n$end"
            else:
                constraint_freeze = ""

            coord_name = os.path.basename(coord_file)
            pdb_name = os.path.basename(pdb)
            structure_name = os.path.basename(structure)
            job_name = f"{pdb_name}{structure_name}"

            qmscript = write_qmscript(minimization, coord_name, basis, method, total_charge, multiplicity, guess, constraint_freeze)
            jobscript = write_jobscript(job_name, gpus, memory)

            with open('qmscript.in', 'w') as f:
                f.write(qmscript)
            with open('jobscript.sh', 'w') as f:
                f.write(jobscript)

            os.system('qsub jobscript.sh')

            os.chdir(base_dir)  # Return to base directory

            job_count -= 1
            if job_count <= 0:
                return


if __name__ == '__main__':

    job_count = int(input("Enter the number of jobs to be submitted simultaneously: "))
    method = input("What functional would you like to use (e.g. uwpbeh)? ").lower()
    spheres = int(input("How many spheres would you like to run (e.g. 0-first, 1-second, 2-third)? "))
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

    url = "https://docs.google.com/spreadsheets/d/1St_4YEKcWzrs7yS1GTehfAKabtGCJeuM0sC0u1JG8ZE/gviz/tq?tqx=out:csv&sheet=Sheet1"
    master_list_path = get_master_list(url)
    submit_jobs(job_count, spheres, master_list_path, minimization, basis, method, guess, gpus, memory)
