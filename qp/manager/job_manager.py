"""Manager for submitting DFT single point calcultions."""

import os
import sys
import glob
import time
import shutil
import getpass
import datetime
import requests
import subprocess
import pandas as pd
from itertools import groupby
from operator import itemgetter
from qp.manager import failure_checkup
from qp.manager import job_scripts

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

def residue_exists(sphere_path, res_name, chain, res_id):
    """Returns True if ligand_name exists in the pdb file."""
    with open(sphere_path, 'r') as f:
        for line in f:
            if line.startswith("HETATM") and res_name == line[17:20].strip() \
                and chain == line[21] and res_id == int(line[22:26]):
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


def get_charge(structure_dir=None):
    """Extract the charge values from charge.csv"""
    if structure_dir is None:
        current_dir = os.getcwd()
        charge_dir = os.path.abspath(os.path.join(current_dir, "../../"))
        structure_dir = os.path.abspath(os.path.join(current_dir, "../"))
    else:
        charge_dir = os.path.abspath(os.path.join(structure_dir, "../"))
    charge_csv_path = os.path.join(charge_dir, "charge.csv")
    
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
                res_name, res_id_full = ligand.split('_')
                chain = res_id_full[0]
                res_id = int(res_id_full[1:])
                for i in range(num_sphere):
                    sphere_path = os.path.join(structure_dir, f"{i + 1}.pdb")
                    if residue_exists(sphere_path, res_name, chain, res_id):
                        charge += int(value)
    
    return charge


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


def create_submission_marker(submission_marker_path, job_name, submission_command, submission_output):
    """Creates a comprehensive summary of submission information that doubles as a tracker."""

    with open(submission_marker_path, 'w') as marker:
        marker.write(f"Jobname: {job_name}\n")
        marker.write(f"Author: {getpass.getuser()}\n")
        marker.write(f"Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        marker.write(f"Command: {submission_command}\n")
        marker.write(f"Output: {submission_output}\n")


def submit_jobs(job_count, master_list_path, minimization, basis, method, guess, gpus, memory):
    """Generate and submit jobs to queueing system, returning True if any jobs were submitted."""
    base_dir = os.getcwd()
    submitted_jobs = 0
    
    pdb_dirs = sorted([d for d in os.listdir() if os.path.isdir(d) and not d == 'Protoss'])
    for pdb in pdb_dirs:
        pdb_dir_path = os.path.join(base_dir, pdb)
        structure_dirs = sorted(glob.glob(os.path.join(pdb_dir_path, '[A-Z][0-9]*')))
        
        if len(structure_dirs) == 0:
            continue

        for structure in structure_dirs:
            qm_path = os.path.join(structure, method)
            os.makedirs(qm_path, exist_ok=True)

            submission_marker_path = os.path.join(qm_path, '.submit_record')
            qm_output_file = os.path.join(qm_path, 'qmscript.out')
            log_files = os.path.join(qm_path, 'Z*')
            
            # Remove previous slurm log files
            for file in glob.glob(log_files):
                os.remove(file)

            # Placeholder for submitted jobs to prevent resubmission
            if os.path.exists(submission_marker_path) or os.path.exists(qm_output_file):
                continue
            
            xyz_files = glob.glob(os.path.join(structure, '*.xyz'))
            if len(xyz_files) != 1:
                print(f"Error: Expected 1 xyz file in {structure}, found {len(xyz_files)}")
                continue
            
            coord_file = os.path.join(qm_path, os.path.basename(xyz_files[0]))
            shutil.copy(xyz_files[0], coord_file)
            os.chdir(qm_path)
            
            oxidation, multiplicity = get_electronic(pdb.lower(), master_list_path)
            charge = get_charge()
            total_charge = charge + oxidation
            
            # Get heavy atoms to fix if geometry optimization was requested
            heavy_list = find_heavy() if minimization else ""
            constraint_freeze = f"$constraint_freeze\n{heavy_list}\n$end" if minimization else ""
            
            coord_name = os.path.basename(coord_file)
            pdb_name = os.path.basename(pdb)
            structure_name = os.path.basename(structure)
            job_name = f"{pdb_name}{structure_name}"
            
            qmscript = job_scripts.write_qmscript(minimization, coord_name, basis, method, total_charge, multiplicity, guess, constraint_freeze)
            jobscript = job_scripts.write_jobscript(job_name, gpus, memory)
            
            with open('qmscript.in', 'w') as f:
                f.write(qmscript)
            with open('jobscript.sh', 'w') as f:
                f.write(jobscript)
            
            time.sleep(.25) # Gives user time to abort with ctrl + C
            
            # Execute the job submission command and capture its output
            submission_command = 'qsub jobscript.sh'
            result = subprocess.run(submission_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            submission_output = result.stdout.strip() if result.returncode == 0 else result.stderr.strip()
            print(f"      > {submission_output}")
            create_submission_marker(submission_marker_path, job_name, submission_command, submission_output)

            os.chdir(base_dir)

            submitted_jobs += 1

            
            # Existing condition to break early if job_count is reached
            if submitted_jobs == job_count:
                print(f"   > Submitted {job_count} jobs")
                return submitted_jobs
            
    return 0


def count_running_jobs():
    """Counts jobs submitted by the user that are currently running or in the queue."""

    try:
        user_name = getpass.getuser()
        cmd = f"qstat -u {user_name} | grep {user_name} | wc -l"
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        job_count = int(result.stdout.strip())
    except subprocess.CalledProcessError:
        job_count = 0
    return job_count


def manage_jobs(target_job_count, master_list_path, minimization, basis, method, guess, gpus, memory, check_interval=300):
    while True:
        current_job_count = count_running_jobs()
        sleep_time_seconds = 300
        print(f"   > Currently, there are {current_job_count} jobs running or queued.")

        if current_job_count < target_job_count:
            jobs_needed = target_job_count - current_job_count
            print(f"   > Attempting to submit {jobs_needed} jobs to reach the target of {target_job_count}.")
            submitted_jobs = submit_jobs(jobs_needed, master_list_path, minimization, basis, method, guess, gpus, memory)
            
            # Check if any new jobs were submitted
            if submitted_jobs == 0:
                print("\033[1;31mComplete.\033[0m\n") # Bold red
                sys.exit(0)

            # Wait a moment to let the system update
            print(f"   > Sleeping for {int(sleep_time_seconds / 60)} minutes before checking queue")
            time.sleep(sleep_time_seconds)
        
        else:
            # Max requested jobs are already running
            print(f"   > Sleeping for {int(sleep_time_seconds / 60)} minutes before checking queue")
            time.sleep(sleep_time_seconds)


if __name__ == '__main__':

    job_count = int(input("Enter the request batch size: "))
    method = input("What functional would you like to use (e.g. uwpbeh)? ").lower()
    minimization = False

    if method == "uwpbeh":
        basis = "lacvps_ecp"
        guess = "generate"
        gpus = 1
        memory = "8G"
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
    manage_jobs(job_count, master_list_path, minimization, basis, method, guess, gpus, memory)
