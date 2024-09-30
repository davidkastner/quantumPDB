import os
import sys
import glob
import time
import getpass
import datetime
import subprocess


def create_marker(submission_marker_path, job_name, submission_command, submission_output):
    """Creates a comprehensive summary of submission information that doubles as a tracker."""
    with open(submission_marker_path, 'w') as marker:
        marker.write(f"Jobname: {job_name}\n")
        marker.write(f"Author: {getpass.getuser()}\n")
        marker.write(f"Queue Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        marker.write(f"Command: {submission_command}\n")
        marker.write(f"Output: {submission_output}\n")


def prepare_jobs(job_count, method):
    """Prepare the jobs but do not submit them yet."""
    base_dir = os.getcwd()
    prepared_jobs = []
    
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
                print(f"ERROR: Expected 1 xyz file in {structure}, found {len(xyz_files)}")
                continue
            
            prepared_jobs.append((qm_path, pdb, structure))

            # Break early if job_count is reached
            if len(prepared_jobs) == job_count:
                return prepared_jobs
            
    return prepared_jobs


def submit_jobs(prepared_jobs, scheduler):
    """Submit the prepared jobs to the scheduler."""
    base_dir = os.getcwd()
    submitted_jobs = 0

    for qm_path, pdb, structure in prepared_jobs:
        os.chdir(qm_path)
        
        pdb_name = os.path.basename(pdb)
        structure_name = os.path.basename(structure)
        job_name = f"{pdb_name}{structure_name}"
        
        time.sleep(.2)  # Gives user time to abort with ctrl + C

        # Execute the job submission command and capture its output
        if scheduler == "sge":
            submission_command = 'qsub jobscript.sh'
        elif scheduler == "slurm":
            submission_command = 'sbatch jobscript.sh'

        result = subprocess.run(submission_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        submission_output = result.stdout.strip() if result.returncode == 0 else result.stderr.strip()
        print(f"      > {submission_output}")

        submission_marker_path = os.path.join(qm_path, '.submit_record')
        create_marker(submission_marker_path, job_name, submission_command, submission_output)
        
        os.chdir(base_dir)
        submitted_jobs += 1
        
    return submitted_jobs

def manage_jobs(output, job_count, method, scheduler):
    """Main function for managing QM jobs."""
    # Change into the directory of the generated cluster models
    os.chdir(output)
    print(f"> Attempting to submit {job_count} jobs.")
    
    # Prepare the jobs
    prepared_jobs = prepare_jobs(job_count, method)
    
    # Submit the jobs
    if prepared_jobs:
        submitted_jobs = submit_jobs(prepared_jobs, scheduler)
        print(f"> Successfully submitted {submitted_jobs} jobs.")
    else:
        print("> No jobs were prepared.")
    
    # Done, exit the script
    print("Done.\n")
    sys.exit(0)

