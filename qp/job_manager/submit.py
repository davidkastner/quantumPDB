"""Manager for submitting DFT single point calcultions."""

import os
import sys
import glob
import time
import getpass
import datetime
import subprocess

def create_submission_marker(submission_marker_path, job_name, submission_command, submission_output):
    """Creates a comprehensive summary of submission information that doubles as a tracker."""

    with open(submission_marker_path, 'w') as marker:
        marker.write(f"Jobname: {job_name}\n")
        marker.write(f"Author: {getpass.getuser()}\n")
        marker.write(f"Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        marker.write(f"Command: {submission_command}\n")
        marker.write(f"Output: {submission_output}\n")


def prepare_submission(job_count, method, scheduler):
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
                print(f"ERROR: Expected 1 xyz file in {structure}, found {len(xyz_files)}")
                continue
            
            os.chdir(qm_path)
            
            pdb_name = os.path.basename(pdb)
            structure_name = os.path.basename(structure)
            job_name = f"{pdb_name}{structure_name}"
            
            time.sleep(.25) # Gives user time to abort with ctrl + C
            
            # Execute the job submission command and capture its output
            if scheduler == "sge":
                submission_command = 'qsub jobscript.sh'
            if scheduler == "slurm":
                submission_command = 'sbatch jobscript.sh'

            result = subprocess.run(submission_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            submission_output = result.stdout.strip() if result.returncode == 0 else result.stderr.strip()
            print(f"      > {submission_output}")
            create_submission_marker(submission_marker_path, job_name, submission_command, submission_output)

            os.chdir(base_dir)

            submitted_jobs += 1

            
            # Existing condition to break early if job_count is reached
            if submitted_jobs == job_count:
                print(f"> Submitted {job_count} jobs")
                return submitted_jobs
            
    return 0


def count_running_jobs(scheduler):
    """Counts jobs submitted by the user that are currently running or in the queue."""
    try:
        user_name = getpass.getuser()
        if scheduler == "slurm":
            cmd = f"squeue -u {user_name} | grep {user_name} | wc -l"
        if scheduler == "sge":
            cmd = f"qstat -u {user_name} | grep {user_name} | wc -l"
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        job_count = int(result.stdout.strip())
    except subprocess.CalledProcessError:
        job_count = 0
    return job_count


def manage_jobs(output, target_job_count, method, scheduler):
    """Main function for managing QM jobs."""
    # Change into the directory of the generated cluster models
    os.chdir(output)
    
    while True:
        current_job_count = count_running_jobs(scheduler)
        sleep_time_seconds = 300 # seconds
        print(f"> Currently, there are {current_job_count} jobs running or queued.")

        if current_job_count < target_job_count:
            jobs_needed = target_job_count - current_job_count
            print(f"> Attempting to submit {jobs_needed} jobs to reach the target of {target_job_count}.")
            submitted_jobs = prepare_submission(jobs_needed, method, scheduler)
            
            # Check if any new jobs were submitted
            if submitted_jobs == 0:
                print("Done.")
                sys.exit(0)

            # Wait a moment to let the system update
            print(f"> Sleeping for {int(sleep_time_seconds / 60)} minutes before checking queue")
            time.sleep(sleep_time_seconds)
        
        else:
            # Max requested jobs are already running
            print(f"> Sleeping for {int(sleep_time_seconds / 60)} minutes before checking queue")
            time.sleep(sleep_time_seconds)

