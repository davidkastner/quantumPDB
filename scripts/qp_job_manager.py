import os
import csv
import sys
import time

PDB_PATH = "/home/kastner/projects/quantumPDB/dataset/holo"

def read_charge(pdb_dir):
    charge_file = os.path.join(pdb_dir, 'charge.csv')
    charges = {}
    with open(charge_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        for row in reader:
            # Stop processing at empty line or rows with fewer than 2 columns
            if len(row) < 2:
                break
            charges[row[0]] = row[1:]
    return charges

def check_and_submit_jobs(pdb_dir, max_jobs, current_count):
    charges = read_charge(pdb_dir)
    
    for structure_dir in sorted([d for d in os.listdir(pdb_dir) if os.path.isdir(os.path.join(pdb_dir, d)) and not d == 'Protoss']):
        structure_path = os.path.join(pdb_dir, structure_dir)
        if len([f for f in os.listdir(structure_path) if os.path.isfile(os.path.join(structure_path, f))]) != 4:
            print(f"> Directory {pdb_dir.split('/')[-1]} is in progress with quantumPDB.")
            continue
        
        xyz_file = next((f for f in os.listdir(structure_path) if f.endswith('.xyz')), None)
        if not xyz_file:
            continue
        with open(os.path.join(structure_path, xyz_file), 'r') as file:
            coords = file.read().replace('FE', 'Fe')
        
        for folder in ['Fe2', 'Fe3']:
            folder_path = os.path.join(structure_path, 'QM', folder)
            os.makedirs(folder_path, exist_ok=True)
            with open(os.path.join(folder_path, xyz_file), 'w') as file:
                file.write(coords)

            qmscript_content = f"""coordinates {xyz_file}
basis lacvps_ecp
method ub3lyp
charge {charges[structure_dir][0 if folder == 'Fe2' else 1]}
spinmult {2 if folder == 'Fe2' else 3}
maxit 1000 
dftd d3
scrdir ./scr
pcm cosmo
epsilon 10
scf diis+a
pcm_radii read
pcm_radii_file /home/kastner/reference/pcm_radii
ml_prop yes
end
"""
            with open(os.path.join(folder_path, 'qmscript.in'), 'w') as file:
                file.write(qmscript_content)

        jobscript_content = f"""#! /bin/bash
#$ -N {structure_dir}_{pdb_dir.split('/')[-1]}
#$ -cwd
#$ -l h_rt=300:00:00
#$ -l h_rss=8G
#$ -q (gpusnew|gpus|gpusbig)
#$ -l gpus=1
#$ -pe smp 2
# -fin qmscript.in
# -fin {xyz_file}
# -fout scr/

module load cuda
module load terachem
export OMP_NUM_THREADS=2

terachem qmscript.in > $SGE_O_WORKDIR/qmscript.out
"""
        with open(os.path.join(folder_path, 'jobscript.sh'), 'w') as file:
            file.write(jobscript_content)
            
        # Change to the directory of the job file
        os.chdir(folder_path)

        # Submit the job to the scheduler
        os.system(f"qsub jobscript.sh")

        # Change back to the original directory
        os.chdir(PDB_PATH)

        current_count[0] += 1
        if current_count[0] >= max_jobs:
            return True  # indicate that the limit has been reached

        # Check if there's a stop signal file to stop the script
        if os.path.exists('stop_script.signal'):
            print("Stop signal detected. Exiting script.")
            sys.exit(0)
        
        # Introducing a delay of 2 seconds between submissions
        time.sleep(2)
    return False  # indicate that the limit hasn't been reached yet


def cleanup(pdb_dir):
    for structure_dir in [d for d in os.listdir(pdb_dir) if os.path.isdir(os.path.join(pdb_dir, d)) and not d == 'Protoss']:
        for folder in ['Fe2', 'Fe3']:
            folder_path = os.path.join(pdb_dir, structure_dir, 'QM', folder)
            if os.path.exists(folder_path):
                for file in os.listdir(folder_path):
                    if file.startswith('PDB'):
                        os.remove(os.path.join(folder_path, file))

if __name__ == "__main__":
    # Get the number of jobs to submit from the user
    max_jobs = int(input("Enter the number of jobs to submit: "))

    job_count = [0]  # using a mutable list to pass the count by reference
    limit_reached = False

    for pdb_dir in sorted([d for d in os.listdir(PDB_PATH) if os.path.isdir(os.path.join(PDB_PATH, d))]):
        if limit_reached:
            break
        cleanup(os.path.join(PDB_PATH, pdb_dir))
        limit_reached = check_and_submit_jobs(os.path.join(PDB_PATH, pdb_dir), max_jobs, job_count)