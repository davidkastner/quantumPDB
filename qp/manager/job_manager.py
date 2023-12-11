"""Submits DFT single point calcultions."""

import os
import glob
import shutil
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

def get_charge():
    """Extract the charge values from charge.csv"""
    current_dir = os.getcwd()
    
    # Calculate relative paths
    charge_dir = os.path.abspath(os.path.join(current_dir, "../../../"))
    structure_dir = os.path.abspath(os.path.join(current_dir, "../../"))
    charge_csv_path = os.path.join(charge_dir, "charge.csv")
    metal_path = os.path.join(structure_dir, "0.pdb")
    first_sphere_path = os.path.join(structure_dir, "1.pdb")
    fe_dir = current_dir.split("/")[-1] # Determine if it is Fe2 or Fe3
    
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
    
    # Add the Fe directory value
    charge += int(fe_dir[-1])
    
    return charge

def write_scripts(coordinate_file,
                  basis, 
                  method, 
                  charge, 
                  multiplicity,
                  guess, 
                  constraint_freeze, 
                  job_name,
                  gpus,
                  memory):
    """Generate TeraChem submission scripts."""

    qmscript_content = f"""levelshift yes
levelshiftvala 1.0
levelshiftvalb 0.1
new_minimizer yes
run minimize
coordinates {coordinate_file}
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

    # jobscript.sh
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

    return jobscript_content, qmscript_content


def submit_jobs(job_count, basis, method, guess, ions, constraint_freeze, gpus, memory):
    """Generate and submit jobs to queueing system."""

    pdb_dirs = sorted([d for d in os.listdir() if os.path.isdir(d) and not d == 'Protoss'])

    for pdb in pdb_dirs:
        structure_dirs = sorted(glob.glob(os.path.join(pdb, '[A-Z][0-9]*')))

        for structure in structure_dirs:
            qm_path = os.path.join(structure, method)
            os.makedirs(qm_path, exist_ok=True)

            for ion in ions:
                ion_path = os.path.join(qm_path, ion)
                os.makedirs(ion_path, exist_ok=True)
                if os.path.exists(os.path.join(ion_path, 'qmscript.out')):
                    continue
                
                # Remove the log files
                for file in glob.glob(os.path.join(ion_path, 'Z*')):
                    os.remove(file)

                # Move the xyz file into the QM job directory
                xyz_files = glob.glob(os.path.join(structure, '*.xyz'))
                if len(xyz_files) != 1:
                    print(f"Error: Expected 1 xyz file in {structure}, found {len(xyz_files)}")
                    continue
                coord_file = os.path.join(ion_path, os.path.basename(xyz_files[0]))
                shutil.copy(xyz_files[0], coord_file)

                # Get total charge
                current_directory = os.getcwd()
                os.chdir(ion_path)
                charge = get_charge()

                if constraint_freeze:
                    heavy_list = find_heavy()
                    constraint_freeze = f"$constraint_freeze\n{heavy_list}\n$end"

                os.chdir(current_directory)
                
                if method[0] == "u":
                    multiplicity = 5 if ion == 'Fe2' else 6
                else:
                    multiplicity = 1

                coord_name = os.path.basename(coord_file)
                pdb_name = os.path.basename(pdb)
                structure_name = os.path.basename(structure)
                job_name = f"{ion}{pdb_name}{structure_name}"
                jobscript, qmscript = write_scripts(coord_name,
                                                    basis, 
                                                    method, 
                                                    charge, 
                                                    multiplicity,
                                                    guess, 
                                                    constraint_freeze, 
                                                    job_name,
                                                    gpus,
                                                    memory)
  
                with open(os.path.join(ion_path, 'qmscript.in'), 'w') as f:
                    f.write(qmscript)
                with open(os.path.join(ion_path, 'jobscript.sh'), 'w') as f:
                    f.write(jobscript)

                # Submit the job
                os.system(f'cd {ion_path} && qsub jobscript.sh')

                job_count -= 1
                if job_count <= 0:
                    return


if __name__ == '__main__':

    job_count = int(input("Enter the number of jobs to be submitted simultaneously: "))
    method = input("What functional would you like to use (e.g. uwpbeh)? ").lower()

    if method == "uwpbeh":
        basis = "lacvps_ecp"
        guess = "generate"
        ions = ['Fe2','Fe3']
        constraint_freeze = ""
        gpus = 2
        memory = "16G"
    elif method == "ugfn2xtb":
        basis = "gfn2xtb"
        guess = "hcore"
        ions = ['Fe2','Fe3']
        constraint_freeze = ""
        gpus = 1
        memory = "8G"
    elif method == "gfn2xtb":
        basis = "gfn2xtb"
        guess = "hcore"
        ions = ['Fe2']
        constraint_freeze = True
        gpus = 1
        memory = "8G"

    submit_jobs(job_count, basis, method, guess, ions, constraint_freeze, gpus, memory)
