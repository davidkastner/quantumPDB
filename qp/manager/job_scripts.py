"""Stored job submission scripts"""

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
pcm_radii_file /data1/groups/HJKgroup/src/pcm_radii
ml_prop yes
end

{constraint_freeze}
"""

    return qmscript_content


def write_jobscript(job_name, gpus, memory):
    """Generate bash submission scripts with conditional sleep time."""

    jobscript_content = f"""#! /bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:{gpus}
#SBATCH --ntasks-per-node={gpus * 20}

source /etc/profile

#---TC setup---
module load terachem/1.9-2023.11-dev
LD_LIBRARY_PATH=/usr/local/pkg/cuda/cuda-11.8/lib64/stubs:${{LD_LIBRARY_PATH}}

# your command to run terachem
terachem qmscript.in > qmscript.out
"""

    return jobscript_content