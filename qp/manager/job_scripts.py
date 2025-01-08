"""Stored job submission scripts"""

def write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding):
    """Generate TeraChem job submission scripts."""
    
    minimization_keywords = """new_minimizer yes\nrun minimize\n""" if optimization else ""

    if use_charge_embedding:
        pcm_section = ""
        pointcharges_section = """pointcharges ptchrges.xyz
pointcharges_self_interaction true
"""
    else:
        pcm_section = f"""pcm cosmo
epsilon {dielectric}
pcm_radii read
pcm_radii_file {pcm_radii_file}
pcm_matrix no
"""
        pointcharges_section = ""


    qmscript_content = f"""levelshift yes
levelshiftvala 0.25
levelshiftvalb 0.25
{minimization_keywords}coordinates {coord_file}
basis {basis}
method {method}
charge {total_charge}
spinmult {multiplicity}
guess {guess}
maxit 500 
dftd d3
scrdir ./scr
scf diis+a
{pcm_section}{pointcharges_section}ml_prop yes
end

{constraint_freeze}
"""
    return qmscript_content



def write_slurm_job(job_name, gpus, memory):
    """Generate bash submission scripts with conditional sleep time and scratch directory usage."""

    jobscript_content = f"""#! /bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=xeon-g6-volta
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:{gpus}
#SBATCH --ntasks-per-node={gpus * 20}
#SBATCH --output={job_name}.log

source /etc/profile

#---TC setup---
module load terachem/1.9-2023.11-dev
export LD_LIBRARY_PATH=/usr/local/pkg/cuda/cuda-11.8/lib64/stubs:${{LD_LIBRARY_PATH}}

# Step 1: Copy input files to the scratch directory
echo "Copying files to local scratch directory ($TMPDIR)"
cp * $TMPDIR/

# Step 2: Change to the scratch directory
echo "Changing to scratch directory ($TMPDIR)"
cd $TMPDIR

# Define a cleanup function to copy output files back to the submission directory
function clean_up {{
  echo "Copying results back to the original directory ($SLURM_SUBMIT_DIR)"
  cp -r * $SLURM_SUBMIT_DIR/
}}

# Set up trap to run cleanup on exit
trap clean_up EXIT

# Step 3: Run TeraChem in scratch directory
echo "Run Start Time: $(date '+%Y-%m-%d %H:%M:%S')" >> $SLURM_SUBMIT_DIR/.submit_record
terachem qmscript.in > qmscript.out
echo "Run End Time: $(date '+%Y-%m-%d %H:%M:%S')" >> $SLURM_SUBMIT_DIR/.submit_record
"""

    return jobscript_content

def write_sge_job(job_name, gpus, memory):
    """Generate bash submission scripts with conditional sleep time."""

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

# Start time
SECONDS=0

echo "Run Start Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record
terachem qmscript.in > $SGE_O_WORKDIR/qmscript.out
echo "Run End Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record

# Calculate elapsed time in seconds
ELAPSED_TIME=$SECONDS

# If elapsed time is less than 600 seconds (10 minutes), sleep for the remainder
# This is just for Gibraltar which can't handle short jobs
if [ $ELAPSED_TIME -lt 600 ]; then
    SLEEP_TIME=$((600 - ELAPSED_TIME))
    sleep $SLEEP_TIME
fi
"""

    return jobscript_content