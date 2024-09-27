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
pcm_radii_file /home/kastner/reference/pcm_radii
ml_prop yes
end

{constraint_freeze}
"""

    return qmscript_content


def write_jobscript(job_name, gpus, memory):
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

terachem qmscript.in > $SGE_O_WORKDIR/qmscript.out

# Calculate elapsed time in seconds
ELAPSED_TIME=$SECONDS

# If elapsed time is less than 600 seconds (10 minutes), sleep for the remainder
if [ $ELAPSED_TIME -lt 600 ]; then
    SLEEP_TIME=$((600 - ELAPSED_TIME))
    sleep $SLEEP_TIME
fi
"""

    return jobscript_content
