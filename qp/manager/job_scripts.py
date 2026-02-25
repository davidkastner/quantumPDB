"""Stored job submission scripts"""

def write_qm(optimization, coord_file, basis, method, total_charge, multiplicity, guess, pcm_radii_file, constraint_freeze, dielectric, use_charge_embedding, use_implicit_solvent=True):
    """Generate a TeraChem input file (qmscript.in).

    Creates the input file content for TeraChem with the specified
    calculation parameters. Supports both single-point energy and
    geometry optimization calculations. PCM implicit solvent and
    point charge embedding can be enabled independently or together.

    Parameters
    ----------
    optimization : bool
        If True, include geometry optimization keywords.
    coord_file : str
        Name of the XYZ coordinate file.
    basis : str
        Basis set name (e.g., ``'lacvps_ecp'``).
    method : str
        DFT functional (e.g., ``'wpbeh'``, ``'ub3lyp'``).
    total_charge : int
        Total system charge.
    multiplicity : int
        Spin multiplicity (1 = singlet, 2 = doublet, etc.).
    guess : str
        Initial guess method (e.g., ``'generate'``).
    pcm_radii_file : str
        Path to PCM radii file for cavity construction.
    constraint_freeze : str
        TeraChem constraint block for frozen atoms.
    dielectric : float
        Dielectric constant for PCM solvent.
    use_charge_embedding : bool
        If True, include MM point charges from ``ptchrges.xyz``.
    use_implicit_solvent : bool, optional
        If True, include PCM implicit solvent (COSMO) block. Can be
        enabled alongside ``use_charge_embedding`` for combined
        QM/MM + implicit solvent calculations. Default is True.

    Returns
    -------
    str
        Complete TeraChem input file content.
    """

    minimization_keywords = """new_minimizer yes\nrun minimize\n""" if optimization else ""

    if use_implicit_solvent:
        pcm_section = f"""pcm cosmo
epsilon {dielectric}
pcm_radii read
pcm_radii_file {pcm_radii_file}
pcm_matrix no
"""
    else:
        pcm_section = ""

    if use_charge_embedding:
        pointcharges_section = """pointcharges ptchrges.xyz
pointcharges_self_interaction true
"""
    else:
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
    """Generate a SLURM submission script for TeraChem.

    Creates a bash script with SLURM directives for GPU job submission.
    Configured for systems with NVIDIA Volta GPUs and the TeraChem module.

    Parameters
    ----------
    job_name : str
        Name for the SLURM job (used in output filenames).
    gpus : int
        Number of GPUs to request.
    memory : str
        Memory allocation (currently unused but kept for API consistency).

    Returns
    -------
    str
        Complete SLURM submission script content.
    """

    jobscript_content = f"""#! /bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=xeon-g6-volta
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:{gpus}
#SBATCH --cpus-per-task={gpus * 20}
#SBATCH --output={job_name}.log

source /etc/profile

#---TC setup---
module load terachem/1.9-2023.11-dev

# your command to run terachem
echo "Run Start Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record
terachem qmscript.in > qmscript.out
echo "Run End Time: $(date '+%Y-%m-%d %H:%M:%S')" >> .submit_record
"""

    return jobscript_content

def write_sge_job(job_name, gpus, memory):
    """Generate a Sun Grid Engine (SGE) submission script for TeraChem.

    Creates a bash script with SGE directives for GPU job submission.
    Includes a minimum runtime enforcement (10 minutes) to prevent
    scheduler issues with very short jobs.

    Parameters
    ----------
    job_name : str
        Name for the SGE job (prefixed with 'Z' in the script).
    gpus : int
        Number of GPUs/parallel threads to request.
    memory : str
        Memory allocation string (e.g., ``'8G'``).

    Returns
    -------
    str
        Complete SGE submission script content.
    """

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