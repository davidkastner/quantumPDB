import os
import subprocess

def check_files(directory):
    """Check if 0.pdb, 1.pdb, 2.pdb are present in the directory."""
    files = ["0.pdb", "1.pdb", "2.pdb"]
    return all(os.path.isfile(os.path.join(directory, file)) for file in files)

def write_jobscript(job_name, out_dir):
    """Write jobscript.sh file to the given directory."""
    with open(os.path.join(out_dir, "jobscript.sh"), "w") as f:
        f.write(f"""#! /bin/bash
#$ -N PDB_{job_name}
#$ -cwd
#$ -l h_rt=300:00:00
#$ -l h_rss=8G
#$ -q (gpusnew|gpus|gpusbig)
#$ -l gpus=1
#$ -pe smp 2
# -fin qmscript.in
# -fin *.xyz
# -fout scr/

module load cuda
module load terachem
export OMP_NUM_THREADS=2
terachem qmscript.in > $SGE_O_WORKDIR/qmscript.out""")

def write_qmscript(coordinate_file, charge, multiplicity, out_dir):
    """Write qmscript.in file to the given directory."""
    with open(os.path.join(out_dir, "qmscript.in"), "w") as f:
        f.write(f"""coordinates {coordinate_file}
basis lacvps_ecp
method ub3lyp
charge {charge}
spinmult {multiplicity}
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
""")

def main():
    # Use the current working directory as the base
    pdb_dir = "."

    for pdb_folder in sorted(os.listdir(pdb_dir)):
        pdb_path = os.path.join(pdb_dir, pdb_folder)
        
        if not os.path.isdir(pdb_path) or pdb_folder == "Protoss":
            continue  # Skip files or Protoss folder
        
        ready = True
        for struct_folder in os.listdir(pdb_path):
            struct_path = os.path.join(pdb_path, struct_folder)
            if not check_files(struct_path):
                print(f"> Directory {pdb_folder} is missing files")
                ready = False
                break
        
        if not ready:
            continue
        
        qm_dir = os.path.join(pdb_path, "QM")
        if not os.path.exists(qm_dir):
            os.makedirs(qm_dir)

        for folder in ["Fe2", "Fe3"]:
            folder_path = os.path.join(qm_dir, folder)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            
            # Cleanup
            for file in os.listdir(folder_path):
                if file.startswith("PDB"):
                    os.remove(os.path.join(folder_path, file))
            
            # Check for qmscript.out
            if os.path.exists(os.path.join(folder_path, "qmscript.out")):
                continue
            
            # Copy coordinate file
            xyz_files = [f for f in os.listdir(struct_path) if f.endswith(".xyz")]
            if len(xyz_files) != 1:
                continue
            
            coordinate_file = xyz_files[0]
            with open(os.path.join(struct_path, coordinate_file), "r") as f:
                content = f.read().replace("FE", "Fe")
            
            with open(os.path.join(folder_path, coordinate_file), "w") as f:
                f.write(content)
            
            # Get charge
            with open(os.path.join(folder_path, "charge.txt"), "r") as f:
                charge = int(f.read().strip())
            
            # Set multiplicity
            multiplicity = 5 if folder == "Fe2" else 6
            
            # Write submission files
            write_jobscript(pdb_folder, folder_path)
            write_qmscript(coordinate_file, charge, multiplicity, folder_path)
            
            # Submit the job
            subprocess.call(["qsub", "jobscript.sh"], cwd=folder_path)

if __name__ == "__main__":
    main()
