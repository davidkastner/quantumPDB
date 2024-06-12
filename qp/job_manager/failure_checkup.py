"""Troubleshooting utility that checks for common QM failure modes."""

import os
import glob
import matplotlib.pyplot as plt


def format_plot() -> None:
    """General plotting parameters for the Kulik Lab."""

    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams["xtick.major.pad"] = 5
    plt.rcParams["ytick.major.pad"] = 5
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["svg.fonttype"] = "none"


def check_failure_mode(filepath):
    """Checks for specific failure mode keywords in generated output."""
    with open(filepath, 'r') as f:
        content = f.read()
        
        if "Incorrect molecular charge or spin multiplicity" in content:
            return "charge"
        elif "In Alloc2D: malloc failed" in content:
            return "memory"
        elif "Job terminated" in content:
            return "unknown"
        elif "Job finished" in content:
            return "done"
        
    return "running"


def plot_failures(failure_counts):
    """Create a bar plot for the failure modes."""
    format_plot()

    labels = failure_counts.keys()
    counts = failure_counts.values()

    plt.bar(labels, counts, color="silver")
    plt.xlabel('job status', fontsize=12, fontweight='bold')
    plt.ylabel('job count', fontsize=12, fontweight='bold')
    plt.savefig('job_status.png', bbox_inches="tight", dpi=600)


def check_all_jobs(qm_job_dir):
    """Loop over all jobs and check if they failed or are still queued."""
    
    print(f"   > Checking for failed QM jobs in the {qm_job_dir} directory.")
    output_name = "failure_modes.txt"
    failure_counts = {"done": 0, "charge": 0, "memory": 0, "unknown": 0, "running": 0, "queue": 0}

    with open(output_name, "w") as output_file:
        for pdb_dir in sorted(glob.glob('[0-9]*')):  # Loop over PDB directories
            for chain_dir in os.listdir(pdb_dir):  # Loop over chain subdirectories
                if chain_dir == "Protoss":
                    continue
                chain_dir_path = os.path.join(pdb_dir, chain_dir)

                if os.path.isdir(chain_dir_path) and chain_dir_path != "Protoss":
                    # Check each chain sub-subdirectory (e.g., A208, C208)
                    if qm_job_dir in os.listdir(chain_dir_path):
                        qm_dir_path = os.path.join(chain_dir_path, qm_job_dir)
                        qmscript_path = os.path.join(qm_dir_path, "qmscript.out")

                        if os.path.exists(qmscript_path):
                            failure_mode = check_failure_mode(qmscript_path)
                            failure_counts[failure_mode] += 1
                            if failure_mode not in ["done", "running"]:
                                output_file.write(f"{chain_dir_path} - {failure_mode}\n")
                    else:
                        # No QM directory has been generated yet
                        failure_counts["queue"] += 1

    print(f"   > Saving checkup results in {output_name}\n")
    
    return failure_counts


if __name__ == '__main__':
    qm_job_dir = input("What is the name of your QM job directory? ")
    failure_counts = check_all_jobs(qm_job_dir)
    plot_failures(failure_counts)