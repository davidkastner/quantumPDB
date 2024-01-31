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
        elif "Unable to find atom" in content:
            return "atom type"
        elif "Job terminated" in content:
            return "killed"
        elif "Job finished" in content:
            return "done"
        
    return "running"


def plot_failures(failure_counts):
    """Create a bar plot for the failure modes."""
    format_plot()

    labels = failure_counts.keys()
    counts = failure_counts.values()

    plt.bar(labels, counts, color="silver")
    plt.xlabel('failure modes', fontsize=12, fontweight='bold')
    plt.ylabel('job count', fontsize=12, fontweight='bold')
    plt.savefig('failure_modes.png', bbox_inches="tight", dpi=600)


def check_all_jobs():
    """Loop over all jobs and check if they failed."""
    
    print(f"   > Checking for failed QM jobs.")
    output_name = "failure_modes.txt"
    failure_counts = {"done": 0, "charge": 0, "memory": 0, "atom type": 0, "unknown": 0}

    with open(output_name, "w") as output_file:
        pdb_directories = sorted(glob.glob("."))
        for pdb in pdb_directories:
            for dir in sorted(glob.glob('[0-9]*')):
                for dirpath, dirnames, filenames in os.walk(dir):
                    for filename in filenames:
                        if filename == "qmscript.out":
                            filepath = os.path.join(dirpath, filename)
                            failure_mode = check_failure_mode(filepath)
                            failure_counts[failure_mode] += 1
                            if failure_mode != "done":
                                output_file.write(f"{filepath} - {failure_mode}\n")

    print(f"   > Saving checkup results in {output_name}\n")
    
    return failure_counts


if __name__ == '__main__':
    failure_counts = check_all_jobs()
    plot_failures(failure_counts)