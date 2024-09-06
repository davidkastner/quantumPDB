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


def check_submit_record(submit_record_path):
    """Check the .submit_record file for backlog, queue, running, or done status."""
    with open(submit_record_path, 'r') as f:
        content = f.read()

        queue_time = "Queue Time:" in content
        run_start_time = "Run Start Time:" in content
        run_end_time = "Run End Time:" in content

        if queue_time and not run_start_time:
            return "queue"
        elif run_start_time and not run_end_time:
            return "running"
        elif run_end_time:
            return "done"
        
    return "backlog"


def classify_job(qm_dir_path):
    """Classify the job status based on the presence of .submit_record and qmscript.out."""
    submit_record_path = os.path.join(qm_dir_path, ".submit_record")
    qmscript_path = os.path.join(qm_dir_path, "qmscript.out")

    # Check if there's no .submit_record -> backlog
    if not os.path.exists(submit_record_path):
        return "backlog"

    # Use the .submit_record file to classify queue, running, or done
    submit_status = check_submit_record(submit_record_path)

    # If it's classified as done, check for failure modes
    if submit_status == "done" and os.path.exists(qmscript_path):
        return check_failure_mode(qmscript_path)

    return submit_status


def plot_failures(failure_counts):
    """Create a bar plot for the failure modes in a specific order."""
    format_plot()

    # Ensure that the statuses are ordered as desired
    ordered_labels = ["done", "backlog", "queue", "running", "charge", "memory", "unknown"]
    counts = [failure_counts[status] for status in ordered_labels]

    plt.figure(figsize=(7, 4))
    plt.bar(ordered_labels, counts, color="silver")
    plt.xlabel('job status', fontsize=10, fontweight='bold')
    plt.ylabel('job count', fontsize=10, fontweight='bold')
    plt.savefig('job_status.png', bbox_inches="tight", dpi=600)


def check_all_jobs(qm_job_dir, output):
    """Loop over all jobs and check if they failed or are still queued."""
    
    print(f"> Checking for failed QM jobs in the {output} directory.")
    output_name = "failure_modes.txt"
    failure_counts = {"done": 0, "backlog": 0, "queue": 0, "running": 0, 
                      "charge": 0, "memory": 0, "unknown": 0}

    with open(output_name, "w") as output_file:
        base_dir = os.getcwd()
        os.chdir(output)

        all_pdb_dirs = sorted(glob.glob('[0-9]*'))
        for pdb_dir in all_pdb_dirs:  # Loop over PDB directories
            for chain_dir in os.listdir(pdb_dir):  # Loop over chain subdirectories
                if chain_dir == "Protoss":
                    continue
                chain_dir_path = os.path.join(pdb_dir, chain_dir)

                if os.path.isdir(chain_dir_path):
                    qm_dir_path = os.path.join(chain_dir_path, qm_job_dir)

                    if os.path.exists(qm_dir_path):
                        job_status = classify_job(qm_dir_path)
                        failure_counts[job_status] += 1
                        if job_status not in ["done", "running", "queue"]:
                            output_file.write(f"{chain_dir_path} - {job_status}\n")
                    else:
                        # No QM directory has been generated yet
                        failure_counts["queue"] += 1

    os.chdir(base_dir)
    print(f"> Saving checkup results in {output_name}\n")
    
    return failure_counts
