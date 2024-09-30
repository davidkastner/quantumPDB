import os
import glob
import matplotlib.pyplot as plt
from collections import defaultdict

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


def extract_author(content):
    """Extract the author from the .submit_record content."""
    for line in content.splitlines():
        if line.startswith("Author:"):
            return line.split("Author:")[1].strip()
    return "Unknown"


def check_submit_record(submit_record_path, delete_queued):
    """Check the .submit_record file for backlog, queue, running, or done status, and return author."""
    with open(submit_record_path, 'r') as f:
        content = f.read()

        author = extract_author(content)
        queue_time = "Queue Time:" in content
        run_start_time = "Run Start Time:" in content
        run_end_time = "Run End Time:" in content

        # if queue_time and not run_start_time:
        if run_start_time and not run_end_time:
            if delete_queued:
                print(f"Deleting queued job record: {submit_record_path}")
                os.remove(submit_record_path)
            return "queue", author
        elif run_start_time and not run_end_time:
            return "running", author
        elif run_end_time:
            return "done", author
        
    return "backlog", author


def classify_job(qm_dir_path, delete_queued):
    """Classify the job status based on the presence of .submit_record and qmscript.out, and return the author."""
    submit_record_path = os.path.join(qm_dir_path, ".submit_record")
    qmscript_path = os.path.join(qm_dir_path, "qmscript.out")

    # Check if there's no .submit_record -> backlog
    if not os.path.exists(submit_record_path):
        return "backlog", "Unknown"

    # Use the .submit_record file to classify queue, running, or done and get the author
    submit_status, author = check_submit_record(submit_record_path, delete_queued)

    # If it's classified as done, check for failure modes
    if submit_status == "done" and os.path.exists(qmscript_path):
        return check_failure_mode(qmscript_path), author

    return submit_status, author


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


def plot_authors(author_counts):
    """Create a bar plot for author job counts."""
    format_plot()

    authors = list(author_counts.keys())
    counts = list(author_counts.values())

    plt.figure(figsize=(7, 4))
    plt.bar(authors, counts, color="silver")
    plt.xlabel('authors', fontsize=10, fontweight='bold')
    plt.ylabel('job count', fontsize=10, fontweight='bold')
    plt.savefig('author_credit.png', bbox_inches="tight", dpi=600)


def check_all_jobs(method, output, delete_queued):
    """Loop over all jobs and check if they failed or are still queued."""
    
    print(f"> Checking for failed QM jobs in the {output} directory.")
    output_name = "failure_modes.txt"
    failure_counts = {"done": 0, "backlog": 0, "queue": 0, "running": 0, 
                      "charge": 0, "memory": 0, "unknown": 0}
    author_counts = defaultdict(int)  # Track job counts per author

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
                    qm_dir_path = os.path.join(chain_dir_path, method)

                    submit_record_path = os.path.join(qm_dir_path, ".submit_record")

                    # Check if .submit_record exists
                    if os.path.exists(submit_record_path):
                        # Classify the job and count failures or running statuses
                        job_status, author = classify_job(qm_dir_path, delete_queued)
                        failure_counts[job_status] += 1

                        # Count author only if a submit record exists
                        author_counts[author] += 1  # Increment author count
                        
                        if job_status not in ["done", "running", "queue"]:
                            output_file.write(f"{chain_dir_path} - {job_status}\n")
                    else:
                        # If no .submit_record, count it as backlog for job_status
                        failure_counts["backlog"] += 1

    os.chdir(base_dir)
    print(f"> Saving checkup results in {output_name}\n")

    return failure_counts, author_counts


