import os
import csv
import glob
import shutil
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

def write_author_credit_csv(author_counts):
    """Write author job count CSV."""
    with open(os.path.join("checkup", "author_credit.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["author", "job_count"])
        for author, count in author_counts.items():
            writer.writerow([author, count])


def plot_failures(failure_counts):
    format_plot()
    ordered_labels = ["done", "backlog", "queue", "running", "charge", "memory", "unknown"]
    counts = [failure_counts[status] for status in ordered_labels]

    plt.figure(figsize=(7, 4))
    plt.bar(ordered_labels, counts, color="silver")
    plt.xlabel('job status', fontsize=10, fontweight='bold')
    plt.ylabel('job count', fontsize=10, fontweight='bold')
    plt.savefig(os.path.join("checkup", 'job_status.png'), bbox_inches="tight", dpi=600)

def plot_authors(author_counts):
    format_plot()
    authors = list(author_counts.keys())
    counts = list(author_counts.values())

    plt.figure(figsize=(7, 4))
    plt.bar(authors, counts, color="silver")
    plt.xlabel('authors', fontsize=10, fontweight='bold')
    plt.ylabel('job count', fontsize=10, fontweight='bold')
    plt.savefig(os.path.join("checkup", 'author_credit.png'), bbox_inches="tight", dpi=600)

def plot_failure_modes_from_csv(csv_path):
    """Plot counts of each failure mode from failure_modes.csv."""
    format_plot()
    failure_mode_counts = defaultdict(int)

    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            failure_mode_counts[row['error']] += 1

    labels = list(failure_mode_counts.keys())
    counts = [failure_mode_counts[mode] for mode in labels]

    plt.figure(figsize=(7, 4))
    plt.bar(labels, counts, color="silver")
    plt.xlabel('failure mode', fontsize=10, fontweight='bold')
    plt.ylabel('count', fontsize=10, fontweight='bold')
    plt.savefig(os.path.join("checkup", 'failure_modes.png'), bbox_inches="tight", dpi=600)

def check_all_jobs(method, output, delete_queued):
    print(f"> Checking for failed QM jobs in the {output} directory.")

    # Ensure checkup directory exists
    checkup_dir = "checkup"
    if os.path.exists(checkup_dir):
        shutil.rmtree(checkup_dir)
    os.makedirs(checkup_dir)

    failure_counts = {"done": 0, "backlog": 0, "queue": 0, "running": 0,
                      "charge": 0, "memory": 0, "unknown": 0}
    author_counts = defaultdict(int)
    job_status_rows = []  # For job_status.csv

    with open(os.path.join(checkup_dir, "failure_modes.csv"), "w", newline='') as output_file:
        writer = csv.writer(output_file)
        writer.writerow(['pdb', 'chain', 'error'])

        base_dir = os.getcwd()
        os.chdir(output)

        all_pdb_dirs = sorted(glob.glob('[0-9]*'))
        for pdb_dir in all_pdb_dirs:
            for chain_dir in os.listdir(pdb_dir):
                if chain_dir == "Protoss":
                    continue
                chain_dir_path = os.path.join(pdb_dir, chain_dir)

                if os.path.isdir(chain_dir_path):
                    qm_dir_path = os.path.join(chain_dir_path, method)
                    submit_record_path = os.path.join(qm_dir_path, ".submit_record")

                    if os.path.exists(submit_record_path):
                        job_status, author = classify_job(qm_dir_path, delete_queued)
                        failure_counts[job_status] += 1
                        author_counts[author] += 1

                        job_status_rows.append([pdb_dir, chain_dir, job_status])
                        if job_status not in ["done", "running", "queue"]:
                            writer.writerow([pdb_dir, chain_dir, job_status])
                    else:
                        failure_counts["backlog"] += 1
                        job_status_rows.append([pdb_dir, chain_dir, "backlog"])

        os.chdir(base_dir)

    # Write job_status.csv
    with open(os.path.join(checkup_dir, "job_status.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["pdb", "chain", "job_status"])
        writer.writerows(job_status_rows)

    print(f"> Saved checkup results to {checkup_dir}/\n")
    return failure_counts, author_counts

