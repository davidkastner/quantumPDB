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
    """Classify a completed QM job based on its output file content.

    Parses the TeraChem output file to identify the job outcome:
    successful completion, charge/spin error, memory error, or unknown failure.

    Parameters
    ----------
    filepath : str
        Path to the ``qmscript.out`` file.

    Returns
    -------
    str
        Status code: ``'done'``, ``'charge'``, ``'memory'``, ``'unknown'``,
        or ``'running'`` (if output exists but has no termination marker).
    """
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
    """Extract the submitting user's name from a submit record.

    Parameters
    ----------
    content : str
        Contents of the ``.submit_record`` file.

    Returns
    -------
    str
        Username of the person who submitted the job, or ``'Unknown'``.
    """
    for line in content.splitlines():
        if line.startswith("Author:"):
            return line.split("Author:")[1].strip()
    return "Unknown"


def check_submit_record(submit_record_path, delete_queued):
    """Determine job status from the submit record file.

    Parses timestamps in the ``.submit_record`` file to determine whether
    the job is queued, running, or completed. Optionally deletes records
    for queued (but never started) jobs to allow resubmission.

    Parameters
    ----------
    submit_record_path : str
        Path to the ``.submit_record`` file.
    delete_queued : bool
        If True, delete records for jobs that were queued but never started.

    Returns
    -------
    tuple of (str, str)
        ``(status, author)`` where status is one of ``'queue'``, ``'running'``,
        ``'done'``, or ``'backlog'``.
    """
    with open(submit_record_path, 'r') as f:
        content = f.read()

    author = extract_author(content)
    queue_time     = "Queue Time:" in content
    run_start_time = "Run Start Time:" in content
    run_end_time   = "Run End Time:" in content

    # Queued but never started
    if queue_time and not run_start_time:
        if delete_queued:
            print(f"Deleting queued job record: {submit_record_path}")
            os.remove(submit_record_path)
        return "queue", author

    # Started but not finished
    if run_start_time and not run_end_time:
        return "running", author

    # Finished
    if run_end_time:
        return "done", author

    # Fallback (record present but missing expected markers)
    return "backlog", author


def classify_job(qm_dir_path, delete_queued):
    """Classify a single QM job's status.

    Combines information from the submit record and output file to
    determine the overall job status. For completed jobs, also checks
    for specific failure modes.

    Parameters
    ----------
    qm_dir_path : str
        Path to the QM calculation directory (e.g., ``output/1os7/A200/wpbeh``).
    delete_queued : bool
        If True, delete records for queued but never-started jobs.

    Returns
    -------
    tuple of (str, str)
        ``(status, author)`` where status is one of ``'backlog'``, ``'queue'``,
        ``'running'``, ``'done'``, ``'charge'``, ``'memory'``, or ``'unknown'``.
    """
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
    """Write a CSV summarizing job counts per submitting author.

    Parameters
    ----------
    author_counts : dict
        Mapping of author names to the number of jobs they submitted.
    """
    with open(os.path.join("checkup", "author_credit.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["author", "job_count"])
        for author, count in author_counts.items():
            writer.writerow([author, count])


def plot_failures(failure_counts):
    """Generate a bar chart of job status counts and save to ``checkup/job_status.png``.

    Parameters
    ----------
    failure_counts : dict
        Mapping of status labels (e.g., ``'done'``, ``'running'``) to counts.
    """
    format_plot()
    ordered_labels = ["done", "backlog", "queue", "running", "charge", "memory", "unknown"]
    counts = [failure_counts[status] for status in ordered_labels]

    plt.figure(figsize=(7, 4))
    plt.bar(ordered_labels, counts, color="silver")
    plt.xlabel('job status', fontsize=10, fontweight='bold')
    plt.ylabel('job count', fontsize=10, fontweight='bold')
    plt.savefig(os.path.join("checkup", 'job_status.png'), bbox_inches="tight", dpi=600)

def plot_authors(author_counts):
    """Generate a bar chart of job counts per author and save to ``checkup/author_credit.png``.

    Parameters
    ----------
    author_counts : dict
        Mapping of author names to job counts.
    """
    format_plot()
    authors = list(author_counts.keys())
    counts = list(author_counts.values())

    plt.figure(figsize=(7, 4))
    plt.bar(authors, counts, color="silver")
    plt.xlabel('authors', fontsize=10, fontweight='bold')
    plt.ylabel('job count', fontsize=10, fontweight='bold')
    plt.savefig(os.path.join("checkup", 'author_credit.png'), bbox_inches="tight", dpi=600)

def plot_failure_modes_from_csv(csv_path):
    """Generate a bar chart of failure mode counts from CSV and save to file.

    Parameters
    ----------
    csv_path : str
        Path to the ``failure_modes.csv`` file.
    """
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
    """Classify all QM jobs by status and generate summary reports.

    Scans the output directory for submitted jobs, classifies each as
    done, running, queued, backlog, charge error, memory error, or unknown,
    and writes ``checkup/failure_modes.csv`` and ``checkup/job_status.csv``.

    Parameters
    ----------
    method : str
        DFT method name (subdirectory under each chain directory).
    output : str
        Path to the top-level output directory.
    delete_queued : bool
        If True, delete ``.submit_record`` files for unfinished jobs so
        they can be resubmitted.

    Returns
    -------
    tuple of (dict, dict)
        ``(failure_counts, author_counts)`` where ``failure_counts`` maps
        status labels to counts and ``author_counts`` maps author names
        to job counts.
    """
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

