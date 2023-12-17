"""Troubleshooting utility that checks for common QM failure modes."""

import os
import glob

def check_failure_mode(filepath):
    """Checks for specific failure mode keywords in generated output."""

    with open(filepath, 'r') as f:
        content = f.read()
        if "Incorrect molecular charge or spin multiplicity" in content:
            return "- Charge"
        elif "In Alloc2D: malloc failed" in content:
            return "- Memory"
        elif "Unable to find atom" in content:
            return "- Atom type"
        elif "Job terminated" in content:
            return "- Unknown"
    return None

def main():
    """Loop over all jobs and check if they failed."""
    
    print(f"   > Checking for failed QM jobs.")
    output_name = "failure_checkup.txt"
    with open(output_name, "w") as output_file:
        pdb_directories = sorted(glob.glob("."))
        for pdb in pdb_directories:
            for dir in sorted(glob.glob('[0-9]*')):
                for dirpath, dirnames, filenames in os.walk(dir):
                    for filename in filenames:
                        if filename == "qmscript.out":
                            filepath = os.path.join(dirpath, filename)
                            failure_mode = check_failure_mode(filepath)
                            if failure_mode:
                                output_file.write(f"{filepath}  {failure_mode}\n")
    print(f"   > Saving checkup results in {output_name}\n")


if __name__ == '__main__':
    main()