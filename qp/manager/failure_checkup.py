"""Troubleshooting utility that checks for common failure modes."""

import os
import glob

def check_failure_mode(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
        if "Incorrect molecular charge or spin multiplicity" in content:
            return "- Charge: "
        elif "In Alloc2D: malloc failed" in content:
            return "- Memory: "
        elif "Unable to find atom" in content:
            return "- Atom type: "
        elif "Job terminated" in content:
            return "- Unknown: "
    return None

def main():
    with open("failure_checkup.txt", "w") as output_file:
        for dir in sorted(glob.glob('[0-9]*')):
            for dirpath, dirnames, filenames in os.walk(dir):
                for filename in filenames:
                    if filename == "qmscript.out":
                        filepath = os.path.join(dirpath, filename)
                        failure_mode = check_failure_mode(filepath)
                        if failure_mode:
                            output_file.write(f"{filepath}  {failure_mode}\n")

if __name__ == '__main__':
    main()