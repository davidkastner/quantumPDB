"""Troubleshooting utility that checks for common failure modes."""

import os

def check_failure_mode(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
        if "Incorrect molecular charge or spin multiplicity" in content:
            return "charge"
        elif "In Alloc2D: malloc failed" in content:
            return "memory"
        elif "Job terminated" in content:
            return "unknown"
    return None

def main():
    with open("failure_checkup.txt", "w") as output_file:
        for dirpath, dirnames, filenames in os.walk('.'):
            for filename in filenames:
                if filename == "qmscript.out":
                    filepath = os.path.join(dirpath, filename)
                    failure_mode = check_failure_mode(filepath)
                    if failure_mode:
                        output_file.write(f"{filepath}  {failure_mode}\n")

if __name__ == '__main__':
    main()
