import os
import re
import csv

def count_atoms_in_pdb(pdb_file):
    """Counts atoms in a PDB file."""
    count = 0
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                count += 1
    return count

def get_earliest_dir(dirs):
    """Get the earliest directory based on name starting with a letter followed by numbers."""
    valid_dirs = sorted([d for d in dirs if re.match(r"^[a-zA-Z][0-9]+", d)])
    return valid_dirs[0] if valid_dirs else None

def main():
    results = []
    cwd = os.getcwd()

    # Iterate over PDB directories
    for pdb_dir in os.listdir(cwd):
        pdb_path = os.path.join(cwd, pdb_dir)
        if os.path.isdir(pdb_path) and len(pdb_dir) == 4:  # Typically PDB codes are 4 characters
            subdirs = os.listdir(pdb_path)
            # Ignore the "Protoss" directory
            subdirs = [sd for sd in subdirs if sd != "Protoss"]
            
            # Get the earliest directory
            earliest_dir = get_earliest_dir(subdirs)
            if earliest_dir:
                atom_count = 0
                subdir_path = os.path.join(pdb_path, earliest_dir)
                for file in sorted(os.listdir(subdir_path)):
                    if re.match(r"^[0-9]+\.pdb$", file):  # Matches 0.pdb, 1.pdb, etc.
                        pdb_file_path = os.path.join(subdir_path, file)
                        atom_count += count_atoms_in_pdb(pdb_file_path)
                results.append((pdb_dir, atom_count))
    
    # Store results in cluster_sizes.csv
    with open("cluster_sizes.csv", "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDB", "ATOMS"])
        writer.writerows(results)

if __name__ == "__main__":
    main()