import os
import csv

def count_atoms_in_pdb_files(chain_path, num_spheres):
    """ Count atoms in specified PDB files within a chain directory. """
    total_atoms = 0
    for i in range(1, num_spheres + 1):
        pdb_file = os.path.join(chain_path, f"{i}.pdb")
        if os.path.exists(pdb_file):
            with open(pdb_file, 'r') as file:
                for line in file:
                    if line.startswith("ATOM"):
                        total_atoms += 1
    return total_atoms + 1  # Adding one for the metal

def analyze_pdb_structure(root_dir, num_spheres):
    """ Analyze the PDB structure and return a list of tuples (pdb_id, chain, atom_count). """
    pdb_data = []

    for pdb_id in os.listdir(root_dir):
        pdb_path = os.path.join(root_dir, pdb_id)
        if os.path.isdir(pdb_path):
            for chain in os.listdir(pdb_path):
                if chain == "Protoss":  # Skip the Protoss directory
                    continue
                chain_path = os.path.join(pdb_path, chain)
                if os.path.isdir(chain_path):
                    atom_count = count_atoms_in_pdb_files(chain_path, num_spheres)
                    pdb_data.append((pdb_id, chain, atom_count))

    return pdb_data

def write_to_csv(data, output_file):
    """ Write the data to a CSV file. """
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['PDB ID', 'Chain', 'Number of Atoms'])
        for row in data:
            writer.writerow(row)

# Set the number of spheres here
num_spheres = 2  # For example, set to 2 to include 1.pdb and 2.pdb

# Example usage
root_directory = '.' # Current working directory
output_csv = 'pdb_atoms_count.csv'
pdb_data = analyze_pdb_structure(root_directory, num_spheres)
write_to_csv(pdb_data, output_csv)

# Now the script is ready for use in your Python environment.
