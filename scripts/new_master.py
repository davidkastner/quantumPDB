import os

def main():
    # Read PDB ID codes from the file
    with open('protein_master_list.txt', 'r') as f:
        master_list = [line.strip() for line in f.readlines()]

    # Get list of directories in current working directory
    current_dirs = [name for name in os.listdir() if os.path.isdir(name)]

    # Find the PDB IDs which do not have corresponding folders
    missing_dirs = [pdb_id for pdb_id in master_list if pdb_id not in current_dirs]

    # Print and save the list of missing IDs to a new file
    with open('missing_ids.txt', 'w') as out_file:
        for pdb_id in missing_dirs:
            out_file.write(pdb_id + '\n')

if __name__ == "__main__":
    main()