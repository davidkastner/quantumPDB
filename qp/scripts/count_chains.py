import os

def count_chain_directories():
    current_directory = os.getcwd()
    total_chain_directories = 0

    # List all directories in the current working directory
    for pdb_dir in os.listdir(current_directory):
        if os.path.isdir(pdb_dir):
            # If it's a PDB directory, check its contents
            for chain_dir in os.listdir(pdb_dir):
                chain_path = os.path.join(pdb_dir, chain_dir)
                if os.path.isdir(chain_path) and chain_dir != "Protoss":
                    total_chain_directories += 1

    return total_chain_directories

if __name__ == '__main__':
    total = count_chain_directories()
    print(f"There are {total} chain directories.")

