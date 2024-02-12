def read_ids_from_file(file_path):
    with open(file_path, 'r') as file:
        return set(file.read().splitlines())

def write_ids_to_file(file_path, ids):
    with open(file_path, 'w') as file:
        for id in sorted(ids):
            file.write(f"{id}\n")

def main():
    original_ids = read_ids_from_file('original.txt')
    others_ext_ids = read_ids_from_file('others-ext.txt')
    others_ids = read_ids_from_file('others.txt')

    # Combine all sets and remove duplicates
    all_unique_ids = original_ids.union(others_ext_ids, others_ids)

    # Write the combined unique IDs to a new file
    write_ids_to_file('combined.txt', all_unique_ids)

    print(f"Combined file 'combined.txt' created with {len(all_unique_ids)} unique IDs.")

if __name__ == "__main__":
    main()
