import numpy as np

def translate_to_origin(coordinates, center_of_mass):
    """Translates the molecule so the center of mass is at the origin."""
    return coordinates - center_of_mass

def process_molden(input_file, output_file, com):
    """Process the Molden file and modify only the coordinates section."""
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Identify the indices of the sections
    atoms_start_idx = None
    gto_start_idx = None

    for idx, line in enumerate(lines):
        stripped_line = line.strip()
        if stripped_line == "[Atoms] Angs":
            atoms_start_idx = idx
        elif stripped_line.startswith("[") and stripped_line != "[Atoms] Angs" and atoms_start_idx is not None:
            # First section after [Atoms] Angs
            gto_start_idx = idx
            break

    if atoms_start_idx is None:
        raise ValueError("The [Atoms] Angs section was not found in the file.")

    if gto_start_idx is None:
        # If [GTO] or another section is not found, assume the rest of the file is tail_lines
        gto_start_idx = len(lines)

    # Header lines up to and including [Atoms] Angs
    header_lines = lines[:atoms_start_idx + 1]

    # Atom lines between [Atoms] Angs and the next section
    atom_lines_raw = lines[atoms_start_idx + 1 : gto_start_idx]

    # Rest of the file from the next section onwards
    tail_lines = lines[gto_start_idx:]

    # Process atom lines
    atoms = []
    coordinates = []
    atom_lines = []

    for line in atom_lines_raw:
        parts = line.strip().split()
        if len(parts) >= 6:
            try:
                x, y, z = map(float, parts[3:6])
            except ValueError as e:
                print(f"Error parsing coordinates in line: {line.strip()}")
                raise e
            atoms.append(parts[0])
            coordinates.append([x, y, z])
            atom_lines.append(parts)
        else:
            print(f"Warning: Line skipped due to unexpected format: {line.strip()}")

    if not atom_lines:
        raise ValueError("No atom lines found in the [Atoms] Angs section.")

    # Convert coordinates to numpy array
    coordinates = np.array(coordinates)
    original_center_of_mass = np.array(com)

    # Translate the coordinates
    new_coordinates = translate_to_origin(coordinates, original_center_of_mass)

    # Now write everything to the output file
    with open(output_file, 'w') as outfile:
        # Write header lines
        outfile.writelines(header_lines)

        # Write updated atom lines
        for i, parts in enumerate(atom_lines):
            # Update the coordinates with new values, formatted to match the original precision
            formatted_x = f"{new_coordinates[i][0]:.5f}"
            formatted_y = f"{new_coordinates[i][1]:.5f}"
            formatted_z = f"{new_coordinates[i][2]:.5f}"
            parts[3], parts[4], parts[5] = formatted_x, formatted_y, formatted_z
            # Reconstruct the line with original spacing
            updated_line = '     '.join(parts) + '\n'
            outfile.write(updated_line)

        # Write the rest of the file
        outfile.writelines(tail_lines)

def center_molden(molden_file, com):
    """Centers the molecule in the Molden file using the provided center of mass."""
    output_molden = f"centered_{molden_file}"
    process_molden(molden_file, output_molden, com)
    print(f"> Centered coordinates written to: {output_molden}")

    return output_molden
