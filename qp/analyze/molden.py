import numpy as np
import os

ECP_ACTIVE_ELECTRONS = {
    # s-block
    'K': 9, 'CA': 10,
    'RB': 9, 'SR': 10,
    'CS': 9, 'BA': 10,

    # 3d transition metals (LANL, 10e core except Zn)
    'SC': 11, 'TI': 12, 'V': 13, 'CR': 14, 'MN': 15,
    'FE': 16, 'CO': 17, 'NI': 18, 'CU': 19, 'ZN': 12,

    # 4d transition metals (ECP28)
    'Y': 11, 'ZR': 12, 'NB': 13, 'MO': 14, 'TC': 15,
    'RU': 16, 'RH': 17, 'PD': 18, 'AG': 19, 'CD': 12,

    # 5d transition metals
    'LA': 11,
    'HF': 12, 'TA': 13, 'W': 14, 'RE': 15,
    'OS': 16, 'IR': 17, 'PT': 18, 'AU': 19, 'HG': 12,

    # p-block heavier main-group
    'GA': 3, 'GE': 4, 'AS': 5, 'SE': 6, 'BR': 7,
    'IN': 3, 'SN': 4, 'SB': 5, 'TE': 6, 'I': 7,
    'TL': 13, 'PB': 4, 'BI': 5,
}

def correct_ecp_charges_inplace(molden_file):
    """
    Corrects the nuclear charge in a Molden file for elements with ECPs,
    overwriting the original file. Returns True if the file was modified.
    """
    try:
        with open(molden_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"> ERROR: Molden file not found at {molden_file}")
        return False

    modified_lines = []
    in_atoms_section = False
    was_modified = False
    for line in lines:
        if line.strip() == "[Atoms] Angs":
            in_atoms_section = True
            modified_lines.append(line)
            continue
        elif in_atoms_section and line.strip().startswith("["):
            in_atoms_section = False
        
        if in_atoms_section:
            parts = line.strip().split()
            if len(parts) >= 6:
                element_symbol = parts[0].upper()
                # Check if the element needs correction and if it hasn't been corrected already
                if element_symbol in ECP_ACTIVE_ELECTRONS and parts[2] != str(ECP_ACTIVE_ELECTRONS[element_symbol]):
                    parts[2] = str(ECP_ACTIVE_ELECTRONS[element_symbol])
                    line_to_write = '     '.join(parts) + '\n'
                    modified_lines.append(line_to_write)
                    was_modified = True
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)
    
    if was_modified:
        with open(molden_file, 'w') as f:
            f.writelines(modified_lines)
        print(f"> Overwrote {os.path.basename(molden_file)} with ECP-corrected nuclear charges.")
    
    return was_modified


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