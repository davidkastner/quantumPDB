import os
import re
import csv
import glob
import shutil
import subprocess
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
from importlib import resources
from periodictable import elements
from qp.analyze import molden


def iterate_qm_output(pdb_all, method, base_output_dir, multiwfn_path, settings_ini_path, atmrad_path, charge_scheme, task_function):
    """Iterate over QM job folders and apply an analysis function.

    Loops through all completed QM calculations, copying necessary Multiwfn
    configuration files and applying the provided task function to each
    job that has a ``.molden`` file in its ``scr/`` directory.

    Parameters
    ----------
    pdb_all : list of tuple
        List of ``(pdb_id, pdb_path)`` tuples from the input CSV.
    method : str
        QM method name (e.g., ``'wpbeh'``) used to locate job directories.
    base_output_dir : str
        Base directory containing all PDB subdirectories with QM jobs.
    multiwfn_path : str
        Path to the Multiwfn executable.
    settings_ini_path : str
        Path to the Multiwfn ``settings.ini`` configuration file.
    atmrad_path : str
        Path to the ``atmrad/`` directory containing atomic radii data
        for Hirshfeld-I calculations.
    charge_scheme : str
        Charge scheme to use (e.g., ``'Hirshfeld'``, ``'CM5'``).
    task_function : callable
        Function to apply to each valid job. Called as
        ``task_function(scr_dir_path, molden_file, multiwfn_path, charge_scheme)``.

    Notes
    -----
    If the molden file was modified by ECP charge correction, old analysis
    files are deleted to force re-calculation.
    """
    pdb_ids = {pdb[0] for pdb in pdb_all}
    # Loop over all PDB directories in the base_output_dir
    all_pdb_dirs = sorted(glob.glob(os.path.join(base_output_dir, '[0-9]*')))
    for pdb_dir in all_pdb_dirs:  # Loop over PDB directories
        current_pdb_dir = pdb_dir.split("/")[-1]
        if current_pdb_dir in pdb_ids: # Check if this PDB is in the users csv
            for chain_dir in os.listdir(pdb_dir):  # Loop over chain subdirectories
                chain_dir_path = os.path.join(pdb_dir, chain_dir)
                if not os.path.isdir(chain_dir_path):
                    continue
                
                method_dir_path = os.path.join(chain_dir_path, method)
                if not os.path.exists(method_dir_path):
                    continue

                scr_dir_path = os.path.join(method_dir_path, 'scr')

                # Copy the atmrad dir and settings.ini into the current working directory
                atmrad_dest_path = os.path.join(scr_dir_path, 'atmrad/')
                # Make sure scr_dir_path exists (defensive)
                os.makedirs(scr_dir_path, exist_ok=True)

                # Copy atmrad with race-safe semantics
                try:
                    # Works on Python 3.8+ (your env is 3.8)
                    shutil.copytree(atmrad_path, atmrad_dest_path, dirs_exist_ok=True)
                except TypeError:
                    # For safety if running under <3.8 somewhere else:
                    if not os.path.exists(atmrad_dest_path):
                        try:
                            shutil.copytree(atmrad_path, atmrad_dest_path)
                        except FileExistsError:
                            # Someone else created it between exists() and copytree(); ignore
                            pass
                except FileExistsError:
                    # Created by another process between the check and copy; ignore
                    pass

                settings_ini_dest_path = os.path.join(scr_dir_path, 'settings.ini')
                if not os.path.exists(settings_ini_dest_path):
                    shutil.copy(settings_ini_path, settings_ini_dest_path)
                
                cpu_count = get_cpu_count(settings_ini_dest_path)
                print(f"> Using {cpu_count} threads for Multiwfn", flush=True)

                if os.path.exists(scr_dir_path):
                    # Move into the QM scr directory
                    original_dir = os.getcwd()
                    os.chdir(scr_dir_path)
                    scr_dir_path = os.getcwd()

                    molden_files = glob.glob(f'{chain_dir}.molden')
                    if molden_files:
                        molden_file = os.path.basename(molden_files[0])
                        
                        # Correct the Molden file and get modification status
                        was_modified = molden.correct_ecp_charges_inplace(molden_file)
                        
                        # If the molden file was updated, delete old analysis files to force re-calculation.
                        if was_modified:
                            print("> Molden file was updated. Deleting old analysis files to force re-calculation.")
                            molden_base_name = os.path.splitext(molden_file)[0]
                            
                            # Potential output filenames to delete
                            charge_file_to_delete = f"{molden_base_name}_{charge_scheme}.chg"
                            dipole_csv_to_delete = f"{molden_base_name}_dipole_com.csv"
                            dipole_raw_to_delete = f"{molden_base_name}_dipole_com.out"

                            # Delete files if they exist
                            for f in [charge_file_to_delete, dipole_csv_to_delete, dipole_raw_to_delete]:
                                if os.path.exists(f):
                                    os.remove(f)
                                    print(f">   Deleted old file: {f}")

                        try:
                            task_function(scr_dir_path, molden_file, multiwfn_path, charge_scheme)
                        except Exception as e:
                            print(f"> ERROR: {current_pdb_dir} {chain_dir} failed: {e}. Skipping.")             
                        # Commenting out the deletion lines to prevent race conditions
                        # shutil.rmtree('atmrad/')  # Delete the atmrad directory when done
                        # os.remove('settings.ini')  # Delete the settings.ini file when done
                    else:
                        print(f"> WARNING: No molden file in {scr_dir_path}", flush=True)
                    os.chdir(original_dir)
                else:
                    print(f"> WARNING: No scr directory in {method_dir_path}.", flush=True)
        else:
            continue


def get_settings_ini_path():
    """Get the path to the bundled Multiwfn settings.ini file.

    Returns
    -------
    str
        Absolute path to the ``settings.ini`` file in the package resources.
    """
    # Dynamically get path to package settings.ini file
    with resources.path('qp.resources', 'settings.ini') as ini_path:
        settings_ini_path = str(ini_path)

    return settings_ini_path


def get_atmrad_path():
    """Get the path to the bundled atomic radii files for Hirshfeld-I.

    The ``atmrad/`` directory contains free-atom electron density data
    required by Multiwfn for Hirshfeld-I charge calculations.

    Returns
    -------
    str
        Absolute path to the ``atmrad/`` directory in the package resources.
    """
    # Dynamically get path to package atmrad/ directory
    with resources.path('qp.resources', '') as resources_path:
        atmrad_path = os.path.join(str(resources_path), 'atmrad')

    return atmrad_path


def get_com(fragment_atom_indices, scr_dir_path):
    """
    Gets the center of mass of the atoms in the PDB file (0.pdb).
    Calculates the center of mass only for the given fragment atom indices.

    Parameters
    ----------
    fragment_atom_indices : list of int
        A list of atom serial numbers to include in the calculation.

    Returns
    -------
    center_of_mass : list
        List containing the x, y, z coordinates
    """
    chain_dir_path = os.path.abspath(os.path.join(scr_dir_path, '..', '..'))
    substrate_pdb_path = os.path.join(chain_dir_path, "0.pdb")
    parser = PDBParser()
    structure = parser.get_structure('substrate', substrate_pdb_path)
    total_mass = 0.0
    center_of_mass = [0.0, 0.0, 0.0]

    fragment_atom_indices_set = set(fragment_atom_indices)  # Convert to set for faster lookup

    for atom in structure.get_atoms():
        atom_serial_number = atom.get_serial_number()
        if atom_serial_number in fragment_atom_indices_set:
            element = atom.element.strip().capitalize()
            # Handle cases where the element symbol might be missing or incorrect
            if not element or len(element) > 2:
                element = atom.get_name().strip()[0].capitalize()
            mass = getattr(elements, element).mass  # Get the atomic mass from periodictable
            coord = atom.coord

            # Update the weighted average for the center of mass
            total_mass += mass
            center_of_mass[0] += mass * coord[0]
            center_of_mass[1] += mass * coord[1]
            center_of_mass[2] += mass * coord[2]

    if total_mass == 0.0:
        raise ValueError("Total mass is zero. No atoms found in the specified fragment_atom_indices.")

    # Finalize the center of mass by dividing by the total mass
    center_of_mass = [coord / total_mass for coord in center_of_mass]

    return center_of_mass



def get_coc(scr_dir_path, molden_file, fragment_atom_indices, selected_charge_scheme, multiwfn_path):
    """
    Gets the center of partial charges of the atoms from the .chg file.
    Calculates the center of partial charges only for the given fragment atom indices.

    Parameters
    ----------
    scr_dir_path : str
        The directory containing the .chg file and molden file.
    molden_file : str
        The molden file name, used to construct the .chg file name.
    fragment_atom_indices : list of int
        A list of atom serial numbers to include in the calculation.
    selected_charge_scheme : str
        The charge scheme used to generate the .chg file.
    multiwfn_path : str
        Path to the Multiwfn executable, needed to generate the .chg file if missing.

    Returns
    -------
    center_of_charge : list
        List containing the x, y, z coordinates of the center of charge.
    """
    # Build the .chg file name
    molden_base_name = os.path.splitext(os.path.basename(molden_file))[0]
    charge_file = f"{molden_base_name}_{selected_charge_scheme}.chg"

    # Check if the .chg file exists
    if not os.path.exists(charge_file):
        print(f"> Charge file {charge_file} not found. Generating it using charge_scheme function.", flush=True)
        # Generate the .chg file using the charge_scheme function
        charge_scheme(scr_dir_path, molden_file, multiwfn_path, selected_charge_scheme)

    # Verify that the .chg file now exists
    if not os.path.exists(charge_file):
        raise FileNotFoundError(f"Charge file {charge_file} could not be generated.")

    # Read the .chg file
    coords = []
    charges = []
    with open(charge_file, 'r') as f:
        lines = f.readlines()

    # Read the PDB file to get the atom serial numbers in order
    chain_dir_path = os.path.abspath(os.path.join(scr_dir_path, '..', '..'))
    substrate_pdb_path = os.path.join(chain_dir_path, "0.pdb")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('substrate', substrate_pdb_path)

    # Build a list of atom serial numbers in order
    pdb_atom_serial_numbers = [atom.get_serial_number() for atom in structure.get_atoms()]

    if len(lines) != len(pdb_atom_serial_numbers):
        raise ValueError("Mismatch between number of atoms in .chg file and PDB file.")

    fragment_atom_indices_set = set(fragment_atom_indices)
    total_charge = 0.0
    weighted_position = np.array([0.0, 0.0, 0.0])

    for i, line in enumerate(lines):
        tokens = line.strip().split()
        if len(tokens) < 5:
            continue
        element = tokens[0]
        x, y, z = map(float, tokens[1:4])
        charge = float(tokens[4])

        # Get the atom serial number corresponding to this atom
        atom_serial_number = pdb_atom_serial_numbers[i]

        if atom_serial_number in fragment_atom_indices_set:
            total_charge += charge
            weighted_position += np.array([x, y, z]) * charge

    if total_charge == 0.0:
        raise ValueError("Total charge is zero. No atoms found in the specified fragment_atom_indices or charges sum to zero.")

    center_of_charge = weighted_position / total_charge

    return center_of_charge.tolist()



def get_atom_range():
    """
    Gets the range of atoms in the file 0.pdb, which is the substrate.
    Returns the atoms as a range string (e.g., "1-32") and a list of atom serial numbers.
    This is primarily because Multiwfn needs the range string format.

    Returns
    -------
    fragment_atom_range : str
        A range of numbers represented as a string e.g., "1-32"
    fragment_atom_indices : list of int
        A list of atom serial numbers (integers)
    """
    chain_dir_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir))
    substrate_pdb_path = os.path.join(chain_dir_path, "0.pdb")
    parser = PDBParser()
    structure = parser.get_structure('substrate', substrate_pdb_path)
    atom_serial_numbers = [atom.get_serial_number() for atom in structure.get_atoms()]
    if not atom_serial_numbers:
        raise ValueError("No atoms found in the PDB file.")

    min_serial = min(atom_serial_numbers)
    max_serial = max(atom_serial_numbers)
    fragment_atom_range = f"{min_serial}-{max_serial}"
    fragment_atom_indices = atom_serial_numbers  # List of atom serial numbers

    return fragment_atom_range, fragment_atom_indices


def get_cpu_count(settings_ini_dest_path):
    """
    Reads the settings.ini file and returns the number of threads specified by the 'nthreads=' line.

    Parameters
    ----------
    settings_ini_dest_path : str
        The path to the settings.ini file.

    Returns
    -------
    int
        The number of threads specified in the settings.ini file, or 1 if not found or an error occurs.
    """
    try:
        with open(settings_ini_dest_path, 'r') as file:
            for line in file:
                # Look for a line that contains 'nthreads=' and extract the number after it
                match = re.search(r'nthreads=\s*(\d+)', line)
                if match:
                    # Extract the first matching group (the number after 'nthreads=')
                    return int(match.group(1))
    except FileNotFoundError:
        print(f"Error: The file {settings_ini_dest_path} was not found.", flush=True)
    except ValueError:
        print(f"Error: Unable to parse 'nthreads' in {settings_ini_dest_path}.", flush=True)
    
    # Default to 1 if nthreads isn't found or there is an error
    return 1



def is_valid_dipole_file(file_path):
    """Check whether a dipole calculation output file is complete.

    Multiwfn dipole calculations can sometimes fail or be interrupted.
    This function checks for the presence of a completion marker in
    the output file to determine if the calculation finished successfully.

    Parameters
    ----------
    file_path : str
        Path to the dipole output file (``*_dipole_com.out``).

    Returns
    -------
    bool
        True if the file exists and contains the completion marker.
    """
    if not os.path.exists(file_path):
        return False
    with open(file_path, 'r') as f:
        for line in f:
            if "Magnitude of molecular dipole moment (a.u.&Debye):" in line:
                return True
    return False


def charge_scheme(scr_dir_path, molden_file, multiwfn_path, charge_scheme):
    """Calculate atomic partial charges using Multiwfn.

    Runs Multiwfn to compute atomic charges using the specified charge
    partitioning scheme and saves the results to a ``.chg`` file.

    Parameters
    ----------
    scr_dir_path : str
        Directory containing the molden file (used for logging).
    molden_file : str
        Filename of the molden file to process.
    multiwfn_path : str
        Path to the Multiwfn executable.
    charge_scheme : str
        Charge partitioning scheme. Must be one of: ``'Hirshfeld'``,
        ``'Voronoi'``, ``'Mulliken'``, ``'ADCH'``, ``'Hirshfeld-I'``,
        ``'CM5'``.

    Raises
    ------
    ValueError
        If ``charge_scheme`` is not one of the supported schemes.

    Notes
    -----
    Output is written to ``{molden_base}_{charge_scheme}.chg`` in the
    current working directory. If this file already exists, the
    calculation is skipped.
    """

    # Mapping of charge schemes to their corresponding Multiwfn codes
    available_charge_schemes = {
        "Hirshfeld": "1",
        "Voronoi": "2",
        "Mulliken": "5",
        "ADCH": "11",
        "Hirshfeld-I": "15",
        "CM5": "16"
    }

    # Check if the provided charge_scheme is valid
    if charge_scheme not in available_charge_schemes:
        raise ValueError(f"Invalid charge scheme '{charge_scheme}'. Must be one of {list(available_charge_schemes.keys())}.")

    # Proceed with the selected charge scheme
    pdb_name = scr_dir_path.split("/")[-4]
    chain_name = scr_dir_path.split("/")[-3]

    molden_base_name = os.path.splitext(os.path.basename(molden_file))[0]
    charge_command = f"{multiwfn_path} {molden_file}"

    scheme_code = available_charge_schemes[charge_scheme]
    new_charge_file = f"{molden_base_name}_{charge_scheme}.chg"

    # If the charge file already exists, skip
    if os.path.exists(new_charge_file):
        print(f"> Charge file {new_charge_file} for {pdb_name} already exists. Skipping {charge_scheme}...", flush=True)
        return

    # Run Multiwfn with the user-selected charge scheme
    proc = subprocess.Popen(charge_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    commands = ["7", scheme_code, "1", "y", "0", "q"]
    output, _ = proc.communicate("\n".join(commands).encode())  # Send the commands to Multiwfn

    # Check if the code finished and output the charge file
    charge_output_file = Path(f"{molden_base_name}.chg")
    if charge_output_file.exists():
        charge_output_file.rename(new_charge_file)
        print(f"> {pdb_name.upper()} {chain_name} {charge_scheme} scheme computed and saved to {new_charge_file}", flush=True)
    else:
        print(f"> ERROR: Multiwfn failed for {pdb_name} {chain_name} using {charge_scheme} scheme.", flush=True)



def calc_dipole(scr_dir_path, molden_file, multiwfn_path, charge_scheme):
    """Calculate the substrate dipole moment using Multiwfn.

    Computes the molecular dipole moment of the substrate fragment using the
    center of mass as the reference point. The molden file is first
    translated so that the substrate center of mass is at the origin.

    Parameters
    ----------
    scr_dir_path : str
        Directory containing the molden file (typically ``method/scr/``).
    molden_file : str
        Filename of the molden file to process.
    multiwfn_path : str
        Path to the Multiwfn executable.
    charge_scheme : str
        Charge scheme name (passed for API consistency but not used in
        dipole calculation).

    Notes
    -----
    Outputs are written to:
    - ``{molden_base}_dipole_com.out`` - raw Multiwfn output
    - ``{molden_base}_dipole_com.csv`` - parsed dipole components (a.u.)

    The substrate atoms are identified from ``0.pdb`` in the cluster directory.
    """
    pdb_name = scr_dir_path.split("/")[-4]
    chain_name = scr_dir_path.split("/")[-3]
    molden_base_name = os.path.splitext(os.path.basename(molden_file))[0]

    fragment_atom_range, fragment_atom_indices = get_atom_range()
    center_of_mass_list = get_com(fragment_atom_indices, scr_dir_path)
    raw_out_file = f"{molden_base_name}_dipole_com.out"
    dipole_results_file = f"{molden_base_name}_dipole_com.csv"


    output_molden_file = molden.center_molden(molden_file, center_of_mass_list)
    dipole_command = f"{multiwfn_path} {output_molden_file} > {raw_out_file}"
    
    if not is_valid_dipole_file(raw_out_file):
        proc = subprocess.Popen(dipole_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        commands = ["15", "-5", fragment_atom_range, "-1", "3", "2", "0", "q"]
        output, _ = proc.communicate("\n".join(commands).encode())  # Send the commands to Multiwfn

    if is_valid_dipole_file(raw_out_file) and not os.path.exists(dipole_results_file):
        # Parse raw_out_file and generate dipole_results_file
        dipole_a_u = None
        magnitude_a_u = None

        with open(raw_out_file, 'r') as f:
            for line in f:
                if "Molecular dipole moment (a.u.)" in line:
                    dipole_a_u = line.split()[4:7]  # Extract the x, y, z components in a.u.
                if "Magnitude of molecular dipole moment (a.u.&Debye)" in line:
                    magnitude_a_u = line.split()[6]  # Extract the magnitude in a.u.

        if dipole_a_u and magnitude_a_u:
            # Write results to CSV
            with open(dipole_results_file, mode='w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Dipole Component (a.u.)', 'X (a.u.)', 'Y (a.u.)', 'Z (a.u.)', 'Magnitude (a.u.)'])
                csv_writer.writerow(['Molecular dipole moment'] + dipole_a_u + [magnitude_a_u])
    
    print(f"> {pdb_name.upper()} {chain_name} dipole computed using center of mass\n", flush=True)
