import os
import csv
import glob
import shutil
import subprocess
from pathlib import Path
from Bio.PDB import PDBParser
from importlib import resources
from periodictable import elements


def get_settings_ini_path():
    # Dynamically get path to package settings.ini file
    with resources.path('qp.resources', 'settings.ini') as ini_path:
        settings_ini_path = str(ini_path)

    return settings_ini_path

def get_atmrad_path():
    # Dynamically get path to package atmrad/ directory
    with resources.path('qp.resources', '') as resources_path:
        atmrad_path = os.path.join(str(resources_path), 'atmrad')

    return atmrad_path

def get_com():
    """
    Gets the center of mass of the atoms in the PDB file (0.pdb)
    
    Returns
    -------
    center_of_mass: list
        List containing the x, y, z coordinates
    """
    scr_dir_path = os.getcwd()
    chain_dir_path = os.path.abspath(os.path.join(scr_dir_path, '..', '..'))
    substrate_pdb_path = os.path.join(chain_dir_path, "0.pdb")
    parser = PDBParser()
    structure = parser.get_structure('substrate', substrate_pdb_path)
    total_mass = 0.0
    center_of_mass = [0.0, 0.0, 0.0]
    
    for atom in structure.get_atoms():
        element = atom.element.strip().capitalize()
        mass = getattr(elements, element).mass  # Get the atomic mass from periodictable
        coord = atom.coord
        
        # Update the weighted average for the center of mass
        total_mass += mass
        center_of_mass[0] += mass * coord[0]
        center_of_mass[1] += mass * coord[1]
        center_of_mass[2] += mass * coord[2]
    
    # Finalize the center of mass by dividing by the total mass
    center_of_mass = [coord / total_mass for coord in center_of_mass]
    
    return center_of_mass


def get_atom_range():
    """
    Gets the range of atoms in the file 0.pdb, which is the substrate.
    Returns the atoms as a range e.g., "1-32".
    
    Returns
    -------
    fragment_atom_range: str
        A range of numbers represented as a string e.g., "1-32"
    """
    scr_dir_path = os.getcwd()
    chain_dir_path = os.path.abspath(os.path.join(scr_dir_path, '..', '..'))
    substrate_pdb_path = os.path.join(chain_dir_path, "0.pdb")
    parser = PDBParser()
    structure = parser.get_structure('substrate', substrate_pdb_path)
    atom_count = len(list(structure.get_atoms()))  # Get total number of atoms
    fragment_atom_range = f"1-{atom_count}"

    return fragment_atom_range


def iterate_qm_output(pdb_all, method, base_output_dir, multiwfn_path, settings_ini_path, atmrad_path, task_function):
    """
    Loop over all QM job folders and apply the provided task function 
    for successful jobs with a .molden file inside the 'scr' directory.

    Parameters
    ----------
    method: str
        The name of the method (e.g., wpbeh) being used for QM jobs.
    
    base_output_dir: str
        The base directory containing all QM job folders.

    settings_ini_path: str
        The path to the settings.ini file to be copied into each QM job folder.
    
    atmrad_path: str
        The path to the atmrad directory to be copied into each QM job folder.

    task_function: callable
        A function that will be applied to each valid QM job directory.
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
                if os.path.exists(atmrad_dest_path):
                    shutil.rmtree(atmrad_dest_path)
                shutil.copytree(atmrad_path, atmrad_dest_path)

                settings_ini_dest_path = os.path.join(scr_dir_path, 'settings.ini')
                shutil.copy(settings_ini_path, settings_ini_dest_path)

                if os.path.exists(scr_dir_path):
                    molden_files = glob.glob(os.path.join(scr_dir_path, '*.molden'))
                    if molden_files:
                        # Call the task function with the directory and molden file
                        molden_file = os.path.basename(molden_files[0])
                        task_function(scr_dir_path, molden_file, multiwfn_path)
                        # Clean up by deleting the atmrad directory and settings.ini file after the task
                        shutil.rmtree(atmrad_dest_path)  # Delete the atmrad directory when done
                        os.remove(settings_ini_dest_path)  # Delete the settings.ini file when done
                    else:
                        print(f"> WARNING: No molden file in {scr_dir_path}")
                else:
                    print(f"> WARNING: No scr directory in {method_dir_path}.")
        else:
            continue


def charge_scheme(scr_dir_path, molden_file,  multiwfn_path):
    """
    Calculate a variety of charge schemes with Multiwfn for a given .molden file and
    process the output to extract relevant information, saving it in a readable format.

    Parameters
    ----------
    scr_dir_path: str
        The directory containing the molden file.
    
    molden_file: str
        The full path to the molden file to be processed with Multiwfn.
    """

    original_dir = os.getcwd()
    os.chdir(scr_dir_path)
    pdb_name = scr_dir_path.split("/")[-4]
    chain_name = scr_dir_path.split("/")[-3]

    molden_base_name = os.path.splitext(os.path.basename(molden_file))[0]
    # charge_schemes = {"Hirshfeld": "1", "Voronoi": "2", "Mulliken": "5", "ADCH": "11", "Hirshfeld-I": "15", "CM5": "16"}
    charge_schemes = {"Hirshfeld-I": "15"}
    charge_command = f"{multiwfn_path} {molden_file}"
    fragment_atom_range = get_atom_range() # Get the atoms to define a Multiwfn fragment

    for scheme_name, scheme_code in charge_schemes.items():
        new_charge_file = f"{molden_base_name}_{scheme_name}.txt"
        if os.path.exists(new_charge_file):
            continue
        proc = subprocess.Popen(charge_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        commands = ["7", "-1", fragment_atom_range, scheme_code, "1", "y", "0", "0", "q"]
        output, _ = proc.communicate("\n".join(commands).encode())  # Send the commands to Multiwfn

        # Check if the code finished
        charge_output_file = Path(f"{molden_base_name}.chg")
        if charge_output_file.exists():
            charge_output_file.rename(new_charge_file)
        else:
            print(f"> ERROR: Multiwfn failed for {pdb_name} {chain_name}")
    
    print(f"> {pdb_name.upper()} {chain_name} partial charge schemes computed", flush=True)
    os.chdir(original_dir)


def is_valid_dipole_file(file_path):
    """
    Check if the dipole output file exists and contains the expected completion text.
    
    We do this instead of just checking if the file exists because sometimes a job will not finish.
    This jobs can take a long time and might not always finish correctly.
    If one does not finish correctly, we rerun it.
    
    """
    if not os.path.exists(file_path):
        return False
    with open(file_path, 'r') as f:
        for line in f:
            if "Magnitude of molecular dipole moment (a.u.&Debye):" in line:
                return True
    return False


def calc_dipole(scr_dir_path, molden_file, multiwfn_path):
    """
    Calculate dipole of the substrate with Multiwfn for a given .molden file and
    process the output to extract relevant information, saving it in a readable format.

    Parameters
    ----------
    scr_dir_path: str
        The directory containing the molden file.
    
    molden_file: str
        The full path to the molden file to be processed with Multiwfn.
    """
    original_dir = os.getcwd()
    os.chdir(scr_dir_path)
    pdb_name = scr_dir_path.split("/")[-4]
    chain_name = scr_dir_path.split("/")[-3]
    molden_base_name = os.path.splitext(os.path.basename(molden_file))[0]
    raw_out_file = f"{molden_base_name}_dipole.out"
    dipole_results_file = f"{molden_base_name}_dipole.csv"

    dipole_command = f"{multiwfn_path} {molden_file} > {raw_out_file}"
    fragment_atom_range = get_atom_range()
    center_of_mass_list = get_com()
    center_of_mass_str = ",".join(map(str, center_of_mass_list))
    
    if not is_valid_dipole_file(raw_out_file):
        proc = subprocess.Popen(dipole_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        commands = ["15", "-10", center_of_mass_str, "2", "-5", fragment_atom_range, "-1", "3", "2", "0", "q"]
        output, _ = proc.communicate("\n".join(commands).encode())  # Send the commands to Multiwfn
    
    if not is_valid_dipole_file(raw_out_file) and not os.path.exists(dipole_results_file):
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
    
    print(f"> {pdb_name.upper()} {chain_name} dipole computed", flush=True)
    os.chdir(original_dir)