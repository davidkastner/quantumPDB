import os
import time
import math
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from Bio.PDB import PDBParser, calc_dihedral
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

def format_plot() -> None:
    """General plotting parameters for the Kulik Lab."""

    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams["xtick.major.pad"] = 5
    plt.rcParams["ytick.major.pad"] = 5
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["svg.fonttype"] = "none"

def get_distance_dihedral(nitrogen, prev_carbon, prev_oxygen, iron_atoms):
    """Return the distance and dihedral angle if they are valid."""
    valid_distances_dihedrals = [(nitrogen - iron, abs(math.degrees(calc_dihedral(iron.get_vector(), nitrogen.get_vector(), prev_carbon.get_vector(), prev_oxygen.get_vector()))))
                                 for iron in iron_atoms
                                 if 1.0 < (nitrogen - iron) < 2.4 and 120 < abs(math.degrees(calc_dihedral(iron.get_vector(), nitrogen.get_vector(), prev_carbon.get_vector(), prev_oxygen.get_vector()))) < 180]

    return valid_distances_dihedrals[0] if valid_distances_dihedrals else (None, None)

def check_motif(structure, filename):
    distance_dihedral_data = []
    hits = []
    iron_atoms = [atom for atom in structure.get_atoms() if atom.element in ['FE', 'FE2', 'FE3']]
    for model in structure:
        for chain in model:
            prev_residue = None
            for residue in chain:
                if 'N' in residue and prev_residue and all(atom in prev_residue for atom in ['C', 'O']):
                    distance, dihedral = get_distance_dihedral(residue['N'], prev_residue['C'], prev_residue['O'], iron_atoms)
                    if distance is not None and dihedral is not None:
                        residue_info = f"{chain.id}_{residue.get_resname()}{residue.id[1]}"
                        distance_dihedral_data.append((distance, dihedral, filename, residue_info))
                        hits.append(residue_info)
                prev_residue = residue

    if hits:
        protein_name = filename.upper().split('.')[0]
        print(f'{protein_name} -> {", ".join(hits)}')

    return distance_dihedral_data

def print_extremes(all_data, start_time, plot_name):
    df = pd.DataFrame(all_data, columns=['Distance', 'Dihedral', 'Filename', 'Residue'])
    highest_dihedral = df.loc[df['Dihedral'].idxmax()]
    lowest_dihedral = df.loc[df['Dihedral'].idxmin()]
    highest_distance = df.loc[df['Distance'].idxmax()]
    lowest_distance = df.loc[df['Distance'].idxmin()]

    total_time = round(time.time() - start_time, 3)
    print(
        f"""
        ------------------------------RESULTS-------------------------------
        RESULT:  Highest distance: {round(highest_distance["Distance"], 2)} Å from {highest_distance["Filename"]} -> {highest_distance["Residue"]}
                 Lowest distance: {round(lowest_distance["Distance"], 2)} Å from {lowest_distance["Filename"]} -> {lowest_distance["Residue"]}
                 Highest dihedral: {round(highest_dihedral["Dihedral"], 2)}° from {highest_dihedral["Filename"]} -> {highest_dihedral["Residue"]}
                 Lowest dihedral: {round(lowest_dihedral["Dihedral"], 2)}° from {lowest_dihedral["Filename"]} -> {lowest_dihedral["Residue"]}
        OUTPUT:  Generated a plot {plot_name}.
        TIME:    Total execution time: {total_time} seconds.
        --------------------------------------------------------------------\n
        """
    )

def get_kde_plot(all_data, plot_name):
    """Generate and save a KDE plot of the data."""
    df = pd.DataFrame(all_data, columns=['Distance', 'Dihedral', 'Filename', 'Residue'])
    plt.figure(figsize=(5, 4))
    cmap = LinearSegmentedColormap.from_list('custom_cmap', ["#F1F8FF", "#002349"])
    sns.kdeplot(data=df, x='Distance', y='Dihedral', cmap=cmap, fill=True, bw_adjust=.5, cbar=True)
    plt.xlabel('distance (Å)', fontsize=12, fontweight='bold')
    plt.ylabel('dihedral angle (°)', fontsize=12, fontweight='bold')
    plt.savefig(plot_name, bbox_inches="tight", format="png", dpi=600)
    plt.close()


def main():
    start_time = time.time()
    parser = PDBParser(QUIET=True)
    all_data = []
    for filename in sorted(os.listdir('.')):
        if filename.lower().endswith('.pdb'):
            structure = parser.get_structure(filename, filename)
            data = check_motif(structure, filename)
            all_data.extend(data)
    all_data = [(float(distance), float(dihedral), str(filename), str(residue_info)) for distance, dihedral, filename, residue_info in all_data]


    # Visualize the data using KDE plot
    plot_name = 'kde.png'
    format_plot()
    get_kde_plot(all_data, plot_name)
    print_extremes(all_data, start_time, plot_name)

if __name__ == '__main__':
    main()

