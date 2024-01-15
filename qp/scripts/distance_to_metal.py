import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def format_figure():
    """
    Sets formatting for matplotlib.

    """
    font = {"family": "sans-serif", "weight": "bold", "size": 14}
    plt.rc("font", **font)
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["axes.linewidth"] = 2.5
    plt.rcParams["xtick.major.size"] = 10
    plt.rcParams["xtick.major.width"] = 2.5
    plt.rcParams["ytick.major.size"] = 10
    plt.rcParams["ytick.major.width"] = 2.5
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["mathtext.default"] = "regular"

def parse_pdb(file_path):
    """Parse a PDB file to extract atom coordinates."""
    coords = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append(np.array([x, y, z]))
                except ValueError as e:
                    print(f"Error parsing line: {line}")
                    print(e)
    return coords

def compute_distances(coord1, coords):
    """Compute distances between a reference coordinate and a list of coordinates."""
    return [np.linalg.norm(coord1 - coord) for coord in coords]

def process_directory(root_dir):
    """Process the directories and compute distances."""
    data = []

    for pdb_id in os.listdir(root_dir):
        pdb_dir = os.path.join(root_dir, pdb_id)
        if os.path.isdir(pdb_dir):
            for chain in os.listdir(pdb_dir):
                chain_dir = os.path.join(pdb_dir, chain)
                if os.path.isdir(chain_dir):
                    files = os.listdir(chain_dir)
                    if '0.pdb' in files:
                        metal_coords = parse_pdb(os.path.join(chain_dir, '0.pdb'))
                        if metal_coords:
                            metal_coord = metal_coords[0]  # Assuming only one metal atom
                            for sphere_num, sphere in enumerate(['1.pdb', '2.pdb', '3.pdb']):
                                sphere_path = os.path.join(chain_dir, sphere)
                                if os.path.exists(sphere_path):
                                    sphere_coords = parse_pdb(sphere_path)
                                    distances = compute_distances(metal_coord, sphere_coords)
                                    for dist in distances:
                                        data.append((sphere_num+1, dist))

    return pd.DataFrame(data, columns=['coordination_sphere', 'distance'])

def plot_data(df):
    """Plot the data using an improved strip plot with short, thick horizontal median lines and vertical capped lines."""
    format_figure()
    plt.figure(figsize=(5, 4))

    colors = ['#f77f00', '#d62828', '#003049']
    for i in range(1, 4):
        sns.stripplot(x='coordination_sphere', y='distance', data=df[df['coordination_sphere'] == i], 
                      jitter=True, alpha=0.0075, color=colors[i-1])

    # Adding short and thick horizontal lines for the median
    medians = df.groupby(['coordination_sphere'])['distance'].median()
    max_values = df.groupby(['coordination_sphere'])['distance'].max()
    min_values = df.groupby(['coordination_sphere'])['distance'].min()

    for i in range(1, 4):
        plt.hlines(medians.loc[i], xmin=i-1-0.12, xmax=i-1+0.12, color='black', lw=4, zorder=5, alpha=1)
        # Adding a vertical line with caps to mark the maximum and minimum values
        range_val = (max_values.loc[i] - min_values.loc[i]) / 2
        plt.errorbar(i-1, min_values.loc[i] + range_val, yerr=range_val, color='black', capsize=5, zorder=6)

    plt.xlabel('coordination sphere', weight="bold")
    plt.ylabel('distance (Å)', weight="bold")
    plt.xticks([0, 1, 2], ['first', 'second', 'third'])
    plt.savefig('strip_plot.png', bbox_inches="tight", dpi=600)

# Directory where the script is executed and filename for the CSV
root_dir = '.'
csv_filename = 'atom_metal_distances.csv'

# Check if the CSV file already exists
if os.path.exists(csv_filename):
    print(f"   > The file {csv_filename} already exists")
    distance_data = pd.read_csv(csv_filename)
else:
    # Process the directories and compute distances
    distance_data = process_directory(root_dir)

    # Check if any data was added
    if not distance_data.empty:
        # Save the distance data to a CSV file
        distance_data.to_csv(csv_filename, index=False)
        print(f"Data saved to {csv_filename}")
    else:
        print("No data was added to the DataFrame. Please check the directory structure and file formats.")

# Plot the data
if not distance_data.empty:
    plot_data(distance_data)
else:
    print("No data available for plotting.")