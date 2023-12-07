import os
import csv
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
from adjustText import adjust_text
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def compute_general_charge_features(pos_centroid, neg_centroid, metal_centroid, feature_types):
    features = {}
    if pos_centroid is not None and neg_centroid is not None:
        if 'general_charge_dist' in feature_types:
            features['neg_pos_dist'] = compute_distance(neg_centroid, pos_centroid)
        if 'general_charge_angles' in feature_types:
            features['neg_pos_angle'] = compute_angle(neg_centroid, metal_centroid, pos_centroid)
    return features

def compute_angle(A, B, C):
    """Compute the angle at B formed by lines BA and BC."""
    BA = A - B
    BC = C - B
    cosine_angle = np.dot(BA, BC) / (np.linalg.norm(BA) * np.linalg.norm(BC))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    return np.degrees(angle)

def compute_centroid(atoms):
    """Compute the geometric centroid of a set of atoms."""
    coords = np.array([atom.get_coord() for atom in atoms])
    return np.mean(coords, axis=0)

def compute_distance(point1, point2):
    """Compute the Euclidean distance between two points."""
    return np.linalg.norm(point1 - point2)

def compute_residue_features(residues, pos_amino_acids, neg_amino_acids):
    """Compute the features corresponding to the groups of residues."""
    aa_counts = {aa: 0 for aa in pos_amino_acids.keys() | neg_amino_acids.keys()}
    aa_atoms = {aa: [] for aa in pos_amino_acids.keys() | neg_amino_acids.keys()}
    pos_charge_count = neg_charge_count = 0
    pos_atoms = []
    neg_atoms = []

    for res in residues:
        resname = res.get_resname()
        res_atoms = [atom.get_id() for atom in res]

        # Count and collect atoms for positively charged amino acids
        if resname in pos_amino_acids:
            if all(atom in res_atoms for atom in pos_amino_acids[resname]):
                pos_charge_count += 1
                aa_counts[resname] += 1
                aa_atoms[resname].extend(atom for atom in res if atom.get_id() in pos_amino_acids[resname])

        # Count and collect atoms for negatively charged amino acids
        elif resname in neg_amino_acids:
            if all(atom not in res_atoms for atom in neg_amino_acids[resname]):
                neg_charge_count += 1
                aa_counts[resname] += 1
                aa_atoms[resname].extend(atom for atom in res if atom.get_id() in neg_amino_acids[resname])

        # Calculate for the general charge features
        if resname in pos_amino_acids and all(atom in res_atoms for atom in pos_amino_acids[resname]):
            pos_atoms.extend(res)
        elif resname in neg_amino_acids and all(atom not in res_atoms for atom in neg_amino_acids[resname]):
            neg_atoms.extend(res)

    aa_centroids = {aa: compute_centroid(atoms) if atoms else None for aa, atoms in aa_atoms.items()}
    pos_centroid = compute_centroid(pos_atoms) if pos_atoms else None
    neg_centroid = compute_centroid(neg_atoms) if neg_atoms else None

    return aa_counts, aa_centroids, pos_charge_count, neg_charge_count, pos_centroid, neg_centroid


def compute_geometric_features(aa_centroids, metal_centroid, feature_types):
    features = {}
    for aa1, centroid1 in aa_centroids.items():
        if centroid1 is not None:
            for aa2, centroid2 in aa_centroids.items():
                if centroid2 is not None and aa1 != aa2:
                    if 'residue_dist' in feature_types:
                        features[f'{aa1}_{aa2}_dist'] = compute_distance(centroid1, centroid2)
                    if 'residue_angles' in feature_types:
                        features[f'{aa1}_{aa2}_angle'] = compute_angle(centroid1, metal_centroid, centroid2)
            if 'residue_dist' in feature_types:
                features[f'{aa1}_metal_dist'] = compute_distance(centroid1, metal_centroid)
    return features



def process_pdb_files(pdb_id, chain_dir, pos_amino_acids, neg_amino_acids, feature_types):
    """Compute the features for the different spheres."""
    results = {}
    pdb_parser = PDBParser(QUIET=True)
    metal_structure = pdb_parser.get_structure('metal', os.path.join(pdb_id, chain_dir, '0.pdb'))
    metal_centroid = compute_centroid([atom for atom in metal_structure.get_atoms()])

    for sphere in range(1, 4):
        pdb_file = f'{sphere}.pdb'
        try:
            structure = pdb_parser.get_structure('temp', os.path.join(pdb_id, chain_dir, pdb_file))
            aa_counts, aa_centroids, pos_count, neg_count, pos_centroid, neg_centroid = compute_residue_features(structure.get_residues(), pos_amino_acids, neg_amino_acids)
            # General charge features
            if 'general_charge_counts' in feature_types:
                results[f'{sphere}_pos_charge_count'] = pos_count
                results[f'{sphere}_neg_charge_count'] = neg_count

            # General charge angles and distances
            if 'general_charge_angles' in feature_types or 'general_charge_dist' in feature_types:
                general_charge_features = compute_general_charge_features(pos_centroid, neg_centroid, metal_centroid, feature_types)
                results.update({f'{sphere}_{k}': v for k, v in general_charge_features.items()})

            # Residue-specific features
            if 'residue_charge_count' in feature_types:
                results.update({f'{sphere}_{k}_count': v for k, v in aa_counts.items()})
            if 'residue_angles' in feature_types or 'residue_dist' in feature_types:
                geom_features = compute_geometric_features(aa_centroids, metal_centroid, feature_types)
                results.update({f'{sphere}_{k}': v for k, v in geom_features.items()})
            
            # Inside process_pdb_files function
            if 'residue_angles' in feature_types or 'residue_dist' in feature_types:
                geom_features = compute_geometric_features(aa_centroids, metal_centroid, feature_types)
                results.update({f'{sphere}_{k}': v for k, v in geom_features.items()})

        except Exception as e:
            print(f"Error processing {pdb_file} in {chain_dir}: {e}")

    return results


def main(feature_types, csv_name):
    """Main handler function for parsing, calculating, and analyzing data."""
    charge_results = []
    neg_amino_acids = {"ASP": ["HOD1", "HOD2"],
                       "GLU": ["HOE1", "HOE2"],
                       "CYS": ["HG"],
                       "TYR": ["HH"],
                       }
    pos_amino_acids = {"ARG": ["HE", "HH11", "HH12", "HH21", "HH22"],
                       "LYS": ["HZ1", "HZ2", "HZ3"],
                       "HIS": ["HD1", "HE2"],
                       }

    pdb_directories = [d for d in os.listdir('.') if os.path.isdir(d)]
    for pdb_id in tqdm(pdb_directories, desc="Processing PDBs"):
        if os.path.isdir(pdb_id):
            for chain_dir in os.listdir(pdb_id):
                if os.path.isdir(os.path.join(pdb_id, chain_dir)) and chain_dir != 'Protoss':
                    charge_data = process_pdb_files(pdb_id, chain_dir, pos_amino_acids, neg_amino_acids, feature_types)
                    charge_results.append({"pdb_id": pdb_id, "chain": chain_dir, **charge_data})

    with open(csv_name, 'w', newline='') as csvfile:
        if charge_results:
            all_fieldnames = set()
            for result in charge_results:
                all_fieldnames.update(result.keys())

            fieldnames = ['pdb_id', 'chain'] + sorted(all_fieldnames - {'pdb_id', 'chain'})
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in charge_results:
                writer.writerow(row)
        else:
            print("No data to write to CSV.")

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


def caclulate_pca(csv_name):
    """Perform PCA analysis."""
    # Load your dataset
    data = pd.read_csv(csv_name)
    data.fillna(9999, inplace=True)

    # Selecting all features except the first two columns
    features = data.columns[2:]
    x = data.loc[:, features].values
    pdb_ids = data.iloc[:, 0].values  # Assumed first column
    chains = data.iloc[:, 1].values  # Assumed second column

    # Standardizing the features
    x = StandardScaler().fit_transform(x)

    # Performing PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    
    # Save out the PCA analysis as a csv
    principalDf.insert(0, 'Chain', chains)
    principalDf.insert(0, 'PDB_ID', pdb_ids)
    principalDf.to_csv('charge_pca.csv', index=False)


    # Print the composition of principal components with feature names
    print("PCA Component 1:")
    for weight, feature in sorted(zip(pca.components_[0], features), reverse=True):
        print(f"{feature}: {weight}")
    print("\nPCA Component 2:")
    for weight, feature in sorted(zip(pca.components_[1], features), reverse=True):
        print(f"{feature}: {weight}")

    return principalDf, pdb_ids


def plot_pca(principalDf, pdb_ids, out_plot_name):
    """Create the plot from the PCA analysis."""
    # Plotting the results
    format_plot()
    plt.figure(figsize=(4, 4))
    plt.axhline(y=0, color='k', linestyle='-', linewidth=2, alpha=.2, zorder=1)
    plt.axvline(x=0, color='k', linestyle='-', linewidth=2, alpha=.2, zorder=1)
    plt.xlabel('Principal Component 1', fontsize=12, fontweight='bold')
    plt.ylabel('Principal Component 2', fontsize=12, fontweight='bold')
    plt.scatter(principalDf['PC1'], principalDf['PC2'], alpha=0.25, color="blue")

    # Labeling outliers
    distances = np.sqrt(principalDf['PC1']**2 + principalDf['PC2']**2)
    outliers = distances.nlargest(10).index

    texts = []
    for outlier in outliers:
        x, y = principalDf.iloc[outlier]['PC1'], principalDf.iloc[outlier]['PC2']
        texts.append(plt.text(x, y, pdb_ids[outlier], ha='right', va='bottom'))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1))

    plt.savefig(out_plot_name, bbox_inches="tight", format="png", dpi=600)
    plt.close()


if __name__ == "__main__":
    # Select which feature categories who would like to extract from quantumPDB
    feature_types = ['general_charge_counts',
                     'general_charge_dist',
                     'general_charge_angles',
                     'residue_charge_count',
                     'residue_dist',
                     'residue_angles',
                     ]
    
    # Generate features from quantumPDB
    csv_name = "charge_results.csv"
    if not os.path.exists(csv_name):
        main(feature_types, csv_name)

    # Generate PCA plots
    out_plot_name = "pca_6.png"
    principalDf, pdb_ids = caclulate_pca(csv_name)
    plot_pca(principalDf, pdb_ids, out_plot_name)
