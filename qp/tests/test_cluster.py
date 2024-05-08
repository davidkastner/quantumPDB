import pytest
import numpy as np
import os
import glob
import filecmp

from qp.cluster import coordination_spheres


def check_clusters(path, out, metal_ids):
    for metal in metal_ids:
        assert os.path.isdir(os.path.join(path, metal)), f"Cluster center {metal} not found"
        sphere_count = len(glob.glob(os.path.join(path, metal, "?.pdb")))
        for i in range(sphere_count):
            expected_pdb = os.path.join(path, metal, f"{i}.pdb")
            output_pdb = os.path.join(out, metal, f"{i}.pdb")
            assert filecmp.cmp(expected_pdb, output_pdb), f"Sphere {i} PDB does not match expected"

    # Expected charge.csv includes additional ligands, which must be excluded
    expected_charge = os.path.join(path, "charge.csv")
    output_charge = os.path.join(out, "charge.csv")
    with open(expected_charge, "r") as e, open(output_charge, "r") as o:
        expected_lines = e.readlines()
        output_lines = o.readlines()
        assert expected_lines[:len(output_lines)] == output_lines, "Charge does not match expected"

    expected_count = os.path.join(path, "count.csv")
    output_count = os.path.join(out, "count.csv")
    assert filecmp.cmp(expected_count, output_count), "Residue count does not match expected"


@pytest.mark.parametrize("sample_cluster", [
    ("1sp9", ("A446", "B446")),
    ("2q4a", ("A901", "B902")),
    ("3a8g", ("A301",))
], indirect=True)
def test_extract_clusters(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"],
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("4ilv", ("A301", "B301"))], indirect=True)
def test_cap_heavy(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], capping=2, 
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("1lm6", ("A204",))], indirect=True)
def test_box_plot(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"],
        smooth_method="box_plot"
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2r6s", ("A501",))], indirect=True)
def test_dbscan(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 
        smooth_method="dbscan", eps=6, min_samples=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [
    ("2chb", (
        "A1_A2_A3_A4",
        "B1_B2_B3_B4_B5",
        "C1_C2_C3_C4",
        "I1_I2_I3_I4_I5",
        "J1_J2_J3_J4_J5"
    )),
    ("4z42", ("C601_C602", "F601_F602", "I601_I602", "L601_L602"))
], indirect=True)
def test_merge_centers(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["BGC", "GAL", "NGA", "SIA", "NI"],
        merge_cutoff=4.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2fd8", ("B501_B502_B503",))], indirect=True)
def test_prune_atoms(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["DT", "MA7"],
        max_atom_count=102, merge_cutoff=2.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)
