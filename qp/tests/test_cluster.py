import pytest
import numpy as np
import os
import filecmp

from qp.cluster import coordination_spheres


def check_clusters(path, out, metal_ids):
    for metal in metal_ids:
        for i in range(4):
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
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=1, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("4ilv", ("A301", "B301"))], indirect=True)
def test_cap_heavy(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=2, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("1lm6", ("A204",))], indirect=True)
def test_box_plot(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=1, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="box_plot"
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2r6s", ("A501",))], indirect=True)
def test_dbscan(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=1, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="dbscan", eps=6, min_samples=3
    )
    check_clusters(path, tmpdir, metal_ids)
