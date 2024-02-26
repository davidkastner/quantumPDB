import pytest
import numpy as np
import os
import filecmp

from qp.cluster import coordination_spheres


@pytest.fixture
def sample_pdb(request):
    np.random.seed(66265)
    pdb, metal_ids = request.param
    return pdb, metal_ids, os.path.join(os.path.dirname(__file__), "samples", pdb)


def check_clusters(path, out, metal_ids):
    for metal in metal_ids:
        expected_xyz = os.path.join(path, metal, f"{metal}.xyz")
        output_xyz = os.path.join(out, metal, f"{metal}.xyz")
        assert filecmp.cmp(expected_xyz, output_xyz), "Cluster XYZ does not match expected"

        expected_pdb = os.path.join(path, metal, f"{metal}.pdb")
        output_pdb = os.path.join(out, metal, f"{metal}.pdb")
        assert filecmp.cmp(expected_pdb, output_pdb), "Cluster PDB does not match expected"

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


@pytest.mark.parametrize("sample_pdb", [
    ("1sp9", ("A446", "B446")),
    ("2q4a", ("A901", "B902")),
    ("3a8g", ("A301",))
], indirect=True)
def test_extract_clusters(tmpdir, sample_pdb):
    pdb, metal_ids, path = sample_pdb
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=1, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_pdb", [("4qm8", ("A201", "B201"))], indirect=True)
def test_cap_heavy(tmpdir, sample_pdb):
    pdb, metal_ids, path = sample_pdb
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=2, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_pdb", [("1lm6", ("A204",))], indirect=True)
def test_box_plot(tmpdir, sample_pdb):
    pdb, metal_ids, path = sample_pdb
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=1, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="box_plot"
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_pdb", [("2r6s", ("A501",))], indirect=True)
def test_dbscan(tmpdir, sample_pdb):
    pdb, metal_ids, path = sample_pdb
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    coordination_spheres.extract_clusters(
        pdb_path, tmpdir, ["FE", "FE2"], 3, 
        capping=1, charge=True, count=True, xyz=True, 
        first_sphere_radius=4.0,
        smooth_method="dbscan", eps=6, min_samples=3
    )
    check_clusters(path, tmpdir, metal_ids)
