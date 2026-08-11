import pytest
import numpy as np
import os
import glob
import filecmp

from qp.cluster import spheres
from qp.cluster.spheres import CenterResidue


FE_CENTER = CenterResidue("FE_FE2")
SUGAR_CENTER = CenterResidue("BGC_GAL_NGA_SIA_NI")
DNA_CENTER = CenterResidue("DT_B501-MA7_B502-DT_B503")


@pytest.mark.parametrize(
    "spec, mode, residues",
    [
        ("FE", "fuzzy", ["FE"]),
        ("FE_FE2", "fuzzy", ["FE", "FE2"]),
        ("FE_A199", "fuzzy", ["FE", "A199"]),
        ("exact:FE_A199", "strict", ["FE_A199"]),
        ("EXACT:FE_A199", "strict", ["FE_A199"]),
        ("CU_A357-CU_A358", "strict", ["CU_A357", "CU_A358"]),
        ("exact:FE_A155-HIS_A93", "strict", ["FE_A155", "HIS_A93"]),
        ("DT_B501-MA7_B502-DT_B503", "strict", ["DT_B501", "MA7_B502", "DT_B503"]),
    ],
)
def test_center_residue_parsing(spec, mode, residues):
    center = CenterResidue(spec)
    assert center.mode == mode
    assert center.residue_list == residues


def test_center_residue_exact_prefix_requires_key():
    with pytest.raises(ValueError, match="exact:"):
        CenterResidue("exact:")


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
    spheres.extract_clusters(
        pdb_path, tmpdir, FE_CENTER,
        smooth_method="dummy_atom", mean_distance=3, capping=0
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("4ilv", ("A301", "B301"))], indirect=True)
def test_cap_heavy(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, FE_CENTER, capping=2,
        smooth_method="dummy_atom", mean_distance=3
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("1lm6", ("A204",))], indirect=True)
def test_box_plot(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, FE_CENTER,
        smooth_method="box_plot", capping=0
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2r6s", ("A501",))], indirect=True)
def test_dbscan(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, FE_CENTER,
        smooth_method="dbscan", eps=6, min_samples=3, capping=0
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
    spheres.extract_clusters(
        pdb_path, tmpdir, SUGAR_CENTER,
        merge_cutoff=4.0,
        smooth_method="dummy_atom", mean_distance=3, capping=0
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [("2fd8", ("B501_B502_B503",))], indirect=True)
def test_prune_atoms(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, DNA_CENTER,
        max_atom_count=102, merge_cutoff=2.0,
        smooth_method="dummy_atom", mean_distance=3, capping=0
    )
    check_clusters(path, tmpdir, metal_ids)


@pytest.mark.parametrize("sample_cluster", [
    ("2chb", (
        "A1_A2_A3_A4",
        "B1_B2_B3_B4_B5",
        "C1_C2_C3_C4",
        "I1_I2_I3_I4_I5",
        "J1_J2_J3_J4_J5"
    ))
], indirect=True)
def test_cluster_name_template(tmpdir, sample_cluster):
    """cluster_name_template should produce short, collision-free names
    (avoiding the long residue-concatenated names in the default case,
    which can exceed path length limits when many centers are merged),
    while the underlying sphere geometry stays identical."""
    pdb, metal_ids, path = sample_cluster
    pdb_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    spheres.extract_clusters(
        pdb_path, tmpdir, SUGAR_CENTER,
        merge_cutoff=4.0,
        smooth_method="dummy_atom", mean_distance=3,
        first_sphere_radius=4.0,
        cluster_name_template="A_{radius}",
        capping=0,
    )

    expected_names = ["A_4", "A_4_1", "A_4_2", "A_4_3", "A_4_4"]
    for name in expected_names:
        assert os.path.isdir(os.path.join(tmpdir, name)), f"Expected cluster dir {name} not found"

    # No unexpected directories, and no collisions/overwrites occurred
    generated_dirs = sorted(
        d for d in os.listdir(tmpdir)
        if os.path.isdir(os.path.join(tmpdir, d))
    )
    assert generated_dirs == sorted(expected_names)

    # Sphere geometry should be identical to the default-naming output, just
    # under new directory names. Centers aren't guaranteed to be processed
    # in a fixed order (get_center_residues iterates over a set), so match
    # each generated cluster to its expected counterpart by content rather
    # than by position.
    remaining_expected = list(metal_ids)
    for new_name in expected_names:
        output_pdb_0 = os.path.join(tmpdir, new_name, "0.pdb")
        match = None
        for old_name in remaining_expected:
            expected_pdb_0 = os.path.join(path, old_name, "0.pdb")
            if filecmp.cmp(expected_pdb_0, output_pdb_0):
                match = old_name
                break
        assert match is not None, f"Generated cluster {new_name} does not match any expected cluster"
        remaining_expected.remove(match)

        sphere_count = len(glob.glob(os.path.join(path, match, "?.pdb")))
        for i in range(sphere_count):
            expected_pdb = os.path.join(path, match, f"{i}.pdb")
            output_pdb = os.path.join(tmpdir, new_name, f"{i}.pdb")
            assert filecmp.cmp(expected_pdb, output_pdb), (
                f"Sphere {i} PDB for {match} -> {new_name} does not match expected"
            )
    assert not remaining_expected, "Not all expected clusters were matched"


def test_cluster_name_template_bad_field(tmpdir):
    """An unknown field in the template should raise a clear error rather
    than failing silently or with a confusing traceback."""
    with pytest.raises(ValueError, match="cluster_name_template"):
        spheres.extract_clusters(
            os.path.join(os.path.dirname(__file__), "samples", "1lm6", "Protoss", "1lm6_protoss.pdb"),
            tmpdir, FE_CENTER,
            smooth_method="box_plot",
            cluster_name_template="A_{not_a_real_field}"
        )
