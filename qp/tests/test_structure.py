import pytest
import os
import glob
import shutil
import filecmp

from qp.protonate import get_protoss, fix
from qp.structure import struct_to_file

# Skip Modeller tests if in Github actions
MISSING_LICENSE = False
try:
    import modeller
except:
    MISSING_LICENSE = True

if not MISSING_LICENSE:
    from qp.structure import missing_loops


# ========== protonate ==========

@pytest.mark.parametrize("sample_pdb", ["1lm6"], indirect=True)
def test_protoss(tmpdir, sample_pdb):
    pdb, path = sample_pdb
    pdb_path = os.path.join(path, f"{pdb}.pdb")
    out = os.path.join(tmpdir, f"{pdb}_protoss.pdb")

    pid = get_protoss.upload(pdb_path)
    job = get_protoss.submit(pid)
    get_protoss.download(job, out)
    assert os.path.getsize(out) > 0, "Found empty PDB file"


@pytest.mark.parametrize("sample_pdb", ["2fd8"], indirect=True)
def test_repair_ligands(tmpdir, sample_pdb):
    pdb, path = sample_pdb
    pdb_path = os.path.join(path, f"{pdb}_modeller.pdb")
    expected_prot = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    input_prot = os.path.join(path, f"{pdb}_protoss_raw.pdb")
    output_prot = os.path.join(tmpdir, f"{pdb}_protoss.pdb")

    shutil.copy(input_prot, output_prot)
    get_protoss.repair_ligands(output_prot, pdb_path)
    assert filecmp.cmp(expected_prot, output_prot), "Repaired Protoss PDB does not match expected"


@pytest.mark.parametrize("sample_pdb", ["1sp9", "2q4a", "3a8g", "3x20", "6f2a"], indirect=True)
def test_adjust_activesites(tmpdir, sample_pdb):
    pdb, path = sample_pdb
    expected_prot = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")
    input_prot = os.path.join(path, "Protoss", f"{pdb}_protoss_old.pdb")
    output_prot = os.path.join(tmpdir, f"{pdb}_protoss.pdb")

    shutil.copy(input_prot, output_prot)
    fix.adjust_activesites(output_prot, ["FE", "FE2"])
    assert filecmp.cmp(expected_prot, output_prot), "Adjusted Protoss PDB does not match expected"


@pytest.mark.parametrize("sample_pdb", ["1lm6", "1sp9", "2q4a", "2r6s", "3a8g", "4ilv"], indirect=True)
def test_compute_charge(tmpdir, sample_pdb):
    pdb, path = sample_pdb
    sdf_path = os.path.join(path, "Protoss", f"{pdb}_ligands.sdf")
    prot_path = os.path.join(path, "Protoss", f"{pdb}_protoss.pdb")

    expected_charge = {}
    with open(os.path.join(path, f"charge.csv"), "r") as f:
        for l in f.readlines()[::-1]:
            if l == "\n":
                break
            ligand, charge = l.split(",")
            expected_charge[ligand] = int(charge)
    output_charge = get_protoss.compute_charge(sdf_path, prot_path)
    assert expected_charge == output_charge, "Ligand charge does not match expected"


# ========== missing_loops ==========

@pytest.mark.skipif(MISSING_LICENSE, reason="Modeller license not found")
@pytest.mark.parametrize("sample_pdb", ["1lm6", "1sp9", "2q4a", "2r6s", "3a8g", "4ilv"], indirect=True)
def test_write_alignment(tmpdir, sample_pdb):
    pdb, path = sample_pdb
    pdb_path = os.path.join(path, f"{pdb}.pdb")
    AA = missing_loops.define_residues()
    residues = missing_loops.get_residues(pdb_path, AA)
    residues = missing_loops.clean_termini(residues)

    expected_ali = os.path.join(path, f"{pdb}.ali")
    output_ali = os.path.join(tmpdir, f"{pdb}.ali")
    missing_loops.write_alignment(residues, pdb, pdb_path, output_ali)
    assert filecmp.cmp(expected_ali, output_ali), "Alignment file does not match expected"


@pytest.mark.skipif(MISSING_LICENSE, reason="Modeller license not found")
@pytest.mark.parametrize("sample_pdb", ["2r6s"], indirect=True)
def test_build_model(tmpdir, sample_pdb):
    pdb, path = sample_pdb
    pdb_path = os.path.join(path, f"{pdb}.pdb")
    ali_path = os.path.join(path, f"{pdb}.ali")

    AA = missing_loops.define_residues()
    residues = missing_loops.get_residues(pdb_path, AA)
    residues = missing_loops.clean_termini(residues)

    expected_modeller = os.path.join(path, f"{pdb}_modeller.pdb")
    output_modeller = os.path.join(tmpdir, f"{pdb}_modeller.pdb")
    missing_loops.build_model(residues, pdb, pdb_path, ali_path, output_modeller)

    # First line contains timestamp, ignore when comparing
    with open(expected_modeller, "r") as e, open(output_modeller, "r") as o:
        expected_lines = e.readlines()
        output_lines = o.readlines()
        assert expected_lines[1:] == output_lines[1:], "Modeller output does not match expected"


# ========== struct_to_file ==========

@pytest.mark.parametrize("sample_cluster", [
    ("1lm6", ("A204",)),
    ("1sp9", ("A446", "B446")),
    ("2q4a", ("A901", "B902")),
    ("2r6s", ("A501",)),
    ("3a8g", ("A301",)),
    ("4ilv", ("A301", "B301"))
], indirect=True)
def test_to_xyz(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    for metal in metal_ids:
        sphere_paths = sorted(glob.glob(os.path.join(path, metal, "?.pdb")))
        expected_xyz = os.path.join(path, metal, f"{metal}.xyz")
        output_xyz = os.path.join(tmpdir, f"{metal}.xyz")
        struct_to_file.to_xyz(output_xyz, *sphere_paths)
        assert filecmp.cmp(expected_xyz, output_xyz), f"XYZ file does not match expected"


@pytest.mark.parametrize("sample_cluster", [
    ("1lm6", ("A204",)),
    ("1sp9", ("A446", "B446")),
    ("2q4a", ("A901", "B902")),
    ("2r6s", ("A501",)),
    ("3a8g", ("A301",)),
    ("4ilv", ("A301", "B301"))
], indirect=True)
def test_combine_pdbs(tmpdir, sample_cluster):
    pdb, metal_ids, path = sample_cluster
    for metal in metal_ids:
        sphere_paths = sorted(glob.glob(os.path.join(path, metal, "?.pdb")))
        expected_pdb = os.path.join(path, metal, f"{metal}.pdb")
        output_pdb = os.path.join(tmpdir, f"{metal}.pdb")
        struct_to_file.combine_pdbs(output_pdb, ["FE", "FE2"], *sphere_paths)
        assert filecmp.cmp(expected_pdb, output_pdb), f"Combined PDB does not match expected"
