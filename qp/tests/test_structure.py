import pytest
import os
import filecmp

# Skip tests if in Github actions
MISSING_LICENSE = False
try:
    import modeller
except:
    MISSING_LICENSE = True

if not MISSING_LICENSE:
    from qp.structure import missing_loops


@pytest.fixture
def sample_pdb(request):
    pdb = request.param
    return pdb, os.path.join(os.path.dirname(__file__), "samples", pdb)


@pytest.mark.skipif(MISSING_LICENSE, reason="MODELLER license not found")
@pytest.mark.parametrize("sample_pdb", ["1lm6", "1sp9", "2q4a", "3a8g", "4qm8"], indirect=True)
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


@pytest.mark.skipif(MISSING_LICENSE, reason="MODELLER license not found")
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
        assert expected_lines[1:] == output_lines[1:], "MODELLER output does not match expected"