import pytest
import os

from qp.structure import setup


def test_fetch_pdb(tmpdir):
    pdb = "1lm6"
    out = os.path.join(tmpdir, f"{pdb}.pdb")

    setup.fetch_pdb(pdb, out)
    assert os.path.getsize(out) > 0, "Found empty PDB file"

    with pytest.raises(ValueError):
        setup.fetch_pdb("XXXX", out)


def test_parse_input(tmpdir):
    pdbs = ["1lm6", "1sp9", "2q4a", "2r6s", "3a8g", "4ilv"]
    batch = os.path.join(tmpdir, f"pdbs.txt")
    with open(batch, "w") as f:
        f.write("\n".join(pdbs[:4]))
    pdb_path = os.path.join(tmpdir, "4ilv", "4ilv.pdb")
    os.makedirs(os.path.join(tmpdir, "4ilv"))
    with open(pdb_path, "w") as f:
        pass
    input_pdbs = [batch, "3a8g", pdb_path]
    center_yaml_residues = ["FE", "FE2"]

    expected_pdbs = [(p, os.path.join(tmpdir, p, f"{p}.pdb"), None) for p in pdbs]
    # Local .pdb path keeps the absolute input path rather than the download target
    expected_pdbs[-1] = ("4ilv", pdb_path, None)
    output_pdbs, output_centers = setup.parse_input(
        input_pdbs, tmpdir, center_yaml_residues
    )
    assert expected_pdbs == output_pdbs, "Parsed input does not match expected"
    assert output_centers == center_yaml_residues, "Center residues do not match expected"
