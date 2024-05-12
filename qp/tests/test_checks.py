import pytest
import os

from qp.checks import fetch_pdb


def test_fetch_pdb(tmpdir):
    pdb = "1lm6"
    out = os.path.join(tmpdir, f"{pdb}.pdb")

    fetch_pdb.fetch_pdb(pdb, out)
    assert os.path.getsize(out) > 0, "Found empty PDB file"

    with pytest.raises(ValueError):
        fetch_pdb.fetch_pdb("XXXX", out)


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
    
    expected_pdbs = [(p, os.path.join(tmpdir, p, f"{p}.pdb")) for p in pdbs]
    output_pdbs = fetch_pdb.parse_input(input_pdbs, tmpdir)
    assert expected_pdbs == output_pdbs, "Parsed input does not match expected"
