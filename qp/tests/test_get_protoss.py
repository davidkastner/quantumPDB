"""Live ProteinsPlus API v2 regression for custom PDB Protoss uploads.

Uses a small RCSB structure (1cbs: CRABP with retinoic acid) so the test
exercises real upload / Protoss / download without mocking HTTP.
"""

import os

import pytest
import requests

from qp.protonate import get_protoss
from qp.tests.protoss_helpers import (
    PROTOSS_NETWORK_ERRORS,
    proteins_plus_reachable,
    skip_if_protoss_unavailable,
)

# Cellular retinoic-acid-binding protein I with all-trans retinoic acid (REA).
# Small (~1.2k atoms) and has a clear non-water ligand for retention checks.
TEST_PDB_CODE = "1cbs"
TEST_LIGAND = "REA"
RCSB_PDB_URL = f"https://files.rcsb.org/download/{TEST_PDB_CODE.upper()}.pdb"


pytestmark = pytest.mark.skipif(
    not proteins_plus_reachable(),
    reason="ProteinsPlus API unreachable; skipping live Protoss integration test",
)


def _hetatm_names(pdb_text):
    return {
        line[17:20].strip()
        for line in pdb_text.splitlines()
        if line.startswith("HETATM")
    }


def _has_hydrogens(pdb_text):
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) > 76 and line[76:78].strip() == "H":
            return True
        if line.split()[-1] == "H":
            return True
    return False


def test_custom_pdb_upload_and_protoss_download(tmp_path):
    """Upload a custom RCSB PDB and download Protoss protein + ligand SDF."""
    pdb_path = tmp_path / f"{TEST_PDB_CODE}.pdb"
    protoss_path = tmp_path / f"{TEST_PDB_CODE}_protoss.pdb"
    ligands_path = tmp_path / f"{TEST_PDB_CODE}_ligands.sdf"

    try:
        r = requests.get(RCSB_PDB_URL, timeout=60)
        r.raise_for_status()
    except requests.RequestException as e:
        pytest.skip(
            f"RCSB unreachable ({type(e).__name__}); skipping live Protoss integration test"
        )
    pdb_path.write_text(r.text)
    assert TEST_LIGAND in _hetatm_names(r.text)

    try:
        pid = get_protoss.upload(str(pdb_path))
        job = get_protoss.submit(pid)
        get_protoss.download(job, str(protoss_path), "protein")
        get_protoss.download(job, str(ligands_path), "ligands")
    except PROTOSS_NETWORK_ERRORS as exc:
        skip_if_protoss_unavailable(exc)

    assert set(job) == {"protein", "ligands"}
    assert os.path.getsize(protoss_path) > 0
    assert os.path.getsize(ligands_path) > 0

    protoss_text = protoss_path.read_text()
    assert TEST_LIGAND in _hetatm_names(protoss_text), (
        f"{TEST_LIGAND} missing from Protoss PDB; ligand retention failed"
    )
    assert _has_hydrogens(protoss_text), "Protoss PDB has no hydrogen atoms"

    sdf = ligands_path.read_text()
    assert TEST_LIGAND in sdf, f"{TEST_LIGAND} missing from ligands SDF"
    assert "$$$$\n" in sdf, "ligands SDF missing newline after $$$$ delimiters"
    assert "$$$$" + TEST_LIGAND not in sdf.replace("$$$$\n", ""), (
        "ligands SDF has concatenated molecules without newline after $$$$"
    )
