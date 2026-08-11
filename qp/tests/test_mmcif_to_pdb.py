"""Tests for temporary mmCIF → PDB conversion."""

import json
import os
from unittest.mock import MagicMock

import pytest
import requests
from Bio.PDB import MMCIFParser, PDBParser

from qp.structure import setup
from qp.structure.mmcif_to_pdb import (
    OversizedStructureError,
    convert_mmcif_to_pdb,
    remap_sidecar_path,
)

# Small in-repo fixture (includes a 5-letter CCD residue for remap coverage)
CIF_MINI = os.path.join(os.path.dirname(__file__), "data", "mini_mmcif.cif")


def _cif_atom_count(cif_path):
    structure = MMCIFParser(QUIET=True).get_structure("cif", cif_path)
    return sum(1 for _ in structure.get_atoms())


def test_convert_mini_mmcif(tmpdir):
    out = os.path.join(tmpdir, "mini.pdb")
    info = convert_mmcif_to_pdb(CIF_MINI, out)

    assert os.path.isfile(out)
    assert os.path.getsize(out) > 0
    assert "A1E3R" in info["resname_map"]
    assert len(info["resname_map"]["A1E3R"]) == 3

    sidecar = remap_sidecar_path(out)
    assert os.path.isfile(sidecar)
    with open(sidecar) as handle:
        sidecar_data = json.load(handle)
    assert sidecar_data["resname_map"]["A1E3R"] == info["resname_map"]["A1E3R"]

    structure = PDBParser(QUIET=True).get_structure("mini", out)
    n_atoms = sum(1 for _ in structure.get_atoms())
    assert n_atoms == _cif_atom_count(CIF_MINI)
    assert n_atoms > 0

    for residue in structure.get_residues():
        assert len(residue.resname) <= 3
    for model in structure:
        for chain in model:
            assert len(chain.id) == 1

    # Occupancy / B-factor columns must not be glued together
    with open(out) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                assert len(line) >= 66
                assert line[54:66].count(".") <= 2


def test_get_pdbs_accepts_cif(tmpdir):
    cif_path = os.path.join(tmpdir, "demo.cif")
    with open(cif_path, "w") as handle:
        handle.write("data_demo\n")
    output = os.path.join(tmpdir, "out")
    pdb_all = setup.get_pdbs([cif_path], output)
    assert len(pdb_all) == 1
    pdb_id, pdb_path, source_cif = pdb_all[0]
    assert pdb_id == "demo"
    assert pdb_path == os.path.join(output, "demo", "demo.pdb")
    assert source_cif == os.path.abspath(cif_path)


def test_ensure_structure_pdb_converts_local_cif(tmpdir):
    pdb_path = os.path.join(tmpdir, "mini", "mini.pdb")
    status = setup.ensure_structure_pdb("mini", pdb_path, source_cif=CIF_MINI)
    assert status == "converted"
    assert os.path.isfile(pdb_path)
    assert setup.ensure_structure_pdb("mini", pdb_path, source_cif=CIF_MINI) == "exists"


def test_fetch_pdb_mmcif_fallback(tmpdir, monkeypatch):
    """When classic PDB 404s, fetch mmCIF and convert."""
    with open(CIF_MINI) as handle:
        cif_text = handle.read()

    def fake_get(url, timeout=30):
        response = MagicMock()
        if url.endswith(".pdb"):
            response.status_code = 404
            response.text = ""
        elif url.endswith(".cif"):
            response.status_code = 200
            response.text = cif_text
        else:
            response.status_code = 500
            response.text = ""
        return response

    monkeypatch.setattr(setup.requests, "get", fake_get)
    out = os.path.join(tmpdir, "mini", "mini.pdb")
    status = setup.fetch_pdb("mini", out)
    assert status == "mmcif"
    assert os.path.isfile(out)
    assert os.path.isfile(os.path.join(tmpdir, "mini", "mini.cif"))
    structure = PDBParser(QUIET=True).get_structure("mini", out)
    assert sum(1 for _ in structure.get_atoms()) == _cif_atom_count(CIF_MINI)


def test_fetch_pdb_invalid_id(tmpdir, monkeypatch):
    def fake_get(url, timeout=30):
        response = MagicMock()
        response.status_code = 404
        response.text = ""
        return response

    monkeypatch.setattr(setup.requests, "get", fake_get)
    out = os.path.join(tmpdir, "XXXX.pdb")
    with pytest.raises(ValueError, match="XXXX"):
        setup.fetch_pdb("XXXX", out)


def test_oversized_structure_raises(tmpdir, monkeypatch):
    """Oversized conversions raise OversizedStructureError for batch skip handling."""
    monkeypatch.setattr("qp.structure.mmcif_to_pdb.MAX_PDB_ATOMS", 1)
    out = os.path.join(tmpdir, "too_big.pdb")
    with pytest.raises(OversizedStructureError, match="atoms"):
        convert_mmcif_to_pdb(CIF_MINI, out)
    assert not os.path.isfile(out)


def test_center_resname_aliases_from_remap():
    """Original 5-letter CCD codes expand to mapped names for center matching."""
    from qp.structure.mmcif_to_pdb import (
        expand_resnames_for_matching,
        normalize_center_key,
    )

    resname_map = {"A1E3R": "A1E"}
    assert expand_resnames_for_matching(["A1E3R"], resname_map) == {"A1E3R", "A1E"}
    assert expand_resnames_for_matching(["A1E"], resname_map) == {"A1E"}
    assert normalize_center_key("A1E3R_A302", resname_map) == "A1E_A302"
    assert normalize_center_key("A1E_A302", resname_map) == "A1E_A302"
    assert normalize_center_key("A1E3R_A302", {}) == "A1E3R_A302"


def _rcsb_available():
    try:
        r = requests.head(
            "https://files.rcsb.org/download/21ZQ.cif",
            timeout=15,
            allow_redirects=True,
        )
        return r.status_code == 200
    except requests.RequestException:
        return False


@pytest.mark.skipif(not _rcsb_available(), reason="RCSB unreachable; skipping 21ZQ integration test")
def test_fetch_21zq_and_cluster_with_a1e3r(tmpdir):
    """Live RCSB fetch of mmCIF-only 21ZQ, then cluster around original CCD A1E3R."""
    import glob

    from qp.cluster.spheres import CenterResidue, extract_clusters
    from qp.structure.mmcif_to_pdb import load_remap_sidecar

    pdb_path = os.path.join(tmpdir, "21ZQ", "21ZQ.pdb")
    try:
        status = setup.fetch_pdb("21ZQ", pdb_path)
    except (IOError, ValueError) as exc:
        pytest.skip(f"RCSB fetch failed ({exc}); skipping 21ZQ integration test")

    assert status == "mmcif"
    assert os.path.isfile(pdb_path)
    cif_path = os.path.join(tmpdir, "21ZQ", "21ZQ.cif")
    assert os.path.isfile(cif_path)

    # Compare conversion atom counts rather than a hard-coded golden number
    # (RCSB entries can be revised over time).
    n_cif = _cif_atom_count(cif_path)
    n_pdb = sum(
        1 for _ in PDBParser(QUIET=True).get_structure("21ZQ", pdb_path).get_atoms()
    )
    assert n_pdb == n_cif
    assert n_pdb > 0

    remap = load_remap_sidecar(pdb_path)
    assert remap["resname_map"].get("A1E3R") == "A1E"

    cluster_out = os.path.join(tmpdir, "clusters")
    os.makedirs(cluster_out, exist_ok=True)
    center = CenterResidue("A1E3R", resname_map=remap["resname_map"])
    cluster_paths = extract_clusters(
        pdb_path,
        cluster_out,
        center,
        sphere_count=2,
        first_sphere_radius=4.0,
        capping=0,
        charge=False,
        count=True,
        xyz=False,
        smooth_method="dummy_atom",
        mean_distance=3,
    )

    assert cluster_paths, "Expected at least one cluster centered on A1E3R/A1E"
    # A1E3R is residue 302 on chain A after conversion
    assert any(os.path.basename(p) == "A302" for p in cluster_paths)

    sphere0 = os.path.join(cluster_out, "A302", "0.pdb")
    assert os.path.isfile(sphere0) and os.path.getsize(sphere0) > 0
    with open(sphere0) as handle:
        sphere0_text = handle.read()
    assert "A1E" in sphere0_text
    assert glob.glob(os.path.join(cluster_out, "A302", "*.pdb"))
