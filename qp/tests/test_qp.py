"""
Unit and regression test for the qp package.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import json
import tempfile
from unittest import mock

import pytest

import qp
from qp.manager.charge_embedding import load_custom_charges, get_charges, parse_pdb_to_xyz


def test_qp_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "qp" in sys.modules


def test_load_custom_charges_valid_json():
    """Test that load_custom_charges correctly loads a valid JSON file."""
    charges = {
        "ALA": {"N": -0.4157, "H": 0.2719},
        "GLY": {"N": -0.4157}
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(charges, f)
        tmppath = f.name
    try:
        result = load_custom_charges(tmppath)
        assert result == charges
        assert result["ALA"]["N"] == -0.4157
        assert result["ALA"]["H"] == 0.2719
        assert result["GLY"]["N"] == -0.4157
    finally:
        os.unlink(tmppath)


def test_load_custom_charges_file_not_found():
    """Test that load_custom_charges raises FileNotFoundError for missing files."""
    with pytest.raises(FileNotFoundError):
        load_custom_charges("/nonexistent/path.json")


def test_load_custom_charges_invalid_json():
    """Test that load_custom_charges raises an error for invalid JSON."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        f.write("not valid json {{")
        tmppath = f.name
    try:
        with pytest.raises(json.JSONDecodeError):
            load_custom_charges(tmppath)
    finally:
        os.unlink(tmppath)


def test_get_charges_default_uses_ff14sb():
    """Test that get_charges uses ff14SB dict when no custom file is provided."""
    with mock.patch('qp.manager.charge_embedding.ff14SB_dict.get_ff14SB_dict') as mock_ff, \
         mock.patch('qp.manager.charge_embedding.rename_and_clean_resnames'), \
         mock.patch('qp.manager.charge_embedding.parse_pdb'), \
         mock.patch('qp.manager.charge_embedding.remove_qm_atoms'), \
         mock.patch('qp.manager.charge_embedding.read_xyz', return_value=__import__('numpy').array([[0, 0, 0]])), \
         mock.patch('qp.manager.charge_embedding.parse_pdb_to_xyz'), \
         mock.patch('os.path.exists', return_value=False), \
         mock.patch('os.mkdir'), \
         mock.patch('os.getcwd', return_value='/fake/output/pdb1/A200/method'), \
         mock.patch('shutil.rmtree'):
        mock_ff.return_value = {"ALA": {"N": -0.4157}}
        get_charges(20)
        mock_ff.assert_called_once()


def test_get_charges_custom_file_used():
    """Test that get_charges uses a custom charges file when provided."""
    charges = {"ALA": {"N": -0.4157, "H": 0.2719}}
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(charges, f)
        tmppath = f.name

    _real_exists = os.path.exists

    def _exists_side_effect(path):
        """Allow the custom charges file to exist, mock others as absent."""
        if path == tmppath:
            return _real_exists(path)
        return False

    try:
        with mock.patch('qp.manager.charge_embedding.ff14SB_dict.get_ff14SB_dict') as mock_ff, \
             mock.patch('qp.manager.charge_embedding.rename_and_clean_resnames'), \
             mock.patch('qp.manager.charge_embedding.parse_pdb'), \
             mock.patch('qp.manager.charge_embedding.remove_qm_atoms'), \
             mock.patch('qp.manager.charge_embedding.read_xyz', return_value=__import__('numpy').array([[0, 0, 0]])), \
             mock.patch('qp.manager.charge_embedding.parse_pdb_to_xyz'), \
             mock.patch('os.path.exists', side_effect=_exists_side_effect), \
             mock.patch('os.mkdir'), \
             mock.patch('os.getcwd', return_value='/fake/output/pdb1/A200/method'), \
             mock.patch('shutil.rmtree'):
            get_charges(20, charge_embedding_charges=tmppath)
            mock_ff.assert_not_called()
    finally:
        os.unlink(tmppath)


def test_parse_pdb_to_xyz_residue_based_selection():
    """Test that parse_pdb_to_xyz includes complete residues.

    If any atom of a residue is within the cutoff distance, all atoms
    of that residue must be included in the output.  Residues entirely
    outside the cutoff must be excluded.
    """
    import numpy as np

    def make_pdb_line(serial, atom_name, res_name, chain, res_seq, x, y, z, charge):
        """Build a fixed-width PDB ATOM line."""
        return (
            f"ATOM  {serial:>5d} {atom_name:<4s} {res_name:>3s} {chain}{res_seq:>4d}"
            f"    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"  1.00{charge:>6.2f}           \n"
        )

    # ALA atom 1 (N)  at (0,0,0)  -> distance 0  <= 10 -> triggers residue inclusion
    # ALA atom 2 (CA) at (12,0,0) -> distance 12 > 10  -> excluded by old code, included by new
    # GLY atom   (N)  at (50,0,0) -> distance 50 > 10  -> residue excluded entirely
    pdb_content = (
        make_pdb_line(1, "N",  "ALA", "A", 1,  0.0, 0.0, 0.0, -0.42)
        + make_pdb_line(2, "CA", "ALA", "A", 1, 12.0, 0.0, 0.0,  0.42)
        + make_pdb_line(3, "N",  "GLY", "A", 2, 50.0, 0.0, 0.0, -0.30)
    )

    qm_centroid = np.array([0.0, 0.0, 0.0])
    cutoff = 10.0

    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as pdb_f:
        pdb_f.write(pdb_content)
        pdb_path = pdb_f.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as out_f:
        out_path = out_f.name

    try:
        parse_pdb_to_xyz(pdb_path, out_path, qm_centroid, cutoff)

        with open(out_path, 'r') as f:
            result_lines = f.readlines()

        # Header: atom count and comment
        assert result_lines[0].strip() == "2", (
            "Expected 2 atoms (full ALA residue), got: " + result_lines[0].strip()
        )
        assert result_lines[1].strip() == "Generated from PDB file"

        # Both ALA atoms present, GLY excluded
        data_lines = result_lines[2:]
        assert len(data_lines) == 2, f"Expected 2 data lines, got {len(data_lines)}"

        # First atom: charge=-0.42, coords=(0,0,0)
        parts_0 = data_lines[0].split()
        assert float(parts_0[0]) == pytest.approx(-0.42)
        assert float(parts_0[1]) == pytest.approx(0.0)

        # Second atom: charge=0.42, coords=(12,0,0)
        parts_1 = data_lines[1].split()
        assert float(parts_1[0]) == pytest.approx(0.42)
        assert float(parts_1[1]) == pytest.approx(12.0)
    finally:
        os.unlink(pdb_path)
        os.unlink(out_path)
