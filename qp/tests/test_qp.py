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
from qp.manager.charge_embedding import load_custom_charges, get_charges


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
