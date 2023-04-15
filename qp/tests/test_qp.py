"""
Unit and regression test for the qp package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import qp


def test_qp_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "qp" in sys.modules
