"""
Unit and regression test for the quantumpdb package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import quantumpdb


def test_quantumpdb_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "quantumpdb" in sys.modules
