import pytest
import os


@pytest.fixture
def sample_cluster(request):
    pdb, metal_ids = request.param
    return pdb, metal_ids, os.path.join(os.path.dirname(__file__), "samples", pdb)


@pytest.fixture
def sample_pdb(request):
    pdb = request.param
    return pdb, os.path.join(os.path.dirname(__file__), "samples", pdb)
