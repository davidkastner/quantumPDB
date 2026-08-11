"""Helpers for live ProteinsPlus / Protoss integration tests."""

import pytest
import requests


def proteins_plus_reachable(timeout=15):
    """Return True if the ProteinsPlus API responds with a non-5xx status."""
    try:
        r = requests.get("https://proteins.plus/api/v2/", timeout=timeout)
        return r.status_code < 500
    except requests.RequestException:
        return False


# Exceptions from get_protoss.upload / submit / download that indicate the
# ProteinsPlus service is unreachable or the request failed after retries.
# Do not include ValueError: that is also raised for Protoss job/business
# failures (e.g. non-success status), which should fail the test.
PROTOSS_NETWORK_ERRORS = (
    KeyError,
    requests.RequestException,
    TimeoutError,
    OSError,
)


def skip_if_protoss_unavailable(exc):
    """Convert a Protoss upload/API failure into pytest.skip."""
    pytest.skip(
        f"Protoss unavailable ({type(exc).__name__}: {exc}); skipping live Protoss test"
    )
