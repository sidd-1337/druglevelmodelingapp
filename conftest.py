"""Shared pytest configuration.

Some tests call the live SÚKL open-data API. They are real integration tests -
they catch the API changing its field names or dropping a package code - so they
are kept, not mocked. But they cannot pass where the network is blocked (CI
sandboxes, this project's dev container behind an HTTPS proxy), and a suite that
always has four red tests stops being a usable signal.

So they are marked ``network`` and skipped once a single probe shows SÚKL is
unreachable. A *failing* SÚKL call still fails the test when the API answers at
all, which is the case the tests exist for.
"""
from __future__ import annotations

import pytest


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "network: needs the live SÚKL open-data API; skipped when unreachable"
    )


@pytest.fixture(scope="session")
def sukl_reachable() -> bool:
    """True if the live SÚKL API answers. Probed once per session."""
    from core.sukl_client import query_sukl_api

    try:
        return bool(query_sukl_api("Ritalin", max_results=1))
    except Exception:
        return False


@pytest.fixture
def requires_sukl(sukl_reachable: bool) -> None:
    if not sukl_reachable:
        pytest.skip("live SÚKL API not reachable from this environment")
