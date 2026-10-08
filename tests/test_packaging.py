"""Packaged browser assets are complete and contain no stale compiled bundles."""

import re
from importlib.resources import files

import pytest


def test_browser_package_matches_exactly_the_current_asset_manifest():
    web = files("dna_design").joinpath("web")
    if not web.joinpath("index.html").is_file():
        pytest.skip("Browser assets have not been built in this source checkout")
    page = web.joinpath("index.html").read_text()
    referenced = set(re.findall(r'(?:src|href)="(/assets/[^"]+)"', page))
    packaged = {
        "/assets/" + path.name for path in web.joinpath("assets").iterdir() if path.is_file()
    }
    assert referenced == packaged, "Package contains missing or obsolete browser assets"
    assert web.joinpath("THIRD_PARTY_NOTICES.txt").is_file()
