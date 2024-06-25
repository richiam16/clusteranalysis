"""
Unit and regression test for the clusterkit package.
"""

# Import package, test suite, and other packages as needed
import clusterkit
import pytest
import sys


def test_clusterkit_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "clusterkit" in sys.modules


def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"
