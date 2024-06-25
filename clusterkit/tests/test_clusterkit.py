"""
Unit and regression test for the clusterkit package.
"""

# Import package, test suite, and other packages as needed
import clusterkit
import pytest
import sys
from MDAnalysis.tests.datafiles import PSF, DCD
import MDAnalysis as mda


def test_clusterkit_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "clusterkit" in sys.modules


def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"


def test_get_features_ramachandran():
    from clusterkit.ClusterAnalysis import ClusterAnalysis
    u = mda.Universe(PSF, DCD)
    C = ClusterAnalysis(u)
    C.featurize('ramachandran')
    assert C.data.shape == (98, 848)
