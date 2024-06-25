from MDAnalysis.analysis.base import AnalysisBase
from sklearn import cluster
import numpy as np
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.analysis import distances


class Featurizer(object):
    def __init__(self, feature):
        """Featurizer class for extracting MD features with MDAnalysis
        Parameters:
        -----------
        feature : String
            name of feature
        """
        if feature == "ramachandran":
            self.get_features = self._get_features_ramachandran
        elif feature == "distance matrix":
            self.get_features = self._get_features_distance_matrix
        elif feature == "coordinates":
            self.get_features = self._get_features_coordinates
        elif feature == "custom":
            self.get_features = self._get_features_custom
        else:
            raise Exception(
                "Features extraction method %s not recognised" % feature)

        self.feature = feature

    def _get_features_coordinates(self, universe):
        """
        alpha carbons coordinates
        """
        crds = []
        ca = universe.select_atoms("name CA")
        for ts in universe.trajectory:
            crds.append(ca.positions.flatten())

        return np.array(crds)

    def _get_features_ramachandran(self, universe):
        """
        dihedral angles
        """
        r = Ramachandran(universe.select_atoms('protein')).run()
        r_sin = np.sin(np.deg2rad(r.angles))
        r_cos = np.cos(np.deg2rad(r.angles))

        r_sin = r_sin.reshape((r_sin.shape[0], np.prod(r_sin.shape[1:])))
        r_cos = r_cos.reshape((r_cos.shape[0], np.prod(r_cos.shape[1:])))
        return np.concatenate((r_sin, r_cos), axis=1)

    def _get_features_distance_matrix(self, universe):
        """
        returns the distance matrix of each conformation (lower diagonal, flattened)
        see MDAnalysis.analysis.distances.self_distance_array
        """
        pass

    def _get_features_custom(self, universe):
        pass


class ClusterAnalysis(AnalysisBase):
    def __init__(self, data, method):
        super(ClusterAnalysis, self).__init__()
        self.data = data
        self.method = method

    def cluster(self):
        Cluster = getattr(cluster, self.method)

