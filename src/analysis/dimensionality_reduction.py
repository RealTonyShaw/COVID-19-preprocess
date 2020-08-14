import pandas as pd
import scipy as sc
from sklearn.decomposition import PCA
from scipy.stats import rankdata
import numpy as np


def reduce_dimension(file_path: str) -> pd.DataFrame:
    """
    Apply PCA to original redundant data.
    :param file_path: Path of materialized embedded sequences.
    :return: Simplified sequences.
    """
    covid_seq = pd.read_csv(file_path)
    pca = PCA(n_components=50) # the value of n_components is the same as the implementation of Tempel
    new_covid_seq = pca.fit_transform(covid_seq)
    return new_covid_seq


def connect_cluster_of_two_consecutive_years(prev_clusters_path: str, subs_clusters_path: str) -> dict:
    """
    Get the evolution path.
    :param prev_clusters_path: Path of materialized embedded sequences.
    :param subs_clusters_path: Path of materialized embedded sequences.
    :return: Evolution path.
    """
    prev_clusters = pd.read_csv(prev_clusters_path)
    subs_clusters = pd.read_csv(subs_clusters_path)
    distance_mat = sc.spatial.distance.cdist(prev_clusters, subs_clusters, metric='euclidean')
    evo_path = dict()
    for _ in range(distance_mat.shape[0]):
        ranked_data = rankdata(distance_mat[_], method='min') - 1
        evo_path[_] = np.where(ranked_data == 0)[0].tolist()
    return evo_path
