import pandas as pd
from sklearn.decomposition import PCA


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

