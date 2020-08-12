import pandas as pd
import src.data_operation.get_data as gd
import src.amino_acids.make_amino_acids_sequence as ma
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt


def save_trigrams(amino_acids_sequence: str, protvec_path: str, file_path: str):
    """
    Materialize embedding results to file.
    :param amino_acids_sequence: Path of 3-grams sequences of amino acids.
    :param protvec_path: Path of protVec file.
    :param file_path: Output path.
    :return:
    """
    # gd.embed_all_sequences_by_concat(ma.covid_amino_acids_path, gd.covid_protvec_path)
    embedded_trigrams = gd.embed_all_sequences_by_concat(amino_acids_sequence, protvec_path)
    embedded_trigrams.to_csv(file_path, index=False)


def kmeans_of_trigrams_sequences(file_path: str):
    """
    Visualize kmeans result
    :param file_path: Path of materialized embedded sequences.
    :return:
    """
    covid_seq = pd.read_csv(file_path)
    sse = {}
    for k in range(1, 5):
        kmeans = KMeans(n_clusters=k, max_iter=1000).fit(covid_seq)
        # Sum of distances of samples to their closest cluster center
        sse[k] = kmeans.inertia_
    plt.figure()
    plt.plot(list(sse.keys()), list(sse.values()))
    plt.xlabel("Number of cluster")
    plt.ylabel("SSE")
    plt.show()


def aggregate_sequences(file_path: str, aggregation_num: int):
    """
    Print kmeans centers
    :param file_path: Path of materialized embedded sequences.
    :param aggregation_num: The number of aggregations.
    :return:
    """
    covid_seq = pd.read_csv(file_path)
    kmeans = KMeans(n_clusters=aggregation_num, max_iter=1000).fit(covid_seq)
    print(kmeans.cluster_centers_)


