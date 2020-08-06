import pandas as pd
import src.data_operation.get_data as gd
import src.amino_acids.make_amino_acids_sequence as ma
from sklearn.cluster import KMeans

covid_overlapping_trigrams_path = "../../data/overlapping.csv"

def save_trigrams(amino_acids_sequence: str, protvec_path: str, file_path: str):
    # gd.embed_all_sequences_by_concat(ma.covid_amino_acids_path, gd.covid_protvec_path)
    embedded_trigrams = gd.embed_all_sequences_by_concat(amino_acids_sequence, protvec_path)
    embedded_trigrams.to_csv(file_path, index=False)


def kmeans_of_trigrams_sequences(file_path: str):
    covid_seq = pd.read_csv(file_path)
    kmeans = KMeans(n_clusters=3).fit(covid_seq)
    centroids = kmeans.cluster_centers_
    print(centroids)


# save_trigrams(ma.covid_amino_acids_path, gd.covid_protvec_path, covid_overlapping_trigrams_path)
kmeans_of_trigrams_sequences(covid_overlapping_trigrams_path)