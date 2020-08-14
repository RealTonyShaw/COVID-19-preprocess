from src.amino_acids import make_amino_acids_sequence
from src.analysis import aggregation
from src.analysis import dimensionality_reduction
from src.data_operation import get_data
from src.paths import file_paths

make_amino_acids_sequence.convert_dna_to_amino_acids(file_paths.covid_dna_sequences_path,
                                                     file_paths.covid_amino_acids_path,
                                                     21562,
                                                     25384)
aggregation.save_trigrams(file_paths.covid_amino_acids_path,
                          file_paths.covid_protvec_path,
                          file_paths.covid_overlapping_trigrams_path)
aggregation.kmeans_of_trigrams_sequences(file_paths.covid_overlapping_trigrams_path)
aggregation.aggregate_sequences(file_paths.covid_overlapping_trigrams_path, 3)
dimensionality_reduction.reduce_dimension(file_paths.covid_overlapping_trigrams_path)
# TODO: Add dimensionality_reduction.connect_cluster_of_two_consecutive_years
