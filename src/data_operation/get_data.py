from Bio import SeqIO
import pandas as pd

# Paths
covid_metadata_path = "../../data/metadata.tsv"
covid_sequences_path = "../../data/alignment.fasta"


def get_metadata(metadata_path: str):
    """
    Get the metadata for different
    :param metadata_path: The path of a metadata file.
    :return: A dataframe object of metadata.
    """
    metadata = pd.read_csv(metadata_path, "	")
    return metadata


def get_sequences(sequences_path: str):
    """
    Get a list of SeqRecords of a given fasta file.
    :param sequences_path: The path of a fasta file.
    :return: A list of SeqRecords.
    """
    records = list(SeqIO.parse(sequences_path, "fasta"))
    return records
