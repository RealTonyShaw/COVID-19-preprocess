from Bio import SeqIO
import pandas as pd

# Paths
covid_metadata_path = "../../data/metadata.tsv"
covid_sequences_path = "../../data/alignment.fasta"
covid_protvec_path = "../../data/protVec_100d_3grams.csv"


def get_metadata(metadata_path: str):
    """
    Get the metadata for different strains
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


def get_protvec(protvec_path: str):
    """
    Get protvec
    :param protvec_path: The path of a protvec file.
    :return: A dataframe object of protvec.
    """
    protvec = pd.read_csv(protvec_path, "	")
    return protvec
