from Bio import SeqIO
import pandas as pd
from src.amino_acids.make_amino_acids_sequence import get_overlapping_trigrams_of_sequence

# Paths
covid_metadata_path = "../../data/metadata.tsv"
covid_dna_sequences_path = "../../data/alignment.fasta"
covid_protvec_path = "../../data/protVec_100d_3grams.csv"


def get_metadata(metadata_path: str) -> pd.DataFrame:
    """
    Get the metadata for different strains
    :param metadata_path: The path of a metadata file.
    :return: A dataframe object of metadata.
    """
    metadata = pd.read_csv(metadata_path, "	")
    return metadata


def get_sequences(sequences_path: str) -> list:
    """
    Get a list of SeqRecords of a given fasta file.
    :param sequences_path: The path of a fasta file.
    :return: A list of SeqRecords.
    """
    sequences = list(SeqIO.parse(sequences_path, "fasta"))
    return sequences


def get_protvec(protvec_path: str) -> pd.DataFrame:
    """
    Get protVec
    :param protvec_path: The path of a protVec file.
    :return: A dataframe object of protVec.
    """
    protvec = pd.read_csv(protvec_path, "	")
    return protvec


def convert_sequence_to_protvec(amino_acids_sequence: str, protvec_path: str) -> pd.Series:
    """
    Embed a 3-grams sequence into the protVec.
    :param amino_acids_sequence: A 3-grams sequence of amino acids.
    :param protvec_path: The path of a protVec file.
    :return: A series object of the embedding protVec.
    """
    trigrams = get_overlapping_trigrams_of_sequence(amino_acids_sequence)
    protvec = get_protvec(protvec_path)
    embedded_trigrams = pd.DataFrame()
    for trigram in trigrams:
        if protvec[protvec['words'] == trigram].empty:
            # if embedded_trigrams.empty:
            #     embedded_trigrams = protvec[protvec['words'] == '<unk>'].iloc[:, 1:]
            # else:
            #     embedded_trigrams = embedded_trigrams + protvec[protvec['words'] == '<unk>'].iloc[:, 1:]
            embedded_trigrams = pd.concat([embedded_trigrams, protvec[protvec['words'] == '<unk>']])
        else:
            # if embedded_trigrams.empty:
            #     embedded_trigrams = protvec[protvec['words'] == trigram].iloc[:, 1:]
            # else:
            #     embedded_trigrams = embedded_trigrams + protvec[protvec['words'] == trigram].iloc[:, 1:]
            embedded_trigrams = pd.concat([embedded_trigrams, protvec[protvec['words'] == trigram]])
    return embedded_trigrams.iloc[:, 1:][:].sum().drop("words")


def convert_sequence_to_protvec_by_concat(amino_acids_sequence: str, protvec_path: str) -> pd.Series:
    """
    Embed a 3-grams sequence into the protVec by concat.
    :param amino_acids_sequence: A 3-grams sequence of amino acids.
    :param protvec_path: The path of a protVec file.
    :return: A series object of the embedding protVec.
    """
    trigrams = get_overlapping_trigrams_of_sequence(amino_acids_sequence)
    protvec = get_protvec(protvec_path)
    embedded_trigrams = pd.Series()
    for trigram in trigrams:
        if protvec[protvec['words'] == trigram].empty:
            embedded_trigrams = pd.concat([embedded_trigrams, protvec[protvec['words'] == '<unk>'].iloc[0]])
        else:
            embedded_trigrams = pd.concat([embedded_trigrams, protvec[protvec['words'] == trigram].iloc[0]])
    return embedded_trigrams.drop("words")


def embed_all_sequences(sequences_path: str, protvec_path: str) -> pd.DataFrame:
    """
    Embed all sequences into protVec.
    :param sequences_path: The path of a fasta file.
    :param protvec_path: The path of a protVec file.
    :return: A dataframe object of the embedding sequences.
    """
    sequences = get_sequences(sequences_path)
    # protvec = get_protvec(protvec_path)
    embedded_trigrams = pd.DataFrame()
    for sequence in sequences:
        # TODO: Fix concat bugs
        embedded_trigrams = \
            pd.concat([embedded_trigrams, convert_sequence_to_protvec(sequence.seq, protvec_path).to_frame().T])
            # pd.concat([embedded_trigrams, convert_sequence_to_protvec(sequence.seq, protvec_path).to_frame().T.iloc[:, 1:]])
        print(embedded_trigrams)
    return embedded_trigrams


def embed_all_sequences_by_concat(sequences_path: str, protvec_path: str) -> pd.DataFrame:
    """
    Embed all sequences into protVec by concat.
    :param sequences_path: The path of a fasta file.
    :param protvec_path: The path of a protVec file.
    :return: A dataframe object of the embedding sequences.
    """
    sequences = get_sequences(sequences_path)
    # protvec = get_protvec(protvec_path)
    embedded_trigrams = pd.DataFrame()
    # TODO: Delete _, it's only for test purpose.
    _ = 0
    for sequence in sequences:
        _ += 1
        if _ == 10:
            break
        # TODO: Fix concat bugs
        embedded_trigrams = \
            pd.concat([embedded_trigrams,
                       convert_sequence_to_protvec_by_concat(sequence.seq, protvec_path).to_frame().T])
        # pd.concat([embedded_trigrams, convert_sequence_to_protvec_by_concat(sequence.seq, protvec_path).to_frame().T.iloc[:, 1:]])
        print(embedded_trigrams)
    return embedded_trigrams
