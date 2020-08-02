from Bio.Seq import Seq
from src.data_operation import get_data

# Paths
covid_amino_acids_path = "../../data/amino_acids.fasta"


def convert_dna_to_amino_acids(sequences_path: str, amino_acids_path: str):
    """
    Convert original DNA fasta file to amino acids fasta file.
    :param amino_acids_path: The path of amino acids fasta file to be generated.
    :param sequences_path: The path of a DNA fasta file.
    :return:
    """
    full_records = get_data.get_sequences(sequences_path)
    f = open(amino_acids_path, 'w')
    for record in full_records:
        f.write(">")
        f.write(record.id)
        f.write('\n')
        f.write(str(record.seq.translate()))
        f.write('\n')
    f.close()


def get_overlapping_trigrams_of_sequence(amino_acids_sequence: str) -> list:
    """
    Get the overlapping 3-grams according to the form proposed by Tempel.
    :param amino_acids_sequence: An amino acids sequence with length >= 3.
    :return: a list that contains all the 3-grams in the original sequence.
    """
    trigrams = []
    for _ in range(len(amino_acids_sequence) - 2):
        trigrams.append(amino_acids_sequence[_:_ + 3])
    return trigrams
