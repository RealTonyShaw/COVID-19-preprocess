import random
from src.data_operation import get_data
from Bio import SeqIO


def convert_dna_to_amino_acids(sequences_path: str, amino_acids_path: str, start_point: int, end_point: int):
    """
    Convert original DNA fasta file to amino acids fasta file.
    :param amino_acids_path: The path of amino acids fasta file to be generated.
    :param sequences_path: The path of a DNA fasta file.
    :param start_point: Indicate the start point of subsequence.
    :param end_point: Indicate the end point of subsequence.
    :return:
    """
    full_records = get_data.get_sequences(sequences_path)
    f = open(amino_acids_path, 'w')
    for record in full_records:
        # TODO: Let date be a param.
        if get_date(record, get_data.covid_metadata_path) == 2:
            f.write(">")
            f.write(record.id)
            f.write('\n')
            f.write(random_replace_ambiguous_acids(str(record.seq[start_point:end_point].translate()))[0:-1])
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


def random_replace_ambiguous_acids(amino_acids_sequence: str) -> str:
    """
    Randomly selects replacements for all uncertain amino acids.
    :param amino_acids_sequence: An amino acids sequence which contains uncertain amino acids.
    :return: A certain amino acids sequence.
    """
    replacements = {'B': 'DN',
                    'J': 'IL',
                    'Z': 'EQ',
                    'X': 'ACDEFGHIKLMNPQRSTVWY'}
    # '*': 'ACDEFGHIKLMNPQRSTVWY'} '*' is the stop symbol.
    certain_sequence = amino_acids_sequence
    for uncertain in replacements.keys():
        for _ in range(len(amino_acids_sequence)):
            if amino_acids_sequence[_] in replacements.keys():
                certain_sequence = certain_sequence[:_] + random.choice(replacements[amino_acids_sequence[_]]) + \
                                   certain_sequence[_ + 1:]
        # amino_acids_sequence.replace() would replace each uncertain amino acid to the same amino acid,
        # which is implemented below.
        # amino_acids_sequence = amino_acids_sequence.replace(uncertain, random.choice(replacements[uncertain]))
    return certain_sequence


def get_date(sequence_record: SeqIO.SeqRecord, metadata_path: str) -> int:
    """
    Get the collected month of a given virus sequence.
    :param sequence_record: A sequence record.
    :param metadata_path: The path of a metadata file.
    :return: Month.
    """
    metadata = get_data.get_metadata(metadata_path)
    # TODO: Specify except: Date may not exist in metadata.
    try:
        date = metadata[metadata['strain'] == str(sequence_record.id)].iloc[0, 2]
    except:
        return 0
    date_month = date.split("-")
    # Test if the month data is not accurate.
    if date_month[1] == "XX":
        return 0
    return int(date_month[1])
