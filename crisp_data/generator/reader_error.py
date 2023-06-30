import enum
import random

import numpy
import typing

TRIMERS = ["GGT", "CGT", "AGT", "CGA", "CCA", "CCT", "GCT", "ACT", "GGA", "AAT"]
LETTERS = ["A", "T", "G", "C"]


class Reader:
    """
    A Class to represent Readers
    """

    def __init__(
        self,
        name: str,
        error_distribution_parameters: typing.Tuple[float, float],
        post_homopolymer_error_frequencies: typing.List[typing.List],
        error_prone_trimer_frequencies: typing.List,
    ):
        self.name = name
        self.error_distribution_parameters = error_distribution_parameters

        # Creating fault-prone trimer error dictionary
        self.error_prone_trimer_frequencies = {
            TRIMERS[i]: error_prone_trimer_frequencies[i] for i in range(len(TRIMERS))
        }

        # Creating post-homopolymer error dictionary - min homopolymer length = 2
        max_homopolymer_length = len(post_homopolymer_error_frequencies[0])
        error_frequencies = {}
        for i in range(len(LETTERS)):
            for j in range(max_homopolymer_length):
                error_frequencies[
                    "".join([LETTERS[i]] * (j + 2))
                ] = post_homopolymer_error_frequencies[i][j]
        self.post_homopolymer_error_frequencies = error_frequencies

    def random_error(self) -> float:
        """
        Generates a random error specific for a Reader as a rate of fault bp to the read overlap bp.
        """
        return abs(
            numpy.random.normal(
                loc=self.error_distribution_parameters[0],
                scale=self.error_distribution_parameters[1],
            )
        )


class Readers:
    """
    A class containing information on the error distributions of Particular Illumina Readers;
    based on https://doi.org/10.1093/nargab/lqab019
    """

    MiSeq: Reader = Reader(
        name="MiSeq",
        error_distribution_parameters=(0.00473, 0.00938),
        error_prone_trimer_frequencies=[
            1.46,
            1.33,
            1.31,
            1.35,
            1.29,
            1.23,
            1.27,
            1.24,
            1.30,
            1.80,
        ],
        post_homopolymer_error_frequencies=[
            [1.4, 1.3, 1.5, 1.4, 1.4, 1.5, 2.0],
            [1.6, 1.65, 1.8, 1.7, 1.9, 2.3, 2.6],
            [2.5, 3.1, 3.6, 3.3, 3.5, 1, 1],
            [2.2, 2.4, 2.5, 5.5, 4.2, 1, 1],
        ],
    )

    MiniSeq: Reader = Reader(
        name="MiniSeq",
        error_distribution_parameters=(0.00613, 0.00459),
        error_prone_trimer_frequencies=[
            1.18,
            1.19,
            1.13,
            1.20,
            1.21,
            1.18,
            1.20,
            1.17,
            1.14,
            1.08,
        ],
        post_homopolymer_error_frequencies=[
            [1.5, 1.4, 1.6, 1.6, 1.4, 2, 2],
            [1.6, 1.6, 1.7, 1.8, 1.8, 2, 2],
            [2.6, 2.7, 2.6, 2.65, 2.7, 1, 1],
            [2.5, 2.3, 3.0, 3.3, 4.0, 1, 1],
        ],
    )

    NextSeq500: Reader = Reader(
        name="NextSeq500",
        error_distribution_parameters=(0.00429, 0.00827),
        error_prone_trimer_frequencies=[
            1.45,
            1.32,
            1.20,
            1.18,
            1.12,
            1.05,
            1.10,
            0.75,
            0.80,
            0.90,
        ],
        post_homopolymer_error_frequencies=[
            [1, 1.1, 1.1, 1, 1.1, 1, 0.9],
            [1.8, 1.85, 1.85, 1.60, 3.2, 1.5, 1],
            [1, 2.3, 1.65, 2.2, 2.3, 1, 1],
            [1, 1.1, 1.2, 1, 5.4, 1.1, 1, 1],
        ],
    )
    NextSeq550: Reader = Reader(
        name="NextSeq550",
        error_distribution_parameters=(0.00593, 0.00435),
        error_prone_trimer_frequencies=[
            1.58,
            1.34,
            1.30,
            1.23,
            1.16,
            1.05,
            1.10,
            1.05,
            1.10,
            0.80,
        ],
        post_homopolymer_error_frequencies=[
            [1, 1.1, 1.1, 1, 1.1, 1, 0.9],
            [1.8, 1.85, 1.85, 1.60, 3.2, 1.5, 1],
            [1, 2.3, 1.65, 2.2, 2.3, 1, 1],
            [1, 1.1, 1.2, 1, 5.4, 1.1, 1, 1],
        ],
    )
    HiSeq2500: Reader = Reader(
        name="HiSeq2500",
        error_distribution_parameters=(0.00112, 0.00544),
        error_prone_trimer_frequencies=[
            1.75,
            1.39,
            1.45,
            1.33,
            1.32,
            1.17,
            1.21,
            1.16,
            1.30,
            1.10,
        ],
        post_homopolymer_error_frequencies=[
            [1.55, 2.1, 1.55, 2.0, 1.80, 1.90, 2.50],
            [1.5, 1.7, 1.9, 1.80, 3.2, 2.5, 3.0],
            [3.2, 8.5, 6.2, 12.1, 1, 1, 1],
            [2.1, 2.6, 3.6, 8.8, 6.4, 1, 1],
        ],
    )
    NovaSeq6000: Reader = Reader(
        name="NovaSeq6000",
        error_distribution_parameters=(0.00109, 0.00350),
        error_prone_trimer_frequencies=[
            1.07,
            1.30,
            1.18,
            1.06,
            1.16,
            1.23,
            1.05,
            1.16,
            0.50,
            1.18,
        ],
        post_homopolymer_error_frequencies=[
            [2.1, 2.0, 2.2, 2.1, 2.1, 2.5, 2.7],
            [1.7, 1.8, 1.9, 2.0, 2.1, 2.5, 4.3],
            [0.8, 0.8, 0.9, 0.8, 0.9, 0.9, 1],
            [1, 1, 1.3, 1.6, 1.6, 2.5, 1],
        ],
    )
    HiSeqXTen: Reader = Reader(
        name="HiSeqXTen",
        error_distribution_parameters=(0.00087, 0.00126),
        error_prone_trimer_frequencies=[
            1.43,
            1.30,
            1.30,
            1.18,
            1.18,
            1.12,
            1.14,
            1.10,
            1.10,
            1.16,
        ],
        post_homopolymer_error_frequencies=[
            [1.8, 2.1, 1.8, 1.8, 2.0, 1.9, 1.2],
            [2.0, 2.2, 2.5, 1.8, 1.9, 1.8, 1.8],
            [2.5, 3.55, 5.1, 3.9, 8.4, 2, 2],
            [2.2, 2.2, 4.9, 2.0, 2.4, 2, 2],
        ],
    )

    @classmethod
    def list(cls) -> list:
        readers_list = []
        for reader in dir(cls):
            if reader[0:2] != "__":
                readers_list.append(reader)
            else:
                break
        return readers_list


class ReadersEnum(enum.Enum):

    """
    This is the enumerator class to parse the reader field in the config.json.
    """

    MiSeq = Readers.MiSeq
    MiniSeq = Readers.MiniSeq
    HiSeq2500 = Readers.HiSeq2500
    HiSeqXTen = Readers.HiSeqXTen
    NextSeq500 = Readers.NextSeq500
    NextSeq550 = Readers.NextSeq550
    NovaSeq6000 = Readers.NovaSeq6000


def introduce_reader_error(
    sequence: str,
    reader: typing.Optional[Reader] = Readers.MiSeq,
    read_length: typing.Optional[int] = None,
    direction: typing.Optional[bool] = True,
) -> (str, typing.List):
    """
    Introduces a Reader-specific error to a single Read of a particular sequence
    based on https://doi.org/10.1093/nargab/lqab019;
    Each error is independent.
    Includes errors in characteristic trimers, random and post-homopolymer errors
    param sequence: Input sequence as a string
    param read_length: Read length of the NGS, is int and frequently equals 600, 300, 150 or 75,
                       set to length of the sequence by default.
    param reader: The Model of the Reader to emulate
    return: An emulated Read of the sequence as a str with a Reader-specific error and list of tuples,
            characterising each error:
            reported_errors [i] = (position:int, original_value:str, new_value:str, type:str)
    """
    sequence_length = len(sequence)

    if read_length is None:
        read_length = sequence_length

    if read_length > sequence_length:
        raise ValueError("Read Length can't be more than Sequence length.")
    if sequence_length >= 2 * read_length:
        raise ValueError("Doubled Read Length should be more than Sequence length.")

    read_overlap = sequence_length - 2 * (sequence_length - read_length)

    if not direction:
        sequence = sequence[::-1]

    total_trimers = sequence_length - 2
    trimer_map = {}
    trimer_probabilities = {}
    homopolymer_map = {}
    post_homopolymer_probabilities = {}

    reported_errors = []

    # Calculating number of fault bases from the given Reader's parameters
    fault_bp_number = int(reader.random_error() * read_overlap)

    # Locate start of error-prone trimer sequences and adjust error location probability
    for trimer in TRIMERS:
        trimer_map[trimer] = [
            i for i in range(sequence_length) if sequence.startswith(trimer, i)
        ]
        trimer_probabilities[trimer] = (
            len(trimer_map[trimer])
            / total_trimers
            * reader.error_prone_trimer_frequencies[trimer]
        )

    # Locate start of homopolymer sequences and adjust error location probability
    for homopolymer in list(reader.post_homopolymer_error_frequencies.keys()):
        homopolymer_map[homopolymer] = [
            i
            for i in range(len(sequence))
            if (
                sequence.startswith(homopolymer, i)
                and (not sequence.endswith(homopolymer[0], i, i + len(homopolymer) + 1))
                and sequence[i - 1] != homopolymer[0]
            )
        ]

        post_homopolymer_probabilities[homopolymer] = (
            len(homopolymer_map[homopolymer])
            / (sequence_length - len(homopolymer))
            * reader.post_homopolymer_error_frequencies[homopolymer]
        )

    # Unify trimer and post-homopolymer error probabilities
    trimers_index = len(list(trimer_probabilities.keys()))
    total_probabilities = {**trimer_probabilities, **post_homopolymer_probabilities}
    sequence_map = {**trimer_map, **homopolymer_map}

    # Add Random option as error location
    total_probabilities["RANDOM"] = 1 - sum(list(total_probabilities.values()))

    # Locate error sites and introduce Reader-specific errors to the sequence

    for i in range(fault_bp_number):
        error_key_position = list(
            numpy.random.multinomial(n=1, pvals=list(total_probabilities.values()))
        ).index(1)
        error_key = list(total_probabilities.keys())[error_key_position]
        if error_key != "RANDOM":
            base_position = random.choice(sequence_map[error_key])

            if error_key_position >= trimers_index:
                position = base_position + len(error_key)
                reported_errors.append(
                    (position, sequence[position], error_key[0], "POST-HOMOPOLYMER")
                )

                # Replace base at defined position with a base of a corresponding homopolymer
                sequence = sequence[:position] + error_key[0] + sequence[position + 1:]

            else:
                position = random.randint(base_position, base_position + 3)
                new_value = sequence[position - 1]
                if new_value == sequence[position]:
                    new_value = random.choice(LETTERS)
                reported_errors.append(
                    (position, sequence[position], new_value, "TRIMER")
                )

                # Replace base at defined position according to substitutions' characteristic to a Reader
                sequence = sequence[:position] + new_value + sequence[position + 1:]

        else:
            position = random.randint(0, sequence_length - 1)
            new_value = sequence[position - 1]
            if new_value == sequence[position]:
                new_value = random.choice(LETTERS)
            reported_errors.append((position, sequence[position], new_value, "RANDOM"))
            # Replace base at defined position according to substitutions' characteristic to a Reader
            sequence = sequence[:position] + new_value + sequence[position + 1:]

    return sequence, reported_errors
