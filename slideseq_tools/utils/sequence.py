"""
Functions to manage DNA sequences.
"""


def hamming(seq1: str, seq2: str) -> int:
    """\
    Returns Hamming distance between to sequences.

    Function raises:

        * `TypeError` if sequences are not `str`.
        * `ValueError` if sequences don't have same length.

    Parameters
    ----------
    seq1
        First sequence as `str`.
    seq2
        Second sequence as `str`.
    """
    if not isinstance(seq1, str):
        raise TypeError(f"{seq1} is not a str.")

    if not isinstance(seq2, str):
        raise TypeError(f"{seq2} is not a str.")

    if len(seq1) != len(seq2):
        raise ValueError(f"{seq1} and {seq2} don't have same length")

    distance = 0

    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            distance += 1

    return distance
