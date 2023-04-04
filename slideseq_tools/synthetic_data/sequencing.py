"""
Creates synthetic Slide-seq sequencing data for testing.
"""

from typing import List, Tuple
from pathlib import Path
from random import randrange
import numpy as np
from Bio import SeqIO

from slideseq_tools.utils.constants import BASES, MUTATIONS
from slideseq_tools.gff import GFF


class Sequencing:
    """Synthetic Slide-seq sequencing data."""

    gff_path: Path = None
    record_dict = None
    features = None

    @classmethod
    def __init__(cls, gff_path: str, fasta_path: str, length: int = 50) -> None:
        """\
        Constructor for Sequencing Slide-seq class.

        Constructor raises:

            * `FileNoFoundError` if  `GFF` file doesn't exist
            * `FileNoFoundError` if  `FASTA` file doesn't exist

        Parameters
        ----------
        gff_path 
            Path of the `GFF` file.
        fasta_path 
            Path of the `FASTA` file.
        length
            Length of transcripts.
        """
        gff_path = Path(gff_path)
        if not gff_path.exists():
            raise FileNotFoundError(f"GFF image {gff_path} doesn't exist.")
        cls.gff_path = gff_path

        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA image {fasta_path} doesn't exist.")
        cls.fasta_path = fasta_path

        cls.length = length

    @classmethod
    def random_sequence(cls, length: int = 14) -> str:
        """\
        Returns a random DNA sequence as `str`.

        Parameters
        ----------
        length
            Sequence length.
        """
        min_val = min(BASES.keys())
        max_val = max(BASES.keys())

        bases = []

        for _ in range(length):
            bases.append(BASES[randrange(min_val, max_val + 1)])

        return "".join(bases)

    def _is_dna(self, sequence: str) -> bool:
        """Returns if a `str` is DNA."""
        bases = set(BASES.values())
        for base in sequence:
            if base not in bases:
                return False
        return True

    @classmethod
    def mutate(cls, sequence: str, n_bases: int = None) -> str:
        """\
        Randomly mutate a DNA sequence.

        Method raises a `ValueError` if sequence not made of DNA base.
        It mutates everything if `n_bases` is greater than sequence length.
        If `n_bases` not specified it mutates everything

        Parameters
        ----------
        sequence
            Sequence to mutate as `str`.
        n_bases
            Number of base to mutate.
        """
        if not cls._is_dna(cls, sequence):
            raise ValueError(f"{sequence} is not DNA.")

        # number of bases to mutate
        if n_bases is None:
            size = len(sequence)
        else:
            size = min(n_bases, len(sequence))

        positions = np.random.choice(a=len(sequence), size=size, replace=False)
        sequence = list(sequence)
        for pos in positions:
            sequence[pos] = MUTATIONS[sequence[pos]]

        return "".join(sequence)

    def get_transcripts(self, n_transcripts: int = 3) -> Tuple:
        """\
        Returns a list of transcripts as a tuples (`seqid`, `start`, `end`, `sequence`).

            * `seqid` is returned as `str`.
            * `start` and `end` are returned as `int`.
            * `sequence` is returned as an upper case `str`.

        Parameters
        ----------
        n_transcripts
            Number of transcripts to return.
        """
        if not self.features:
            gff = GFF(self.gff_path)
            self.features = gff.get_features(min_length=self.length)

        if not self.record_dict:
            self.record_dict = SeqIO.index(str(self.fasta_path.absolute()), "fasta")

        size = min(len(self.features), n_transcripts)
        indexes = np.random.choice(a=len(self.features), size=size, replace=False)

        transcripts = []

        for i in indexes:
            seqid, start, end = self.features[i]
            end = start + self.length
            transcript = self.record_dict[seqid][start:end]
            transcripts.append((seqid, start, end, str(transcript.seq).upper()))

        return transcripts

    @classmethod
    def generate_q_score_string(cls, n_bases: int = 50) -> str:
        """\
        Returns Q-Score string.

        Parameters
        ----------
        n_bases
            Number of bases.
        """
        ascii_codes = list(range(33, 74))
        codes = np.random.choice(a=ascii_codes[-10:], size=n_bases)
        return "".join(list(map(chr, codes)))

    @classmethod
    def get_phred_scores(cls, n_bases: int = 50) -> List:
        """\
        Returns Solexa quality score list.

        Parameters
        ----------
        n_bases
            Number of bases.
        """
        scores = list(range(41))
        return np.random.choice(a=scores[-10:], size=n_bases)
