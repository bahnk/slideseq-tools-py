"""
Creates synthetic Slide-seq data for testing.
"""

import gzip
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from slideseq_tools.utils.constants import UP_PRIMER
from slideseq_tools.synthetic_data.spatial import Puck
from slideseq_tools.synthetic_data.sequencing import Sequencing


# pylint: disable=too-many-arguments
class SlideSeq:
    """Synthetic Slide-seq data."""

    tiff_path: Path = None
    gff_path: Path = None
    fasta_path: Path = None
    length: int = None
    n_beads: int = None
    puck = None

    def __init__(
        self,
        tiff_path: str,
        gff_path: str,
        fasta_path: str,
        length: int = 50,
        n_beads: int = int(8 * 1e4),
    ) -> None:
        """\
        Constructor for Slide-seq class.

        Constructor raises:

            * `FileNoFoundError` if  `TIFF` file doesn't exist
            * `FileNoFoundError` if  `GFF` file doesn't exist
            * `FileNoFoundError` if  `FASTA` file doesn't exist

        Parameters
        ----------
        tiff_path:
            Path of the `TIFF` file.
        gff_path 
            Path of the `GFF` file.
        fasta_path 
            Path of the `FASTA` file.
        length
            Length of transcripts.
        n_beads
            Number of beads.
        """
        tiff_path = Path(tiff_path)
        if not tiff_path.exists():
            raise FileNotFoundError(f"tiff image {tiff_path} doesn't exist.")
        self.tiff_path = tiff_path

        gff_path = Path(gff_path)
        if not gff_path.exists():
            raise FileNotFoundError(f"GFF image {gff_path} doesn't exist.")
        self.gff_path = gff_path

        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA image {fasta_path} doesn't exist.")
        self.fasta_path = fasta_path

        self.length = length
        self.n_beads = n_beads
        self.seq = Sequencing(gff_path=gff_path, fasta_path=fasta_path, length=length)

    def generate_puck(self, barcode_length=14) -> None:
        """\
        Generates bead barcodes and coordinates.

        Parameters
        ----------
        barcode_length
            Barcodes sequence lenght.
        """
        dframe = Puck.coordinates(self.tiff_path)
        dframe = dframe.sample(min(dframe.shape[0], self.n_beads))

        barcodes = []

        for _ in range(self.n_beads):
            barcodes.append(Sequencing.random_sequence(length=barcode_length))

        dframe.index = pd.Index(data=barcodes, name="Barcode")
        dframe = dframe.reset_index()

        self.puck = dframe

    def save_coordinates(self, path: str) -> None:
        """\
        Saves beads coordinates.

        Parameters
        ----------
        path
            `CSV` file path.
        """
        self.puck.to_csv(path, header=False, index=False, float_format="%.4f")

    def _randint(self, max_value: int = 10):
        """\
        Returns a random integer between 0 and `max_value`.

        Parameters
        ----------
        max_value
            Upper bound of the interval.
        """
        return np.random.randint(max_value, size=1)[0]

    def generate_reads(self, prefix: str = "sample", n_reads: int = 10) -> Tuple:
        """\
        Returns a tuple of Read 1 and Read 2 as `SeqRecord`.

        Parameters
        ----------
        n_reads
            Number of reads to return.
        """
        if not self.puck is None:
            self.generate_puck()

        reads1 = []
        reads2 = []

        transcripts = self.seq.get_transcripts(n_transcripts=n_reads)

        for counter in range(n_reads):

            # read 1
            barcode = self.puck.sample(1).Barcode.values[0]
            barcode = Sequencing.mutate(barcode, self._randint(3))
            umi = Sequencing.random_sequence(9)
            up_primer = Sequencing.mutate(UP_PRIMER, self._randint(3))
            read1 = barcode[:8] + up_primer + barcode[8:] + "TC" + umi
            if 0 == np.random.randint(21, size=1)[0]:
                read1 = read1[: np.random.randint(len(read1), size=1)[0]]
            read1 = SeqRecord(
                seq=Seq(read1),
                id=f"{prefix}-read{counter}",
                name=f"{prefix}-read{counter}",
                description=f"Synthetic read 1 {prefix}-{counter}",
                letter_annotations={
                    "phred_quality": Sequencing.get_phred_scores(len(read1))
                },
            )
            reads1.append(read1)

            # read 2
            _, _, _, transcript = transcripts[counter]
            transcript = Sequencing.mutate(transcript, self._randint(6))
            read2 = SeqRecord(
                seq=Seq(transcript),
                id=f"{prefix}-{counter}",
                name=f"{prefix}-{counter}",
                description=f"Synthetic read 2 {prefix}-{counter}",
                letter_annotations={
                    "phred_quality": Sequencing.get_phred_scores(len(transcript))
                },
            )
            reads2.append(read2)

        return reads1, reads2

    @classmethod
    def write_fastq(cls, reads1: List, reads2: List, path_prefix: str) -> Tuple:
        """\
        Writes 2 `FASTQ` files for Read 1 and Read 2.
        Returns the `FASTQ` files paths.

        Parameters
        ----------
        reads1
            List of Read 1 as `SeqRecord`.
        reads2
            List of Read 2 as `SeqRecord`.
        path_prefix
            Path prefix for `FASTQ`.`gz`.
        """
        fastq1 = f"{path_prefix}.R1.fastq.gz"
        fastq2 = f"{path_prefix}.R2.fastq.gz"

        with gzip.open(fastq1, "wt") as file_obj:
            SeqIO.write(sequences=reads1, handle=file_obj, format="fastq")

        with gzip.open(fastq2, "wt") as file_obj:
            SeqIO.write(sequences=reads2, handle=file_obj, format="fastq")

        return fastq1, fastq2
