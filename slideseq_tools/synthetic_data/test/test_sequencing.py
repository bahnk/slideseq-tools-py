"""
Testing module for the slideseq_tools.synthethic_data.sequencing module.
"""

import os
import pytest

from ..sequencing import Sequencing


class TestSequencing:
    """The test class associated with the Sequencing class."""

    gff_subpath = "Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.gtf"
    fasta_subpath = (
        "Bacillus_subtilis_168/Ensembl/EB2/Sequence/WholeGenomeFasta/genome.fa"
    )

    def test_constructor_not_existing_gff_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when `GFF` path
        doesn't exist.
        """
        gff_path = tmp_path / "file"
        fasta_path = os.path.join(os.getenv("AWS_IGENOMES"), self.fasta_subpath)
        with pytest.raises(FileNotFoundError):
            Sequencing(gff_path=gff_path, fasta_path=fasta_path)

    def test_constructor_not_existing_fasta_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when `FASTA` path
        doesn't exist.
        """
        gff_path = os.path.join(os.getenv("AWS_IGENOMES"), self.gff_subpath)
        fasta_path = tmp_path / "file"
        with pytest.raises(FileNotFoundError):
            Sequencing(gff_path=gff_path, fasta_path=fasta_path)

    def test_get_transcripts_returns_proper_length(self):
        """Tests if `get_transcripts` returns proper lenght transcripts."""
        gff_path = os.path.join(os.getenv("AWS_IGENOMES"), self.gff_subpath)
        fasta_path = os.path.join(os.getenv("AWS_IGENOMES"), self.fasta_subpath)
        length = 50
        sequencing = Sequencing(gff_path=gff_path, fasta_path=fasta_path, length=length)
        for _, _, _, transcript in sequencing.get_transcripts():
            assert len(transcript) == length
