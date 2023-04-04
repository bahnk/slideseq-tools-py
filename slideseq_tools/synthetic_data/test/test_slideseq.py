"""
Testing module for the slideseq_tools.synthethic_data.slideseq module.
"""

import os
from pathlib import Path
import pytest

import slideseq_tools
from ..slideseq import SlideSeq


class TestSlideSeq:
    """The test class associated with the SlideSeq class."""

    tiff_path = str(Path(slideseq_tools.__file__).parent / "assets/puck/puck.tif")

    gff_subpath = "Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.gtf"
    fasta_subpath = (
        "Bacillus_subtilis_168/Ensembl/EB2/Sequence/WholeGenomeFasta/genome.fa"
    )

    def test_constructor_not_existing_tiff_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when `TIFF` path
        doesn't exist.
        """
        tiff_path = tmp_path / "file"
        gff_path = os.path.join(os.getenv("AWS_IGENOMES"), self.gff_subpath)
        fasta_path = os.path.join(os.getenv("AWS_IGENOMES"), self.fasta_subpath)
        with pytest.raises(FileNotFoundError):
            SlideSeq(tiff_path=tiff_path, gff_path=gff_path, fasta_path=fasta_path)

    def test_constructor_not_existing_gff_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when `GFF` path
        doesn't exist.
        """
        tiff_path = self.tiff_path
        gff_path = tmp_path / "file"
        fasta_path = os.path.join(os.getenv("AWS_IGENOMES"), self.fasta_subpath)
        with pytest.raises(FileNotFoundError):
            SlideSeq(tiff_path=tiff_path, gff_path=gff_path, fasta_path=fasta_path)

    def test_constructor_not_existing_fasta_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when `FASTA` path
        doesn't exist.
        """
        tiff_path = self.tiff_path
        gff_path = os.path.join(os.getenv("AWS_IGENOMES"), self.gff_subpath)
        fasta_path = tmp_path / "file"
        with pytest.raises(FileNotFoundError):
            SlideSeq(tiff_path=tiff_path, gff_path=gff_path, fasta_path=fasta_path)

    def test_write_fastq(self, tmp_path):
        """Tests if `write_fastq` method ouputs something."""
        tiff_path = self.tiff_path
        gff_path = os.path.join(os.getenv("AWS_IGENOMES"), self.gff_subpath)
        fasta_path = os.path.join(os.getenv("AWS_IGENOMES"), self.fasta_subpath)
        slideseq = SlideSeq(
            tiff_path=tiff_path, gff_path=gff_path, fasta_path=fasta_path
        )
        slideseq.generate_puck()
        reads1, reads2 = slideseq.generate_reads()
        path_prefix = tmp_path / "file"
        SlideSeq.write_fastq(reads1=reads1, reads2=reads2, path_prefix=path_prefix)
        assert os.path.exists(str(path_prefix) + ".R1.fastq.gz")
        assert os.path.exists(str(path_prefix) + ".R2.fastq.gz")
