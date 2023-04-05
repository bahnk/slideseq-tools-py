"""
Creates synthetic Slide-seq data.
"""

# coding: utf-8

import os
import sys
import logging
from pathlib import Path
import click
import pandas as pd

from slideseq_tools.synthetic_data.slideseq import SlideSeq

# pylint: disable=no-value-for-parameter
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
@click.command()
@click.option("--n-samples", default=2, help="number of samples")
@click.option("--n-files", default=5, help="number of files per sample")
@click.option("--n-reads", default=int(2 * 1e4), help="number of reads per file")
@click.option("--read-structure", default="8C18U6C2X9M", help="read 1 structure")
@click.option("--out-dir", default="data", help="number of reads per file")
@click.argument("tiff_path")
@click.argument("genome_path")
def main(n_samples, n_files, n_reads, read_structure, out_dir, tiff_path, genome_path):
    """
    Create synthetic Slide-seq data.
    """
    # tiff file
    if not os.path.exists(tiff_path):
        raise FileNotFoundError(f"{tiff_path} doesn't exist.")

    # genome
    genome_path = Path(genome_path)
    if not genome_path.exists():
        raise FileNotFoundError(f"{genome_path} doesn't exist.")
    gff_path = genome_path / "Annotation/Genes/genes.gtf"
    if not gff_path.exists():
        raise FileNotFoundError(f"{gff_path} doesn't exist.")
    fasta_path = genome_path / "Sequence/WholeGenomeFasta/genome.fa"
    if not fasta_path.exists():
        raise FileNotFoundError(f"{fasta_path} doesn't exist.")

    # output directory
    out_dir = Path(out_dir)
    if out_dir.exists() and not out_dir.is_dir():
        raise FileExistsError(f"{out_dir} already exists and not a directory.")
    if not out_dir.exists():
        try:
            os.makedirs(out_dir)
        except OSError as _:
            pass

    slideseq = SlideSeq(tiff_path=tiff_path, gff_path=gff_path, fasta_path=fasta_path)

    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

    rows = []

    for sample_num in range(1, n_samples + 1):

        sample = f"sample{sample_num}"

        # puck
        slideseq.generate_puck()
        puck_path = str(out_dir / f"{sample}.csv")
        slideseq.save_coordinates(puck_path)

        # reads
        for file_num in range(1, n_files + 1):

            prefix = f"{sample}-file{file_num}"
            path_prefix = str(out_dir / f"{sample}_L{file_num:03d}")

            logging.info("Creating %s", path_prefix)

            reads1, reads2 = slideseq.generate_reads(prefix=prefix, n_reads=n_reads)

            fastq1, fastq2 = SlideSeq.write_fastq(
                reads1=reads1, reads2=reads2, path_prefix=path_prefix
            )

            row = {
                "sample": sample,
                "fastq_1": fastq1,
                "fastq_2": fastq2,
                "puck": puck_path,
                "read_structure": read_structure,
                "genome": genome_path,
            }

            rows.append(row)

    samplesheet = pd.DataFrame.from_records(rows)
    samplesheet.to_csv(out_dir / "samplesheet.csv", index=False)


if __name__ == "__main__":
    main()
