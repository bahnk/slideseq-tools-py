"""
Checks Slide-seq sample sheet and new columns for downstream processing.
"""

# coding: utf-8

import click

from slideseq_tools.config.samplesheet import SampleSheet

# pylint: disable=no-value-for-parameter
@click.command()
@click.option("--launch-dir", default="./", help="nextflow launch directory")
@click.argument("in_samplesheet")
@click.argument("out_samplesheet")
def main(launch_dir, in_samplesheet, out_samplesheet):
    """
    Opens the sample sheet as `CSV`, checks it and add required columns for
    downstream processing. Finally, save the new sample sheet as `CSV`.
    """
    samplesheet = SampleSheet(in_samplesheet, launch_dir)
    samplesheet.create_samplesheet()
    samplesheet.save(out_samplesheet)


if __name__ == "__main__":
    main()
