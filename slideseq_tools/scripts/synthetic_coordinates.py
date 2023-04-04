"""
Creates puck bead coordinates from a TIFF image.
"""

# coding: utf-8

import click

from slideseq_tools.synthetic_data.spatial import Puck

# pylint: disable=no-value-for-parameter
@click.command()
@click.option("--n-beads", default=int(8 * 1e4), help="number of beads")
@click.argument("tiff_path")
@click.argument("csv_path")
def main(n_beads, tiff_path, csv_path):
    """
    Opens TIFF image and creates coordinates, then subsamples beads and
    saves coordinates in a `CSV` file.
    """
    dframe = Puck.coordinates(tiff_path)
    dframe = dframe.sample(min(dframe.shape[0], n_beads))
    dframe.to_csv(csv_path, header=False, index=False, float_format="%.15f")


if __name__ == "__main__":
    main()
