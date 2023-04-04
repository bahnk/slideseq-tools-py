"""
Creates puck coordinates for testing.
"""

from pathlib import Path
from skimage.io import imread

import numpy as np
import pandas as pd


# pylint: disable=too-few-public-methods
class Puck:
    """Testing puck."""

    path: Path

    @classmethod
    def coordinates(cls, tiff_path: str) -> pd.DataFrame:
        """\
        Creates bead coordinates from a `TIFF` image.

        Raises a `FileNoFoundError` if  `TIFF` image doesn't exist.

        Parameters
        ----------
        tiff_path
            Path of the `TIFF` file.

        Returns
        -------
        pd.DataFrame
            Data frame with two columns `x` and `y` representing coordinates.
        """
        path = Path(tiff_path)

        if not path.exists():
            raise FileNotFoundError(f"TIFF image {path} doesn't exist.")

        img = imread(path)
        zeros = np.where(img == 0)

        rotation = np.array(
            [
                np.cos(3 * np.pi / 2),
                -np.sin(3 * np.pi / 2),
                np.sin(3 * np.pi / 2),
                np.cos(3 * np.pi / 2),
            ]
        ).reshape((2, 2))

        coords = np.matmul(rotation, np.array([zeros[0], zeros[1]]))
        dframe = pd.DataFrame(coords.T)
        dframe.columns = ["x", "y"]

        return dframe
