"""
For ample sheet validity.
"""

import re
from pathlib import Path
from typing import Dict
import pandas as pd

# pylint: disable=no-name-in-module
from pydantic import BaseModel, FilePath, DirectoryPath

from slideseq_tools.config.read_structure import ReadStructure


class SampleSheetRow(BaseModel):
    """Data class for for Slide-seq sample sheet row."""

    sample: str
    fastq_1: FilePath
    fastq_2: FilePath
    puck: FilePath
    read_structure: str
    genome: DirectoryPath

    def min_length(self) -> int:
        """ "Returns minimum read length considering the structure."""
        structure = ReadStructure(self.read_structure)
        return structure.min_length()

    def umi_tools_regex(self) -> str:
        """ "Returns UMI-tools regex considering the structure."""
        structure = ReadStructure(self.read_structure)
        return structure.umi_tools_regex()

    def gff_path(self) -> str:
        """ "Returns GTF file path."""
        path = Path(self.genome) / "Annotation/Genes/genes.gtf"
        return str(path)

    def star_index(self) -> str:
        """ "Returns STAR index directory path."""
        path = Path(self.genome) / "Sequence/STARIndex"
        return str(path)

    def puck_name(self) -> str:
        """Returns puck name from path."""
        return re.sub(r"\.csv$", "", Path(self.puck).name)

    def dict(self) -> Dict:
        """Returns `dict` containing original and additional info."""
        return {
            **super().dict(),
            "min_length": self.min_length(),
            "umi_tools_regex": self.umi_tools_regex(),
            "gff": self.gff_path(),
            "star_index": self.star_index(),
            "puck_name": self.puck_name(),
        }


class SampleSheet:
    """Slide-seq experiment sample sheet."""

    path: Path
    original_dframe: None
    dframe: None

    def __init__(self, path: str) -> None:
        """\
        Constructor taking sample sheet path.

        Constructor raises:

            * `FileNoFoundError` if sample sheet file doesn't exist
            * `IOError` if sample sheet can't be read

        Parameters
        ----------
        path
            Path of the sample sheet `CSV` file.
        """
        path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"Sample sheet {path} doesn't exist.")

        try:
            dframe = pd.read_csv(path)
        except Exception as exc:
            raise IOError(f"Sample sheet {path} can't be read.") from exc

        self.original_dframe = dframe

    def create_samplesheet(self) -> None:
        """
        Creates sample sheet with additional required columns for downstream
        processing.
        """
        sheet_rows = []

        for _, row in self.original_dframe.iterrows():
            sheet_row = SampleSheetRow(**row.to_dict())
            sheet_rows.append(sheet_row.dict())

        self.dframe = pd.DataFrame.from_records(sheet_rows)

    def save(self, path: str) -> None:
        """
        Saves modified sample sheet as `CSV`.

        Parameters
        ----------
        path
            `CSV` path to use.
        """
        if not hasattr(self, "dframe"):
            self.create_samplesheet()
        self.dframe.to_csv(path, index=False)
