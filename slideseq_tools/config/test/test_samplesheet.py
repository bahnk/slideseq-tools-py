"""
Testing module for the slideseq_tools.config.samplesheet module.
"""

from pathlib import Path
import pytest

import slideseq_tools
from ..samplesheet import SampleSheet


class TestSampleSheet:
    """The test class associated with the SampleSheet class."""

    def test_constructor_not_existing_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when samplesheet path
        doesn't exist.
        """
        path = tmp_path / "file"
        with pytest.raises(FileNotFoundError):
            SampleSheet(path)

    def test_save(self, tmp_path):
        """Tests if `save` method ouputs something."""
        path = tmp_path / "file.csv"
        root = Path(slideseq_tools.__file__).parent
        samplesheet_path = root / "config/test/data/samplesheet.csv"
        samplesheet = SampleSheet(path=samplesheet_path)
        samplesheet.save(path)
        assert path.exists()
