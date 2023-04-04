"""
Testing module for the slideseq_tools.gff module.
"""

import pytest

from ..gff import GFF

# pylint: disable=too-few-public-methods
class TestGFF:
    """The test class associated with the GFF class."""

    def test_constructor_not_existing_path(self, tmp_path):
        """
        Tests if constructor returns `FileNotFoundError` when `GFF` path
        doesn't exist.
        """
        path = tmp_path / "file"
        with pytest.raises(FileNotFoundError):
            GFF(path)
