"""
Testing module for the slideseq_tools.config.read_structure module.
"""

import pytest

from ..read_structure import ReadStructure


class TestReadStructure:
    """The test class associated with the ReadStructure class."""

    def test_returns_valuerror_with_improper_structure(self):
        """
        Tests if constructor returns `ValueError` when given an improper
        structure doesn't exist.
        """
        structures = [
            "8C18U",
            "8C9M",
            "18U9M",  # missing C, U, M
            "8C18U5B",
            "8CB18",  # invalid letters
        ]
        with pytest.raises(ValueError):
            for structure in structures:
                ReadStructure(structure)

    def test_returns_valid_min_length(self):
        """Tests if `min_length` method returns a proper value."""
        definitions = [
            ("8C18U6C2X9M", 43),
            ("8C18U6C2X7M", 41),
            ("8C18U6C9M", 41),
            ("8C18U6C9M2X", 41),
        ]
        for struct_def, length in definitions:
            structure = ReadStructure(struct_def)
            assert length == structure.min_length()
