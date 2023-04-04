"""
Manages Slide-seq read 1 structure containing bead barcode and UMI.
"""

import re
from typing import List, Tuple

# pylint: disable=too-few-public-methods
class ReadStructure:
    """Structure of read 1 containing bead barcode and UMI."""

    definition: str
    sequence: str
    segments: List

    def __init__(self, structure: str) -> None:
        """\
        Constructor taking structure definition.

        Raises a `ValueError` if the structure is not valid.

        Parameters
        ----------
        structure
            A `str` specifying the read structure, for example `8C18U7C2X8M`.
            `C` is the bead barcode, `M` is the UMI, `U` is the UP primer and
            `X` is ignored.
        """
        structure = structure.upper()
        self._check(structure)
        self.sequence, self.segments = self._parse(structure)
        self.structure = structure

    def _check(self, structure: str) -> str:
        """\
        Checks structure `str` validity.

        Parameters
        ----------
        structure
            Structure `str` to validate.
        """
        letters = set(list(re.sub(r"\d", "", structure)))

        # missing barcode
        if "C" not in letters:
            raise ValueError(
                f"Read structure {structure} doesn't have "
                "at least a C to locate the bead barcode."
            )

        # missing UP primer
        if "U" not in letters:
            raise ValueError(
                f"Read structure {structure} doesn't have "
                "a U to locate the UP primer."
            )

        # missing UMI
        if "M" not in letters:
            raise ValueError(
                f"Read structure {structure} doesn't have "
                "at least a M to locate the UMI."
            )

        # invalid letters
        for letter in ["C", "M", "U", "X"]:
            try:
                letters.remove(letter)
            except KeyError as _:
                pass

        if len(letters) > 0:
            raise ValueError(
                f"Read structure {structure} shouldn't contain " f"{','.join(letters)}."
            )

    def _parse(self, structure: str) -> Tuple[List]:
        """\
        Returns sequence and segments from read structure definition.

        The method raises a `ValueError` if the structure definition is not
        well formed.

        Parameters
        ----------
        structure
            A `str` specifying the read structure, for example `8C18U7C2X8M`.
            `C` is the bead barcode, `M` is the UMI, `U` is the UP primer and
            `X` is ignored.
        """
        symbols = re.split(r"\d+", structure)
        lengths = re.split("C|M|U|X", structure)

        # letter-number pairs
        if len(symbols) != len(lengths):
            raise ValueError(
                f"Read structure {structure} "
                "doesn't have as many numbers as letters."
            )

        sequence = ""
        segments = []
        counter = []

        for symbol, length in zip(symbols[1:], lengths[:-1]):
            counter.append(symbol)
            segments.append((symbol, length, counter.count(symbol)))
            sequence += "".join([symbol] * int(length))

        return sequence, segments

    def min_length(self) -> int:
        """\
        Returns minimum read length considering the structure.
        """
        return len(re.sub("X+$", "", self.sequence))

    def umi_tools_regex(self) -> str:
        """\
        Returns a UMI tools `bc-pattern` regex.
        """
        names = {"C": "cell", "M": "umi", "X": "discard"}
        regex = ""
        counter = 0
        for symbol, length, num in self.segments:
            if symbol == "U":
                regex += f".{{1.{length}}}"
            else:
                regex += f"(?P<{names[symbol]}_{num}>.{{1,{length}}})"
                if symbol == "X":
                    counter += 1

        return "^" + regex + f"(?P<discard_{counter+1}>.*)$"
