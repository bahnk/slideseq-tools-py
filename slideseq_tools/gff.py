"""
Manages GFF files.
"""

from pathlib import Path
from typing import Dict

import re


class GFF:
    """GFF file."""

    path: str

    def __init__(self, gff_path: str) -> None:
        """\
        Constructor of GFF class.

        Raises a `FileNoFoundError` if  `GFF` file doesn't exist.

        Parameters
        ----------
        path 
            Path of the `GFF` file.
        """
        path = Path(gff_path)

        if not path.exists():
            raise FileNotFoundError(f"GFF file {gff_path} doesn't exist.")

        self.path = path

    def parse_attributes(self, attribute_str: str) -> Dict:
        """Returns a dictionary containg attributes."""
        attributes = {}

        regex = re.compile(r'^\s*(?P<tag>\S+)\s"(?P<value>\S+)"\s*$')

        for pair in attribute_str.split(";"):

            if not pair:
                continue

            mtch = regex.match(pair)

            if mtch:
                attributes[mtch.group("tag")] = mtch.group("value")
            else:
                raise ValueError(f"Cannot parse field: {pair}")

        return attributes

    def parse_record(self, record_str: str) -> Dict:
        """Returns dictionary containing columns of a record. The method raises a
        `ValueError` if the record cannot be parsed."""
        colnames = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]

        cols = record_str.split("\t")
        cols[-1] = self.parse_attributes(cols[-1])

        if len(colnames) != len(cols):
            raise ValueError(
                f"Record {record_str} doesn't contain {len(colnames)} columns."
            )

        record = dict(list(zip(colnames, cols)))

        return record

    def get_features(self, min_length: int = 50) -> Dict:
        """Returns a list of features (`seqname`, `start`, `end`)."""
        features = []

        with open(self.path, "r", encoding="utf-8") as file_obj:
            for line in file_obj:
                rec = self.parse_record(line.rstrip())
                if (int(rec["end"]) - int(rec["end"]) + 1) < min_length:
                    feature = (rec["seqname"], int(rec["start"]), int(rec["end"]))
                    features.append(feature)

        return features
