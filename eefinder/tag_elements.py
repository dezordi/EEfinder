"""Flag overlapping elements and compute per-element average identity."""

from __future__ import annotations
import numpy as np
import pandas as pd
from eefinder.get_taxonomy import PROTEIN_ID_SEPARATOR
from eefinder.lineage import merge_key

OVERLAP_MARGIN = 100


def _list_to_string(overlaped_elements: list) -> str:
    """Join a list of element IDs into a comma-separated string."""
    return ",".join(map(str, overlaped_elements))


def _average_pident(protein_ids: str) -> float:
    """Return the mean percent identity encoded in a ``protein_ids`` value.

    The field looks like ``"ACC1|30.0;ACC2|45.0"``; each entry carries the
    hit's identity after its ``|``.

    Returns
    -------
    float
        The rounded mean identity, or ``nan`` when no identity is present.
    """
    entries = str(protein_ids).split(PROTEIN_ID_SEPARATOR)
    values = [float(entry.split("|")[1]) for entry in entries if "|" in entry]
    return round(np.mean(values), 1) if values else np.nan


class TagElements:
    """Annotate a taxonomy table with overlap tags and average identity.

    Adds two columns in a single pass over the table:

    * ``overlaped_element_id`` / ``tag`` — elements on the same contig within
      :data:`OVERLAP_MARGIN` of each other but assigned to a *different* family
      are cross-referenced and tagged ``"overlaped"``; the rest are ``"unique"``.
    * ``average_pident`` — the mean percent identity of the element's hits.

    Overwrites ``tax_file`` in place. Runs on instantiation.

    Parameters
    ----------
    tax_file : str
        Path to the taxonomy TSV to annotate.
    """

    def __init__(self, tax_file: str) -> None:
        self.tax_file = tax_file

        self.tag_elements()

    def _coordinates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract contig/start/end/taxon from the table.

        The taxon is the family when the lineage has one, otherwise the lineage
        itself, so overlaps between unclassified viruses are only called within
        one lineage.
        """
        element_id = df["element_id"]
        return pd.DataFrame(
            {
                "id": element_id,
                "contig": element_id.str.replace(":.*", "", regex=True),
                "start": element_id.str.replace("-.*", "", regex=True)
                .str.replace(".*:", "", regex=True)
                .astype(int),
                "end": element_id.str.replace(".*-", "", regex=True).astype(int),
                "family": df["taxonomy"].map(lambda t: merge_key(t, "family")),
            }
        )

    def tag_elements(self) -> None:
        """Compute overlap tags and average identity, then rewrite the file."""
        df = pd.read_csv(self.tax_file, sep="\t")
        df["element_id"] = df["element_id"].str.replace(".*/", "", regex=True)

        coords = self._coordinates(df)
        overlaps = []
        for i, row in coords.iterrows():
            same_contig = coords[coords["contig"] == row["contig"]]
            matched = [
                other["id"]
                for j, other in same_contig.iterrows()
                if i != j
                and other["start"] <= row["end"] + OVERLAP_MARGIN
                and other["end"] >= row["start"] - OVERLAP_MARGIN
                and other["family"] != row["family"]
            ]
            overlaps.append(_list_to_string(matched))

        df["overlaped_element_id"] = overlaps
        df["tag"] = np.where(df["overlaped_element_id"] == "", "unique", "overlaped")
        df["average_pident"] = df["protein_ids"].apply(_average_pident)

        df.to_csv(self.tax_file, sep="\t", index=False, header=True)
