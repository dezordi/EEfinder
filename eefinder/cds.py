"""Coding regions of the elements, from the viral protein hits that built them.

Each validated hit of an element yields the ORF (ATG to stop codon, standard
genetic code) of its reading frame that overlaps its aligned region most. Where
hits of one element overlap, only the longest one is used (the best-scoring one
on a tie). A hit with a stop codon in its aligned region is flagged ``pseudo``:
its ORF covers only part of the homology.
"""

from __future__ import annotations
from dataclasses import dataclass
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from eefinder.get_taxonomy import PROTEIN_ID_SEPARATOR


@dataclass
class Cds:
    """One coding region, in GFF3 coordinates (1-based, inclusive)."""

    id: str
    start: int
    end: int
    strand: str
    protein: str
    name: str
    ref_start: int
    ref_end: int
    pident: float
    evalue: float
    bitscore: float
    sequence: str
    aligned_start: int
    aligned_end: int
    stop_codons: int
    partial: bool


def _orf(seq: str, start: int, end: int, sense: str):
    """Return the ORF of a hit's reading frame that overlaps the hit most.

    ``seq`` is the element, ``[start, end)`` the hit in it (0-based, whole
    codons). The frame is cut at its stop codons; in each piece the ORF runs
    from the first ATG to the stop, or to the element's end when the frame is
    still open there (``partial``). Returns 0-based ``(start, end, protein,
    stops_in_hit, partial)`` in ``seq``, the end including the stop codon, or
    ``None`` when no ORF overlaps the hit.
    """
    length = len(seq)
    if sense == "neg":
        seq = str(Seq(seq).reverse_complement())
        start, end = length - end, length - start
    offset = start % 3
    protein = str(Seq(seq[offset : offset + (length - offset) // 3 * 3]).translate())
    first, last = (start - offset) // 3, (end - offset) // 3

    best, piece_start = None, 0
    for piece_end in [i for i, aa in enumerate(protein) if aa == "*"] + [len(protein)]:
        atg = protein.find("M", piece_start, piece_end)
        overlap = min(piece_end, last) - max(atg, first)
        if atg != -1 and overlap > 0 and (best is None or overlap > best[0]):
            best = (overlap, atg, piece_end)
        piece_start = piece_end + 1
    if best is None:
        return None

    _, atg, stop = best
    partial = stop == len(protein)
    orf_start, orf_end = offset + 3 * atg, offset + 3 * stop + (0 if partial else 3)
    if sense == "neg":
        orf_start, orf_end = length - orf_end, length - orf_start
    return (
        orf_start,
        orf_end,
        protein[atg:stop],
        protein[first:last].count("*"),
        partial,
    )


def _fragment_names(products: "list[str]", ref_starts: "list[int]") -> "list[str]":
    """Number the products an element's CDSs share, along the reference protein.

    ``polyprotein fragment 1`` is the one aligned nearest the protein's
    N-terminus, whatever the strand.
    """
    names = list(products)
    for product in set(products):
        shared = [i for i, other in enumerate(products) if other == product]
        if product and len(shared) > 1:
            by_protein = sorted(shared, key=lambda i: ref_starts[i])
            for n, i in enumerate(by_protein, start=1):
                names[i] = f"{product} fragment {n}"
    return names


def _accessions(protein_ids: str) -> "set[str]":
    return {
        entry.split("|")[0].strip()
        for entry in str(protein_ids).split(PROTEIN_ID_SEPARATOR)
    }


def element_cds(
    elements_fasta: str, tax_file: str, taxonomy_signature: str, prefix: str = ""
) -> "dict[str, list[Cds]]":
    """Return the CDSs of each element, keyed by its FASTA header.

    A hit belongs to an element when it lies on the element's contig and
    strand, within it, and its protein is one of the element's
    ``protein_ids`` (elements of different taxa may overlap).

    Parameters
    ----------
    elements_fasta : str
        Element sequences, headed ``{prefix}/{contig}:{start}-{end}`` (BED-style).
    tax_file : str
        Their taxonomy table (``element_id``, ``sense``, ``protein_ids``).
    taxonomy_signature : str
        The per-hit table of the taxonomy stage (genome coordinates, protein,
        product, scores).
    prefix : str
        Run prefix, as in the FASTA headers.
    """
    sequences = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(elements_fasta, "fasta")
    }
    tax = pd.read_csv(tax_file, sep="\t")
    hits = pd.read_csv(taxonomy_signature, sep=",")
    hits["contig"] = hits["qseqid"].str.replace(r"\:.*", "", regex=True)

    cds = {}
    for _, element in tax.iterrows():
        header = (
            f"{prefix}/{element['element_id']}" if prefix else element["element_id"]
        )
        if header not in sequences:
            continue
        contig, _, coords = header.rpartition(":")
        lo, hi = (int(x) for x in coords.split("-"))
        candidates = hits[
            (hits["contig"] == contig)
            & (hits["sense"] == element["sense"])
            & (hits["qstart"] - 1 >= lo)
            & (hits["qend"] <= hi)
            & hits["sseqid"].isin(_accessions(element["protein_ids"]))
        ]
        candidates = candidates.assign(
            span=candidates["qend"] - candidates["qstart"]
        ).sort_values(["span", "bitscore"], ascending=False, kind="stable")

        kept, orfs = [], set()
        for _, hit in candidates.iterrows():
            region = (int(hit["qstart"]) - 1, int(hit["qend"]))
            if any(region[0] < end and start < region[1] for start, end in kept):
                continue
            kept.append(region)
            orf = _orf(sequences[header], region[0] - lo, region[1] - lo, hit["sense"])
            # hits in one open frame share their ORF: the longest hit names it
            if orf is None or orf[:2] in orfs:
                continue
            orfs.add(orf[:2])
            genome_orf = (lo + orf[0], lo + orf[1], *orf[2:])
            cds.setdefault(header, []).append((region, hit, genome_orf))

    result = {}
    for header, found in cds.items():
        found.sort(key=lambda item: item[2][:2])
        names = _fragment_names(
            [_product(hit) for _, hit, _ in found],
            [int(hit["sstart"]) for _, hit, _ in found],
        )
        result[header] = [
            Cds(
                id=f"{header}.cds{n}",
                start=orf[0] + 1,
                end=orf[1],
                strand="+" if hit["sense"] == "pos" else "-",
                protein=hit["sseqid"],
                name=name,
                ref_start=int(hit["sstart"]),
                ref_end=int(hit["send"]),
                pident=float(hit["pident"]),
                evalue=float(hit["evalue"]),
                bitscore=float(hit["bitscore"]),
                sequence=orf[2],
                aligned_start=region[0] + 1,
                aligned_end=region[1],
                stop_codons=orf[3],
                partial=orf[4],
            )
            for n, ((region, hit, orf), name) in enumerate(zip(found, names), start=1)
        ]
    return result


def _product(hit) -> str:
    """The hit's protein name from the database metadata, whitespace collapsed."""
    return "" if pd.isna(hit["Protein"]) else " ".join(str(hit["Protein"]).split())


def _element_order(header: str):
    contig, _, coords = header.rpartition(":")
    start, _, end = coords.partition("-")
    return contig, int(start), int(end)


def write_faa(cds: "dict[str, list[Cds]]", output_file: str) -> None:
    """Write one protein per CDS, in the GFF3's order (contig, start)."""
    with open(output_file, "w") as faa:
        for header in sorted(cds, key=_element_order):
            for region in cds[header]:
                faa.write(
                    f">{region.id} hit={region.protein} "
                    f"ref_range={region.ref_start}-{region.ref_end} "
                    f"pident={region.pident} "
                    f"pseudo={'true' if region.stop_codons else 'false'} "
                    f"partial={'true' if region.partial else 'false'} "
                    f"name={region.name}\n{region.sequence}\n"
                )
