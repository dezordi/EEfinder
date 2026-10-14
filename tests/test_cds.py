"""Unit tests for eefinder.cds: the elements' CDSs and their proteins."""

from __future__ import annotations

import pandas as pd

from eefinder.cds import _fragment_names, _orf, element_cds, write_faa
from eefinder.gff import WriteGFF3

SIGNATURE_COLUMNS = (
    "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,"
    "evalue,bitscore,sense,bed_name,tag,Taxonomy,Molecule_type,Protein,Host"
).split(",")


def _hit(sseqid, qstart, qend, sense, bitscore):
    return [
        f"P/c1:{qstart}-{qend}", sseqid, 50.0, 3, 0, 0, qstart, qend, 1, 3,
        1e-5, bitscore, sense, "", "EE", "f__Fam", "ssRNA", f"{sseqid} product", "h",
    ]  # fmt: skip


def _inputs(tmp_path):
    # c1:10-40 (+), frame "MA*GGFPKPP": hit PB (14..25, "A*GG") lies in the
    # ORF ATG GCC TAA (11..19) and past its stop
    # c1:100-130 (-), reverse complement GGG ATG TTT GGG AAA TAA ...: hit PC
    # (116..124, "FGK") lies in the ORF ATG..TAA (113..127, "MFGK")
    fasta = tmp_path / "elements.fa"
    fasta.write_text(
        ">P/c1:10-40\nATGGCCTAAGGGGGGTTTCCCAAACCCCCC\n"
        ">P/c1:100-130\nGGGGGGGGGGGGTTATTTCCCAAACATCCC\n"
    )
    tax = tmp_path / "elements.tax.tsv"
    pd.DataFrame(
        {
            "element_id": ["c1:10-40", "c1:100-130"],
            "sense": ["pos", "neg"],
            "protein_ids": ["PA|50.0;PB|50.0", "PC|50.0"],
        }
    ).to_csv(tax, sep="\t", index=False)
    signature = tmp_path / "taxonomy_signature.csv"
    pd.DataFrame(
        [
            _hit("PA", 11, 19, "pos", 50.0),
            _hit("PB", 14, 25, "pos", 40.0),  # overlaps PA: longer, scores lower
            _hit("PX", 11, 19, "pos", 90.0),  # not one of the element's proteins
            _hit("PC", 116, 124, "neg", 30.0),
        ],
        columns=SIGNATURE_COLUMNS,
    ).assign(Protein=["pa", "polyprotein, partial", "px", "capsid"]).to_csv(
        signature, index=False
    )
    return fasta, tax, signature


def test_element_cds_is_the_orf_of_the_longest_hit(tmp_path):
    fasta, tax, signature = _inputs(tmp_path)
    cds = element_cds(str(fasta), str(tax), str(signature), prefix="P")

    (plus,) = cds["P/c1:10-40"]
    assert (plus.id, plus.start, plus.end, plus.strand) == (
        "P/c1:10-40.cds1",
        11,
        19,
        "+",
    )
    assert (plus.protein, plus.sequence) == ("PB", "MA")
    assert (plus.aligned_start, plus.aligned_end, plus.stop_codons) == (14, 25, 1)
    (minus,) = cds["P/c1:100-130"]
    assert (minus.start, minus.end, minus.strand, minus.sequence) == (
        113,
        127,
        "-",
        "MFGK",
    )
    assert not plus.partial and not minus.partial


def test_orf_open_to_the_element_end_is_partial_and_needs_an_atg():
    assert _orf("ATGAAAAAA", 3, 9, "pos") == (0, 9, "MKK", 0, True)
    assert _orf("AAAAAAAAA", 0, 9, "pos") is None


def test_gff3_cds_children_and_faa(tmp_path):
    fasta, tax, signature = _inputs(tmp_path)
    cds = element_cds(str(fasta), str(tax), str(signature), prefix="P")
    tax_full = tmp_path / "full.tax.tsv"
    full = pd.read_csv(tax, sep="\t")
    for column in ("protein_products", "molecule_type", "taxonomy", "host", "tag"):
        full[column] = ""
    full["average_pident"] = 50.0
    full.to_csv(tax_full, sep="\t", index=False)

    gff = tmp_path / "out.gff3"
    WriteGFF3(str(tax_full), str(gff), prefix="P", cds=cds)
    lines = [line.split("\t") for line in gff.read_text().splitlines()[1:]]
    assert [line[2] for line in lines] == [
        "endogenous_viral_element",
        "CDS",
        "endogenous_viral_element",
        "CDS",
    ]
    assert lines[1][3:8] == ["11", "19", "50.0", "+", "0"]
    assert "Parent=P/c1:10-40;" in lines[1][8]
    assert "Name=polyprotein%2C partial;hit=PB;" in lines[1][8]
    assert "aligned_range=14-25;" in lines[1][8]
    assert "product=" not in lines[1][8]
    assert lines[1][8].endswith("pseudo=true;stop_codons=1")
    assert "pseudo" not in lines[3][8]

    faa = tmp_path / "out.faa"
    write_faa(cds, str(faa))
    assert faa.read_text().splitlines() == [
        ">P/c1:10-40.cds1 hit=PB ref_range=1-3 pident=50.0 pseudo=true partial=false"
        " name=polyprotein, partial",
        "MA",
        ">P/c1:100-130.cds1 hit=PC ref_range=1-3 pident=50.0 pseudo=false"
        " partial=false name=capsid",
        "MFGK",
    ]


def test_cds_sharing_a_product_are_numbered_along_the_protein():
    # in genome order; on the minus strand the protein runs the other way
    names = _fragment_names(
        ["capsid", "capsid", "polyprotein", "capsid", ""],
        [392, 267, 10, 145, 1],
    )
    assert names == [
        "capsid fragment 3",
        "capsid fragment 2",
        "polyprotein",
        "capsid fragment 1",
        "",
    ]
