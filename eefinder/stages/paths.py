"""Canonical file names for a ``screening`` run.

Every file a stage reads or writes is named here, once. The step classes still
write their own historical suffixes; the stages rename those outputs to the
names below, so the stage interface does not depend on them.

The ``PREFIX.EEs.*`` names are the published outputs and must not change.
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class ScreeningPaths:
    """Every path a ``screening`` run uses, derived from the output directory.

    Parameters
    ----------
    outdir : str
        Run output directory (no trailing slash).
    prefix : str
        Run prefix, also the stem of every file below.

    Example
    -------
    >>> paths = ScreeningPaths(outdir="results", prefix="Aedes")
    >>> paths.cleaned_genome
    'results/Aedes.cleaned_genome.fa'
    >>> paths.ee_elements_tax
    'results/Aedes.elements.tax.tsv'
    """

    outdir: str
    prefix: str

    # -- helpers ----------------------------------------------------------
    def _p(self, suffix: str) -> str:
        return f"{self.outdir}/{self.prefix}.{suffix}"

    @property
    def cleaned_genome(self) -> str:
        """Return the prefixed, length-filtered genome."""
        return self._p("cleaned_genome.fa")

    @property
    def genome_lengths(self) -> str:
        """Return the ``<contig>\\t<length>`` index ``bedtools slop`` needs."""
        return self._p("cleaned_genome.lengths.tsv")

    @property
    def ee_hits(self) -> str:
        """Return the raw hits of the genome against the EE database."""
        return self._p("ee_hits.tsv")

    @property
    def ee_hits_filtered(self) -> str:
        """Return the EE hits with redundant alignments collapsed."""
        return self._p("ee_hits.filtered.tsv")

    @property
    def ee_candidates_bed(self) -> str:
        """Return the candidate element intervals."""
        return self._p("ee_candidates.bed")

    @property
    def ee_candidates(self) -> str:
        """Return the candidate element sequences."""
        return self._p("ee_candidates.fa")

    @property
    def host_hits(self) -> str:
        """Return the raw hits of the candidates against the baits."""
        return self._p("host_hits.tsv")

    @property
    def host_hits_filtered(self) -> str:
        """Return the host-bait hits with redundant alignments collapsed."""
        return self._p("host_hits.filtered.tsv")

    @property
    def ee_hits_validated(self) -> str:
        """Return the hits that outscored their best host-bait hit."""
        return self._p("ee_hits.validated.tsv")

    @property
    def taxonomy_signature(self) -> str:
        """Return the per-hit taxonomy signature."""
        return self._p("taxonomy_signature.csv")

    @property
    def merge_bed(self) -> str:
        """Return the per-hit intervals annotated with the merge key."""
        return self._p("merge.annotated.bed")

    @property
    def merge_bed_merged(self) -> str:
        """Return the result of ``bedtools merge`` over :attr:`merge_bed`."""
        return self._p("merge.merged.bed")

    @property
    def merge_elements_bed(self) -> str:
        """Return the merged intervals with the annotation stripped off."""
        return self._p("merge.elements.bed")

    @property
    def ee_elements(self) -> str:
        """Return the merged element sequences."""
        return self._p("elements.fa")

    @property
    def ee_elements_tax(self) -> str:
        """Return the taxonomy table, one row per merged element."""
        return self._p("elements.tax.tsv")

    # -- post-processing ---------------------------------------------------
    @property
    def ee_elements_cleaned(self) -> str:
        """Return the elements surviving the soft-mask filter."""
        return self._p("elements.cleaned.fa")

    @property
    def ee_elements_cleaned_tax(self) -> str:
        """Return the taxonomy table of the mask-cleaned elements."""
        return self._p("elements.cleaned.tax.tsv")

    @property
    def ee_elements_gff3(self) -> str:
        """Return the GFF3 annotation of the merged elements."""
        return self._p("elements.gff3")

    @property
    def ee_elements_cleaned_gff3(self) -> str:
        """Return the GFF3 annotation of the mask-cleaned elements."""
        return self._p("elements.cleaned.gff3")

    @property
    def flanks_bed(self) -> str:
        """Return the element intervals."""
        return self._p("flanks.bed")

    @property
    def flanks_slop_bed(self) -> str:
        """Return the element intervals widened by ``--flank`` nt."""
        return self._p("flanks.slop.bed")

    @property
    def flanks(self) -> str:
        """Return the element sequences plus their flanking regions."""
        return self._p("flanks.fa")

    # -- published outputs (do not rename) --------------------------------
    @property
    def final_fasta(self) -> str:
        return self._p("EEs.fa")

    @property
    def final_tax(self) -> str:
        return self._p("EEs.tax.tsv")

    @property
    def final_gff3(self) -> str:
        return self._p("EEs.gff3")

    @property
    def final_flanks(self) -> str:
        return self._p("EEs.flanks.fa")

    @property
    def final_cleaned_fasta(self) -> str:
        return self._p("EEs.cleaned.fa")

    @property
    def final_cleaned_tax(self) -> str:
        return self._p("EEs.cleaned.tax.tsv")

    @property
    def final_cleaned_gff3(self) -> str:
        return self._p("EEs.cleaned.gff3")

    # -- housekeeping ------------------------------------------------------
    @property
    def tmp_dir(self) -> str:
        """Return the directory intermediates are archived in."""
        return f"{self.outdir}/tmp_files"

    @property
    def removed_dir(self) -> str:
        """Return the directory overlap filtering preserves drops in."""
        return f"{self.outdir}/tmp_outputs"

    def final_outputs(self) -> "list[str]":
        """Return the published outputs of a run, plus the run log.

        Everything else directly in :attr:`outdir` is an intermediate.
        """
        return [
            self.final_fasta,
            self.final_tax,
            self.final_gff3,
            self.final_flanks,
            self.final_cleaned_fasta,
            self.final_cleaned_tax,
            self.final_cleaned_gff3,
            f"{self.outdir}/eefinder.log",
        ]

    def intermediates(self) -> "list[str]":
        """Return the named intermediates, in pipeline order."""
        return [
            self.cleaned_genome,
            self.genome_lengths,
            self.ee_hits,
            self.ee_hits_filtered,
            self.ee_candidates_bed,
            self.ee_candidates,
            self.host_hits,
            self.host_hits_filtered,
            self.ee_hits_validated,
            self.taxonomy_signature,
            self.merge_bed,
            self.merge_bed_merged,
            self.merge_elements_bed,
            self.flanks_bed,
            self.flanks_slop_bed,
        ]
