"""The ``eefinder screening`` command group.

``screening`` runs the pipeline. Invoked without a subcommand it behaves as it
always has and runs ``all``; each stage is also available on its own, taking
explicit inputs and defaulting them to the canonical paths of ``-od``/``-pr``.
"""

from __future__ import annotations
import json
import sys
import time
from dataclasses import asdict
import click
from eefinder import __version__
from eefinder.log import enable_debug, logger
from eefinder.run_message import PaperInfo
from eefinder.stages import (
    DataCleaning,
    FlanksExtraction,
    MergeFragmentedElements,
    PostProcessing,
    PrepareInputs,
    PutativeElementsFilter,
    ScreeningPaths,
    ScreeningPipeline,
    SequenceAlignment,
    TaxonomyAssignment,
)
from eefinder.translation import TRANSLATION_METHODS
from eefinder.utils import RunArguments, RunInfo, check_outdir
from eefinder.versions import (
    collect_dependency_versions,
    collect_system_info,
    find_env_yml,
    report_run_context,
    write_versions_yml,
)


class DefaultGroup(click.Group):
    """Group that runs :attr:`default_command` when none is named.

    Keeps ``eefinder screening <options>`` working now that ``screening`` is a
    group.
    """

    default_command = "all"

    def parse_args(self, ctx, args):
        """Insert the default command when the first argument is not one."""
        if args and args[0] not in self.commands and args[0] not in ("-h", "--help"):
            args = [self.default_command] + list(args)
        return super().parse_args(ctx, args)


# -- option groups -----------------------------------------------------------
def _outdir_options(func):
    """Attach the -od/-pr options."""
    func = click.option(
        "-pr",
        "--prefix",
        help="Write the prefix name for output files. This prefix will be used to create the EEname (The Endogenous Element name will be formated as PREFIX|CONTIG/SCAFFOLD:START-END) default = input file name.",
    )(func)
    func = click.option(
        "-od",
        "--outdir",
        help="Path and dir to store output results.",
        required=True,
    )(func)
    return func


def _database_options(func):
    """Attach the -db/-mt/-bt options."""
    func = click.option(
        "-bt",
        "--hostgenesbaits",
        help="Host genes baits proteins, used to filter putative EEs .fasta file.",
        required=True,
    )(func)
    func = click.option(
        "-mt",
        "--dbmetadata",
        help="Proteins from viruses or bacterias metadata .csv file.",
        required=True,
    )(func)
    func = click.option(
        "-db",
        "--database",
        help="Proteins from viruses or bacterias database .fasta file.",
        required=True,
    )(func)
    return func


def _search_options(func):
    """Attach the -md/-tm/-p/-rj options."""
    func = click.option(
        "-rj",
        "--range_junction",
        help="Sets the range for junction of BLAST/DIAMOND redudant hits, default=100",
        type=int,
        default=100,
    )(func)
    func = click.option(
        "-p",
        "--threads",
        help="Threads for multi-thread analysis, default = 1.",
        type=int,
        default=1,
    )(func)
    func = click.option(
        "-tm",
        "--translation_method",
        help="How proteins are obtained for the similarity searches (applied to "
        "BOTH the main and the host-bait search): 'default' = six-frame "
        "blastx/diamond blastx; 'gv' = pyrodigal-gv prediction; 'rv' = pyrodigal-rv "
        "prediction; 'gv-rv' = both predictions clustered with cd-hit (100%/100%). "
        "Prediction modes align with blastp and map coordinates back to nucleotides. "
        "default = default",
        default="default",
        type=click.Choice(list(TRANSLATION_METHODS)),
    )(func)
    func = click.option(
        "-md",
        "--mode",
        help="Choose between BLAST or the DIAMOND strategies (fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensisitve) to run analysis, default = blastx.",
        default="blastx",
        type=click.Choice(
            [
                "blastx",
                "fast",
                "mid-sensitive",
                "sensitive",
                "more-sensitive",
                "very-sensitive",
                "ultra-sensitive",
            ]
        ),
    )(func)
    return func


def _postprocess_options(func):
    """Attach the -cm/-mp/-ov/-tf/-ntf/-an options."""
    func = click.option(
        "-an",
        "--analysis",
        help="Type of endogenous elements being screened; sets the GFF3 feature "
        "type ('virus' -> endogenous_viral_element, 'bacteria' -> "
        "endogenous_bacterial_element). default = virus",
        default="virus",
        type=click.Choice(["virus", "bacteria"]),
    )(func)
    func = click.option(
        "-ntf",
        "--non_target_families",
        help="Family to DROP when --overlap targets is used. Repeat for multiple "
        "families (e.g. -ntf Retroviridae -ntf Metaviridae). Mutually exclusive "
        "with --target_families.",
        multiple=True,
    )(func)
    func = click.option(
        "-tf",
        "--target_families",
        help="Family to KEEP when --overlap targets is used. Repeat for multiple "
        "families (e.g. -tf Flaviviridae -tf Caulimoviridae). Mutually exclusive "
        "with --non_target_families.",
        multiple=True,
    )(func)
    func = click.option(
        "-ov",
        "--overlap",
        help="What to do with elements tagged as overlaped: 'keep' (default, keep "
        "all), 'longest' (keep the longest element of each overlap), or 'targets' "
        "(keep only elements from --target_families). Filtered-out elements are "
        "saved to tmp_outputs/. default = keep",
        default="keep",
        type=click.Choice(["keep", "longest", "targets"]),
    )(func)
    func = click.option(
        "-cm",
        "--clean_masked",
        help="Remove EEs in regions considered repetitive?",
        is_flag=True,
    )(func)
    func = click.option(
        "-mp",
        "--mask_per",
        help="Limit of lowercase letters in percentage to consider a putative Endogenous Elements as a repetitive region, default = 50.",
        type=int,
        default=50,
    )(func)
    return func


def _common_options(func):
    """Attach the --debug and versions-file options."""
    func = click.option(
        "--debug",
        help="Emit verbose debug logging (intermediate file paths, per-step "
        "details). default = off",
        is_flag=True,
    )(func)
    func = click.option(
        "--versions-key",
        help="Top-level key used in --versions-yml. default = EEFINDER",
        default="EEFINDER",
    )(func)
    func = click.option(
        "--versions-yml",
        help="Write the detected dependency versions as YAML to this path.",
        default=None,
    )(func)
    return func


# -- helpers -----------------------------------------------------------------
def _resolve_prefix(prefix, genome_file=None):
    """Return the run prefix, deriving it from the genome file when omitted."""
    if prefix:
        return prefix
    if not genome_file:
        raise click.UsageError(
            "-pr/--prefix is required for this stage, since there is no "
            "-in/--genome_file to derive it from."
        )
    import re

    stem = re.sub(r"\..*", "", genome_file)
    return re.sub(r".*/", "", stem).rstrip("\n")


def _report(outputs, versions_yml, versions_key):
    """Print a stage's outputs and optionally write the versions file."""
    for name, path in outputs.existing().items():
        click.echo(f"  {name}: {path}")
    logger.info(
        f"{outputs.step_info.step} finished in "
        f"{outputs.step_info.total_time_minutes} min"
    )
    if versions_yml:
        write_versions_yml(versions_yml, versions_key)
        click.echo(f"  versions: {versions_yml}")


def _run_stage(stage, versions_yml=None, versions_key="EEFINDER"):
    """Run one stage, reporting its outputs and exiting non-zero on failure."""
    try:
        outputs = stage.run()
    except Exception as err:
        click.secho(f"Failed to run stage '{stage.name}': {err}", err=True, fg="red")
        sys.exit(1)
    _report(outputs, versions_yml, versions_key)
    return outputs


@click.group(cls=DefaultGroup, invoke_without_command=False)
@click.version_option(__version__)
def screening():
    """Run the EEfinder screening pipeline, whole or one stage at a time."""


# -- stage 0 -----------------------------------------------------------------
@screening.command(name="prepare")
@_database_options
@_outdir_options
@click.option(
    "-md",
    "--mode",
    help="Index for BLAST ('blastx') or DIAMOND (any sensitivity). default = blastx.",
    default="blastx",
)
@click.option(
    "-p",
    "--threads",
    help="Threads for multi-thread analysis, default = 1.",
    type=int,
    default=1,
)
@click.option(
    "-id",
    "--index_databases",
    help="Index databases?",
    is_flag=True,
)
@_common_options
def prepare_cmd(
    database,
    dbmetadata,
    hostgenesbaits,
    outdir,
    prefix,
    mode,
    threads,
    index_databases,
    versions_yml,
    versions_key,
    debug,
):
    """Validate the reference inputs and build the search indexes."""
    if debug:
        enable_debug()
    _run_stage(
        PrepareInputs(
            database=database,
            dbmetadata=dbmetadata,
            hostgenesbaits=hostgenesbaits,
            outdir=check_outdir(outdir),
            prefix=_resolve_prefix(prefix) if prefix else "eefinder",
            mode=mode,
            threads=threads,
            index_databases=index_databases,
        ),
        versions_yml,
        versions_key,
    )


# -- stage I -----------------------------------------------------------------
@screening.command(name="clean")
@click.option(
    "-in",
    "--genome_file",
    help="Input genome fasta file (nucleotides).",
    required=True,
)
@_outdir_options
@click.option(
    "-ln",
    "--length",
    help="Minimum length of contigs used for BLAST or DIAMOND, default = 10000.",
    type=int,
    default=10000,
)
@_common_options
def clean_cmd(genome_file, outdir, prefix, length, versions_yml, versions_key, debug):
    """Drop short contigs and prefix the FASTA headers."""
    if debug:
        enable_debug()
    _run_stage(
        DataCleaning(
            genome_file=genome_file,
            outdir=check_outdir(outdir),
            prefix=_resolve_prefix(prefix, genome_file),
            length=length,
        ),
        versions_yml,
        versions_key,
    )


# -- stage II ----------------------------------------------------------------
@screening.command(name="align")
@click.option(
    "--genome",
    help="Cleaned genome from the 'clean' stage. "
    "default = <outdir>/<prefix>.cleaned_genome.fa",
    default=None,
)
@click.option(
    "-db",
    "--database",
    help="Proteins from viruses or bacterias database .fasta file.",
    required=True,
)
@_outdir_options
@_search_options
@_common_options
def align_cmd(
    genome,
    database,
    outdir,
    prefix,
    mode,
    translation_method,
    threads,
    range_junction,
    versions_yml,
    versions_key,
    debug,
):
    """Search the genome against the reference proteins."""
    if debug:
        enable_debug()
    outdir = check_outdir(outdir)
    prefix = _resolve_prefix(prefix, genome)
    paths = ScreeningPaths(outdir=outdir, prefix=prefix)
    _run_stage(
        SequenceAlignment(
            genome=genome or paths.cleaned_genome,
            database=database,
            outdir=outdir,
            prefix=prefix,
            mode=mode,
            threads=threads,
            translation_method=translation_method,
            range_junction=range_junction,
        ),
        versions_yml,
        versions_key,
    )


# -- stage III ---------------------------------------------------------------
@screening.command(name="filter")
@click.option(
    "--candidates",
    help="Candidate sequences from the 'align' stage. "
    "default = <outdir>/<prefix>.ee_candidates.fa",
    default=None,
)
@click.option(
    "--ee-hits-filtered",
    help="Filtered database hits from the 'align' stage. "
    "default = <outdir>/<prefix>.ee_hits.filtered.tsv",
    default=None,
)
@click.option(
    "-bt",
    "--hostgenesbaits",
    help="Host genes baits proteins, used to filter putative EEs .fasta file.",
    required=True,
)
@_outdir_options
@_search_options
@_common_options
def filter_cmd(
    candidates,
    ee_hits_filtered,
    hostgenesbaits,
    outdir,
    prefix,
    mode,
    translation_method,
    threads,
    range_junction,
    versions_yml,
    versions_key,
    debug,
):
    """Drop candidates that match the host better than the database."""
    if debug:
        enable_debug()
    outdir = check_outdir(outdir)
    prefix = _resolve_prefix(prefix)
    paths = ScreeningPaths(outdir=outdir, prefix=prefix)
    _run_stage(
        PutativeElementsFilter(
            candidates=candidates or paths.ee_candidates,
            ee_hits_filtered=ee_hits_filtered or paths.ee_hits_filtered,
            hostgenesbaits=hostgenesbaits,
            outdir=outdir,
            prefix=prefix,
            mode=mode,
            threads=threads,
            translation_method=translation_method,
            range_junction=range_junction,
        ),
        versions_yml,
        versions_key,
    )


# -- stage IV ----------------------------------------------------------------
@screening.command(name="taxonomy")
@click.option(
    "--ee-hits-validated",
    help="Surviving hits from the 'filter' stage. "
    "default = <outdir>/<prefix>.ee_hits.validated.tsv",
    default=None,
)
@click.option(
    "-mt",
    "--dbmetadata",
    help="Proteins from viruses or bacterias metadata .csv file.",
    required=True,
)
@_outdir_options
@_common_options
def taxonomy_cmd(
    ee_hits_validated,
    dbmetadata,
    outdir,
    prefix,
    versions_yml,
    versions_key,
    debug,
):
    """Join the surviving hits to the reference metadata."""
    if debug:
        enable_debug()
    outdir = check_outdir(outdir)
    prefix = _resolve_prefix(prefix)
    paths = ScreeningPaths(outdir=outdir, prefix=prefix)
    _run_stage(
        TaxonomyAssignment(
            ee_hits_validated=ee_hits_validated or paths.ee_hits_validated,
            dbmetadata=dbmetadata,
            outdir=outdir,
            prefix=prefix,
        ),
        versions_yml,
        versions_key,
    )


# -- stage V -----------------------------------------------------------------
@screening.command(name="merge")
@click.option(
    "--taxonomy-signature",
    help="Per-hit signature from the 'taxonomy' stage. "
    "default = <outdir>/<prefix>.taxonomy_signature.csv",
    default=None,
)
@click.option(
    "--genome",
    help="Cleaned genome. default = <outdir>/<prefix>.cleaned_genome.fa",
    default=None,
)
@_outdir_options
@click.option(
    "-lm",
    "--limit",
    help="Limit of bases used to merge regions on bedtools merge, default = 1.",
    type=int,
    default=1,
)
@click.option(
    "-ml",
    "--merge_level",
    help="Taxonomy level to merge elements by genus or family, default = family",
    default="family",
    type=click.Choice(["family", "genus"]),
)
@_common_options
def merge_cmd(
    taxonomy_signature,
    genome,
    outdir,
    prefix,
    limit,
    merge_level,
    versions_yml,
    versions_key,
    debug,
):
    """Merge neighbouring same-taxon fragments into elements."""
    if debug:
        enable_debug()
    outdir = check_outdir(outdir)
    prefix = _resolve_prefix(prefix)
    paths = ScreeningPaths(outdir=outdir, prefix=prefix)
    _run_stage(
        MergeFragmentedElements(
            taxonomy_signature=taxonomy_signature or paths.taxonomy_signature,
            genome=genome or paths.cleaned_genome,
            outdir=outdir,
            prefix=prefix,
            limit=limit,
            merge_level=merge_level,
        ),
        versions_yml,
        versions_key,
    )


# -- post-processing ---------------------------------------------------------
@screening.command(name="postprocess")
@click.option(
    "--ee-elements",
    help="Merged elements from the 'merge' stage. "
    "default = <outdir>/<prefix>.elements.fa",
    default=None,
)
@click.option(
    "--ee-elements-tax",
    help="Per-element taxonomy table from the 'merge' stage. "
    "default = <outdir>/<prefix>.elements.tax.tsv",
    default=None,
)
@_outdir_options
@_postprocess_options
@_common_options
def postprocess_cmd(
    ee_elements,
    ee_elements_tax,
    outdir,
    prefix,
    mask_per,
    clean_masked,
    overlap,
    target_families,
    non_target_families,
    analysis,
    versions_yml,
    versions_key,
    debug,
):
    """Apply the repeat filter, resolve overlaps and write the GFF3."""
    if debug:
        enable_debug()
    _check_overlap_options(overlap, target_families, non_target_families)
    outdir = check_outdir(outdir)
    prefix = _resolve_prefix(prefix)
    paths = ScreeningPaths(outdir=outdir, prefix=prefix)
    _run_stage(
        PostProcessing(
            ee_elements=ee_elements or paths.ee_elements,
            ee_elements_tax=ee_elements_tax or paths.ee_elements_tax,
            outdir=outdir,
            prefix=prefix,
            clean_masked=clean_masked,
            mask_per=mask_per,
            overlap=overlap,
            target_families=target_families,
            non_target_families=non_target_families,
            analysis=analysis,
        ),
        versions_yml,
        versions_key,
    )


# -- stage VI ----------------------------------------------------------------
@screening.command(name="flanks")
@click.option(
    "--ee-elements",
    help="Merged elements. default = <outdir>/<prefix>.elements.fa",
    default=None,
)
@click.option(
    "--genome",
    help="Cleaned genome. default = <outdir>/<prefix>.cleaned_genome.fa",
    default=None,
)
@_outdir_options
@click.option(
    "-fl",
    "--flank",
    help="Length of flanking regions of Endogenous Elements to be extracted, default = 10000.",
    type=int,
    default=10000,
)
@_common_options
def flanks_cmd(
    ee_elements, genome, outdir, prefix, flank, versions_yml, versions_key, debug
):
    """Extract each element plus its flanking regions."""
    if debug:
        enable_debug()
    outdir = check_outdir(outdir)
    prefix = _resolve_prefix(prefix)
    paths = ScreeningPaths(outdir=outdir, prefix=prefix)
    _run_stage(
        FlanksExtraction(
            ee_elements=ee_elements or paths.ee_elements,
            genome=genome or paths.cleaned_genome,
            outdir=outdir,
            prefix=prefix,
            flank=flank,
        ),
        versions_yml,
        versions_key,
    )


def _check_overlap_options(overlap, target_families, non_target_families):
    """Reject an --overlap targets call without exactly one family list."""
    if overlap == "targets" and bool(target_families) == bool(non_target_families):
        click.secho(
            "--overlap targets requires exactly one of --target_families or "
            "--non_target_families.",
            err=True,
            fg="red",
        )
        sys.exit(1)


# -- every stage -------------------------------------------------------------
@screening.command(name="all")
@click.option(
    "-in",
    "--genome_file",
    help="Input genome fasta file (nucleotides).",
    required=True,
)
@_database_options
@_outdir_options
@_search_options
@_postprocess_options
@click.option(
    "-ln",
    "--length",
    help="Minimum length of contigs used for BLAST or DIAMOND, default = 10000.",
    type=int,
    default=10000,
)
@click.option(
    "-fl",
    "--flank",
    help="Length of flanking regions of Endogenous Elements to be extracted, default = 10000.",
    type=int,
    default=10000,
)
@click.option(
    "-lm",
    "--limit",
    help="Limit of bases used to merge regions on bedtools merge, default = 1.",
    type=int,
    default=1,
)
@click.option(
    "-ml",
    "--merge_level",
    help="Taxonomy level to merge elements by genus or family, default = family",
    default="family",
    type=click.Choice(["family", "genus"]),
)
@click.option(
    "-rm",
    "--removetmp",
    help="Remove temporary files generated through analysis?",
    is_flag=True,
)
@click.option(
    "-id",
    "--index_databases",
    help="Index databases?",
    is_flag=True,
)
@_common_options
def all_cmd(
    genome_file,
    database,
    dbmetadata,
    hostgenesbaits,
    outdir,
    prefix,
    mode,
    translation_method,
    threads,
    range_junction,
    mask_per,
    clean_masked,
    overlap,
    target_families,
    non_target_families,
    analysis,
    length,
    flank,
    limit,
    merge_level,
    removetmp,
    index_databases,
    versions_yml,
    versions_key,
    debug,
):
    """Run every stage of the screening pipeline."""
    if debug:
        enable_debug()
    _check_overlap_options(overlap, target_families, non_target_families)

    start_running_time = time.time()
    print_info = PaperInfo()
    print_info.print_start(__version__)

    system_info = collect_system_info()
    try:
        dependencies = collect_dependency_versions(find_env_yml())
    except Exception as err:  # pragma: no cover - defensive
        logger.warning(f"Could not collect dependency versions: {err}")
        dependencies = []
    report_run_context(system_info, dependencies)

    prefix = _resolve_prefix(prefix, genome_file)
    outdir = check_outdir(outdir)

    try:
        outputs = ScreeningPipeline(
            genome_file=genome_file,
            database=database,
            dbmetadata=dbmetadata,
            hostgenesbaits=hostgenesbaits,
            outdir=outdir,
            prefix=prefix,
            mode=mode,
            threads=threads,
            translation_method=translation_method,
            range_junction=range_junction,
            length=length,
            limit=limit,
            merge_level=merge_level,
            flank=flank,
            clean_masked=clean_masked,
            mask_per=mask_per,
            overlap=overlap,
            target_families=target_families,
            non_target_families=non_target_families,
            analysis=analysis,
            index_databases=index_databases,
            removetmp=removetmp,
        ).run()
    except Exception as err:
        click.secho(f"Failed to run screening: {err}", err=True, fg="red")
        sys.exit(1)

    _print_outputs(outputs, flank)
    if versions_yml:
        write_versions_yml(versions_yml, versions_key)

    run_arguments = RunArguments(
        genome_file=genome_file,
        prefix=prefix,
        outdir=outdir,
        database=database,
        dbmetadata=dbmetadata,
        baits=hostgenesbaits,
        mode=mode,
        length=length,
        flank=flank,
        limit=limit,
        range_junction=range_junction,
        mask_per=mask_per,
        clean_masked=clean_masked,
        threads=threads,
        removetmp=removetmp,
        index_databases=index_databases,
        merge_level=merge_level,
        analysis=analysis,
        overlap=overlap,
        target_families=list(target_families),
        non_target_families=list(non_target_families),
        translation_method=translation_method,
    )
    run_info = RunInfo.from_run(
        __version__,
        system_info,
        run_arguments,
        dependencies,
        start_running_time,
        time.time(),
        outputs.step_infos,
    )
    logger.debug(f"Writing run summary to {outdir}/eefinder.log")
    with open(f"{outdir}/eefinder.log", "w") as json_out:
        json.dump(asdict(run_info), json_out, indent=4)
    print_info.print_finish()


def _print_outputs(outputs, flank):
    """Print the published outputs of a full run."""
    print("")
    print("Output files:\n")
    print(
        f"{outputs.final_fasta} ----------------------------- Fasta file with Endogenous Elements nucleotide sequences."
    )
    print(
        f"{outputs.final_tax} ------------------------ TSV file with Endogenous Elements taxonomy."
    )
    print(
        f"{outputs.final_gff3} --------------------------- GFF3 annotation of Endogenous Elements."
    )
    print(
        f"{outputs.final_flanks} ---------------------- Fasta file with Endogenous Elements plus {flank}nt in each flanking regions."
    )
    if outputs.final_cleaned_fasta:
        print(
            f"{outputs.final_cleaned_fasta} --------------------- Fasta file with Cleaned Endogenous Elements."
        )
        print(
            f"{outputs.final_cleaned_tax} ---------------- TSV file with Cleaned Endogenous Elements."
        )
        print(
            f"{outputs.final_cleaned_gff3} ------------------- GFF3 annotation of Cleaned Endogenous Elements."
        )
    print("")
