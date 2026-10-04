"""Unit tests for the eefinder.stages package (no external binaries)."""

from __future__ import annotations
from dataclasses import dataclass
import pytest
from eefinder.models import GENOME_FASTA, ModelError
from eefinder.stages import STAGES, STAGES_BY_NAME, DataCleaning, ScreeningPaths
from eefinder.stages.base import Stage, StageOutputs
from eefinder.stages.prepare import PrepareInputs, check_database_pairing


def test_every_stage_has_a_unique_subcommand_name():
    names = [stage.name for stage in STAGES]

    assert len(names) == len(set(names))
    assert set(names) == set(STAGES_BY_NAME)


def test_the_expected_stages_are_exposed():
    assert set(STAGES_BY_NAME) == {
        "prepare",
        "clean",
        "align",
        "filter",
        "taxonomy",
        "merge",
        "postprocess",
        "flanks",
        "all",
    }


def test_paths_are_derived_from_outdir_and_prefix():
    paths = ScreeningPaths(outdir="out", prefix="run")

    assert paths.cleaned_genome == "out/run.cleaned_genome.fa"
    assert paths.ee_hits_filtered == "out/run.ee_hits.filtered.tsv"
    assert paths.ee_elements_tax == "out/run.elements.tax.tsv"
    assert paths.final_tax == "out/run.EEs.tax.tsv"


def test_final_outputs_and_intermediates_do_not_overlap():
    paths = ScreeningPaths(outdir="out", prefix="run")

    assert not set(paths.final_outputs()) & set(paths.intermediates())


# -- the base class ----------------------------------------------------------
@dataclass
class _DemoOutputs(StageOutputs):
    genome: str


class _DemoStage(Stage):
    """Stage that writes whatever it was told to, for testing the base class."""

    name = "demo"
    title = "Demo stage"
    stage_id = "X"
    models = {"genome": GENOME_FASTA}

    def __init__(self, outdir, prefix, content=">a\nACGT\n"):
        super().__init__(outdir, prefix)
        self.content = content

    def _execute(self):
        with open(self.paths.cleaned_genome, "w") as handle:
            handle.write(self.content)
        return _DemoOutputs(step_info=None, genome=self.paths.cleaned_genome), "did it"


def test_run_fills_in_step_info(tmp_path):
    result = _DemoStage(str(tmp_path / "out"), "run").run()

    assert result.step_info.step == "Demo stage"
    assert result.step_info.message == "did it"
    assert float(result.step_info.total_time_minutes) >= 0


def test_run_creates_the_output_directory(tmp_path):
    outdir = tmp_path / "nested" / "out"

    _DemoStage(str(outdir), "run").run()

    assert outdir.is_dir()


def test_run_enforces_the_declared_models(tmp_path):
    # Not a FASTA, so the declared model rejects it.
    with pytest.raises(ModelError, match="cleaned genome"):
        _DemoStage(str(tmp_path / "out"), "run", content="ACGT\n").run()


def test_outputs_report_their_files(tmp_path):
    result = _DemoStage(str(tmp_path / "out"), "run").run()

    assert list(result.files()) == ["genome"]
    assert list(result.existing()) == ["genome"]


def test_require_names_every_missing_input(tmp_path):
    stage = _DemoStage(str(tmp_path / "out"), "run")

    with pytest.raises(FileNotFoundError, match="absent_a.*absent_b"):
        stage._require("absent_a", "absent_b")


def test_rename_moves_a_step_output(tmp_path):
    stage = _DemoStage(str(tmp_path / "out"), "run")
    produced = tmp_path / "produced.txt"
    produced.write_text("x")
    canonical = tmp_path / "canonical.txt"

    assert stage._rename(str(produced), str(canonical)) == str(canonical)
    assert canonical.read_text() == "x"
    assert not produced.exists()


def test_rename_is_a_noop_for_the_same_path(tmp_path):
    stage = _DemoStage(str(tmp_path / "out"), "run")
    path = tmp_path / "same.txt"

    # Not created on disk: an identical path must not be touched at all.
    assert stage._rename(str(path), str(path)) == str(path)


def test_rename_reports_a_missing_step_output(tmp_path):
    stage = _DemoStage(str(tmp_path / "out"), "run")

    with pytest.raises(FileNotFoundError, match="to have been produced"):
        stage._rename(str(tmp_path / "absent"), str(tmp_path / "target"))


# -- data cleaning (pure Python, so it runs without binaries) ----------------
def test_data_cleaning_prefixes_and_filters(tmp_path, fasta_factory):
    genome = fasta_factory(
        "genome.fa", {"long_contig": "A" * 100, "short_contig": "A" * 10}
    )

    result = DataCleaning(
        genome_file=str(genome), outdir=str(tmp_path / "out"), prefix="run", length=50
    ).run()

    assert result.cleaned_genome == f"{tmp_path}/out/run.cleaned_genome.fa"
    assert result.kept == 1
    assert result.total == 2
    content = open(result.cleaned_genome).read()
    assert ">run/long_contig" in content
    assert "short_contig" not in content


def test_data_cleaning_reports_a_missing_genome(tmp_path):
    stage = DataCleaning(
        genome_file=str(tmp_path / "absent.fa"),
        outdir=str(tmp_path / "out"),
        prefix="run",
    )

    with pytest.raises(FileNotFoundError, match="missing required input"):
        stage.run()


# -- input preparation -------------------------------------------------------
def test_prepare_inputs_validates_without_indexing(
    tmp_path, fasta_factory, taxonomy_csv
):
    database = fasta_factory("db.fa", {"PROT_A": "MKV"})
    baits = fasta_factory("baits.fa", {"HOST_1": "MKV"})

    result = PrepareInputs(
        database=str(database),
        dbmetadata=str(taxonomy_csv),
        hostgenesbaits=str(baits),
        outdir=str(tmp_path / "out"),
        prefix="run",
        index_databases=False,
    ).run()

    assert result.indexed is False
    assert result.database == str(database)


def test_prepare_inputs_rejects_a_non_fasta_database(tmp_path, taxonomy_csv):
    database = tmp_path / "db.fa"
    database.write_text("MKV\n")
    baits = tmp_path / "baits.fa"
    baits.write_text(">HOST_1\nMKV\n")

    stage = PrepareInputs(
        database=str(database),
        dbmetadata=str(taxonomy_csv),
        hostgenesbaits=str(baits),
        outdir=str(tmp_path / "out"),
        prefix="run",
    )

    with pytest.raises(ValueError, match="does not look like a FASTA"):
        stage.run()


def test_prepare_inputs_rejects_an_empty_database(tmp_path, taxonomy_csv):
    database = tmp_path / "db.fa"
    database.touch()
    baits = tmp_path / "baits.fa"
    baits.write_text(">HOST_1\nMKV\n")

    stage = PrepareInputs(
        database=str(database),
        dbmetadata=str(taxonomy_csv),
        hostgenesbaits=str(baits),
        outdir=str(tmp_path / "out"),
        prefix="run",
    )

    with pytest.raises(ValueError, match="is empty"):
        stage.run()


def test_database_pairing_accepts_a_fully_described_database(
    tmp_path, fasta_factory, taxonomy_csv
):
    # Both accessions are in the metadata fixture.
    database = fasta_factory("db.fa", {"PROT_A": "MKV", "PROT_C": "MKV"})

    assert check_database_pairing(str(database), str(taxonomy_csv)) == 2


def test_database_pairing_allows_extra_metadata_rows(
    tmp_path, fasta_factory, taxonomy_csv
):
    # The fixture also describes PROT_C; a database may be a subset.
    database = fasta_factory("db.fa", {"PROT_A": "MKV"})

    assert check_database_pairing(str(database), str(taxonomy_csv)) == 1


def test_database_pairing_rejects_a_protein_without_metadata(
    tmp_path, fasta_factory, taxonomy_csv
):
    database = fasta_factory("db.fa", {"PROT_A": "MKV", "PROT_ZZ": "MKV"})

    with pytest.raises(ValueError, match="PROT_ZZ"):
        check_database_pairing(str(database), str(taxonomy_csv))


def test_database_pairing_reports_how_many_are_missing(
    tmp_path, fasta_factory, taxonomy_csv
):
    database = fasta_factory("db.fa", {f"PROT_MISSING_{i}": "MKV" for i in range(9)})

    with pytest.raises(ValueError, match=r"9 of the 9 protein\(s\)"):
        check_database_pairing(str(database), str(taxonomy_csv))


def test_database_pairing_truncates_the_example_list(
    tmp_path, fasta_factory, taxonomy_csv
):
    database = fasta_factory("db.fa", {f"PROT_MISSING_{i}": "MKV" for i in range(9)})

    with pytest.raises(ValueError, match=r"\(4 more\)"):
        check_database_pairing(str(database), str(taxonomy_csv))


def test_prepare_inputs_rejects_an_undescribed_database(
    tmp_path, fasta_factory, taxonomy_csv
):
    database = fasta_factory("db.fa", {"PROT_ZZ": "MKV"})
    baits = fasta_factory("baits.fa", {"HOST_1": "MKV"})

    stage = PrepareInputs(
        database=str(database),
        dbmetadata=str(taxonomy_csv),
        hostgenesbaits=str(baits),
        outdir=str(tmp_path / "out"),
        prefix="run",
    )

    with pytest.raises(ValueError, match="have no row in"):
        stage.run()
