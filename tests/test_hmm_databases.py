"""Unit tests for the profile-database assembly helpers.

The builds themselves need cd-hit/mafft/CIAlign/hmmbuild and are covered by the
integration tests; what is tested here is the logic that decides which proteins
reach a profile and what metadata the profile carries.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from eefinder.hmm_builder import BuilderSettings
from eefinder.hmm_databases import (
    HMM_METADATA_COLUMNS,
    HOST_METADATA_COLUMNS,
    SOURCE_NAMESPACES,
    is_unclassified,
    is_uninformative,
    plurality,
    read_builder_settings,
    read_sequences,
)
from eefinder.get_databases import METADATA_COLUMNS


class TestMetadataColumns:
    def test_screening_columns_come_first_and_unchanged(self):
        """GetFinalTaxonomy reads the taxonomy table by column position."""
        assert HMM_METADATA_COLUMNS[: len(METADATA_COLUMNS)] == METADATA_COLUMNS

    def test_extra_columns_are_appended_after_host(self):
        assert HMM_METADATA_COLUMNS[len(METADATA_COLUMNS) :] == [
            "Protein_votes",
            "Protein_agreement",
            "Profile_seqs",
            "Profile_length",
            "Source",
            "LCA_rank",
        ]

    def test_host_table_has_no_viral_taxonomy(self):
        assert "Family" not in HOST_METADATA_COLUMNS
        assert HOST_METADATA_COLUMNS[0] == "Accession"

    def test_every_source_has_an_id_namespace(self):
        assert SOURCE_NAMESPACES["ncbi-refseq"] == "NCBIREFSEQ"
        assert SOURCE_NAMESPACES["host"] == "HOST"


class TestPlurality:
    def test_returns_the_dominant_value(self):
        assert plurality(["Aedes"] * 9 + ["Culex"]) == "Aedes"

    def test_falls_back_when_no_value_dominates(self):
        assert (
            plurality(["Aedes"] * 5 + ["Culex"] * 5, fallback="Unclassified")
            == "Unclassified"
        )

    def test_ignores_unclassified_entries(self):
        assert plurality(["Unclassified", "", "Aedes", "Aedes"]) == "Aedes"

    def test_empty_input_returns_the_fallback(self):
        assert plurality([], fallback="Unknown") == "Unknown"

    def test_lower_threshold_accepts_a_simple_majority(self):
        assert plurality(["Aedes"] * 6 + ["Culex"] * 4, minimum=0.5) == "Aedes"


class TestUninformativeAndUnclassified:
    @pytest.mark.parametrize(
        "product",
        [
            "hypothetical protein",
            "Hypothetical protein LOC123",
            "uncharacterized protein",
            "uncharacterised protein",
            "unnamed protein product",
            "predicted protein",
            "putative uncharacterized protein",
            "Unknown",
        ],
    )
    def test_placeholder_products_are_uninformative(self, product):
        assert is_uninformative(product)

    @pytest.mark.parametrize("product", ["Glycoprotein", "actin", "RdRp"])
    def test_real_products_are_informative(self, product):
        assert not is_uninformative(product)

    @pytest.mark.parametrize("value", ["", "Unclassified", "unknown", "nan", None])
    def test_missing_taxonomy_is_unclassified(self, value):
        assert is_unclassified(value)

    def test_a_real_family_is_classified(self):
        assert not is_unclassified("Chuviridae")


def test_read_sequences_indexes_by_the_first_header_token(tmp_path):
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">PROT_A some product [organism=X]\nMKALL\n>PROT_B\nKKTTT\n")
    assert read_sequences(str(fasta)) == {"PROT_A": "MKALL", "PROT_B": "KKTTT"}


class TestReadBuilderSettings:
    def _write_log(self, tmp_path, settings):
        (tmp_path / "db.hmm").write_text("HMMER3/f\n")
        (tmp_path / "db.log").write_text(
            json.dumps({"builder_settings": settings.to_dict()})
        )
        return str(tmp_path / "db.hmm")

    def test_settings_round_trip_through_the_log(self, tmp_path):
        settings = BuilderSettings(cluster_identity=0.6, min_profile_seqs=4)
        loaded = read_builder_settings(self._write_log(tmp_path, settings))
        assert loaded == settings
        assert loaded.differences(settings) == []

    def test_missing_log_returns_none_rather_than_a_mismatch(self, tmp_path):
        (tmp_path / "db.hmm").write_text("HMMER3/f\n")
        assert read_builder_settings(str(tmp_path / "db.hmm")) is None

    def test_unknown_keys_in_the_log_are_ignored(self, tmp_path):
        (tmp_path / "db.hmm").write_text("HMMER3/f\n")
        (tmp_path / "db.log").write_text(
            json.dumps({"builder_settings": {"cluster_identity": 0.8, "future": 1}})
        )
        loaded = read_builder_settings(str(tmp_path / "db.hmm"))
        assert loaded.cluster_identity == 0.8

    def test_corrupt_log_returns_none(self, tmp_path):
        (tmp_path / "db.hmm").write_text("HMMER3/f\n")
        (tmp_path / "db.log").write_text("not json")
        assert read_builder_settings(str(tmp_path / "db.hmm")) is None


class TestVogdbSourceIntegration:
    """The VOGDB source mapped onto EEfinder metadata, without any download."""

    @pytest.fixture
    def vogdb_build(self, tmp_path, monkeypatch, vogdb_tables):
        from eefinder import hmm_databases
        from eefinder.hmm_databases import GetViralHmmDatabase
        from eefinder.hmm_sources import VogdbFilters

        monkeypatch.setattr(
            hmm_databases,
            "download_vogdb",
            lambda workdir, release="latest", metadata_only=True: vogdb_tables,
        )
        monkeypatch.setattr(hmm_databases, "vogdb_release", lambda workdir: "236")

        outdir = tmp_path / "db"
        GetViralHmmDatabase(
            outdir=str(outdir),
            prefix="vogdb_hmm",
            sources=["vogdb"],
            metadata_only=True,
            vogdb_filters=VogdbFilters(min_profile_seqs=3),
        )
        return outdir

    def test_metadata_only_writes_no_profile_database(self, vogdb_build):
        assert not (vogdb_build / "vogdb_hmm.hmm").exists()
        assert (vogdb_build / "vogdb_hmm.csv").is_file()
        assert (vogdb_build / "vogdb_hmm.log").is_file()

    def test_rows_carry_namespaced_ids_and_exact_columns(self, vogdb_build):
        frame = pd.read_csv(vogdb_build / "vogdb_hmm.csv")
        assert list(frame.columns) == HMM_METADATA_COLUMNS
        assert set(frame.Accession) == {"VOGDB__VOG00002", "VOGDB__VOG00005"}
        assert set(frame.Source) == {"vogdb"}
        assert set(frame.LCA_rank) == {"genus"}

    def test_taxonomy_and_protein_are_mapped(self, vogdb_build):
        frame = pd.read_csv(vogdb_build / "vogdb_hmm.csv").set_index("Accession")
        row = frame.loc["VOGDB__VOG00002"]
        assert row.Family == "Chuviridae"
        assert row.Genus == "Mivirus"
        assert row.Protein == "Glycoprotein"  # canonicalised by the bundled map
        assert row.Host == "Aedes aegypti"

    def test_host_is_taken_from_the_members_that_have_one(self, vogdb_build):
        """`vogdb.host.txt` fills the host column sparsely; that is not an error.

        VOG00005 has two member species and only one of them records a host, so
        the known host is reported rather than discarded. ``Unknown`` is reserved
        for groups where no member has one at all -- ``Host`` is an optional
        field (proposal.md §1.1).
        """
        frame = pd.read_csv(vogdb_build / "vogdb_hmm.csv").set_index("Accession")
        assert frame.loc["VOGDB__VOG00005"].Host == "Aedes aegypti"

    def test_log_records_the_filters_and_counts(self, vogdb_build):
        log = json.loads((vogdb_build / "vogdb_hmm.log").read_text())
        assert log["arguments"]["sources"] == ["vogdb"]
        assert log["arguments"]["metadata_only"] is True
        assert log["arguments"]["vogdb_filters"]["min_lca_rank"] == "genus"
        assert log["profile_counts"]["profiles_built"] == 2

    def test_metadata_only_is_rejected_for_the_building_source(self, tmp_path):
        from eefinder.hmm_databases import GetViralHmmDatabase

        with pytest.raises(ValueError, match="metadata-only"):
            GetViralHmmDatabase(
                outdir=str(tmp_path / "x"),
                sources=["ncbi-refseq"],
                metadata_only=True,
            )

    def test_unknown_source_is_rejected(self, tmp_path):
        from eefinder.hmm_databases import GetViralHmmDatabase

        with pytest.raises(ValueError, match="Unknown source"):
            GetViralHmmDatabase(outdir=str(tmp_path / "x"), sources=["efam"])
