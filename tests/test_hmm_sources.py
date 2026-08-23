"""Unit tests for the VOGDB adapter — fixtures only, nothing is downloaded.

The fixture tables reproduce the real headers and row shapes of a VOGDB release
(the ``#GroupName``-style headers are quoted verbatim), so the parsers are
tested against the format they will meet without pulling a 554 MB archive onto a
storage-constrained machine.
"""

from __future__ import annotations

import io
import json
import subprocess
import tarfile

import pytest

from eefinder import hmm_sources
from eefinder.hmm_sources import (
    FilterCounts,
    clean_vogdb_description,
    extract_profiles_from_hmm,
    load_neordrp,
    load_rvdb,
    protein_from_keywords,
    resolve_taxids,
    VogRecord,
    VogdbFilters,
    download_vogdb,
    extract_profiles,
    filter_vogs,
    load_vogdb,
    parse_annotations,
    parse_hosts,
    parse_lca,
    parse_members,
    parse_species,
    parse_virusonly,
    profile_length,
    vogdb_release,
)


class TestParsers:
    def test_annotations(self, vogdb_tables):
        records = parse_annotations(vogdb_tables["vog.annotations.tsv.gz"])
        assert set(records) == {f"VOG0000{n}" for n in range(1, 6)}
        assert records["VOG00004"].protein_count == 287
        assert records["VOG00002"].description == "glycoprotein"

    def test_lca_lineage_is_split_on_semicolons(self, vogdb_tables):
        lineages = parse_lca(vogdb_tables["vog.lca.tsv.gz"])
        names, taxid = lineages["VOG00001"]
        assert names == ["Viruses", "Riboviria", "Chuviridae", "Mivirus"]
        assert taxid == "111111"

    def test_members_keep_only_the_taxid(self, vogdb_tables):
        members = parse_members(vogdb_tables["vog.members.tsv.gz"])
        assert members["VOG00001"] == ["1001", "1002"]

    def test_virusonly_flags_per_stringency(self, vogdb_tables):
        flags = parse_virusonly(vogdb_tables["vog.virusonly.tsv.gz"])
        assert flags["VOG00005"] == {"high": False, "medium": True, "low": True}

    def test_species_are_keyed_by_taxid(self, vogdb_tables):
        assert parse_species(vogdb_tables["vogdb.species.txt"])["1003"] == "Zika virus"

    def test_hosts_carry_the_phage_flag(self, vogdb_tables):
        hosts = parse_hosts(vogdb_tables["vogdb.host.txt"])
        assert hosts["1001"] == (False, "Aedes aegypti")
        assert hosts["2001"] == (True, "")

    def test_load_merges_every_table(self, vogdb_tables):
        records = load_vogdb(vogdb_tables)
        record = records["VOG00002"]
        assert record.lineage[-1] == "Mivirus"
        assert record.member_taxids == ["1001", "1002"]
        assert record.virus_only["medium"] is True

    def test_release_is_read_from_the_tree(self, vogdb_tables, tmp_path):
        assert vogdb_release(str(tmp_path)) == "236"


class TestVogRecord:
    def test_genus_and_family_come_from_the_lineage(self):
        record = VogRecord(
            group="VOG1", lineage=["Viruses", "Riboviria", "Chuviridae", "Mivirus"]
        )
        assert record.genus_family == ("Mivirus", "Chuviridae")
        assert record.lca_rank == "genus"

    def test_family_only_lineage_resolves_to_family(self):
        record = VogRecord(group="VOG1", lineage=["Viruses", "Flaviviridae"])
        assert record.lca_rank == "family"

    def test_broad_lineage_resolves_to_nothing(self):
        record = VogRecord(group="VOG1", lineage=["Viruses", "Duplodnaviria"])
        assert record.lca_rank == "none"

    def test_phage_clades_are_recognised(self):
        assert VogRecord(
            group="V", lineage=["Viruses", "Caudoviricetes"]
        ).is_phage_clade()
        assert not VogRecord(
            group="V", lineage=["Viruses", "Chuviridae"]
        ).is_phage_clade()

    def test_uninformative_by_category_or_description(self):
        assert VogRecord(group="V", category="Xu").is_uninformative()
        assert VogRecord(
            group="V", category="Xr", description="REFSEQ hypothetical protein"
        ).is_uninformative()
        assert not VogRecord(
            group="V", category="Xr", description="glycoprotein"
        ).is_uninformative()


class TestFilters:
    def _filter(self, vogdb_tables, **overrides):
        records = load_vogdb(vogdb_tables)
        hosts = parse_hosts(vogdb_tables["vogdb.host.txt"])
        return filter_vogs(records, hosts, VogdbFilters(**overrides))

    def test_default_filters_keep_only_usable_groups(self, vogdb_tables):
        kept, counts = self._filter(vogdb_tables, min_profile_seqs=3)
        # VOG00001 uninformative, VOG00003 family-level LCA and too small,
        # VOG00004 phage and not virus-only.
        assert {record.group for record in kept} == {"VOG00002", "VOG00005"}
        assert counts.total == 5
        assert counts.kept == 2

    def test_phage_groups_are_dropped_by_the_eukaryotic_filter(self, vogdb_tables):
        _, counts = self._filter(vogdb_tables, virus_only_stringency="none")
        assert counts.dropped_not_eukaryotic == 1

    def test_virus_only_stringency_is_applied(self, vogdb_tables):
        kept_medium, _ = self._filter(
            vogdb_tables,
            virus_only_stringency="medium",
            min_lca_rank="none",
            min_profile_seqs=1,
            exclude_uninformative=False,
        )
        kept_high, _ = self._filter(
            vogdb_tables,
            virus_only_stringency="high",
            min_lca_rank="none",
            min_profile_seqs=1,
            exclude_uninformative=False,
        )
        # VOG00005 is virus-only at medium but not at high stringency.
        assert "VOG00005" in {r.group for r in kept_medium}
        assert "VOG00005" not in {r.group for r in kept_high}

    def test_min_lca_rank_enforces_the_mandatory_taxonomy_field(self, vogdb_tables):
        genus, _ = self._filter(
            vogdb_tables,
            min_lca_rank="genus",
            min_profile_seqs=1,
            exclude_uninformative=False,
        )
        family, _ = self._filter(
            vogdb_tables,
            min_lca_rank="family",
            min_profile_seqs=1,
            exclude_uninformative=False,
        )
        # VOG00003 resolves only to Flaviviridae: kept at family, dropped at genus.
        assert "VOG00003" in {r.group for r in family}
        assert "VOG00003" not in {r.group for r in genus}

    def test_taxon_filter_matches_the_lineage(self, vogdb_tables):
        kept, counts = self._filter(
            vogdb_tables,
            taxon="Flaviviridae",
            min_lca_rank="family",
            min_profile_seqs=1,
            exclude_uninformative=False,
        )
        assert {r.group for r in kept} == {"VOG00003"}
        assert counts.dropped_taxon == 4

    def test_small_groups_are_dropped(self, vogdb_tables):
        _, counts = self._filter(
            vogdb_tables,
            min_profile_seqs=10,
            min_lca_rank="none",
            exclude_uninformative=False,
        )
        assert counts.dropped_too_small >= 1

    def test_counts_are_reported_as_a_dataclass(self, vogdb_tables):
        _, counts = self._filter(vogdb_tables)
        assert isinstance(counts, FilterCounts)
        assert (
            counts.kept
            + sum(
                (
                    counts.dropped_taxon,
                    counts.dropped_not_eukaryotic,
                    counts.dropped_not_virus_only,
                    counts.dropped_lca_rank,
                    counts.dropped_too_small,
                    counts.dropped_uninformative,
                )
            )
            == counts.total
        )


class TestProfileExtraction:
    def _archive(self, tmp_path, groups):
        path = tmp_path / "vog.hmm.tar.gz"
        with tarfile.open(path, "w:gz") as archive:
            for group in groups:
                body = (
                    f"HMMER3/f [3.4 | Aug 2023]\nNAME  {group}\nLENG  145\n"
                    "ALPH  amino\nHMM   A C\n//\n"
                ).encode()
                info = tarfile.TarInfo(f"hmm/{group}.hmm")
                info.size = len(body)
                archive.addfile(info, io.BytesIO(body))
        return str(path)

    def test_only_the_kept_groups_are_extracted_and_renamed(self, tmp_path):
        archive = self._archive(tmp_path, ["VOG00001", "VOG00002", "VOG00003"])
        extracted = extract_profiles(
            archive, {"VOG00001", "VOG00003"}, str(tmp_path / "out")
        )
        assert set(extracted) == {"VOG00001", "VOG00003"}
        text = open(extracted["VOG00001"]).read()
        assert "NAME  VOGDB__VOG00001\n" in text
        assert "NAME  VOG00001\n" not in text

    def test_profile_length_is_read_back(self, tmp_path):
        archive = self._archive(tmp_path, ["VOG00001"])
        extracted = extract_profiles(archive, {"VOG00001"}, str(tmp_path / "out"))
        assert profile_length(extracted["VOG00001"]) == 145


class TestDownload:
    def test_a_failing_mirror_falls_back_to_the_next(self, tmp_path, monkeypatch):
        import urllib.error

        attempts = []

        def fake_fetch(url, dest):
            attempts.append(url)
            if url.startswith("https://broken"):
                raise urllib.error.URLError("mirror down")
            if url.endswith(".md5"):
                raise urllib.error.URLError("no sidecar")
            open(dest, "w").write("236\n")

        monkeypatch.setattr(hmm_sources, "_fetch", fake_fetch)
        paths = download_vogdb(
            str(tmp_path),
            metadata_only=True,
            mirrors=("https://broken.example", "https://working.example"),
        )
        assert set(paths) >= set(hmm_sources.VOGDB_METADATA_FILES)
        assert any(url.startswith("https://working") for url in attempts)

    def test_metadata_only_never_requests_the_profile_archive(
        self, tmp_path, monkeypatch
    ):
        import urllib.error

        requested = []

        def fake_fetch(url, dest):
            requested.append(url)
            if url.endswith(".md5"):
                raise urllib.error.URLError("no sidecar")
            open(dest, "w").write("236\n")

        monkeypatch.setattr(hmm_sources, "_fetch", fake_fetch)
        download_vogdb(str(tmp_path), metadata_only=True, mirrors=("https://m",))
        assert not any(hmm_sources.VOGDB_PROFILE_ARCHIVE in url for url in requested)

    def test_every_mirror_failing_is_an_error(self, tmp_path, monkeypatch):
        import urllib.error

        monkeypatch.setattr(
            hmm_sources,
            "_fetch",
            lambda url, dest: (_ for _ in ()).throw(urllib.error.URLError("down")),
        )
        with pytest.raises(RuntimeError, match="could not download VOGDB"):
            download_vogdb(str(tmp_path), mirrors=("https://a", "https://b"))


class TestDescriptionCleaning:
    def test_uniprot_provenance_prefix_is_stripped(self):
        assert (
            clean_vogdb_description(
                "sp|P84400|MB43_EHV1V Membrane protein UL43 homolog"
            )
            == "Membrane protein UL43 homolog"
        )

    def test_case_variants_are_handled(self):
        assert (
            clean_vogdb_description("Sp|B8XTP8|POLG_COSAA Genome polyprotein")
            == "Genome polyprotein"
        )

    def test_plain_descriptions_are_untouched(self):
        assert clean_vogdb_description("terminase large subunit") == (
            "terminase large subunit"
        )

    def test_no_pipe_survives(self):
        """`|` is a field separator in bed.py, so it must not reach the table."""
        assert "|" not in clean_vogdb_description("tr|X|Y some|protein")


# --------------------------------------------------------------------------- #
# RVDB-prot
# --------------------------------------------------------------------------- #
class TestRvdbKeywordNames:
    """RVDB annotates a family with a bag of scored tokens, not a phrase."""

    def test_a_token_set_rule_yields_the_canonical_name(self):
        keywords = [
            ("rna", 2311),
            ("polymerase", 1000),
            ("viral", 1000),
            ("dependent", 863),
        ]
        assert protein_from_keywords(keywords) == "RdRp"

    def test_reverse_transcriptase_is_recognised(self):
        assert (
            protein_from_keywords([("reverse", 90), ("transcriptase", 88)])
            == "Reverse Transcriptase"
        )

    def test_unmatched_keywords_are_joined_by_frequency(self):
        keywords = [("tail", 40), ("fiber", 38), ("viral", 90), ("protein", 200)]
        assert protein_from_keywords(keywords) == "Tail fiber"

    def test_only_stopwords_is_unknown(self):
        assert protein_from_keywords([("viral", 90), ("protein", 200)]) == "Unknown"

    def test_no_keywords_is_unknown(self):
        assert protein_from_keywords([]) == "Unknown"

    def test_the_result_is_independent_of_input_order(self):
        keywords = [("capsid", 50), ("major", 30)]
        assert protein_from_keywords(keywords) == protein_from_keywords(
            list(reversed(keywords))
        )


class TestLoadRvdb:
    @pytest.fixture
    def rvdb_sqlite(self, tmp_path):
        """A miniature copy of the RVDB-prot annotation schema."""
        import sqlite3

        path = tmp_path / "rvdb.sqlite"
        connection = sqlite3.connect(path)
        connection.executescript("""
            create table family (id integer, size integer, nbseq integer, LCAtaxid integer);
            create table keyword (id integer, str text);
            create table fam_kw_ref (famId integer, kwId integer, freq integer);
            create table fam_kw_seqnames (famId integer, kwId integer, freq integer);
            insert into family values (1, 1824, 881, 11292), (2, 12, 4, 10239);
            insert into keyword values (1,'rna'), (2,'polymerase'), (3,'dependent'), (4,'capsid');
            insert into fam_kw_ref values (1,1,2311), (1,2,1000), (1,3,863);
            insert into fam_kw_seqnames values (2,4,50);
            """)
        connection.commit()
        connection.close()
        return str(path)

    def test_families_become_records_with_composed_names(
        self, rvdb_sqlite, monkeypatch
    ):
        monkeypatch.setattr(
            hmm_sources,
            "resolve_taxids",
            lambda taxids, cache_path=None: {
                "11292": {
                    "species": "Rabies lyssavirus",
                    "genus": "Lyssavirus",
                    "family": "Rhabdoviridae",
                    "mol_type": "ssRNA(-)",
                }
            },
        )
        records = load_rvdb(rvdb_sqlite)
        assert set(records) == {"FAM000001", "FAM000002"}

        first = records["FAM000001"]
        assert first.description == "RdRp"
        assert first.genus_family == ("Lyssavirus", "Rhabdoviridae")
        assert first.species == "Rabies lyssavirus"
        assert first.mol_type == "ssRNA(-)"
        assert first.protein_count == 881

    def test_families_whose_lca_is_root_get_no_taxonomy(self, rvdb_sqlite, monkeypatch):
        """LCA = Viruses is 954 of RVDB 32.0's 13,679 families."""
        monkeypatch.setattr(
            hmm_sources, "resolve_taxids", lambda taxids, cache_path=None: {}
        )
        records = load_rvdb(rvdb_sqlite)
        assert records["FAM000002"].lca_rank == "none"

    def test_sequence_name_keywords_are_the_fallback(self, rvdb_sqlite, monkeypatch):
        monkeypatch.setattr(
            hmm_sources, "resolve_taxids", lambda taxids, cache_path=None: {}
        )
        records = load_rvdb(rvdb_sqlite)
        assert records["FAM000002"].description == "Capsid Protein"


class TestResolveTaxids:
    def test_the_cache_is_reused_and_extended(self, tmp_path, monkeypatch):
        cache = tmp_path / "taxonomy.json"
        cache.write_text(json.dumps({"11292": {"genus": "Lyssavirus"}}))

        calls = []

        def fake_run(command, **kwargs):
            calls.append(command)
            payload = {
                "taxonomy": {
                    "query": ["10244"],
                    "classification": {
                        "family": {"name": "Poxviridae"},
                        "genus": {"name": "Orthopoxvirus"},
                        "species": {"name": "Orthopoxvirus monkeypox"},
                    },
                    "genomic_moltype": "dsDNA",
                }
            }
            return subprocess.CompletedProcess(
                command, 0, stdout=json.dumps(payload) + "\n", stderr=""
            )

        monkeypatch.setattr(hmm_sources.subprocess, "run", fake_run)
        resolved = resolve_taxids(["11292", "10244"], cache_path=str(cache))

        assert resolved["11292"]["genus"] == "Lyssavirus"  # from the cache
        assert resolved["10244"]["family"] == "Poxviridae"  # newly resolved
        assert resolved["10244"]["mol_type"] == "dsDNA"
        # Only the missing taxid was queried.
        assert calls and "11292" not in calls[0]
        assert json.loads(cache.read_text())["10244"]["genus"] == "Orthopoxvirus"

    def test_invalid_taxids_are_not_queried(self, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(
            hmm_sources.subprocess,
            "run",
            lambda command, **kwargs: calls.append(command)
            or subprocess.CompletedProcess(command, 0, stdout="", stderr=""),
        )
        assert resolve_taxids(["0", "", "not-a-taxid"]) == {}
        assert not calls


# --------------------------------------------------------------------------- #
# NeoRdRp
# --------------------------------------------------------------------------- #
NEORDRP_ANNOTATION = (
    "\t".join(
        [
            "Seed_RdRp_Sequence Name",
            "Seed_RdRp_HMM Name",
            "InterProScan_CDD_Signature Description",
            "hmmsearch_RdRp-scan_description_of_target",
        ]
    )
    + "\n"
    + "\n".join(
        [
            "AAA43036.1\t810.aln_RDRP_0.25_1351-2305\tCardiovirus_RdRp\tfamily:Picornaviridae",
            "AAA43037.1\t810.aln_RDRP_0.25_1351-2305\tCardiovirus_RdRp\t-",
            "AAA46565.1\t-\tps-ssRNAv_Solemoviridae_RdRp\t-",
            "AAA46566.1\t2505.aln_RDRP_0.25_282-734\tps-ssRNAv_Bromoviridae_RdRp\t-",
            "AAA46567.1\t999.aln_RDRP_0.25_1-100\t-\t-",
        ]
    )
)


class TestLoadNeoRdRp:
    @pytest.fixture
    def annotation(self, tmp_path):
        path = tmp_path / "neordrp.annotation.tsv"
        path.write_text(NEORDRP_ANNOTATION + "\n")
        return str(path)

    def test_one_record_per_profile_not_per_seed(self, annotation):
        records = load_neordrp(annotation)
        # The seed with no profile ("-") is skipped.
        assert set(records) == {
            "810.aln_RDRP_0.25_1351-2305",
            "2505.aln_RDRP_0.25_282-734",
            "999.aln_RDRP_0.25_1-100",
        }

    def test_taxonomy_is_inferred_from_signature_names(self, annotation):
        records = load_neordrp(annotation)
        assert records["810.aln_RDRP_0.25_1351-2305"].genus_family == (
            "Cardiovirus",
            "Picornaviridae",
        )
        assert records["2505.aln_RDRP_0.25_282-734"].genus_family == (
            "",
            "Bromoviridae",
        )

    def test_profiles_without_a_signature_have_no_taxonomy(self, annotation):
        """Most of NeoRdRp is in this state; --min-lca-rank then drops them."""
        assert (
            records_rank(load_neordrp(annotation), "999.aln_RDRP_0.25_1-100") == "none"
        )

    def test_every_profile_is_an_rdrp_of_an_rna_virus(self, annotation):
        records = load_neordrp(annotation)
        assert {record.description for record in records.values()} == {"RdRp"}
        assert {record.mol_type for record in records.values()} == {"RNA"}

    def test_a_changed_layout_is_reported(self, tmp_path):
        path = tmp_path / "broken.tsv"
        path.write_text("some\tother\tcolumns\n1\t2\t3\n")
        with pytest.raises(ValueError, match="Seed_RdRp_HMM Name"):
            load_neordrp(str(path))


def records_rank(records, group):
    return records[group].lca_rank


class TestExtractFromConcatenatedHmm:
    def test_only_the_kept_profiles_are_written_and_renamed(self, tmp_path):
        database = tmp_path / "all.hmm"
        database.write_text(
            "".join(
                f"HMMER3/f [3.4 | Aug 2023]\nNAME  {name}\nLENG  10\nHMM  A C\n//\n"
                for name in ("FAM1", "FAM2", "FAM3")
            )
        )
        extracted = extract_profiles_from_hmm(
            str(database),
            {"FAM1": "RVDB__FAM000001", "FAM3": "RVDB__FAM000003"},
            str(tmp_path / "out"),
        )
        assert set(extracted) == {"FAM1", "FAM3"}
        assert "NAME  RVDB__FAM000001\n" in open(extracted["FAM1"]).read()

    def test_xz_archives_are_read_directly(self, tmp_path):
        import lzma

        database = tmp_path / "all.hmm.xz"
        with lzma.open(database, "wt") as handle:
            handle.write("HMMER3/f\nNAME  FAM1\nLENG  10\nHMM  A C\n//\n")
        extracted = extract_profiles_from_hmm(
            str(database), {"FAM1": "RVDB__FAM000001"}, str(tmp_path / "out")
        )
        assert set(extracted) == {"FAM1"}
