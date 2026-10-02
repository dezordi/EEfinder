# EEfinder bundled data

## `viral_proteins.tsv` — viral protein-name standardization map

A starting reference to normalise the free-text `Protein` names found in the
viral database metadata (built by `eefinder get-databases virus`) into a small
set of canonical names.

This is intentionally **not** exhaustive — it covers the most frequent /
biologically relevant proteins and is meant to be extended.

The map is loaded by `eefinder/normalization.py` (`standardize_protein`, target
`virus`).

### Columns

```
suggested_name <TAB> current_name <TAB> match_type <TAB> molecule_type_scope <TAB> notes
```

- **suggested_name** — the canonical name emitted when the rule matches.
- **current_name** — the already-normalised key to match against (lower-case,
  separators collapsed; see below).
- **match_type** — `exact` or `contains` (see step 2).
- **molecule_type_scope** — restricts the rule to a genome group (see shorthand).
- **notes** — free-text documentation for the rule (optional).

### How the table is applied

1. **Normalise** the raw `Protein` string before matching:
   - remove leaked NCBI `[key=value]` tags (e.g. `[organism=...]`);
   - strip a leading `CDS:` / `ORF:` naming directive;
   - normalise molecular-weight tokens (`33 kDa`, `33-kDa`, `33K-like protein`
     → `33 kDa protein`);
   - fix common misspellings (e.g. `membran` → `membrane`);
   - lowercase;
   - replace hyphens/underscores AND compound separators `/\();,` with a single
     space (so `CP/RdRp fusion`, `...; RdRp` become matchable tokens);
   - collapse runs of whitespace to one space; strip ends;
   - strip leading qualifiers: `putative`, `predicted`, `probable`, `possible`,
     `presumed`, `presumptive` (also from the emitted name);
   - strip trailing `, partial` / ` precursor`.

   e.g. `Putative RNA-dependent RNA Polymerase` → `rna dependent rna polymerase`.

2. **Match** against `current_name` (already normalised in this file):
   - `match_type=exact` — the normalised string equals `current_name`;
   - `match_type=contains` — `current_name` occurs as a whole-word substring.

   Prefer `exact` for short/ambiguous names (e.g. `l`, `n`, `cp`) so that, for
   example, `rna polymerase sigma factor` (a phage enzyme) does **not** become
   RdRp. On the **output** name, quotes and the characters `:,/\?!` are removed,
   the first letter is capitalised, and a name that is only a directive becomes
   `Unknown`.

3. **Scope** — only apply the row when the record's `Molecule_type` is within
   `molecule_type_scope`; rows scoped to a group must not rewrite proteins from
   other groups.

### `molecule_type_scope` shorthand

| Token    | Meaning                                                        |
|----------|----------------------------------------------------------------|
| `RNA`    | any RNA genome, RT excluded [ssRNA(+), ssRNA(-), ssRNA(+/-), dsRNA] |
| `+ssRNA` | ssRNA(+)                                                       |
| `-ssRNA` | ssRNA(-)                                                       |
| `dsRNA`  | dsRNA                                                          |
| `RT`     | reverse-transcribing [ssRNA-RT, dsDNA-RT]                     |
| `dsDNA`  | dsDNA                                                          |
| `ssDNA`  | ssDNA                                                          |
| `any`    | no restriction                                                |

### Row groupings (order in the file)

- **RNA viruses** — RdRp (the worked example: the viral replicase / L protein;
  `contains` catches compound names like `P2-RdRp`, `RdRp protein`, `CP/RdRp
  fusion`, `...; RdRp`, `RdRp-like`), Nucleocapsid Protein, Phosphoprotein (P),
  Matrix protein (M), Glycoprotein (envelope glycoprotein; covers RNA **and**
  DNA viruses), Fusion (F) — kept **separate** from the generic glycoprotein,
  Capsid / coat protein (CP), and other (+)RNA proteins.
- **RT viruses** — Reverse Transcriptase, Gag, Env, integrase.
- **dsDNA phages / large DNA viruses** — dominant in a whole-virome download;
  structural naming is a large separate domain, so this is a representative
  starter set to be expanded. Note: in dsDNA phages the major capsid protein
  (MCP) is a different context from the (+)RNA/dsRNA capsid protein.
- **Unresolved** — leading `CDS:`/`ORF:` naming directives are stripped in code
  before matching; a name that is only such a directive (e.g. `CDS:`, `ORF`)
  becomes `Unknown`. Names beginning with `hypothetical` (any spelling) are also
  flagged `Unknown` and dropped from the FASTA and CSV.

## `ictv_genome_composition.tsv` — ICTV genome-composition table

Maps a virus **family** to its genome composition (`dsDNA`, `ssRNA(+)`,
`ssRNA-RT`, …). This is the only source of the `Molecule_type` column in the
metadata CSV: the NCBI datasets taxonomy report does not carry it.

Loaded by `eefinder/get_databases.py` (`_load_genome_composition`, exposed as
`molecule_type_for_family`).

### Format

```
Family <TAB> Genome
```

- One header row, skipped on load. Only the first two tab-separated fields are
  read, so trailing columns are ignored.
- The lookup is **exact and case-sensitive**. A family that is absent — or
  present under a different spelling — yields an **empty** `Molecule_type`, with
  no warning. Bacterial families are absent by design.
- A family with more than one composition carries them `"; "`-joined, e.g.
  `Pleolipoviridae → ssDNA; ssDNA(+/-); dsDNA`.

> **Watch the token order.** For multi-value families `_in_scope` in
> `normalization.py` tests the DNA/dsRNA molecule-type scopes with `startswith`,
> so only the *first* token decides them: a family stored as `dsDNA; ssDNA` does
> **not** satisfy an `ssDNA`-scoped protein rule, while `ssDNA; dsDNA` does. Keep
> the order the VMR delivers rather than sorting the tokens.

### Source

The ICTV **Virus Metadata Resource (VMR)**: <https://ictv.global/vmr>, which
publishes one `.xlsx` per Master Species List release
(`VMR_MSL<n>.v<n>.<YYYYMMDD>.xlsx`). Its data sheet has `Family` and `Genome`
columns that correspond one-to-one to this file's two columns.

Note that <https://ictv.global/virus-properties> — the ICTV Report chapter that
describes these categories — is **an HTML page with no downloadable table**. It
explains the vocabulary; the VMR is what you regenerate from.

### Regenerating it

Use [`accessory_scripts/update_ictv_table.py`](../../accessory_scripts/update_ictv_table.py).
It downloads the newest VMR, rebuilds the table, and reports what changed. It
parses the spreadsheet with the standard library only, so it needs no extra
dependency:

```bash
# see what a new MSL release would change, without writing anything
python accessory_scripts/update_ictv_table.py --dry-run

# write the regenerated table over the bundled one
python accessory_scripts/update_ictv_table.py

# work from a VMR file you already downloaded
python accessory_scripts/update_ictv_table.py -in VMR_MSL41.v1.20260729.xlsx
```

The grouping rule the script implements, if you ever need to redo it by hand:
collect the distinct genome compositions per family, in the order the VMR lists
them, and join them with `"; "`.

> **The one trap.** A VMR `Genome` cell **may itself already be a
> `;`-separated list**, so the cells must be split into atomic tokens before
> deduplicating. Joining the distinct *cell values* instead emits a composite
> value next to its own parts — `Arenaviridae` comes out as
> `ssRNA(-); ssRNA(+/-); ssRNA(-); ssRNA(+/-)`.

### What to check after regenerating

The script prints each of these; all three are silent-failure modes:

1. **Removed families.** A family that vanished was more likely *renamed* than
   retired — the ICTV renames taxa every release. Every record still carrying
   the old name would get an empty `Molecule_type`.
2. **Families with no composition** in the VMR, which also produce an empty
   `Molecule_type`.
3. **Multi-value families**, whose first token drives the scope rules above.
   Review these when protein standardisation output changes unexpectedly.

Then run `pytest -m "not integration" tests/test_get_databases.py` — it asserts
`Molecule_type` values taken from this table (e.g. `Flaviviridae → ssRNA(+)`).
