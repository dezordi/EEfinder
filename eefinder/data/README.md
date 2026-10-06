# EEfinder bundled data

## `protein_rules.yaml` — protein-name normalisation rules

Normalises the free-text `Protein` names of a RefSeq download (`eefinder
get-databases`) into a small set of canonical names. Loaded by
`eefinder/normalization.py` (`standardize_protein`).

This is intentionally **not** exhaustive — it covers the most frequent /
biologically relevant proteins and is meant to be extended. Everything is data:
adding a term never requires a code change.

### The three kinds of adjustment

The file has one section per kind, applied in this order:

| Section | What it does | Shape |
|---------|--------------|-------|
| `typos` | whole-word spelling corrections, everywhere | correct spelling -> list of misspellings |
| `rewrites` | ordered regex edits to the *shape* of a name | list of `pattern` / `replacement` / `notes` |
| `proteins` | synonym -> canonical name | nested under the scopes the rule applies in |

A **typo** says a word is misspelled; a **rewrite** says a name is
*shaped* wrong (a redundant suffix, a reversed designation, a fused word); a
**protein** rule says one name *means* another. Reach for the narrowest one.

### `typos`

```yaml
typos:
  glycoprotein:
    - glyocprotein   # transposition
    - gylcoprotein   # transposition
```

Keyed by the **correct** spelling, so every misspelling of one word sits
together. Matching is whole-word and case-insensitive, and the correction is
applied both to the match key (so it reaches a `proteins` rule) and to the
emitted name (so typo variants of an *unmapped* protein still converge). Because
it is whole-word, a typo that is a prefix of the correct spelling (`membran` ->
`membrane`) does not corrupt the already-correct word.

### `rewrites`

```yaml
rewrites:
  - pattern: '\b(orf\s*\d+[a-z]?)\s+protein\b'
    replacement: \1
    notes: an ORF designation already says it is a protein
```

Python regexes, applied **in file order** with `re.IGNORECASE`; `replacement`
takes `\1`-style group references and may be empty (a deletion). They run after
the typo corrections, so a pattern may rely on the corrected spelling. Since
`re.IGNORECASE` is on, use an un-ignored group — `(?-i:[A-Z])` — when a rule has
to be case-sensitive.

### `proteins`

```yaml
proteins:
  Chuviridae:           # taxon scope: a taxon name at any rank, or "any"
    "-ssRNA":           # molecule_type scope (table below), or "any"
      Glycoprotein:     # the canonical name emitted
        exact:          # match type
          - s protein
        contains:
          - g protein
        regex:
          - ^gp\d+$
```

The nesting **is** the scope, so it is stated once for a whole group of terms
instead of being repeated per term.

Keys are matched against the name after normalisation: lower-cased, hyphens,
underscores and the compound separators `/\();,` replaced by single spaces,
whitespace collapsed, leading `putative`/`predicted`/… and trailing `, partial`
/ ` precursor` stripped. So `Putative RNA-dependent RNA Polymerase` is matched as
`rna dependent rna polymerase` — write the keys in that form.

| Match type | Meaning |
|------------|---------|
| `exact` | the normalised name equals the key |
| `contains` | the key occurs in it as a whole-word substring |
| `regex` | a Python regex, searched against it (case already folded) |

Prefer `exact` for short or ambiguous keys (`l`, `n`, `cp`, `g`) so that, say,
`rna polymerase sigma factor` (a phage enzyme) does not become RdRp.

#### Precedence

Explicit, not a consequence of the layout — regrouping the file cannot silently
change results:

1. rules restricted to a **taxon** beat rules scoped to `any`, which is how a
   family that names a protein against its genome group's convention wins
   (`coat protein` is the capsid in (+)RNA viruses but the nucleocapsid in
   *Chuviridae*);
2. `exact` beats `contains`/`regex`;
3. the **earliest match** in the name, since a product string leads with the
   protein's name (`nucleocapsid phosphoprotein` is a nucleocapsid protein);
4. file order.

#### `molecule_type` scope shorthand

The record's `Molecule_type` comes from `ictv_genome_composition.tsv` via its
family, so a rule scoped to a genome group never rewrites proteins of another.

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

Several tokens may be `;`-joined (`+ssRNA; dsRNA`). Quote `"-ssRNA"` in YAML so
it cannot be read as a list item.

#### `taxon` scope

A taxon name at **any** rank — family, order, genus — matched case-insensitively
against the record's lineage, or `any` for no restriction. Several names may be
`;`-joined. Scoping by family is the usual case; an order (`Jingchuvirales`)
covers its families in one rule.

> An NCBI placeholder node (`unclassified Chuviridae`) is not a taxon: the
> lineage records the real family above it, and `molecule_type_for_family`
> falls back to the stripped name. Write rules against the real name.

### What is still in code

The generic pipeline around the tables: leaked NCBI `[key=value]` tags are
removed, a leading `CDS:`/`ORF:` directive is stripped, molecular-weight tokens
are normalised (`33 kDa`, `33-kDa`, `33K-like protein` -> `33 kDa protein`),
`NSxx protein` is reduced to `NSxx`, quotes and `:,/\?!` are removed from the
output, the first letter is capitalised, a name that is only a directive becomes
`Unknown`, and a name beginning with `hypothetical` (any spelling, via `typos`)
is flagged `Unknown` and dropped from the FASTA and CSV.

### Row groupings (order in the file)

- **`any` taxon** — the group-wide vocabulary: RdRp (the worked example: the
  viral replicase / L protein; `contains` catches compound names like `P2-RdRp`,
  `RdRp protein`, `CP/RdRp fusion`), Nucleocapsid Protein, Phosphoprotein (P),
  Matrix protein (M), Glycoprotein, Fusion (F) — kept **separate** from the
  generic glycoprotein — Capsid / coat protein (CP), the other (+)RNA proteins,
  the RT set (Reverse Transcriptase, Gag, Env, integrase), and a starter set for
  dsDNA phages and large DNA viruses, where structural naming is a large
  separate domain. Note that in dsDNA phages the major capsid protein (MCP) is a
  different context from the (+)RNA/dsRNA capsid protein.
- **Per-family sections** — only where a family disagrees with the above.

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
