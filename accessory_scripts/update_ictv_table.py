#!/usr/bin/python3

# libraries
import argparse
import re
import shutil
import sys
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

# arguments
parser = argparse.ArgumentParser(
    description=(
        "Regenerate eefinder/data/ictv_genome_composition.tsv (family -> genome\n"
        "composition) from the ICTV Virus Metadata Resource (VMR) spreadsheet.\n"
        "\n"
        "The VMR is parsed with the standard library only (an .xlsx is a zip of\n"
        "XML), so this script needs no extra dependency.\n"
        "\n"
        "Reports what changed against the table currently on disk and warns about\n"
        "the three cases that silently degrade EEfinder output: families that\n"
        "disappeared (a rename leaves Molecule_type empty), families the VMR has\n"
        "no composition for, and families with more than one composition (their\n"
        "token order drives the protein standardisation scopes)."
    ),
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "-in",
    "--input",
    help="Local VMR .xlsx to read instead of downloading one.",
    default=None,
)
parser.add_argument(
    "-u",
    "--url",
    help="VMR .xlsx URL to download (default: the newest one listed on\n"
    "https://ictv.global/vmr).",
    default=None,
)
parser.add_argument(
    "-o",
    "--output",
    help="Where to write the TSV (default: the bundled table).",
    default=str(
        Path(__file__).resolve().parent.parent
        / "eefinder"
        / "data"
        / "ictv_genome_composition.tsv"
    ),
)
parser.add_argument(
    "-n",
    "--dry-run",
    help="Report the differences without writing the output file.",
    action="store_true",
)
args = parser.parse_args()

VMR_PAGE = "https://ictv.global/vmr"
ICTV_ROOT = "https://ictv.global"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELS_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKG_RELS_NS = "{http://schemas.openxmlformats.org/package/2006/relationships}"

#: ictv.global answers 403 to the default urllib User-Agent, so send a real one.
HEADERS = {"User-Agent": f"EEfinder-update-ictv-table/{__version__}"}


def fetch(url: str, destination: str = None) -> bytes:
    """GET a URL, either returning the body or streaming it to a file."""
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=300) as response:
        if destination is None:
            return response.read()
        with open(destination, "wb") as handle:
            shutil.copyfileobj(response, handle)
    return b""


def latest_vmr_url() -> str:
    """Return the newest VMR .xlsx URL listed on the ICTV VMR page.

    The files are named ``VMR_MSL<n>.v<n>.<YYYYMMDD>.xlsx``; the one with the
    highest embedded date wins, so a new MSL release is picked up automatically.
    """
    page = fetch(VMR_PAGE).decode("utf-8", "replace")

    links = re.findall(r'href="([^"]*/VMR/[^"]*\.xlsx)"', page)
    if not links:
        sys.exit(f"No VMR .xlsx link found on {VMR_PAGE} -- pass one with --url.")

    def sort_key(link: str) -> str:
        found = re.search(r"(\d{8})", link)
        return found.group(1) if found else ""

    best = max(links, key=sort_key)
    return best if best.startswith("http") else f"{ICTV_ROOT}{best}"


def shared_strings(book: zipfile.ZipFile) -> list:
    """Read the workbook's shared-string table (cells reference it by index)."""
    if "xl/sharedStrings.xml" not in book.namelist():
        return []
    root = ET.fromstring(book.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.iter(f"{NS}t"))
        for item in root.iter(f"{NS}si")
    ]


def data_sheet_path(book: zipfile.ZipFile) -> str:
    """Path of the worksheet holding the VMR records.

    The sheet is named after the release (``VMR MSL41``), so it is matched by
    prefix rather than by name, and resolved through the workbook relationships
    because sheet order does not have to match the ``sheetN.xml`` numbering.
    """
    workbook = ET.fromstring(book.read("xl/workbook.xml"))
    rels = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
    targets = {
        rel.get("Id"): rel.get("Target")
        for rel in rels.iter(f"{PKG_RELS_NS}Relationship")
    }

    for sheet in workbook.iter(f"{NS}sheet"):
        name = sheet.get("name") or ""
        if name.upper().startswith("VMR"):
            target = targets.get(sheet.get(f"{RELS_NS}id"), "")
            return f"xl/{target.lstrip('/')}" if target else ""
    return ""


def row_cells(row, strings: list) -> dict:
    """Map a row's column letters to their text values."""
    values = {}
    for cell in row:
        column = "".join(char for char in (cell.get("r") or "") if char.isalpha())
        if cell.get("t") == "inlineStr":
            node = cell.find(f"{NS}is")
            text = "".join(t.text or "" for t in node.iter(f"{NS}t")) if node else ""
        else:
            node = cell.find(f"{NS}v")
            if node is None:
                text = ""
            elif cell.get("t") == "s":
                text = strings[int(node.text)]
            else:
                text = node.text or ""
        values[column] = text
    return values


def family_compositions(xlsx_path: str) -> dict:
    """Build the ``family -> "comp; comp"`` map from a VMR spreadsheet.

    A VMR ``Genome`` cell may itself already be a ``;``-separated list, so the
    cells are split into atomic tokens before being deduplicated -- joining the
    distinct *cells* instead would emit a composite value next to its own parts
    (e.g. ``ssRNA(-); ssRNA(+/-); ssRNA(-); ssRNA(+/-)`` for Arenaviridae).

    Token order is first-seen in VMR row order, which is what the bundled table
    carries. It is not arbitrary: ``normalization._in_scope`` tests some
    molecule-type scopes with ``startswith``, so the first token is the one that
    decides them.
    """
    with zipfile.ZipFile(xlsx_path) as book:
        strings = shared_strings(book)
        sheet_path = data_sheet_path(book)
        if not sheet_path or sheet_path not in book.namelist():
            sys.exit("Could not find the VMR data sheet in the spreadsheet.")
        rows = list(ET.fromstring(book.read(sheet_path)).iter(f"{NS}row"))

    if not rows:
        sys.exit("The VMR data sheet is empty.")

    header = row_cells(rows[0], strings)
    columns = {text.strip(): letter for letter, text in header.items()}
    missing = [name for name in ("Family", "Genome") if name not in columns]
    if missing:
        sys.exit(
            f"The VMR sheet has no {' and no '.join(missing)} column "
            f"(found: {', '.join(sorted(columns))})."
        )
    family_column, genome_column = columns["Family"], columns["Genome"]

    compositions = {}
    for row in rows[1:]:
        cells = row_cells(row, strings)
        family = cells.get(family_column, "").strip()
        if not family:
            continue
        tokens = compositions.setdefault(family, [])
        for token in (part.strip() for part in cells.get(genome_column, "").split(";")):
            if token and token not in tokens:
                tokens.append(token)

    return {family: "; ".join(tokens) for family, tokens in compositions.items()}


def read_current(path: str) -> dict:
    """Read the table already on disk, so the run can report what changes."""
    table = {}
    if not Path(path).is_file():
        return table
    with open(path) as handle:
        next(handle, None)  # header
        for line in handle:
            family, _, genome = line.rstrip("\n").partition("\t")
            if family:
                table[family] = genome
    return table


def write_table(path: str, compositions: dict) -> None:
    """Write the table in the format eefinder.get_databases expects."""
    with open(path, "w") as handle:
        handle.write("Family\tGenome\n")
        for family in sorted(compositions):
            handle.write(f"{family}\t{compositions[family]}\n")


def report(current: dict, new: dict) -> None:
    """Print the added / changed / removed families and the two warnings."""
    added = sorted(set(new) - set(current))
    removed = sorted(set(current) - set(new))
    changed = sorted(
        family for family in set(current) & set(new) if current[family] != new[family]
    )

    print(f"families on disk : {len(current)}")
    print(f"families in VMR  : {len(new)}")
    print(f"  added   : {len(added)}")
    print(f"  changed : {len(changed)}")
    print(f"  removed : {len(removed)}")

    for family in added:
        print(f"  + {family}\t{new[family]}")
    for family in changed:
        print(f"  ~ {family}\t{current[family]}  ->  {new[family]}")
    for family in removed:
        print(f"  - {family}\t{current[family]}")

    if removed:
        print(
            f"\nWARNING: {len(removed)} family/families are no longer in the VMR. "
            "If they were\nrenamed rather than retired, every record still "
            "carrying the old name gets an\nempty Molecule_type -- the lookup is "
            "exact and failures are silent. Check them\nagainst the current ICTV "
            "taxonomy before shipping the table."
        )

    empty = sorted(family for family, genome in new.items() if not genome)
    if empty:
        print(
            f"\nWARNING: {len(empty)} family/families have no genome composition "
            "in the VMR and\nwill yield an empty Molecule_type: "
            f"{', '.join(empty)}"
        )

    multi = sorted(family for family, genome in new.items() if ";" in genome)
    if multi:
        print(
            f"\nNOTE: {len(multi)} family/families carry more than one composition. "
            "Their first\ntoken decides the startswith-based molecule_type_scope "
            "rules in\nnormalization.py, so review these if protein "
            "standardisation changes:"
        )
        for family in multi:
            print(f"  ! {family}\t{new[family]}")


# main
if args.input:
    spreadsheet = args.input
    print(f"reading {spreadsheet}")
else:
    url = args.url or latest_vmr_url()
    spreadsheet = Path.cwd() / Path(url).name
    print(f"downloading {url}")
    fetch(url, spreadsheet)
    print(f"saved to {spreadsheet}")

compositions = family_compositions(spreadsheet)
current = read_current(args.output)
report(current, compositions)

if args.dry_run:
    print("\n--dry-run: the output file was not written.")
else:
    write_table(args.output, compositions)
    print(f"\nwrote {args.output} ({len(compositions)} families)")
