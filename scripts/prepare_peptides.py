"""
Prepare an FDR-controlled peptide list for PEXMap from a search-engine PSM file.

PEXMap annotates peptides that have already been identified. This script is the
optional step before it: it reads peptide-spectrum matches (PSMs) from a search
engine, controls the false discovery rate (FDR) at the peptide level, removes
decoys and contaminants, converts modified peptides to plain sequences, and
writes one peptide per line for scripts/generate_kmers.py.

Three ways of controlling FDR, chosen automatically from the input:

  tdc       The file contains target and decoy PSMs with a score. Peptide-level
            q-values are estimated by target-decoy competition: each peptide is
            represented by its best-scoring PSM, and
            FDR = (decoy peptides + 1) / target peptides at each score threshold.
  qvalue    The file already contains q-values from a post-processor such as
            Percolator or mokapot. Peptides are kept if their q-value is at or
            below the threshold.
  filtered  The file was already FDR-filtered by the search pipeline and has
            neither decoys nor q-values (e.g. FragPipe psm.tsv). All target
            peptides are kept and the summary states that FDR was not
            re-estimated.

Supported inputs (--format, detected automatically from the header):
  sage        Sage results.sage.tsv (peptide, label, peptide_q)
  maxquant    MaxQuant msms.txt / evidence.txt (Sequence, Score, Reverse)
  pin         Percolator input (SpecId, Label, ..., Peptide, Proteins); needs --score_col
  percolator  Percolator or mokapot results (peptide, q-value columns)
  fragpipe    FragPipe psm.tsv / peptide.tsv (Peptide, Protein)
  generic     any tab- or comma-separated file; give the columns with the options below

Only the Python standard library is used.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


# --------------------------------------------------------------------------- #
# Command line                                                                #
# --------------------------------------------------------------------------- #

parser = argparse.ArgumentParser(
    description="Create an FDR-controlled peptide list for PEXMap from a PSM file.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument("--input", required=True, help="PSM file from a search engine or post-processor")
parser.add_argument("--output", required=True,
                    help="output peptide list (one plain sequence per line), input for generate_kmers.py")
parser.add_argument("--format", default="auto",
                    choices=["auto", "sage", "maxquant", "pin", "percolator", "fragpipe", "generic"],
                    help="input format (default: detect from the header)")
parser.add_argument("--fdr", type=float, default=0.01,
                    help="peptide-level FDR threshold (default: 0.01 = 1%%)")
parser.add_argument("--peptide_col", help="column with the peptide sequence")
parser.add_argument("--score_col", help="column with the PSM score (tdc mode)")
parser.add_argument("--lower_is_better", action="store_true",
                    help="lower scores are better (e.g. e-values); default: higher is better")
parser.add_argument("--qvalue_col", help="column with q-values (qvalue mode)")
parser.add_argument("--protein_col", help="column with protein accessions")
parser.add_argument("--decoy_col", help="column that marks decoys (e.g. MaxQuant 'Reverse', Percolator 'Label')")
parser.add_argument("--decoy_values", default="+,-1,true,decoy",
                    help="comma-separated values in --decoy_col that mean decoy (default: +,-1,true,decoy)")
parser.add_argument("--decoy_prefix", default="REV__,rev_,DECOY_,decoy_,XXX_",
                    help="protein prefixes that mark decoys, used when no decoy column is given "
                         "(a PSM is a decoy only if all its proteins are decoys)")
parser.add_argument("--contaminant_prefix", default="CON__,contam_,Cont_,CONT_",
                    help="protein prefixes that mark contaminants; peptides matching any contaminant "
                         "protein are removed after FDR estimation (empty string to disable)")
parser.add_argument("--assume_filtered", action="store_true",
                    help="treat the input as already FDR-filtered (no decoys or q-values needed)")
parser.add_argument("--min_length", type=int, default=1,
                    help="drop peptides shorter than this (default 1 = keep all; "
                         "generate_kmers.py reports peptides shorter than k)")
args = parser.parse_args()

if not 0 < args.fdr < 1:
    parser.error("--fdr must be between 0 and 1 (e.g. 0.01 for 1%)")


def split_list(text):
    return [t.strip() for t in text.split(",") if t.strip()]


DECOY_VALUES = {v.lower() for v in split_list(args.decoy_values)}
DECOY_PREFIXES = split_list(args.decoy_prefix)
CONTAMINANT_PREFIXES = split_list(args.contaminant_prefix)
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------- #
# Peptide sequence cleaning                                                   #
# --------------------------------------------------------------------------- #

def clean_peptide(raw):
    """Return the plain amino-acid sequence of a (possibly modified) peptide string.

    Handles flanking residues (K.PEPTIDE.R), bracketed or parenthesised
    modifications (PEPT[+79.97]IDE, PEPT(ph)IDE, ProForma [Phospho]),
    MaxQuant underscores (_PEPTIDE_), terminal markers (n[42]PEPTIDE) and
    lower-case modified residues (PEPTmIDE).
    """
    s = raw.strip()
    m = re.match(r"^[A-Za-z\-]\.(.+)\.[A-Za-z\-]$", s)   # K.PEPTIDE.R
    if m:
        s = m.group(1)
    s = re.sub(r"\[[^\]]*\]|\([^)]*\)|\{[^}]*\}", "", s)  # modifications
    s = s.strip("_-.")
    s = re.sub(r"^n(?=[A-Z])", "", s)                    # n-terminal marker
    s = re.sub(r"(?<=[A-Z])c$", "", s)                   # c-terminal marker
    s = re.sub(r"[^A-Za-z]", "", s)
    return s.upper()


# --------------------------------------------------------------------------- #
# Reading the input                                                           #
# --------------------------------------------------------------------------- #

def read_table(path):
    """Read a tab- or comma-separated file. Returns (header, rows as lists)."""
    with open(path, newline="") as f:
        first = f.readline()
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        f.seek(0)
        reader = csv.reader(f, delimiter=delim)
        header = [h.strip() for h in next(reader)]
        rows = [r for r in reader if r and any(c.strip() for c in r)]
    return header, rows


def find_col(header, wanted, required=True, label=""):
    """Find a column by exact name, then case-insensitively. `wanted` may be a list of candidates."""
    candidates = [wanted] if isinstance(wanted, str) else list(wanted)
    for c in candidates:
        if c in header:
            return header.index(c)
    lower = [h.lower() for h in header]
    for c in candidates:
        if c.lower() in lower:
            return lower.index(c.lower())
    if required:
        sys.exit(f"Error: could not find the {label or 'required'} column "
                 f"({' / '.join(candidates)}). Columns in the file: {', '.join(header)}")
    return None


def detect_format(header):
    h = {x.lower() for x in header}
    if "peptide_q" in h and "label" in h and "peptide" in h:
        return "sage"
    if "reverse" in h and "sequence" in h and "score" in h:
        return "maxquant"
    if "specid" in h and "label" in h and "peptide" in h and not any("q-value" in x for x in h):
        return "pin"
    if any(x in h for x in ("q-value", "mokapot q-value", "q_value", "qvalue")) and "peptide" in h:
        return "percolator"
    if "peptide" in h and "protein" in h and ("hyperscore" in h or "probability" in h
                                               or "peptideprophet probability" in h
                                               or "modified peptide" in h):
        return "fragpipe"
    return "generic"


fmt = args.format
header, rows = read_table(args.input)
if fmt == "auto":
    fmt = detect_format(header)

# Column choices per format (options given on the command line always win).
pep_i = score_i = q_i = prot_i = decoy_i = None
protein_tail = False       # Percolator files list proteins in all columns after 'Proteins'
mode = None

if fmt == "sage":
    # Sage reports a peptide-level q-value (peptide_q) on every PSM row, decoys as label = -1
    pep_i = find_col(header, args.peptide_col or "peptide", label="peptide")
    q_i = find_col(header, args.qvalue_col or "peptide_q", label="q-value")
    prot_i = find_col(header, args.protein_col or "proteins", required=False)
    decoy_i = find_col(header, args.decoy_col or "label", label="decoy")
    mode = "qvalue"
elif fmt == "maxquant":
    pep_i = find_col(header, args.peptide_col or "Sequence", label="peptide")
    score_i = find_col(header, args.score_col or "Score", label="score")
    prot_i = find_col(header, args.protein_col or ["Proteins", "Leading razor protein"], required=False)
    decoy_i = find_col(header, args.decoy_col or "Reverse", label="decoy")
    mode = "tdc"
elif fmt == "pin":
    pep_i = find_col(header, args.peptide_col or "Peptide", label="peptide")
    if not args.score_col:
        sys.exit("Error: Percolator input (pin) files contain several scores. "
                 "Choose one with --score_col (and --lower_is_better if lower is better), "
                 "or run Percolator/mokapot first and give their results file instead.")
    score_i = find_col(header, args.score_col, label="score")
    prot_i = find_col(header, args.protein_col or "Proteins", required=False)
    protein_tail = prot_i is not None and prot_i == len(header) - 1
    decoy_i = find_col(header, args.decoy_col or "Label", label="decoy")
    mode = "tdc"
elif fmt == "percolator":
    pep_i = find_col(header, args.peptide_col or ["peptide", "Peptide"], label="peptide")
    q_i = find_col(header, args.qvalue_col or ["q-value", "mokapot q-value", "q_value", "qvalue"], label="q-value")
    prot_i = find_col(header, args.protein_col or ["proteinIds", "Proteins", "protein"], required=False)
    protein_tail = prot_i is not None and prot_i == len(header) - 1
    decoy_i = find_col(header, args.decoy_col or ["Label", "is_decoy"], required=False)
    mode = "qvalue"
elif fmt == "fragpipe":
    pep_i = find_col(header, args.peptide_col or "Peptide", label="peptide")
    prot_i = find_col(header, args.protein_col or "Protein", required=False)
    mode = "filtered"
else:  # generic
    if not args.peptide_col:
        sys.exit("Error: for a generic file give at least --peptide_col, plus either "
                 "--score_col (with decoys) or --qvalue_col, or use --assume_filtered.")
    pep_i = find_col(header, args.peptide_col, label="peptide")
    score_i = find_col(header, args.score_col, label="score") if args.score_col else None
    q_i = find_col(header, args.qvalue_col, label="q-value") if args.qvalue_col else None
    prot_i = find_col(header, args.protein_col, label="protein") if args.protein_col else None
    decoy_i = find_col(header, args.decoy_col, label="decoy") if args.decoy_col else None
    mode = "qvalue" if q_i is not None else ("tdc" if score_i is not None else "filtered")

if args.assume_filtered:
    mode = "filtered"
elif mode == "filtered" and fmt == "generic":
    sys.exit("Error: no --score_col or --qvalue_col given. Add one, or use --assume_filtered "
             "if the file is already FDR-filtered.")


def proteins_of(row):
    if prot_i is None or prot_i >= len(row):
        return []
    cells = row[prot_i:] if protein_tail else [row[prot_i]]
    out = []
    for c in cells:
        out.extend(p.strip() for p in re.split(r"[;,\t]", c) if p.strip())
    return out


def is_decoy(row, proteins):
    if decoy_i is not None:
        v = row[decoy_i].strip().lower() if decoy_i < len(row) else ""
        if header[decoy_i].lower() == "label":   # Percolator: 1 / -1; mokapot: True / False
            return v in ("-1", "false")
        return v in DECOY_VALUES
    if proteins and DECOY_PREFIXES:
        return all(any(p.startswith(d) for d in DECOY_PREFIXES) for p in proteins)
    return False


def is_contaminant(proteins):
    return any(any(p.startswith(c) for c in CONTAMINANT_PREFIXES) for p in proteins)


def to_float(text, what, line_no):
    try:
        return float(text)
    except (TypeError, ValueError):
        sys.exit(f"Error: line {line_no}: {what} '{text}' is not a number.")


# --------------------------------------------------------------------------- #
# Collect PSMs                                                                #
# --------------------------------------------------------------------------- #

psms = []                 # dicts: seq, score, q, decoy, contaminant, proteins
skipped_empty = 0
ambiguous = set()

rank_i = find_col(header, "rank", required=False) if fmt == "sage" else None
skipped_rank = 0

for n, row in enumerate(rows, start=2):
    if rank_i is not None and rank_i < len(row) and row[rank_i].strip() not in ("", "1"):
        skipped_rank += 1
        continue
    if pep_i >= len(row):
        skipped_empty += 1
        continue
    seq = clean_peptide(row[pep_i])
    if not seq:
        skipped_empty += 1
        continue
    proteins = proteins_of(row)
    psm = {
        "seq": seq,
        "decoy": is_decoy(row, proteins),
        "contaminant": is_contaminant(proteins),
        "proteins": proteins,
        "score": None,
        "q": None,
    }
    if mode == "tdc":
        s = to_float(row[score_i], "score", n)
        psm["score"] = -s if args.lower_is_better else s
    elif mode == "qvalue":
        psm["q"] = to_float(row[q_i], "q-value", n)
    psms.append(psm)

if not psms:
    sys.exit("Error: no PSMs with a peptide sequence were found in the input.")

n_target_psms = sum(not p["decoy"] for p in psms)
n_decoy_psms = len(psms) - n_target_psms

if mode == "tdc" and n_decoy_psms == 0:
    sys.exit("Error: no decoy PSMs were found, so the FDR cannot be estimated. Search with a "
             "target-decoy database, check --decoy_col / --decoy_prefix, give --qvalue_col, "
             "or use --assume_filtered if the file is already FDR-filtered.")


# --------------------------------------------------------------------------- #
# Target-decoy competition                                                    #
# --------------------------------------------------------------------------- #

def tdc_qvalues(items):
    """items: list of (score, is_decoy), higher score = better.
    Returns q-values in the same order. FDR at a threshold = (decoys + 1) / targets,
    tied scores are treated as one threshold, and q-values are made monotonic."""
    order = sorted(range(len(items)), key=lambda i: -items[i][0])
    fdr_at = [0.0] * len(items)
    targets = decoys = 0
    k = 0
    while k < len(order):
        j = k
        while j < len(order) and items[order[j]][0] == items[order[k]][0]:
            if items[order[j]][1]:
                decoys += 1
            else:
                targets += 1
            j += 1
        fdr = (decoys + 1) / targets if targets else float("inf")
        for idx in order[k:j]:
            fdr_at[idx] = fdr
        k = j
    q = [0.0] * len(items)
    running = float("inf")
    for idx in reversed(order):
        running = min(running, fdr_at[idx])
        q[idx] = min(running, 1.0)
    return q


# Group PSMs by plain sequence: one entry per peptide.
by_seq = defaultdict(list)
for p in psms:
    by_seq[p["seq"]].append(p)

peptides = {}  # seq -> dict
for seq, group in by_seq.items():
    targets = [p for p in group if not p["decoy"]]
    decoys = [p for p in group if p["decoy"]]
    entry = {
        "seq": seq,
        "n_psms": len(targets),
        "proteins": sorted({pr for p in targets for pr in p["proteins"]}),
        "contaminant": any(p["contaminant"] for p in targets),
        "decoy": not targets,
    }
    if mode == "tdc":
        best_t = max((p["score"] for p in targets), default=None)
        best_d = max((p["score"] for p in decoys), default=None)
        # competition: a sequence seen as both target and decoy counts as whichever scored higher
        if best_t is not None and (best_d is None or best_t >= best_d):
            entry["score"], entry["decoy"] = best_t, False
        else:
            entry["score"], entry["decoy"] = best_d, True
    elif mode == "qvalue":
        entry["q"] = min((p["q"] for p in targets), default=1.0)
    peptides[seq] = entry

if mode == "tdc":
    plist = list(peptides.values())
    qs = tdc_qvalues([(e["score"], e["decoy"]) for e in plist])
    for e, q in zip(plist, qs):
        e["q"] = q
    psm_q = tdc_qvalues([(p["score"], p["decoy"]) for p in psms])
    psms_passing = sum(1 for p, q in zip(psms, psm_q) if not p["decoy"] and q <= args.fdr)

psm_level_qvalues = (mode == "qvalue" and fmt != "sage"
                     and any(len([p for p in g if not p["decoy"]]) > 1 for g in by_seq.values()))


# --------------------------------------------------------------------------- #
# Select peptides                                                             #
# --------------------------------------------------------------------------- #

status = {}
for seq, e in peptides.items():
    if e["decoy"]:
        status[seq] = "decoy"
    elif mode != "filtered" and e["q"] > args.fdr:
        status[seq] = "above_fdr"
    elif e["contaminant"]:
        status[seq] = "contaminant"
    elif set(seq) - STANDARD_AA:
        status[seq] = "nonstandard_residue"
    elif len(seq) < args.min_length:
        status[seq] = "too_short"
    else:
        status[seq] = "kept"

kept = sorted(s for s, st in status.items() if st == "kept")

with open(args.output, "w") as out:
    for seq in kept:
        out.write(seq + "\n")

table_file = args.output + ".table.tsv"
with open(table_file, "w") as out:
    out.write("Peptide\tStatus\tPeptide_qvalue\tTarget_PSMs\tProteins\n")
    for seq in sorted(peptides):
        e = peptides[seq]
        if e["decoy"]:
            continue
        q = "" if mode == "filtered" else f"{e['q']:.6g}"
        out.write(f"{seq}\t{status[seq]}\t{q}\t{e['n_psms']}\t{';'.join(e['proteins'])}\n")


# --------------------------------------------------------------------------- #
# Summary                                                                     #
# --------------------------------------------------------------------------- #

count = defaultdict(int)
for st in status.values():
    count[st] += 1
target_peptides = sum(1 for e in peptides.values() if not e["decoy"])
short = sum(1 for s in kept if len(s) < 8)

method = {
    "tdc": "target-decoy competition at the peptide level, FDR = (decoys + 1) / targets",
    "qvalue": f"q-values taken from the input column '{header[q_i]}'" if q_i is not None else "q-values",
    "filtered": "none (input treated as already FDR-filtered; FDR not re-estimated)",
}[mode]

lines = [
    "PEXMap peptide preparation summary",
    "----------------------------------",
    f"Input file: {os.path.basename(args.input)}",
    f"Detected format: {fmt}",
    f"FDR control: {method}",
]
if mode != "filtered":
    lines.append(f"Peptide-level FDR threshold: {args.fdr:g}")
lines += [
    "",
    f"PSMs read: {len(psms)} ({n_target_psms} target, {n_decoy_psms} decoy)",
]
if skipped_rank:
    lines.append(f"Lower-ranked PSMs ignored (rank > 1): {skipped_rank}")
if skipped_empty:
    lines.append(f"Rows without a peptide sequence (skipped): {skipped_empty}")
if mode == "tdc":
    lines.append(f"Target PSMs at {args.fdr:g} PSM-level FDR (for reference): {psms_passing}")
lines += [
    f"Distinct target peptides (plain sequence): {target_peptides}",
]
if mode != "filtered":
    lines.append(f"  Above the FDR threshold: {count['above_fdr']}")
lines += [
    f"  Contaminants removed: {count['contaminant']}",
    f"  Non-standard residues removed (B, J, O, U, X, Z): {count['nonstandard_residue']}",
]
if args.min_length > 1:
    lines.append(f"  Shorter than {args.min_length} aa removed: {count['too_short']}")
lines += [
    f"  Kept (written to {os.path.basename(args.output)}): {len(kept)}",
    f"    of which shorter than 8 aa (generate_kmers.py will report these): {short}",
]
if mode == "tdc" and n_decoy_psms < 20:
    lines.append("WARNING: very few decoys; FDR estimates are imprecise.")
if psm_level_qvalues:
    lines.append("NOTE: the input lists several PSMs per peptide, so its q-values are PSM-level. "
                 "Each peptide was given its best PSM q-value. For true peptide-level FDR, use the "
                 "post-processor's peptide-level output (e.g. Percolator --results-peptides or "
                 "mokapot.peptides.txt).")
if mode == "filtered":
    lines.append("NOTE: FDR was not re-estimated. State the FDR used by the search pipeline when "
                 "reporting results.")
lines.append(f"Per-peptide details: {os.path.basename(table_file)}")

summary = "\n".join(lines)
print(summary)
with open(args.output + ".summary.txt", "w") as out:
    out.write(summary + "\n")
