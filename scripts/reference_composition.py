#!/usr/bin/env python3
"""
How much of the octamerDB reference is curated (NP_/NM_)
vs. computationally predicted (XP_/XM_) RefSeq accessions, and how many
example peptides receive an assignment involving at least one predicted model.

Usage:
    python3 scripts/reference_composition.py \
        --database data/octamerDB_latest.pkl \
        --kmers example/example_kmers.txt \
        --output results/reference_composition.tsv

Notes:
- This only reads octamerDB.pkl and (optionally) an existing k-mer file; it
  does not need PeptideAtlas or the raw ENACT annotation files.
- Loading an ~1GB pickle can transiently use several GB of RAM (Python's
  per-object overhead is much larger than the on-disk pickle size for a
  dict of many small dicts/lists).
"""
import argparse
import pickle
import re
import sys
from collections import defaultdict, Counter

ACCESSION_PREFIX_RE = re.compile(r"^([A-Z]{1,3})_")

# RefSeq prefix -> curated vs predicted, per NCBI's own scheme
# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
CURATED_PREFIXES = {"NP", "NM", "NR", "NG", "NC"}
PREDICTED_PREFIXES = {"XP", "XM", "XR"}


def classify(accession):
    m = ACCESSION_PREFIX_RE.match(accession)
    if not m:
        return "other/unrecognized"
    prefix = m.group(1)
    if prefix in CURATED_PREFIXES:
        return "curated"
    if prefix in PREDICTED_PREFIXES:
        return "predicted"
    return "other/unrecognized"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--database", required=True, help="Path to octamerDB(.pkl) (merged exon+junction database)")
    ap.add_argument("--kmers", required=False, help="Optional: an existing k-mer file (peptide\\tkmer per line) to also report the per-peptide XP_-involvement statistic")
    ap.add_argument("--output", required=True, help="Output TSV summary path")
    args = ap.parse_args()

    print(f"Loading database: {args.database}", file=sys.stderr)
    with open(args.database, "rb") as f:
        db = pickle.load(f)
    print(f"Total k-mer keys in database: {len(db):,}", file=sys.stderr)

    # --------------------------------------------------------------
    # Part A: distinct accession composition of the whole reference
    # --------------------------------------------------------------
    all_accessions = set()
    for kmer_entries in db.values():
        for entry in kmer_entries:
            for tr in entry.get("transcript_ids", []) or []:
                all_accessions.add(tr)

    class_counts = Counter(classify(a) for a in all_accessions)
    total_acc = len(all_accessions)

    print("\n=== Part A: distinct RefSeq accession composition ===")
    for cls in ("curated", "predicted", "other/unrecognized"):
        n = class_counts.get(cls, 0)
        pct = (n / total_acc * 100) if total_acc else 0.0
        print(f"{cls:>20}: {n:,} ({pct:.1f}%)")
    print(f"{'TOTAL':>20}: {total_acc:,}")

    lines = []
    lines.append("metric\tvalue")
    lines.append(f"total_distinct_accessions\t{total_acc}")
    for cls in ("curated", "predicted", "other/unrecognized"):
        n = class_counts.get(cls, 0)
        pct = (n / total_acc * 100) if total_acc else 0.0
        lines.append(f"accessions_{cls.replace('/', '_')}\t{n}")
        lines.append(f"accessions_{cls.replace('/', '_')}_pct\t{pct:.2f}")

    # --------------------------------------------------------------
    # Part B (optional): for a peptide k-mer file, how many peptides have
    # >=1 matched database entry whose transcript_ids include a predicted
    # (XP_/XM_/XR_) accession -- i.e. how many "unique" calls could be
    # broken specifically by admitting predicted models to the reference.
    # --------------------------------------------------------------
    if args.kmers:
        print(f"\nLoading k-mers: {args.kmers}", file=sys.stderr)
        pep_kmers = defaultdict(set)
        with open(args.kmers) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue
                pep, kmer = parts
                pep_kmers[pep].add(kmer)

        total_peptides = len(pep_kmers)
        involves_predicted = 0
        involves_only_predicted = 0
        involves_only_curated = 0
        unmatched = 0

        for pep, kmers in pep_kmers.items():
            saw_curated = False
            saw_predicted = False
            matched_any = False
            for kmer in kmers:
                entries = db.get(kmer)
                if not entries:
                    continue
                matched_any = True
                for entry in entries:
                    for tr in entry.get("transcript_ids", []) or []:
                        c = classify(tr)
                        if c == "curated":
                            saw_curated = True
                        elif c == "predicted":
                            saw_predicted = True
            if not matched_any:
                unmatched += 1
                continue
            if saw_predicted:
                involves_predicted += 1
            if saw_predicted and not saw_curated:
                involves_only_predicted += 1
            if saw_curated and not saw_predicted:
                involves_only_curated += 1

        print("\n=== Part B: per-peptide predicted-model involvement ===")
        print(f"Total peptides in k-mer file: {total_peptides:,}")
        print(f"Unmatched (no database hit at all): {unmatched:,}")
        print(f"Peptides with >=1 matched entry involving a PREDICTED (XP_/XM_/XR_) accession: {involves_predicted:,} ({involves_predicted/total_peptides*100:.1f}%)")
        print(f"  ...of which involve ONLY predicted accessions (no curated competitor at all): {involves_only_predicted:,}")
        print(f"Peptides with matched entries involving ONLY curated accessions: {involves_only_curated:,}")

        lines.append(f"kmers_total_peptides\t{total_peptides}")
        lines.append(f"kmers_unmatched\t{unmatched}")
        lines.append(f"kmers_involves_predicted\t{involves_predicted}")
        lines.append(f"kmers_involves_predicted_pct\t{involves_predicted/total_peptides*100:.2f}")
        lines.append(f"kmers_involves_only_predicted\t{involves_only_predicted}")
        lines.append(f"kmers_involves_only_curated\t{involves_only_curated}")

    with open(args.output, "w") as out:
        out.write("\n".join(lines) + "\n")
    print(f"\nWrote summary to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
