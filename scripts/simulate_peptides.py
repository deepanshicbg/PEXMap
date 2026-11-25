import argparse
import pickle
import random
import re
from collections import defaultdict


parser = argparse.ArgumentParser(
    description="Simulation benchmark with known ground truth. "
                "Performs in-silico tryptic digestion (cut after K/R, not before P) on real "
                "reconstructed transcript sequences, so every synthetic peptide's true "
                "originating gene/transcript(s) are known by construction, then samples a "
                "peptide set matched to the real MS/MS peptide length distribution."
)
parser.add_argument("--sequences", required=True, help="transcript_sequences.pkl")
parser.add_argument("--database", required=True, help="any octamerDB_k*.pkl, used only to derive transcript->gene mapping")
parser.add_argument("--n", type=int, default=10000, help="number of synthetic peptides to sample")
parser.add_argument("--min_length", type=int, default=8)
parser.add_argument("--max_length", type=int, default=52)
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("--output_peptides", required=True)
parser.add_argument("--output_truth", required=True)
args = parser.parse_args()

random.seed(args.seed)

print("Loading transcript sequences...")
with open(args.sequences, "rb") as f:
    sequences = pickle.load(f)
print(f"Total reference transcripts: {len(sequences):,}")

print("Deriving transcript -> gene mapping from database...")
with open(args.database, "rb") as f:
    db = pickle.load(f)

transcript_to_gene = defaultdict(set)
for kmer_entries in db.values():
    for entry in kmer_entries:
        gene = entry.get("gene_id")
        if gene is None:
            continue
        for tr in entry.get("transcript_ids", []) or []:
            transcript_to_gene[tr].add(gene)
del db


def tryptic_digest(seq):
    """Cut after K or R, except when followed by P (standard trypsin specificity rule)."""
    peptides = []
    start = 0
    for i, aa in enumerate(seq):
        if aa in "KR" and (i + 1 >= len(seq) or seq[i + 1] != "P"):
            peptides.append(seq[start:i + 1])
            start = i + 1
    if start < len(seq):
        peptides.append(seq[start:])
    return peptides


print("Performing in-silico tryptic digestion of all reference transcripts...")
# peptide -> set of transcript_ids it truly originates from (ground truth)
peptide_origin = defaultdict(set)

for tid, seq in sequences.items():
    for pep in tryptic_digest(seq):
        if args.min_length <= len(pep) <= args.max_length:
            peptide_origin[pep].add(tid)

print(f"Total distinct in-silico tryptic peptides in length range [{args.min_length},{args.max_length}]: {len(peptide_origin):,}")

candidates = list(peptide_origin.keys())
random.shuffle(candidates)
sampled = candidates[:args.n]

print(f"Sampled {len(sampled):,} synthetic peptides.")

with open(args.output_peptides, "w") as f:
    for pep in sampled:
        f.write(pep + "\n")

with open(args.output_truth, "w") as f:
    f.write("Peptide\tTrue_transcripts\tTrue_genes\n")
    for pep in sampled:
        true_transcripts = sorted(peptide_origin[pep])
        true_genes = set()
        for tid in true_transcripts:
            true_genes |= transcript_to_gene.get(tid, set())
        f.write(f"{pep}\t{';'.join(true_transcripts)}\t{';'.join(sorted(true_genes))}\n")

lengths = [len(p) for p in sampled]
print(f"Sampled peptide length: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.1f}")
print(f"Saved peptides to {args.output_peptides}")
print(f"Saved ground truth to {args.output_truth}")
