import argparse
import pickle
import time
from collections import defaultdict


parser = argparse.ArgumentParser(
    description="Exact-substring-search "
                "baseline. For each MS/MS peptide, scan every reference transcript sequence "
                "directly for a literal, contiguous occurrence -- no k-mer decomposition, no "
                "hashing, no voting. This is the 'obviously correct' method PEXMap's octamer "
                "approach is benchmarked against, to make the tradeoff (runtime/memory vs. "
                "guaranteed exactness) explicit rather than asserted."
)
parser.add_argument("--peptides", required=True, help="one peptide sequence per line")
parser.add_argument("--sequences", required=True, help="transcript_sequences.pkl (transcript_id -> full protein sequence)")
parser.add_argument("--database", required=True, help="any built octamerDB_k*.pkl, used only to derive transcript->gene mapping")
parser.add_argument("--output", required=True)
args = parser.parse_args()

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
del db  # free the multi-GB k-mer index; only the derived mapping is needed below

peptides = []
with open(args.peptides) as f:
    for line in f:
        pep = line.strip()
        if pep:
            peptides.append(pep)

print(f"Total query peptides: {len(peptides):,}")
print("Running exact-substring scan against every reference transcript (no indexing)...")

items = list(sequences.items())

results = []
unmatched = 0
unique_gene = 0
multi_gene = 0

t0 = time.time()
progress_every = max(1, len(peptides) // 500)

for i, pep in enumerate(peptides, 1):
    matched_transcripts = [tid for tid, seq in items if pep in seq]

    if not matched_transcripts:
        unmatched += 1
        results.append((pep, "", "", 0))
    else:
        genes = set()
        for tid in matched_transcripts:
            genes |= transcript_to_gene.get(tid, set())

        if len(genes) == 1:
            unique_gene += 1
        else:
            multi_gene += 1

        results.append((pep, ";".join(sorted(genes)), ";".join(sorted(matched_transcripts)), len(matched_transcripts)))

    if i % progress_every == 0 or i == len(peptides):
        rate = i / (time.time() - t0)
        eta_s = (len(peptides) - i) / rate if rate > 0 else 0
        print(f"[progress] {i:,}/{len(peptides):,} ({i/len(peptides)*100:.1f}%) "
              f"elapsed={time.time()-t0:.0f}s rate={rate:.1f} pep/s ETA={eta_s/3600:.2f}h", flush=True)

    if i % 50000 == 0:
        with open(args.output, "w") as out:
            out.write("Peptide\tGene_ids\tMatched_transcripts\tN_matched_transcripts\n")
            for r in results:
                out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
        print(f"[checkpoint] wrote {len(results):,} results to {args.output}", flush=True)

elapsed = time.time() - t0

with open(args.output, "w") as out:
    out.write("Peptide\tGene_ids\tMatched_transcripts\tN_matched_transcripts\n")
    for pep, genes, transcripts, n in results:
        out.write(f"{pep}\t{genes}\t{transcripts}\t{n}\n")

print("\n=== Exact-substring baseline results ===")
print(f"Total peptides: {len(peptides):,}")
print(f"Matched >=1 transcript (found somewhere): {len(peptides) - unmatched:,} ({(len(peptides)-unmatched)/len(peptides)*100:.1f}%)")
print(f"  ...of which unique to one gene: {unique_gene:,}")
print(f"  ...of which multi-gene: {multi_gene:,}")
print(f"Unmatched (not found anywhere as exact substring): {unmatched:,}")
print(f"\nTotal wall-clock time: {elapsed:.2f} s ({elapsed/len(peptides)*1000:.2f} ms/peptide)")
print(f"Comparisons performed: {len(peptides):,} peptides x {len(items):,} transcripts = {len(peptides)*len(items):,} substring checks")
