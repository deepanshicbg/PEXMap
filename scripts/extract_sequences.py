import os
import pickle
import argparse


parser = argparse.ArgumentParser(
    description="Reconstruct full per-transcript protein sequences from the ENACT per-gene "
                "annotation folder (needed for the exact-substring "
                "verification step, verify_mapping.py). Lightweight, single pass, independent of k."
)
parser.add_argument("--input_folder", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()


def parse_gene_file(path):
    transcript = None
    ordered = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#Transcript:"):
                if transcript is not None:
                    yield transcript, ordered
                transcript = line.split(":")[1].split(",")[0]
                ordered = []
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            exon_id = parts[0].split(",")[0]
            seq = parts[1].strip()
            ordered.append((exon_id, seq))
    if transcript is not None:
        yield transcript, ordered


transcript_sequences = {}
gene_files = os.listdir(args.input_folder)
print("Total gene files:", len(gene_files))

for gi, gene_file in enumerate(gene_files, 1):
    path = os.path.join(args.input_folder, gene_file)
    for transcript, ordered in parse_gene_file(path):
        full_seq = "".join(seq for _, seq in ordered if seq)
        if full_seq:
            transcript_sequences[transcript] = full_seq
    if gi % 2000 == 0:
        print(f"  processed {gi}/{len(gene_files)}")

print("Total transcripts with reconstructed sequence:", len(transcript_sequences))

with open(args.output, "wb") as f:
    pickle.dump(transcript_sequences, f, protocol=pickle.HIGHEST_PROTOCOL)

print("Saved:", args.output)
