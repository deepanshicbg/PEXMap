import os
import pickle
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description="Build peptide annotation database")

parser.add_argument("--input_folder", required=True)
parser.add_argument("--kmer", type=int, default=8)
parser.add_argument("--organism", required=True)
parser.add_argument("--output", default="peptide_dataset.pkl")

args = parser.parse_args()

input_folder = args.input_folder
k = args.kmer
organism = args.organism
output_file = args.output


########################################
# helper
########################################

def kmers(seq, k):
    L = len(seq)
    for i in range(L - k + 1):
        yield seq[i:i+k]


########################################
# main storage
########################################

# peptide -> (gene,feature_type,feature_id) -> transcripts
peptide_db = defaultdict(lambda: defaultdict(set))


########################################
# parse gene files
########################################

gene_files = os.listdir(input_folder)

print("Total gene files:", len(gene_files))

for gene_file in gene_files:

    gene_id = gene_file
    path = os.path.join(input_folder, gene_file)

    transcript = None
    prev_seq = None
    prev_exon = None

    with open(path) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue


            ################################
            # transcript line
            ################################

            if line.startswith("#Transcript:"):

                transcript = line.split(":")[1].split(",")[0]

                prev_seq = None
                prev_exon = None

                continue


            ################################
            # exon line
            ################################

            parts = line.split("\t")

            if len(parts) != 2:
                continue

            exon_id = parts[0].split(",")[0]
            seq = parts[1].strip()

            if not seq:
                continue


            ################################
            # exon kmers
            ################################

            key_exon = (gene_id, "exon", exon_id)

            for pep in kmers(seq, k):
                peptide_db[pep][key_exon].add(transcript)


            ################################
            # junction kmers
            ################################

            if prev_seq is not None:

                junction_id = prev_exon + "," + exon_id

                # boundary sequence
                junction_seq = prev_seq[-(k-1):] + seq[:(k-1)]

                key_junc = (gene_id, "junction", junction_id)

                for pep in kmers(junction_seq, k):
                    peptide_db[pep][key_junc].add(transcript)


            prev_seq = seq
            prev_exon = exon_id


########################################
# convert structure
########################################

print("Converting database structure...")

final_db = {}

for pep, groups in peptide_db.items():

    entries = []

    for (gene, ftype, feature), transcripts in groups.items():

        entries.append({
            "gene_id": gene,
            "feature_type": ftype,
            "feature_id": feature,
            "transcript_ids": sorted(transcripts),
            "organism": organism
        })

    final_db[pep] = entries


########################################
# save pickle
########################################

print("Saving database...")

with open(output_file, "wb") as f:
    pickle.dump(final_db, f, protocol=pickle.HIGHEST_PROTOCOL)


print("Organism:", organism)
print("k-mer length:", k)
print("Total peptides:", len(final_db))
print("Database saved to:", output_file)