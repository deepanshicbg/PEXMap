import argparse


parser = argparse.ArgumentParser(description="Build peptide annotation database")
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--kmer_length", type=int, default=8)

args = parser.parse_args()
k = args.kmer
output_file = args.output
input_file = args.input


def generate_kmers(peptide, k):
    return [peptide[i:i+k] for i in range(len(peptide) - k + 1)]


total_kmers = 0

with open(input_file) as f, open(output_file, "w") as out:

    for line in f:

        pep = line.strip()

        if not pep:
            continue

        if len(pep) < 8:
            continue

        # Removes duplicate k-mers within same peptide
        kmers = set(generate_kmers(pep))

        for kmer in kmers:
            out.write(f"{pep}\t{kmer}\n")
            total_kmers += 1


print("Total unique k-mers:", total_kmers)