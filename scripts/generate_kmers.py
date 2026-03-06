import sys

input_file = sys.argv[1]
output_file = sys.argv[2]


def generate_kmers(peptide, k=8):
    return [peptide[i:i+k] for i in range(len(peptide) - k + 1)]


total_kmers = 0

with open(input_file) as f, open(output_file, "w") as out:

    for line in f:

        pep = line.strip()

        if not pep:
            continue

        if len(pep) < 8:
            continue

        # remove duplicate kmers within same peptide
        kmers = set(generate_kmers(pep))

        for kmer in kmers:
            out.write(f"{pep}\t{kmer}\n")
            total_kmers += 1


print("Total kmers:", total_kmers)