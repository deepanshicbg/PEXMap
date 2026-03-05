import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

def generate_kmers(peptide, k=8):
    return [peptide[i:i+k] for i in range(len(peptide)-k+1)]

kmers = set()

with open(input_file) as f:

    for line in f:

        pep = line.strip()

        if len(pep) < 8:
            continue

        for kmer in generate_kmers(pep):
            kmers.add(kmer)

with open(output_file,"w") as out:

    for k in sorted(kmers):
        out.write(k+"\n")

print("Total kmers:", len(kmers))