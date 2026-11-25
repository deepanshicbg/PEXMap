import argparse


parser = argparse.ArgumentParser(
    description="Generate overlapping k-mers from MS/MS peptides (step 1 of PEXMap)"
)
parser.add_argument("--input", required=True, help="text file, one peptide sequence per line")
parser.add_argument("--output", required=True, help="output k-mer file (peptide<TAB>k-mer)")
parser.add_argument("--kmer_length", type=int, default=8, help="k-mer length (default: 8)")
parser.add_argument(
    "--excluded_output",
    default=None,
    help="optional file listing peptides excluded for being shorter than --kmer_length",
)

args = parser.parse_args()
k = args.kmer_length
output_file = args.output
input_file = args.input


def generate_kmers(peptide, k):
    return [peptide[i:i+k] for i in range(len(peptide) - k + 1)]


input_peptides = 0
blank_lines = 0
excluded_short = []
seen = set()
duplicates = 0
processed_peptides = 0
total_kmers = 0

with open(input_file) as f, open(output_file, "w") as out:

    for line in f:

        pep = line.strip()

        if not pep:
            blank_lines += 1
            continue

        input_peptides += 1

        if len(pep) < k:
            # Peptides shorter than k cannot yield a k-mer; they are counted and
            # reported rather than silently dropped.
            excluded_short.append(pep)
            continue

        if pep in seen:
            # identical peptide already processed; its k-mers are not written twice
            duplicates += 1
            continue
        seen.add(pep)

        processed_peptides += 1

        # Removes duplicate k-mers within same peptide
        kmers = set(generate_kmers(pep, k))

        for kmer in kmers:
            out.write(f"{pep}\t{kmer}\n")
            total_kmers += 1


if args.excluded_output:
    with open(args.excluded_output, "w") as ex:
        for pep in excluded_short:
            ex.write(pep + "\n")


print("\nk-mer Generation Summary:")
print("-------------------------")
print(f"Input peptides: {input_peptides}")
print(f"Excluded because shorter than {k} aa: {len(excluded_short)}")
if duplicates:
    print(f"Duplicate peptides (counted once): {duplicates}")
print(f"Processed unique peptides (k-mers generated): {processed_peptides}")
print(f"Total unique k-mers written: {total_kmers}")
if blank_lines:
    print(f"(Blank lines ignored: {blank_lines})")
if excluded_short and not args.excluded_output:
    print("Tip: use --excluded_output to save the list of excluded peptides.")
print("-------------------------")
