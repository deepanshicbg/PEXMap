import argparse
import pickle


parser = argparse.ArgumentParser(
    description="Post-hoc exact-substring verification. "
                "Octamer voting confirms that every 8-mer of a peptide occurs SOMEWHERE in the "
                "reference, but not that the peptide occurs CONTIGUOUSLY in the assigned isoform "
                "(e.g. all 4 octamers of MFGHKLTWDYR are present in a "
                "scrambled decoy sequence that does not contain the peptide itself). This script "
                "takes an annotate_peptides.py output and confirms each assignment against "
                "the assigned transcript(s)' real, full-length protein sequence."
)
parser.add_argument("--annotations", required=True, help="output TSV from annotate_peptides.py")
parser.add_argument("--sequences", required=True, help="transcript sequence pickle from extract_sequences.py")
parser.add_argument("--output", required=True, help="annotated TSV with an added Exact_match column")
args = parser.parse_args()

with open(args.sequences, "rb") as f:
    transcript_sequences = pickle.load(f)

total = 0
verified = 0
failed = 0
no_transcript = 0
skipped_ambiguous = 0

with open(args.annotations) as f, open(args.output, "w") as out:

    header = f.readline().rstrip("\n")
    out.write(header + "\tExact_match\n")

    for line in f:
        fields = line.rstrip("\n").split("\t")
        pep = fields[0]
        feature_type = fields[2]
        transcripts_field = fields[4]

        total += 1

        if feature_type in ("multi_gene", "NA") or not transcripts_field:
            skipped_ambiguous += 1
            out.write(line.rstrip("\n") + "\tNA\n")
            continue

        transcripts = transcripts_field.split(";")
        seqs = [transcript_sequences[t] for t in transcripts if t in transcript_sequences]

        if not seqs:
            no_transcript += 1
            out.write(line.rstrip("\n") + "\tNo_reference_sequence\n")
            continue

        is_exact = any(pep in seq for seq in seqs)

        if is_exact:
            verified += 1
            out.write(line.rstrip("\n") + "\tPASS\n")
        else:
            failed += 1
            out.write(line.rstrip("\n") + "\tFAIL\n")

print("Total annotated peptides:", total)
if total == 0:
    raise SystemExit("No annotated peptides found in the input file.")
print(f"  Exact-substring verified (PASS): {verified} ({verified/total*100:.1f}%)")
print(f"  Failed verification (FAIL, false-positive mode): {failed} ({failed/total*100:.1f}%)")
print(f"  Skipped (multi_gene / NA / no transcript call): {skipped_ambiguous}")
print(f"  Assigned transcript missing from sequence reference: {no_transcript}")
