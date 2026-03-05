import pickle
import sys
from collections import Counter

kmer_file = sys.argv[1]
db_file = sys.argv[2]
output_file = sys.argv[3]

with open(db_file, "rb") as f:
    db = pickle.load(f)

gene_counter = Counter()
transcript_counter = Counter()
feature_counter = Counter()

annotations = []

with open(kmer_file) as f:

    for line in f:

        pep = line.strip()

        if pep not in db:
            continue

        for entry in db[pep]:

            gene = entry["gene_id"]
            feature = entry["feature_id"]
            transcripts = entry["transcript_ids"]

            gene_counter[gene] += 1
            feature_counter[feature] += 1

            for tr in transcripts:
                transcript_counter[tr] += 1

            annotations.append((pep, gene, feature, transcripts))


#####################################
# determine dominant annotations
#####################################

top_gene = gene_counter.most_common(1)
top_transcript = transcript_counter.most_common(1)
top_feature = feature_counter.most_common(1)


#####################################
# write detailed output
#####################################

with open(output_file, "w") as out:

    out.write("peptide\tgene\tfeature\ttranscripts\n")

    for pep, gene, feature, transcripts in annotations:

        out.write(
            pep + "\t" +
            gene + "\t" +
            str(feature) + "\t" +
            ",".join(transcripts) + "\n"
        )


#####################################
# print summary
#####################################

print("\n===== SUMMARY =====")

if top_gene:
    print("Top gene:", top_gene[0][0], "hits:", top_gene[0][1])

if top_transcript:
    print("Top transcript:", top_transcript[0][0], "hits:", top_transcript[0][1])

if top_feature:
    print("Top feature:", top_feature[0][0], "hits:", top_feature[0][1])