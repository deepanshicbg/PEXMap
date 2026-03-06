import pickle
import argparse
from collections import defaultdict, Counter


parser = argparse.ArgumentParser(description="Annotate experimental peptides using 8-mer database")

parser.add_argument("--kmers", required=True)
parser.add_argument("--database", required=True)
parser.add_argument("--organism", required=False)
parser.add_argument("--output", required=True)

args = parser.parse_args()

kmer_file = args.kmers
db_file = args.database
output_file = args.output


#####################################
# load database
#####################################

print("Loading database:", db_file)

with open(db_file, "rb") as f:
    db = pickle.load(f)

print("Total peptides in DB:", len(db))


#####################################
# store kmers per experimental peptide
#####################################

pep_kmers = defaultdict(list)

with open(kmer_file) as f:
    for line in f:

        parts = line.strip().split()

        if len(parts) != 2:
            continue

        exp_pep, kmer = parts
        pep_kmers[exp_pep].append(kmer)


#####################################
# annotation aggregation
#####################################

results = []

for exp_pep, kmers in pep_kmers.items():

    gene_counter = Counter()
    transcript_counter = Counter()

    junction_counter = Counter()
    exon_counter = Counter()

    for kmer in kmers:

        if kmer not in db:
            continue

        for entry in db[kmer]:

            gene = entry.get("gene_id")
            feature = entry.get("feature_id")
            transcripts = entry.get("transcript_ids", [])

            gene_counter[gene] += 1

            for tr in transcripts:
                transcript_counter[tr] += 1

            if feature and "," in feature:
                junction_counter[feature] += 1
            else:
                exon_counter[feature] += 1


    if not gene_counter:
        continue


    ################################
    # dominant gene
    ################################

    max_gene = max(gene_counter.values())
    top_genes = [g for g,c in gene_counter.items() if c == max_gene]


    ################################
    # dominant transcript
    ################################

    max_tr = max(transcript_counter.values())
    top_transcripts = [t for t,c in transcript_counter.items() if c == max_tr]


    ################################
    # feature preference
    ################################

    if junction_counter:

        max_feat = max(junction_counter.values())
        top_features = [f for f,c in junction_counter.items() if c == max_feat]

        feature_type = "junction"

    else:

        max_feat = max(exon_counter.values())
        top_features = [f for f,c in exon_counter.items() if c == max_feat]

        feature_type = "exon"


    ################################
    # store results
    ################################

    results.append(
        (
            exp_pep,
            ";".join(sorted(map(str,top_genes))),
            feature_type,
            ";".join(sorted(map(str,top_features))),
            ";".join(sorted(map(str,top_transcripts))),
            max_gene
        )
    )


#####################################
# write output
#####################################

with open(output_file,"w") as out:

    out.write("Experimental_MS_peptide\tGene_ids\tFeature_type\tFeatures\tTranscripts\tKmer_hits\n")

    for r in results:
        out.write("\t".join(map(str,r))+"\n")


print("Annotated experimental peptides:",len(results))