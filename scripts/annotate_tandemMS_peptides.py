import pickle
import argparse
from collections import defaultdict, Counter


parser = argparse.ArgumentParser(
    description="PEXMAp annotations of MS/MS peptides"
)

parser.add_argument("--kmers", required=True) # specify k-mer file generated from MS/MS peptide dataset
parser.add_argument("--database", required=True) # specify octamerDB.pkl or your customized k-mer DB
parser.add_argument("--organism", required=False) # optional: We have used human here
parser.add_argument("--output", required=True) # specify the output file name

args = parser.parse_args()

kmer_file = args.kmers
db_file = args.database
output_file = args.output


####################################################################
# load database (octamerDB.pkl); downloaded from the link provided #
####################################################################

print("Loading database:", db_file)

with open(db_file, "rb") as f:
    db = pickle.load(f)

print("Total 8-mer sequences in octamerDB:", len(db))


####################################
# Build transcript -> gene mapping #
####################################

transcript_to_gene = defaultdict(set)

for kmer_entries in db.values():
    for entry in kmer_entries:

        gene = entry.get("gene_id")
        transcripts = entry.get("transcript_ids", [])

        if gene is None:
            continue

        for tr in transcripts:
            transcript_to_gene[tr].add(gene)

print("Total transcripts mapped:", len(transcript_to_gene))


############################
# Store 8-mers per peptide #
###########################

pep_kmers = defaultdict(list)

with open(kmer_file) as f:
    for line in f:

        parts = line.strip().split()

        if len(parts) != 2:
            continue

        exp_pep, kmer = parts
        pep_kmers[exp_pep].append(kmer)


######################
# PEXMap annotations #
######################

results = []

# Mapping coverage (Number of MS/MS peptide derived unique 8-mers matched / Total number of MS/MS peptide derived unique 8-mers)
cov_100 = 0
cov_80 = 0
cov_50 = 0
cov_30 = 0
cov_below_30 = 0


#For exon-level peptide mapping summary
exon_peptides = 0
junction_peptides = 0


for exp_pep, kmers in pep_kmers.items():

    unique_kmers = set(kmers)

    total_kmers = len(unique_kmers)
    matched_kmers_set = set()

    gene_counter = Counter()
    all_entries = []


    ###############
    # Match kmers #
    ###############

    for kmer in unique_kmers:

        if kmer not in db:
            continue

        matched_kmers_set.add(kmer)

        for entry in db[kmer]:

            gene = entry.get("gene_id")
            feature = entry.get("feature_id")
            transcripts = entry.get("transcript_ids", [])

            if gene is None:
                continue

            all_entries.append((gene, feature, transcripts))
            gene_counter[gene] += 1


    matched_kmers = len(matched_kmers_set)

    coverage = (matched_kmers / total_kmers * 100) if total_kmers > 0 else 0.0


    ####################
    # Coverage binning #
    ####################

    if coverage == 100:
        cov_100 += 1
    elif coverage >= 80:
        cov_80 += 1
    elif coverage >= 50:
        cov_50 += 1
    elif coverage >= 30:
        cov_30 += 1
    else:
        cov_below_30 += 1


    if not gene_counter:
        continue


    ####################################################
    # Select gene ID (Maximal 8-mer matching criteria) #
    ####################################################

    max_gene = max(gene_counter.values())
    top_genes = [g for g, c in gene_counter.items() if c == max_gene]
    selected_gene = sorted(top_genes)[0]


    ####################################
    # Filter entries for selected gene #
    ####################################

    filtered_entries = [e for e in all_entries if e[0] == selected_gene]

    if not filtered_entries:
        continue


    ############################################
    # Annotations at transcript and exon-level #
    ############################################

    transcript_counter = Counter()
    junction_counter = Counter()
    exon_counter = Counter()

    for gene, feature, transcripts in filtered_entries:

        # STRICT transcript filtering
        for tr in transcripts:

            # keep only transcripts belonging to selected gene
            if selected_gene not in transcript_to_gene.get(tr, set()):
                continue

            transcript_counter[tr] += 1

        # feature counting
        if feature:
            if "," in str(feature):
                junction_counter[feature] += 1
            else:
                exon_counter[feature] += 1


    ######################################
    # Transcript with maximum 8-mer hits #
    ######################################

    if transcript_counter:
        max_tr = max(transcript_counter.values())
        top_transcripts = [t for t, c in transcript_counter.items() if c == max_tr]
    else:
        top_transcripts = []


    ############################
    # Exonic-feature selection #
    ############################

    if junction_counter:
        max_feat = max(junction_counter.values())
        top_features = [f for f, c in junction_counter.items() if c == max_feat]
        feature_type = "junction"
        junction_peptides += 1

    elif exon_counter:
        max_feat = max(exon_counter.values())
        top_features = [f for f, c in exon_counter.items() if c == max_feat]
        feature_type = "exon"
        exon_peptides += 1

    else:
        top_features = []
        feature_type = "NA"


    ###########################
    # Save PEXMap annotations #
    ###########################

    results.append(
        (
            exp_pep,
            selected_gene,
            feature_type,
            ";".join(sorted(map(str, top_features))),
            ";".join(sorted(map(str, top_transcripts))),
            max_gene,
            total_kmers,
            matched_kmers,
            round(coverage, 2)
        )
    )


######################################
# Write results into the output file # 
######################################

with open(output_file, "w") as out:

    out.write(
        "Experimental_MS_peptide\tGene_id\tFeature_type\tFeatures\tTranscripts\tKmer_hits\tTotal_unique_kmers\tMatched_kmers\tCoverage_percent\n"
    )

    for r in results:
        out.write("\t".join(map(str, r)) + "\n")


######################
# Summary statistics #
######################

total_peptides = len(pep_kmers)

print("\nCoverage Summary:")
print("-------------------------")
print(f"Total peptides processed: {total_peptides}")
print(f"Coverage = 100% : {cov_100}")
print(f"100% > Coverage > = 80% : {cov_80}")
print(f"80% > Coverage > = 50% : {cov_50}")
print(f"50% > Coverage > = 30% : {cov_30}")
print(f"coverage < 30% : {cov_below_30}")
print("-------------------------")

print("\nFeature Mapping Summary:")
print("-------------------------")
print(f"Peptides mapped to exons: {exon_peptides}")
print(f"Peptides mapped to junctions: {junction_peptides}")
print("-------------------------")


print("Annotated experimental peptides:", len(results))