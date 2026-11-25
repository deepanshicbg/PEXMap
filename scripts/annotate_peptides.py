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
parser.add_argument(
    "--peptides", required=False, default=None,
    help="optional: the original peptide input file given to generate_kmers.py; "
         "lets the summary report input peptides and those excluded for being shorter than k"
)
parser.add_argument(
    "--kmer_length", type=int, default=8,
    help="k-mer length used in generate_kmers.py (only used with --peptides; default: 8)"
)
parser.add_argument(
    "--min_coverage", type=float, default=0.0,
    help="minimum k-mer coverage (%%) for a peptide to be written to the annotation output "
         "(default: 0, i.e. keep every peptide with a gene match). Peptides below the "
         "threshold are counted in the summary and listed in <output>.below_coverage.txt"
)
parser.add_argument(
    "--unmapped_output", required=False, default=None,
    help="file listing processed peptides with no match in the database "
         "(default: <output>.unmapped.txt)"
)

args = parser.parse_args()

kmer_file = args.kmers
db_file = args.database
output_file = args.output
organism_filter = args.organism.lower() if args.organism else None
min_coverage = args.min_coverage
unmapped_file = args.unmapped_output or (output_file + ".unmapped.txt")
below_cov_file = output_file + ".below_coverage.txt"

if not 0 <= min_coverage <= 100:
    parser.error("--min_coverage must be between 0 and 100")


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
multi_gene_peptides = 0

# For QC reporting (every processed peptide ends up in exactly one of these)
unmapped_peptides = []          # no k-mer matched any gene in the database
below_coverage_peptides = []    # matched, but coverage < --min_coverage


for exp_pep, kmers in pep_kmers.items():

    unique_kmers = set(kmers)

    total_kmers = len(unique_kmers)
    matched_kmers_set = set()

    # entries_by_kmer holds, per matched k-mer, the organism-filtered database
    # entries for that k-mer -- kept so gene/transcript/feature support below
    # can be counted by DISTINCT SUPPORTING K-MER rather than by database entry
    # (a single k-mer that appears in more than one EUID entry of the same gene
    # must only count once for that gene).
    entries_by_kmer = {}

    for kmer in unique_kmers:

        if kmer not in db:
            continue

        entries = db[kmer]

        if organism_filter is not None:
            # backward compatible: entries without an "organism" field are not filtered out
            entries = [
                e for e in entries
                if e.get("organism") is None or str(e.get("organism")).lower() == organism_filter
            ]
            if not entries:
                continue

        matched_kmers_set.add(kmer)
        entries_by_kmer[kmer] = entries


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


    if not entries_by_kmer:
        unmapped_peptides.append(exp_pep)
        continue


    ##########################################################
    # Gene-level support = number of DISTINCT k-mers per gene #
    ##########################################################

    gene_kmer_support = defaultdict(set)

    for kmer, entries in entries_by_kmer.items():
        genes_this_kmer = {e.get("gene_id") for e in entries if e.get("gene_id") is not None}
        for gene in genes_this_kmer:
            gene_kmer_support[gene].add(kmer)

    if not gene_kmer_support:
        unmapped_peptides.append(exp_pep)
        continue


    ###############################################
    # User-defined coverage threshold (QC filter) #
    ###############################################

    if coverage < min_coverage:
        below_coverage_peptides.append((exp_pep, round(coverage, 2)))
        continue

    gene_counter = Counter({g: len(s) for g, s in gene_kmer_support.items()})


    ####################################################
    # Select gene ID (Maximal distinct-kmer matching)   #
    ####################################################

    max_gene = max(gene_counter.values())
    top_genes = sorted(g for g, c in gene_counter.items() if c == max_gene)


    if len(top_genes) > 1:

        # Genuine tie at the gene level: report all tied genes explicitly
        # instead of silently picking one by lexicographic string sort
        # Feature/transcript resolution is not
        # meaningful while the gene call itself is ambiguous.
        multi_gene_peptides += 1

        results.append(
            (
                exp_pep,
                ";".join(top_genes),
                "multi_gene",
                "",
                "",
                max_gene,
                total_kmers,
                matched_kmers,
                round(coverage, 2)
            )
        )
        continue


    selected_gene = top_genes[0]


    ##########################################################################
    # Transcript / exon / junction support, restricted to the selected gene, #
    # also counted by DISTINCT supporting k-mer                              #
    ##########################################################################

    transcript_kmer_support = defaultdict(set)
    exon_kmer_support = defaultdict(set)
    junction_kmer_support = defaultdict(set)

    for kmer, entries in entries_by_kmer.items():
        for entry in entries:

            if entry.get("gene_id") != selected_gene:
                continue

            feature = entry.get("feature_id")
            transcripts = entry.get("transcript_ids", []) or []

            # STRICT transcript filtering: keep only transcripts belonging to selected gene
            for tr in transcripts:
                if selected_gene not in transcript_to_gene.get(tr, set()):
                    continue
                transcript_kmer_support[tr].add(kmer)

            if feature:
                if "," in str(feature):
                    junction_kmer_support[feature].add(kmer)
                else:
                    exon_kmer_support[feature].add(kmer)


    ######################################
    # Transcript with maximum k-mer hits #
    ######################################

    if transcript_kmer_support:
        transcript_counter = Counter({t: len(s) for t, s in transcript_kmer_support.items()})
        max_tr = max(transcript_counter.values())
        top_transcripts = sorted(t for t, c in transcript_counter.items() if c == max_tr)
    else:
        top_transcripts = []


    ###########################################################################
    # Exon-vs-junction: genuine maximal-matching comparison across evidence   #
    ###########################################################################

    exon_counter = Counter({f: len(s) for f, s in exon_kmer_support.items()})
    junction_counter = Counter({f: len(s) for f, s in junction_kmer_support.items()})

    best_exon = max(exon_counter.values()) if exon_counter else 0
    best_junction = max(junction_counter.values()) if junction_counter else 0

    if best_exon == 0 and best_junction == 0:
        feature_type = "NA"
        top_features = []

    elif best_junction > best_exon:
        feature_type = "junction"
        top_features = sorted(f for f, c in junction_counter.items() if c == best_junction)
        junction_peptides += 1

    elif best_exon > best_junction:
        feature_type = "exon"
        top_features = sorted(f for f, c in exon_counter.items() if c == best_exon)
        exon_peptides += 1

    else:
        # genuine tie between best exon and best junction support
        feature_type = "exon+junction"
        top_features = (
            sorted(f for f, c in exon_counter.items() if c == best_exon)
            + sorted(f for f, c in junction_counter.items() if c == best_junction)
        )
        exon_peptides += 1
        junction_peptides += 1


    ###########################
    # Save PEXMap annotations #
    ###########################

    results.append(
        (
            exp_pep,
            selected_gene,
            feature_type,
            ";".join(map(str, top_features)),
            ";".join(map(str, top_transcripts)),
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


#################################
# Unmapped / below-coverage lists #
#################################

with open(unmapped_file, "w") as out:
    for pep in unmapped_peptides:
        out.write(pep + "\n")

if min_coverage > 0:
    with open(below_cov_file, "w") as out:
        out.write("Experimental_MS_peptide\tCoverage_percent\n")
        for pep, cov in below_coverage_peptides:
            out.write(f"{pep}\t{cov}\n")


######################
# Summary statistics #
######################

total_peptides = len(pep_kmers)

print("\nInput / QC Summary:")
print("-------------------------")

if args.peptides:
    input_peptides = set()
    with open(args.peptides) as f:
        for line in f:
            pep = line.strip()
            if pep:
                input_peptides.add(pep)
    short_peptides = sum(1 for p in input_peptides if len(p) < args.kmer_length)
    print(f"Input peptides (unique): {len(input_peptides)}")
    print(f"Excluded because shorter than {args.kmer_length} aa: {short_peptides}")
else:
    print("Input peptides: not reported (pass --peptides <input file> to include)")

print(f"Total peptides processed: {total_peptides}")
print(f"  Annotated (written to {output_file}): {len(results)}")
if min_coverage > 0:
    print(f"  Below coverage threshold ({min_coverage}%): {len(below_coverage_peptides)}"
          f" (listed in {below_cov_file})")
print(f"  Unmapped (no database match): {len(unmapped_peptides)} (listed in {unmapped_file})")

check = len(results) + len(below_coverage_peptides) + len(unmapped_peptides)
if check != total_peptides:
    print(f"WARNING: counts do not reconcile ({check} vs {total_peptides} processed)")
print("-------------------------")

print("\nCoverage Summary (all processed peptides):")
print("-------------------------")
print(f"Coverage = 100% : {cov_100}")
print(f"100% > Coverage > = 80% : {cov_80}")
print(f"80% > Coverage > = 50% : {cov_50}")
print(f"50% > Coverage > = 30% : {cov_30}")
print(f"coverage < 30% : {cov_below_30}")
print("-------------------------")

print("\nFeature Mapping Summary (annotated peptides):")
print("-------------------------")
print(f"Peptides mapped to exons: {exon_peptides}")
print(f"Peptides mapped to junctions: {junction_peptides}")
print(f"Peptides with tied (multi_gene) gene-level calls: {multi_gene_peptides}")
print("(exon+junction ties are counted under both exons and junctions)")
print("-------------------------")

print("\nNote: low-coverage mappings can arise from contaminant proteins (e.g. keratins,")
print("trypsin) if they were not removed before running PEXMap; consider --min_coverage.")
print("Next step (recommended): scripts/verify_mapping.py confirms each assignment")
print("by exact substring match against the assigned transcript sequence(s).")
