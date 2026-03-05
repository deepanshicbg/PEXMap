# PEXMap – Peptide Exon Mapping Tool

PEXMap is a proteogenomics tool designed to map experimental MS/MS–derived peptide sequences to their potential genomic origins. The tool uses a customized reference database of **8-mer peptides** that links each peptide fragment to its corresponding **gene IDs**, **transcript IDs**, and **exon or exon-junction features**.

The workflow processes experimentally identified peptides by generating overlapping **8-mer fragments**, which are then queried against the reference annotation database. Matching entries allow the identification of candidate genomic regions that may give rise to the observed peptide sequences, enabling downstream analysis of exon usage and exon-junction–derived peptides in proteogenomic studies.

---

## Workflow

1. Provide experimentally identified peptide sequences.
2. Filter peptides with length ≥ 8 amino acids.
3. Generate overlapping **8-mer fragments** from the input peptides.
4. Search generated 8-mers in the **reference peptide annotation database**.
5. Retrieve associated **gene, transcript, exon, or exon-junction annotations**.
6. Summarize dominant gene or transcript matches based on peptide hits.

---

## Reference Dataset

The peptide annotation database is large and therefore hosted externally.

Download the database from:

https://drive.google.com/uc?id=1jPU8HE6Fcwk4mAU8Fk5m7VJGLtnrrKKF

After downloading, place the file in:

```
data/peptide_dataset.pkl
```

---

## Installation

Clone the repository:

```
git clone https://github.com/deepanshicbg/PEXMap.git
cd PEXMap
```

Install dependencies:

```
pip install -r requirements.txt
```

---

## Usage

### Step 1 — Generate 8-mer peptides

This step filters peptides shorter than 8 amino acids and generates overlapping 8-mer fragments.

```
python scripts/generate_kmers.py input_peptides.txt kmers.txt
```

---

### Step 2 — Annotate peptides

Search generated 8-mers against the reference peptide database.

```
python scripts/annotate_peptides.py \
    --kmers kmers.txt \
    --database data/peptide_dataset.pkl \
    --organism human \
    --output annotations.tsv
```

---

## Arguments

| Argument     | Description                                |
| ------------ | ------------------------------------------ |
| `--kmers`    | File containing generated 8-mer peptides   |
| `--database` | Reference peptide annotation database      |
| `--organism` | Organism name (e.g. human)                 |
| `--output`   | Output file containing peptide annotations |

---

## Input Format

Input peptide file should contain **one peptide sequence per line**:

```
MTEYKLVVVGAG
ADLASRDE
VAVWPTMV
```

---

## Example Run

An example input dataset is provided in the repository.

Example peptide input file:

```
example/example_peptides.txt
```

Generate 8-mer fragments from the example peptides:

```
python scripts/generate_kmers.py example/example_peptides.txt example/example_kmers.txt
```

Annotate the generated peptides:

```
python scripts/annotate_peptides.py \
    --kmers example/example_kmers.txt \
    --database data/peptide_dataset.pkl \
    --organism human \
    --output example/example_output.tsv
```

---

## Output Format

The annotation output reports peptide matches and associated genomic features.

Example output:

```
peptide    gene_id    feature_type    feature_id    transcript_ids
ADLASRDE   92283      junction        T.1.A.7.0.0,T.1.A.9.c.1    NP_694989.2,NP_001309750.1
```

---

## Repository Structure

```
PEXMap
│
├── scripts
│   ├── generate_kmers.py
│   └── annotate_peptides.py
│
├── data
│   └── (place peptide_dataset.pkl here)
│
├── example
│   ├── example_peptides.txt
│   ├── example_kmers.txt
│   └── example_output.tsv
│
├── README.md
├── requirements.txt
└── .gitignore
```

---

## Citation

If you use **PEXMap** in your research, please cite the associated publication (to be added).

---

## Author

Deepanshi Awasthi, PhD Research Scholar, Computational Biology Group, Indian Institute of Science Education and Research (IISER) Mohali, India

