# PEXMap – Peptide Exon Mapping Tool

PEXMap is a proteogenomics tool designed to map experimental MS/MS–derived peptide sequences to their potential genomic origins. The tool uses a customized reference database of **k-mer peptides (default: 8-mers)** that links each peptide fragment to its corresponding **gene IDs**, **transcript IDs**, and **exon or exon-junction features**.

The workflow processes experimentally identified peptides by generating overlapping **k-mer fragments**, which are then queried against the reference annotation database. Matching entries allow the identification of candidate genomic regions that may give rise to the observed peptide sequences, enabling downstream analysis of exon usage and exon-junction–derived peptides in proteogenomic studies.

---

## Workflow

1. Provide experimentally identified peptide sequences.
2. Filter peptides with length ≥ k (default **8 amino acids**).
3. Generate overlapping **k-mer fragments** from the input peptides.
4. Search generated k-mers in the **reference peptide annotation database**.
5. Retrieve associated **gene, transcript, exon, or exon-junction annotations**.
6. Summarize dominant gene or transcript matches based on peptide hits.

---

## Reference Dataset

The peptide annotation database is large and therefore hosted externally.

Download the database from:

https://drive.google.com/uc?id=1jPU8HE6Fcwk4mAU8Fk5m7VJGLtnrrKKF

After downloading, place the file in:


data/peptide_dataset.pkl


---

## Installation

Clone the repository:


git clone https://github.com/deepanshicbg/PEXMap.git

cd PEXMap


Install dependencies:


pip install -r requirements.txt


---

## Usage

### Step 1 — Generate k-mer peptides

This step filters peptides shorter than the selected k-mer length and generates overlapping k-mer fragments.


python scripts/generate_kmers.py input_peptides.txt kmers.txt


---

### Step 2 — Annotate peptides

Search generated k-mers against the reference peptide database.


python scripts/annotate_peptides.py
--kmers kmers.txt
--database data/peptide_dataset.pkl
--organism human
--output annotations.tsv


---

## Arguments

| Argument     | Description                                |
| ------------ | ------------------------------------------ |
| `--kmers`    | File containing generated k-mer peptides   |
| `--database` | Reference peptide annotation database      |
| `--organism` | Organism name (e.g. human)                 |
| `--output`   | Output file containing peptide annotations |

---

## Building Your Own Peptide Database

PEXMap also allows users to generate their own peptide annotation database from organism annotation data.

If you have **ENACT-based transcript–exon annotation files** for an organism, you can generate the peptide database using the provided script:


scripts/build_peptide_database.py


This script reads gene-level annotation files containing:

- transcript IDs  
- exon identifiers  
- amino acid sequences  

and generates overlapping **k-mer peptides** indexed by:

- gene ID  
- transcript ID  
- exon ID  
- exon-junction ID  

The resulting database can then be used directly with the **PEXMap annotation pipeline**.

---

## Generate Database from ENACT Annotation Data

Example command:

```
python scripts/build_peptide_database.py
--input_folder organism_gene_files
--kmer 8
--organism human
--output peptide_dataset.pkl
```

---

## Input Format

Input peptide file should contain **one peptide sequence per line**:


MTEYKLVVVGAG
ADLASRDE
VAVWPTMV


---

## Example Run

Example peptide input file:
```
example/example_peptides.txt
```

Generate k-mer fragments from the example peptides:

```
 python scripts/generate_kmers.py example/example_peptides.txt example/example_kmers.txt
```

Annotate the generated peptides:


python scripts/annotate_peptides.py
--kmers example/example_kmers.txt
--database data/peptide_dataset.pkl
--organism human
--output example/example_output.tsv


---

## Output Format

The annotation output reports peptide matches and associated genomic features.

Example output:


| Experimental_MS_peptide | Gene_ids | Feature_type  |         Features        | Transcripts | Kmer_hits |
|-------------------------|----------|---------------|-------------------------|-------------|-----------|
|       ADLASRDEK         |   92283  | Exon-junction | T.1.A.7.0.0,T.1.A.9.c.1 | NP_694989.2 |     18    |
|     MGTFATLSELHCDK      |    3043  |     Exon      |       T.1.G.2.0.0       | NP_000509.1 |      6    |

---

## Repository Structure
```

PEXMap
│
├── scripts
│ ├── generate_kmers.py
│ ├── annotate_peptides.py
│ └── build_peptide_database.py
│
├── data
│ └── (place peptide_dataset.pkl here)
│
├── example
│ ├── example_peptides.txt
│ ├── example_kmers.txt
│ └── example_output.tsv
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

Deepanshi Awasthi, PhD Research Scholar, Computational Biology Group  
Indian Institute of Science Education and Research (IISER) Mohali, India

