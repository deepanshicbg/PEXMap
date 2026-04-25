# PEXMap – Peptide Exon Mapping Tool

PEXMap **(PeptideEXonMapper)** is an exon-aware proteogenomic framework developed to systematically map experimental **MS/MS-derived peptide sequences** to their genomic and transcriptomic origins. Unlike conventional peptide annotation methods that mainly assign peptides to genes or proteins, PEXMap enables **multi-level mapping** of peptides to **genes, transcript isoforms, exons, and exon–exon junctions**.

The method uses a customized searchable reference database built from human protein-coding transcript isoforms, where sequences are decomposed into overlapping **8-mer subsequences (octamerDB)**. Each 8-mer is indexed with its associated **gene ID, transcript/isoform ID, exon identifier (EUID), or exon-junction context**. A complementary **exon-junction database (ExonjunctionDB)** is also used to improve isoform-specific detection.

For analysis, input MS/MS peptides are filtered (minimum length ≥8 amino acids, excluding low-complexity peptides) and similarly decomposed into overlapping 8-mers. These are matched exactly against the indexed reference databases using **fast dictionary-based lookup**. Peptide assignments are then inferred using maximal matching evidence, allowing reliable identification of shared or uniquely mapped peptides.

PEXMap is particularly useful for detecting **isoform-specific peptide evidence**, resolving peptides originating from alternatively spliced exons, and identifying tissue- or disease-specific transcript usage directly from proteomics datasets.

---

## Workflow

1. Input experimentally identified **MS/MS peptide sequences**.
2. Filter peptides (**≥8 aa**, remove **low-complexity** and **ambiguous sequences**).
3. Generate overlapping **8-mer k-mers**.
4. Query k-mers against the **precomputed**`` .pkl`` **annotation database**.
5. Retrieve mapped **gene, transcript/isoform, exon, and exon–exon junction (EXj)** annotations.
6. Assign peptides using a **maximum k-mer match strategy** to determine dominant mappings.
7. Compute **k-mer coverage statistics**

---


## 📊 Coverage Metric

Coverage (%) = (Matched k-mers / Total unique k-mers) × 100

Where:
- Total unique k-mers = unique 8-mers derived from peptide  
- Matched k-mers = those found in the reference database  

---

## Reference Dataset

The peptide annotation database is large and therefore hosted externally.

Download the database from:

**ENACT v0.5 dataset (used in this study):**
```
https://drive.google.com/uc?id=1jPU8HE6Fcwk4mAU8Fk5m7VJGLtnrrKKF
```
**Latest dataset version (recommended):**
```
https://drive.google.com/file/d/124p7-jL0uxcfmfKGvxBqu2JAOXGH2GqM/view?usp=sharing
```
After downloading, place the file in:

``
data/peptide_dataset.pkl
``

---

## Installation

Clone the repository:

```
git clone https://github.com/deepanshicbg/PEXMap.git
cd PEXMap
```

---

## Usage

### Step 1 — Generate k-mer peptides

This step filters peptides shorter than the selected k-mer length and generates overlapping k-mer fragments.

```
python scripts/generate_kmers.py input_peptides.txt kmers.txt
```

---

### Step 2 — Annotate peptides

Search generated k-mers against the reference peptide database.

```
python scripts/annotate_peptides.py
--kmers kmers.txt
--database data/peptide_dataset.pkl
--organism human
--output annotations.tsv
```

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

``
scripts/build_peptide_database.py
``

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

``
MTEYKLVVVGAG
``

``
ADLASRDE
``

``
VAVWPTMV
``

---

## Example Run

Example peptide input file:
``
example/example_peptides.txt
``

Generate k-mer fragments from the example peptides:

```
 python scripts/generate_kmers.py example/example_MSpeptides.txt example/example_kmers.txt
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

## 📤 Output Columns

| Column                     | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| Experimental_MS_peptide  | Input peptide sequence from MS/MS experiment                               |
| Gene_id                  | Dominant gene selected based on maximum k-mer support                      |
| Feature_type             | Type of feature: `exon` or `junction`                                      |
| Features                 | Exon IDs or exon–exon junction identifiers associated with the peptide     |
| Transcripts              | Dominant transcript(s) belonging to the selected gene                      |
| Kmer_hits                | Total number of k-mer matches supporting the selected gene                 |
| Total_unique_kmers       | Number of unique k-mers derived from the peptide                           |
| Matched_kmers            | Number of k-mers that matched entries in the reference database            |
| Coverage_percent         | Percentage of peptide k-mers matched to the database                       |


## 📌 Example Output

| Experimental_MS_peptide | Gene_id | Feature_type | Features                         | Transcripts                         | Kmer_hits | Total_unique_kmers | Matched_kmers | Coverage_percent |
|------------------------|---------|--------------|----------------------------------|-------------------------------------|-----------|--------------------|---------------|------------------|
| AGSYGAQPVVQTQLNSYGAQA  | 10432   | exon         | T.1.A.2.0.0                      | NP_006319.1                         | 14        | 14                 | 14            | 100.0            |
| NYEENRQVNL             | 1825    | exon         | T.1.G.10.0.0                     | NP_001932.2;NP_077741.2             | 3         | 3                  | 3             | 100.0            |
| GALTGKQPDGSAE          | 9941    | junction     | T.1.A.1.n.2,D.1.A.3.0.0          | NP_005098.2;XP_005265692.1          | 5         | 6                  | 5             | 83.33            |

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


