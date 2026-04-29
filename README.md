# PEXMap – Peptide Exon Mapping Tool

PEXMap **(PeptideEXonMapper)** is an exon-aware proteogenomic framework developed to systematically map experimental **MS/MS-derived peptide sequences** to their genomic and transcriptomic origins. Unlike conventional peptide annotation methods that mainly assign peptides to genes or proteins, PEXMap enables **multi-level mapping** of peptides to **genes, transcript isoforms, exons, and exon–exon junctions**.

The method uses a customized searchable reference database built from human protein-coding transcript isoforms, where sequences are decomposed into overlapping **8-mer subsequences (octamerDB)**. Each 8-mer is indexed with its associated **gene ID, transcript/isoform ID, exon identifier (EUID), or exon-exon junction context**. A complementary **exon-exon junction database (ExonjunctionDB)** is also used to improve isoform-specific detection.

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

Coverage (%) = (Matched MS peptide derived k-mers / Total unique MS peptide derived k-mers) × 100

Where:
- Total unique MS peptide derived 8-mers = unique 8-mers derived from MS peptides  (Experimental MS/MS or query peptides)
- Matched MS peptide derived 8-mers  = derived 8-mers that could match the 8-mers in reference database (exon octamerDB + exonjunctionDB) 

---
## Installation

Clone the repository:

```
git clone https://github.com/deepanshicbg/PEXMap.git
cd PEXMap
```

---


## Reference Dataset

The peptide annotation database is large and therefore hosted externally.

Download the database from:

**ENACT v0.5 dataset: octamerDB (used in this study):**
```
https://drive.google.com/uc?id=1jPU8HE6Fcwk4mAU8Fk5m7VJGLtnrrKKF
```
**Latest dataset version: octamerDB_latest (recommended):**
```
https://drive.google.com/file/d/124p7-jL0uxcfmfKGvxBqu2JAOXGH2GqM/view?usp=sharing
```
After downloading, create a 'data' directory in PEXMap folder and place the file there:

``
PEXMap/data/octamerDB.pkl
``

---

## Usage

### Step 1 — Generate k-mer peptides

This step filters peptides shorter than the selected k-mer length and generates overlapping k-mer fragments. You can customize k-mer length as per your customized kmerDB. Here, the default length is 8. 

```
python scripts/generate_kmers.py --input input_peptides.txt --output kmers.txt --kmer_length 8
```

---

### Step 2 — Annotate peptides

Search generated k-mers against the reference peptide database.

```
python scripts/annotate_tandemMS_peptides.py \
--kmers kmers.txt \
--database data/octamerDB.pkl \
--organism human \
--output PEXMap_annotations.tsv
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
- exon- exon junction ID  

The resulting database can then be used directly with the **PEXMap annotation pipeline**.

---

## Generate Database from ENACT Annotation Data

Example command:

```
python scripts/build_peptide_database.py \
--input_folder organism_gene_files \
--kmer 8 \
--organism human \
--output kmerDB.pkl
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
example/example_tandem_MSpeptides.txt
``

Generate k-mer fragments from the example peptides:

```
 python scripts/generate_kmers.py example/example_tandem_MSpeptides.txt example/example_kmers.txt
```

Annotate the generated peptides:

```
python scripts/annotate_tandemMS_peptides.py \
--kmers example/example_kmers.txt \
--database data/octamerDB.pkl \
--organism human \
--output example/example_PEXMap_annotations.tsv
```

---

## Output Format

The annotation output reports peptide matches and associated genomic features.

**Output Columns**

| Column                     | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| Experimental_MS_peptide  | Input peptide sequence from MS/MS experimental data                               |
| Gene_id                  | Gene ID selected based on maximum k-mer support for mapped MS/MS peptide   |
| Feature_type             | Type of feature: `exon` or `junction`                                      |
| Features                 | Exon IDs or exon–exon junction identifiers associated with the peptide     |
| Transcripts              | Maximally matched (highest k-mer hits) transcript(s) belonging to the selected gene |
| Total_unique_kmers       | Number of unique k-mers derived from the MS/MS peptide                           |
| Matched_kmers            | Number of k-mers that matched entries in the reference database (OctamerDB or kmerDB)   |
| Coverage_percent         | Percentage of MS/MS peptide derived k-mers matched to the database                       |



**Example Output**


| Experimental_MS_peptide | Gene_id | Feature_type | Features                         | Transcripts                              | Total_unique_kmers | Matched_kmers | Coverage_percent |
|------------------------|---------|--------------|----------------------------------|------------------------------------------|--------------------|---------------|------------------|
| TKAIDMCPKNASY          | 7266    | junction     | D.1.G.4.0.0,T.1.G.5.0.0          | NP_003306.3                              | 6                  | 6             | 100.0            |
| EEEDDSALPQEVSI         | 80829   | exon         | T.1.G.3.0.0                      | NP_001183980.1;NP_444251.1               | 7                  | 7             | 100.0            |
| IGKAKTKENRQSIINPDWNFEKM| 474170  | junction     | T.1.A.7.0.0,T.1.A.8.0.0          | XP_047292103.1                           | 16                 | 16            | 100.0            |
| SYAAQQHPQAAASY         | 10432   | exon         | T.1.A.2.0.0                      | NP_006319.1                              | 7                  | 7             | 100.0            |
| GQSEADSDKNATILELR      | 1832    | exon         | T.1.F.23.0.0                     | NP_004406.2                              | 10                 | 8             | 80.0             |

### Interpretation

- Peptides with **100% coverage** indicate that all possible k-mers are supported by the reference database, providing high-confidence mapping.

- **GQSEADSDKNATILELR**  
  → Shows **80% coverage**, meaning some k-mers are not found in the database.  
  → This may indicate:
    - partial sequence support  
    - database incompleteness  
    - biological variation (e.g., mutations or alternative splicing)

- Junction peptides (e.g., *TKAIDMCPKNASY*, *IGKAKTKENRQSIINPDWNFEKM*) provide **splice-aware evidence**, supporting transcript-specific mapping.

- Exon-mapped peptides may be:
  - **unique** to a transcript  
  - or **shared across multiple isoforms** of the same gene


## Summary Statistics

At the end of execution, PEXMap reports:

- Number of peptides with **100% coverage**
- Number of peptides with **≥80%, ≥50%, ≥30%, and <30% coverage**
- Number of peptides mapped to:
  - **exons**
  - **exon–exon junctions**

These metrics provide a quick assessment of mapping quality and splice-aware peptide evidence.
---

## Repository Structure
```

PEXMap
│
├── scripts
│ ├── generate_kmers.py
│ ├── annotate_tandemMS_peptides.py
│ └── build_kmer_database.py
│
├── data
│ └── (place peptide_dataset.pkl here)
│
├── example
│ ├── example_tandem_MSpeptides.txt
│ ├── example_kmers.txt
│ └── example_PEXMap_annotations.tsv
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


