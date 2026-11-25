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
7. Compute **k-mer coverage statistics** and a QC summary (excluded, unmapped and below-threshold peptides).
8. (Recommended) **Verify** each assignment by exact substring match against the assigned transcript sequence(s).

---


## 📊 MS/MS peptide coverage metric

**Coverage (%)** = (Matched MS/MS peptide derived k-mers / Total unique MS/MS peptide derived k-mers) × 100

Where:
- **Total unique MS/MS peptide derived k-mers** = Number of unique MS/MS peptide derived k-mers  
- **Matched MS/MS peptide derived k-mers**  = Number of unique MS/MS peptide derived k-mers that could match the k-mers in the reference database (In our case, k=8, Reference database is: octamerDB) 

---
## Prerequisites

- **Python ≥ 3.8**, accessible from the command line.
- **No third-party Python packages are required.** All scripts use only the Python standard library (`argparse`, `pickle`, `collections`, `os`); `requirements.txt` documents this.
- **Operating system:** developed on Linux/macOS; also runs on Windows (see [Troubleshooting](#troubleshooting) for the command-line syntax).
- **Memory:** the reference database (`octamerDB.pkl`, ~1 GB on disk) is loaded fully into memory, so a machine with ample RAM is needed for full-scale runs (see [Runtime and resources](#runtime-and-resources)).

---
## Installation

1. Clone the repository:

```
git clone https://github.com/deepanshicbg/PEXMap.git
cd PEXMap
```

2. (Recommended) create an isolated environment so the Python version is fixed:

```
# with conda
conda create -n pexmap python=3.10
conda activate pexmap

# or with venv
python3 -m venv pexmap_env
source pexmap_env/bin/activate        # Windows: pexmap_env\Scripts\activate
```

3. Install dependencies (none beyond the standard library; this command simply confirms the environment):

```
pip install -r requirements.txt
```

4. Check the installation:

```
python scripts/generate_kmers.py --help
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
After downloading, create a `data` directory in the PEXMap folder and place the file there:

```
PEXMap/data/octamerDB.pkl
```

If you downloaded the latest version (`octamerDB_latest.pkl`), either rename it to `octamerDB.pkl` or pass its actual file name to `--database`.

> Note: `.pkl` files are Python pickle objects. Only load database files obtained from the links above or built yourself with `scripts/build_database.py`.

---

## Input Format for the MS/MS peptide data (input_peptides.txt)

Input should be **already identified, FDR-controlled peptides** (PEXMap does not perform spectrum matching). The file should contain **one peptide sequence per line**, in uppercase one-letter amino-acid code, without modifications or flanking residues:

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

## Steps to annotate your MS/MS peptide data using PEXMap:

### Step 1 — Generate k-mer peptides

This step generates overlapping k-mer fragments. Peptides shorter than the k-mer length cannot be mapped; they are **counted and reported** (not silently dropped), and can optionally be saved to a file. Duplicate peptides are processed once. You can customize the k-mer length to match your customized kmerDB; the default is 8.

```
python scripts/generate_kmers.py --input input_peptides.txt --output kmers.txt --kmer_length 8 --excluded_output excluded_short_peptides.txt
```

The script prints a summary, for example:

```
k-mer Generation Summary:
-------------------------
Input peptides: 3725
Excluded because shorter than 8 aa: 336
Processed unique peptides (k-mers generated): 3389
Total unique k-mers written: ...
-------------------------
```

---

### Step 2 — Annotate peptides

Search generated k-mers against the reference peptide database.

```
python scripts/annotate_peptides.py \
--kmers kmers.txt \
--database data/octamerDB.pkl \
--organism human \
--peptides input_peptides.txt \
--min_coverage 0 \
--output PEXMap_annotations.tsv
```

This writes:

- `PEXMap_annotations.tsv`: the annotation table (see [Output Format](#output-format))
- `PEXMap_annotations.tsv.unmapped.txt`: processed peptides with no match in the database
- `PEXMap_annotations.tsv.below_coverage.txt`: peptides below `--min_coverage` (only when a threshold > 0 is set)

---

### Step 3 — Verify assignments by exact substring match (recommended)

k-mer voting confirms that every 8-mer of a peptide occurs in the reference, but not that the whole peptide occurs **contiguously** in the assigned isoform. This step checks each assignment against the full protein sequence(s) of the assigned transcript(s) and adds an `Exact_match` column (`PASS`, `FAIL`, `NA` for multi_gene/no-transcript calls, or `No_reference_sequence`). In our benchmarks ~0.9% of assignments failed this check; we recommend discarding `FAIL` rows.

First build the transcript-sequence reference once, from the same ENACT per-gene annotation files used for the database:

```
python scripts/extract_sequences.py --input_folder organism_gene_files --output data/transcript_sequences.pkl
```

Then verify:

```
python scripts/verify_mapping.py \
--annotations PEXMap_annotations.tsv \
--sequences data/transcript_sequences.pkl \
--output PEXMap_annotations_verified.tsv
```

---

## Arguments

**`generate_kmers.py`**

| Argument            | Description                                                         |
| ------------------- | ------------------------------------------------------------------- |
| `--input`           | Input peptide file (one sequence per line)                          |
| `--output`          | Output k-mer file (`peptide<TAB>k-mer`)                             |
| `--kmer_length`     | k-mer length (default: 8; must match the database)                  |
| `--excluded_output` | Optional file listing peptides shorter than `--kmer_length`         |

**`annotate_peptides.py`**

| Argument            | Description                                                                                       |
| ------------------- | ------------------------------------------------------------------------------------------------- |
| `--kmers`           | k-mer file produced by `generate_kmers.py`                                                         |
| `--database`        | Reference annotation database (octamerDB or your own kmerDB)                                      |
| `--organism`        | Optional organism filter (e.g. `human`); database entries from other organisms are ignored        |
| `--output`          | Output annotation table (TSV)                                                                      |
| `--peptides`        | Optional: the original input peptide file, so the summary also reports input and excluded counts  |
| `--kmer_length`     | k-mer length used in Step 1 (only used with `--peptides`; default: 8)                              |
| `--min_coverage`    | Minimum k-mer coverage (%) for a peptide to be written to the output (default: 0 = keep all)      |
| `--unmapped_output` | File for unmapped peptides (default: `<output>.unmapped.txt`)                                       |

**Choosing `--min_coverage`:** 100% coverage means every k-mer of the peptide is supported by the reference. Lower-coverage mappings rest on partial evidence and are more likely to be noise (for example from contaminants or sequence variants). A threshold of 100 gives the most conservative annotation; 80 retains peptides with minor mismatches. The default (0) keeps everything so that users can apply their own cut-off downstream.

---

## Building Your Own Peptide Database

PEXMap also allows users to generate their own peptide annotation database from organism annotation data.

If you have **ENACT-based transcript–exon annotation files** for an organism, then you can generate the peptide database using the provided script:

``
scripts/build_database.py
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

## Generate Database from ENACTDB Annotated Genomes 

Example command:

```
python scripts/build_database.py \
--input_folder organism_gene_files \
--kmer 8 \
--organism human \
--output kmerDB.pkl
```


## Example Run

Example peptide input file:
``
example/example_tandem_MSpeptides.txt
``

The `example/` folder contains:

| File | Description |
| ---- | ----------- |
| `example_tandem_MSpeptides.txt` | 1,000 human MS/MS peptides (input) |
| `example_kmers.txt` | k-mers generated from them in Step 1 (for comparison) |
| `example_output.tsv` | expected annotation output of Step 2 |

Generate k-mer fragments from the example peptides:

```
python scripts/generate_kmers.py --input example/example_tandem_MSpeptides.txt --output example/my_kmers.txt --kmer_length 8
```

Expected summary:

```
Input peptides: 1000
Excluded because shorter than 8 aa: 0
Processed unique peptides (k-mers generated): 1000
Total unique k-mers written: 9946
```

Annotate the generated peptides:

```
python scripts/annotate_peptides.py \
--kmers example/my_kmers.txt \
--database data/octamerDB.pkl \
--organism human \
--peptides example/example_tandem_MSpeptides.txt \
--output example/my_PEXMap_annotations.tsv
```

Expected summary (compare `example/my_PEXMap_annotations.tsv` with `example/example_output.tsv`):

```
Total peptides processed: 1000
  Annotated: 1000
  Unmapped (no database match): 0

Coverage = 100% : 990
100% > Coverage > = 80% : 5
80% > Coverage > = 50% : 4
50% > Coverage > = 30% : 1
coverage < 30% : 0

Peptides mapped to exons: 884
Peptides mapped to junctions: 128
Peptides with tied (multi_gene) gene-level calls: 3
```

(Feature types in the output: 869 `exon`, 113 `junction`, 15 `exon+junction`, 3 `multi_gene`; `exon+junction` ties are counted under both exons and junctions in the summary.) Minor differences are expected if you use a different database version than the one used to generate `example_output.tsv`.

---

## Output Format

The annotation output reports peptide matches and associated genomic features.

**Output Columns**

| Column                     | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| Experimental_MS_peptide  | Input peptide sequence from MS/MS experimental data                               |
| Gene_id                  | Gene ID selected based on maximum distinct-k-mer support for the mapped MS/MS peptide; semicolon-separated list of tied gene IDs when `Feature_type` is `multi_gene` |
| Feature_type             | Type of feature: `exon`, `junction` (exon-exon junction), `exon+junction` (genuine tie between the best exon and best junction support), `multi_gene` (gene-level call is tied and therefore not resolved), or `NA` (no feature-level evidence) |
| Features                 | Exon IDs or exon–exon junction identifiers associated with the peptide; empty when `Feature_type` is `multi_gene` |
| Transcripts              | Maximally matched (highest distinct-k-mer support) transcript(s) belonging to the selected gene; empty when `Feature_type` is `multi_gene` |
| Kmer_hits                | Number of distinct k-mers supporting the selected (or tied) gene call        |
| Total_unique_kmers       | Number of unique k-mers derived from the MS/MS peptide                           |
| Matched_kmers            | Number of k-mers that matched entries in the reference database (OctamerDB or kmerDB)   |
| Coverage_percent         | Percentage of MS/MS peptide derived k-mers matched to the reference database (k-mer hit percentage)                       |



**Example Output**


| Experimental_MS_peptide | Gene_id | Feature_type | Features | Transcripts | Kmer_hits | Total_unique_kmers | Matched_kmers | Coverage_percent |
|------------------------|---------|--------------|----------|-------------|-----------|--------------------|---------------|-------------------|
| AGSYGAQPVVQTQLNSYGAQA | 10432 | exon | T.1.A.2.0.0 | NP_006319.1 | 14 | 14 | 14 | 100.0 |
| CLGVHCDFHQGLLSFYNAR | 79187 | junction | T.1.F.11.0.0,T.1.G.12.0.0 | NP_077309.1 | 12 | 12 | 12 | 100.0 |
| INNPMFRFLFENQTPAHVYYR | 23350 | exon+junction | T.1.G.16.0.0;T.1.G.15.0.0,T.1.G.16.0.0;T.2.G.15.0.0,T.1.G.16.0.0 | NP_001073884.1;NP_001307148.1;NP_001307149.1;NP_001307151.1 | 14 | 14 | 14 | 100.0 |
| IGKAKTKENRQSIINPDWNFEKM | 474170;4905 | multi_gene | | | 16 | 16 | 16 | 100.0 |
| GQSEADSDKNATILELR | 1832 | exon | T.1.F.23.0.0 | NP_004406.2 | 8 | 10 | 8 | 80.0 |

### Interpretation

- Peptides with **100% coverage** indicate that all possible k-mers are supported by the reference database, providing high-confidence mapping.

- **GQSEADSDKNATILELR**
  → Shows **80% coverage**, meaning some k-mers are not found in the database.
  → This may indicate:
    - partial sequence support
    - database incompleteness
    - biological variation (e.g., mutations or alternative splicing)

- Junction peptides (e.g., *CLGVHCDFHQGLLSFYNAR*) provide **splice-aware evidence**, supporting transcript-specific mapping.

- **exon+junction** (e.g., *INNPMFRFLFENQTPAHVYYR*) marks a genuine tie between the best exon-level support and the best junction-level support for the same gene — both feature sets are reported rather than arbitrarily picking one.

- **multi_gene** (e.g., *IGKAKTKENRQSIINPDWNFEKM*) marks a genuine tie at the gene level itself: two or more genes are supported by the same number of distinct k-mers, so all tied gene IDs are listed and no feature/transcript call is made, rather than silently resolving the tie by sorting gene IDs as strings.

- Exon-mapped peptides may be:
  - **unique** to a transcript
  - or **shared across multiple isoforms** of the same gene


## Summary Statistics

At the end of execution, PEXMap reports:

- **Input / QC summary:** input peptides and those excluded for being shorter than k (with `--peptides`), total peptides processed, and how many were **annotated**, **below the coverage threshold**, or **unmapped**. These always add up to the number processed, so no peptide is lost without being reported.
- Number of peptides with **100% coverage** and with **≥80%, ≥50%, ≥30%, and <30% coverage** (all processed peptides)
- Number of annotated peptides mapped to **exons**, **exon–exon junctions**, and with tied **multi_gene** calls

These metrics provide a quick assessment of mapping quality and splice-aware peptide evidence.

### Contaminants and low-coverage mappings

Common contaminants (e.g. keratins, trypsin, BSA from the cRAP list) should be **removed before running PEXMap**, typically during peptide identification. Contaminant peptides can partially match human sequences and appear as low-coverage mappings, i.e. potential false positives. Use `--min_coverage` and the exact-substring verification step (Step 3) to control these.

---

## Runtime and resources

| Task | Approximate runtime |
| ---- | ------------------- |
| Example dataset (1,000 peptides), Steps 1–3 | ~2–3 min, dominated by loading the database |
| Complete human PeptideAtlas build (Step 2) | ~30 min |

Runtime is dominated by loading the ~1 GB pickled database into memory; the k-mer lookup itself scales linearly with the number of peptides. Full-scale runs need a machine with enough RAM to hold the database (several times its on-disk size). Timings are indicative and depend on hardware.

---

## Troubleshooting

| Problem | Cause / solution |
| ------- | ---------------- |
| `FileNotFoundError: data/octamerDB.pkl` | The database is not in the repository (it is too large). Download it (see [Reference Dataset](#reference-dataset)) and place it in `data/`, or give its full path to `--database`. |
| `MemoryError` or the process is killed while loading the database | Not enough RAM to load the database. Use a machine with more memory, close other programs, or build a smaller organism-specific database with `build_database.py`. |
| `_pickle.UnpicklingError` / `EOFError` when loading | The download is incomplete or corrupted (e.g. Google Drive saved an HTML page). Download the file again and check its size (~1 GB). |
| `Total peptides processed: 0` | The k-mer file is empty or not in `peptide<TAB>k-mer` format. Check that Step 1 ran on a file with one peptide per line and that the peptides are ≥ k aa long. |
| Many peptides reported as excluded (< 8 aa) | Expected for short peptides; they cannot form an 8-mer. Use `--excluded_output` to inspect them. |
| Many unmapped peptides | Check that the organism and database match your sample, that `--organism` matches the database, that sequences are uppercase and contain no modification symbols (e.g. `M(ox)`), and that `--kmer_length` matches the database. |
| Many low-coverage mappings | Often contaminants or sequence variants; filter contaminants beforehand and/or set `--min_coverage`. |
| Commands with `\` fail on Windows | `\` line continuation is for Linux/macOS shells. On Windows, put the command on one line, or use `^` (Command Prompt) or `` ` `` (PowerShell) at line ends. |
| `python: command not found` | Use `python3` instead of `python`, or activate the environment created during installation. |

---

## Repository Structure
```

PEXMap
│
├── scripts
│ ├── generate_kmers.py               # Step 1: peptide -> k-mers
│ ├── annotate_peptides.py            # Step 2: k-mer annotation
│ ├── extract_sequences.py            # builds transcript sequences for Step 3
│ ├── verify_mapping.py               # Step 3: exact-substring verification
│ ├── build_database.py               # build your own k-mer database
│ ├── substring_search.py             # benchmark: exact substring search
│ ├── simulate_peptides.py            # benchmark: simulated peptides of known origin
│ └── npxp_split.py                   # curated (NP_) vs predicted (XP_) database composition
│
├── data                              # not in the repository; create it
│ └── octamerDB.pkl                   # downloaded reference database
│
├── example
│ ├── example_tandem_MSpeptides.txt
│ ├── example_kmers.txt
│ └── example_output.tsv
│
├── results
│ ├── k_sensitivity_summary.txt       # k = 6–12 sensitivity benchmark
│ └── npxp_results.tsv                # curated vs predicted accessions in octamerDB
│
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore

```
---

## Citation

If you use **PEXMap** in your research, please cite:

Awasthi D, *et al.* PEXMap: A proteogenomic method for exon and isoform level mapping of mass spectrometry derived peptides. *bioRxiv* (2026). https://doi.org/10.64898/2026.04.29.721330

PEXMap uses the ENACT exon nomenclature: Verma P, *et al.* Exon Nomenclature And Classification of Transcripts (ENACT) provides a systematic framework to annotate exon attributes. *Genome Research* 35(6):1440–1455 (2025). https://doi.org/10.1101/gr.279878.124

---

## License

PEXMap is released under the MIT License; see [LICENSE](LICENSE).

---

## Author

Deepanshi Awasthi, PhD Research Scholar, Computational Biology Group  
Indian Institute of Science Education and Research (IISER) Mohali, India


