# PEXMap – Peptide Exon Mapping Tool

PEXMap (**P**eptide **EX**on **Map**per) maps peptides that have already been identified by MS/MS to the **genes, transcript isoforms, exons and exon–exon junctions** they come from.

Conventional tools usually stop at the gene or protein level. PEXMap adds an exon- and junction-level layer, using the ENACT exon nomenclature (exon unique identifiers, EUIDs), so you can see which exons and splice junctions each peptide supports and whether that evidence is shared across isoforms or specific to some of them.

**How it works, in one paragraph:** every protein-coding isoform in the reference is cut into overlapping 8-amino-acid fragments (8-mers), both within exons and across exon–exon junctions, and each 8-mer is indexed with its gene, transcript and exon/junction. Your peptides are cut into 8-mers in the same way and looked up in this index. Each peptide is assigned to the gene, transcripts and exon/junction with the most distinct matching 8-mers. A final check then confirms that the whole peptide really occurs, uninterrupted, in the assigned isoform.

> **Note on file names (v1.1):** the scripts were renamed to shorter names in v1.1. If you are following an older description: `annotate_tandemMS_peptides.py` → `annotate_peptides.py`, `build_kmer_database.py` → `build_database.py`, `extract_transcript_sequences.py` → `extract_sequences.py`, `verify_exact_substring.py` → `verify_mapping.py`, `exact_substring_baseline.py` → `substring_search.py`, `generate_synthetic_peptides.py` → `simulate_peptides.py`. The code is otherwise unchanged by the renaming.

---

## Contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [Reference database](#reference-database)
- [Input format](#input-format)
- [Running PEXMap](#running-pexmap)
- [Output](#output)
- [Example run with expected results](#example-run-with-expected-results)
- [Building your own database](#building-your-own-database)
- [Benchmark and analysis scripts](#benchmark-and-analysis-scripts)
- [Limitations](#limitations)
- [Runtime and resources](#runtime-and-resources)
- [Troubleshooting](#troubleshooting)
- [Repository structure](#repository-structure)
- [Citation](#citation) · [License](#license) · [Author](#author)

---

## Quick start

```
git clone https://github.com/deepanshicbg/PEXMap.git
cd PEXMap
mkdir data        # put the downloaded octamerDB.pkl here (see "Reference database")

# 1. peptides -> 8-mers
python scripts/generate_kmers.py --input example/example_tandem_MSpeptides.txt --output my_kmers.txt

# 2. 8-mers -> gene / transcript / exon / junction annotation
python scripts/annotate_peptides.py --kmers my_kmers.txt --database data/octamerDB.pkl --organism human --peptides example/example_tandem_MSpeptides.txt --output my_annotations.tsv
```

Then check the results against the expected numbers in [Example run with expected results](#example-run-with-expected-results).

---

## Installation

**Requirements**

- **Python ≥ 3.8**, accessible from the command line.
- **No third-party Python packages.** All scripts use only the Python standard library; `requirements.txt` documents this.
- **Operating system:** Linux, macOS or Windows (on Windows, see [Troubleshooting](#troubleshooting) for multi-line commands).
- **Memory:** the reference database (~1 GB on disk) is loaded fully into memory, so full-scale runs need a machine with ample RAM (see [Runtime and resources](#runtime-and-resources)).

**Steps**

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
3. Install dependencies (there are none beyond the standard library; this simply confirms the environment):
   ```
   pip install -r requirements.txt
   ```
4. Check the installation:
   ```
   python scripts/generate_kmers.py --help
   ```

---

## Reference database

The reference database (`octamerDB.pkl`) is too large for GitHub and is hosted externally. It contains the 8-mers of human protein-coding isoforms from both exons and exon–exon junctions.

| Version | Use | Link |
| ------- | --- | ---- |
| ENACT v0.5 | used in the manuscript | https://drive.google.com/uc?id=1jPU8HE6Fcwk4mAU8Fk5m7VJGLtnrrKKF |
| Latest | recommended for new analyses | https://drive.google.com/file/d/124p7-jL0uxcfmfKGvxBqu2JAOXGH2GqM/view?usp=sharing |

After downloading, create a `data` folder inside PEXMap and place the file there:

```
PEXMap/data/octamerDB.pkl
```

If you downloaded `octamerDB_latest.pkl`, either rename it to `octamerDB.pkl` or give its actual name to `--database`.

> `.pkl` files are Python pickle objects, which can run code when loaded. Only load database files from the links above or ones you built yourself with `scripts/build_database.py`.

---

## Input format

PEXMap does **not** identify peptides from spectra. Its input is a list of peptides that have **already been identified and FDR-controlled** by your search engine or resource (e.g. PeptideAtlas).

- One peptide per line, uppercase one-letter amino-acid code.
- No modification symbols (e.g. `M(ox)`), charges or flanking residues.
- Remove common contaminants (keratins, trypsin, BSA, e.g. the cRAP list) **before** running PEXMap; see [Contaminants](#contaminants-and-low-coverage-mappings).
- Peptides shorter than 8 residues cannot form an 8-mer; they are counted and reported, not silently dropped.

Example (`input_peptides.txt`):

```
MTEYKLVVVGAG
ADLASRDE
VAVWPTMV
```

The length filter is the only filter applied by these scripts. Any other filtering (for example removal of low-complexity peptides) should be done before input.

---

## Running PEXMap

PEXMap runs in three steps. Steps 1 and 2 produce the annotation; Step 3 is a recommended check.

### Step 1 — Generate k-mers

```
python scripts/generate_kmers.py --input input_peptides.txt --output kmers.txt --kmer_length 8 --excluded_output excluded_short_peptides.txt
```

Splits each peptide into overlapping k-mers (default k = 8, which must match the database). Duplicate peptides are processed once. Prints a summary such as:

```
Input peptides: 3725
Excluded because shorter than 8 aa: 336
Processed unique peptides (k-mers generated): 3389
Total unique k-mers written: ...
```

### Step 2 — Annotate peptides

```
python scripts/annotate_peptides.py \
  --kmers kmers.txt \
  --database data/octamerDB.pkl \
  --organism human \
  --peptides input_peptides.txt \
  --min_coverage 0 \
  --output PEXMap_annotations.tsv
```

Writes:

| File | Contents |
| ---- | -------- |
| `PEXMap_annotations.tsv` | the annotation table (see [Output](#output)) |
| `PEXMap_annotations.tsv.unmapped.txt` | peptides with no match in the database |
| `PEXMap_annotations.tsv.below_coverage.txt` | peptides below `--min_coverage` (only when a threshold > 0 is set) |

### Step 3 — Verify each assignment (recommended)

Matching all 8-mers of a peptide does **not** guarantee that the whole peptide occurs, uninterrupted, in the assigned isoform: the 8-mers can be found in the right gene but in a different order or position. Step 3 checks each assignment by looking for the complete peptide sequence in the protein sequence of the assigned transcript(s). On the benchmark, ~0.9% of assignments failed this check at every k tested (`results/k_sensitivity_summary.txt`).

First build the transcript sequence reference once, from the same ENACT per-gene annotation files used to build the database:

```
python scripts/extract_sequences.py --input_folder organism_gene_files --output data/transcript_sequences.pkl
```

Then verify:

```
python scripts/verify_mapping.py \
  --annotations PEXMap_annotations.tsv \
  --sequences data/transcript_sequences.pkl \
  --output PEXMap_annotations_verified.tsv \
  --drop_failed
```

This adds an `Exact_match` column: `PASS`, `FAIL`, `NA` (multi_gene or no transcript call) or `No_reference_sequence`. With `--drop_failed`, `FAIL` rows are removed from the output; without it they are kept and flagged. We recommend using `--drop_failed` (or discarding `FAIL` rows) for any downstream analysis.

### Arguments

**`generate_kmers.py`**

| Argument | Description |
| -------- | ----------- |
| `--input` | input peptide file (one sequence per line) |
| `--output` | output k-mer file (`peptide<TAB>k-mer`) |
| `--kmer_length` | k-mer length (default 8; must match the database) |
| `--excluded_output` | optional file listing peptides shorter than `--kmer_length` |

**`annotate_peptides.py`**

| Argument | Description |
| -------- | ----------- |
| `--kmers` | k-mer file from Step 1 |
| `--database` | reference database (`octamerDB.pkl` or your own) |
| `--output` | output annotation table (TSV) |
| `--organism` | optional organism filter (e.g. `human`); database entries from other organisms are ignored |
| `--peptides` | optional original input file, so the summary also reports input and excluded counts |
| `--kmer_length` | k used in Step 1 (only used with `--peptides`; default 8) |
| `--min_coverage` | minimum k-mer coverage (%) for a peptide to be written to the output (default 0 = keep all) |
| `--unmapped_output` | file for unmapped peptides (default `<output>.unmapped.txt`) |

**`verify_mapping.py`**

| Argument | Description |
| -------- | ----------- |
| `--annotations` | output table from Step 2 |
| `--sequences` | transcript sequence file from `extract_sequences.py` |
| `--output` | output table with the added `Exact_match` column |
| `--drop_failed` | remove `FAIL` rows from the output instead of flagging them |

---

## Output

### Columns of the annotation table

| Column | Description |
| ------ | ----------- |
| `Experimental_MS_peptide` | input peptide |
| `Gene_id` | gene with the most distinct supporting k-mers; for `multi_gene`, all tied gene IDs separated by `;` |
| `Feature_type` | `exon`, `junction` (exon–exon junction), `exon+junction` (best exon and best junction equally supported), `multi_gene` (tie between genes, so no exon/transcript call is made) or `NA` (no exon/junction evidence) |
| `Features` | EUID(s) of the supporting exon(s), or junction(s) written as `EUID1,EUID2`; several features separated by `;` |
| `Transcripts` | transcript(s) of the selected gene with the most distinct supporting k-mers; separated by `;` |
| `Kmer_hits` | number of distinct k-mers supporting the selected (or tied) gene |
| `Total_unique_kmers` | number of unique k-mers in the peptide |
| `Matched_kmers` | number of those k-mers found in the database |
| `Coverage_percent` | `Matched_kmers / Total_unique_kmers × 100` |

### Example rows (from `example/example_output.tsv`)

| Experimental_MS_peptide | Gene_id | Feature_type | Features | Transcripts | Kmer_hits | Total_unique_kmers | Matched_kmers | Coverage_percent |
|---|---|---|---|---|---|---|---|---|
| AGSYGAQPVVQTQLNSYGAQA | 10432 | exon | T.1.A.2.0.0 | NP_006319.1 | 14 | 14 | 14 | 100.0 |
| CLGVHCDFHQGLLSFYNAR | 79187 | junction | T.1.F.11.0.0,T.1.G.12.0.0 | NP_077309.1 | 12 | 12 | 12 | 100.0 |
| INNPMFRFLFENQTPAHVYYR | 23350 | exon+junction | T.1.G.16.0.0;T.1.G.15.0.0,T.1.G.16.0.0;T.2.G.15.0.0,T.1.G.16.0.0 | NP_001073884.1;NP_001307148.1;NP_001307149.1;NP_001307151.1 | 14 | 14 | 14 | 100.0 |
| IGKAKTKENRQSIINPDWNFEKM | 474170;4905 | multi_gene | | | 16 | 16 | 16 | 100.0 |
| GQSEADSDKNATILELR | 1832 | exon | T.1.F.23.0.0 | NP_004406.2 | 8 | 10 | 8 | 80.0 |

### How to read the results

- **Coverage** tells you how many of the peptide's 8-mers are in the database. **100% coverage alone is not proof that the peptide occurs in the assigned isoform**: all 8-mers can be present without the full peptide being present uninterrupted. Treat an assignment as confirmed only after it passes Step 3 (`Exact_match = PASS`).
- **Coverage below 100%** (e.g. *GQSEADSDKNATILELR*, 80%) means some 8-mers are missing from the database. Possible reasons: the reference lacks that sequence, a sequence variant, or a contaminant. Use `--min_coverage` to exclude such peptides.
- **`junction`** (e.g. *CLGVHCDFHQGLLSFYNAR*): the peptide spans an exon–exon boundary, which can distinguish isoforms that share exons but join them differently. A junction call supported by only one or two k-mers rests on very few residues crossing the boundary; check `Kmer_hits` and Step 3 before relying on it.
- **`exon+junction`** (e.g. *INNPMFRFLFENQTPAHVYYR*): the best exon and the best junction are equally supported, so both are reported instead of choosing one arbitrarily.
- **`multi_gene`** (e.g. *IGKAKTKENRQSIINPDWNFEKM*): two or more genes are equally supported, so all are listed and no exon or transcript call is made.
- An exon-mapped peptide may be **specific to one transcript** or **shared by several isoforms** of the gene; the `Transcripts` column shows which.

### Summary printed at the end of Step 2

- **Input / QC summary:** input peptides and those excluded as too short (with `--peptides`), peptides processed, and how many were **annotated**, **below the coverage threshold** or **unmapped**. These always add up to the number processed, so no peptide disappears without being reported.
- **Coverage summary:** peptides with 100%, ≥ 80%, ≥ 50%, ≥ 30% and < 30% coverage.
- **Feature summary:** annotated peptides mapped to exons, to junctions, and tied between genes (`exon+junction` ties are counted under both exons and junctions).

### Contaminants and low-coverage mappings

Contaminant peptides (e.g. keratins, trypsin, BSA) can partly match human proteins and appear as low-coverage mappings, i.e. potential false positives. Remove contaminants before running PEXMap, and use `--min_coverage` and Step 3 to control the rest.

**Choosing `--min_coverage`:** 100 gives the most conservative annotation; 80 also keeps peptides with a few missing 8-mers. The default (0) keeps everything so that you can apply your own cut-off later.

---

## Example run with expected results

The `example/` folder contains:

| File | Description |
| ---- | ----------- |
| `example_tandem_MSpeptides.txt` | 1,000 human MS/MS peptides (input) |
| `example_kmers.txt` | the k-mers Step 1 produces from them (for comparison) |
| `example_output.tsv` | the expected Step 2 output |

**Step 1**

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

**Step 2**

```
python scripts/annotate_peptides.py --kmers example/my_kmers.txt --database data/octamerDB.pkl --organism human --peptides example/example_tandem_MSpeptides.txt --output example/my_annotations.tsv
```

Expected summary (compare `example/my_annotations.tsv` with `example/example_output.tsv`):

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

Feature types in the output: 869 `exon`, 113 `junction`, 15 `exon+junction`, 3 `multi_gene`. Small differences are expected between database versions.

---

## Building your own database

If you have ENACT-format transcript–exon annotation files for an organism (one file per gene, listing each transcript's exons with their EUIDs and amino-acid sequences), you can build your own database:

```
python scripts/build_database.py --input_folder organism_gene_files --kmer 8 --organism human --output kmerDB.pkl
```

The database indexes every k-mer from within exons and across consecutive exon–exon junctions, each with its gene ID, transcript IDs and exon or junction ID. Use it with `--database kmerDB.pkl` in Step 2, and build the matching sequence file for Step 3 with `extract_sequences.py` from the same folder.

---

## Benchmark and analysis scripts

These scripts reproduce the benchmarks reported in the manuscript. They are not needed for routine use.

| Script | Purpose |
| ------ | ------- |
| `substring_search.py` | maps peptides by direct exact substring search of every transcript sequence (no k-mers), as a reference method for accuracy and runtime comparison |
| `simulate_peptides.py` | creates synthetic tryptic peptides of known origin from the reference sequences, to measure gene- and transcript-level precision |
| `npxp_split.py` | reports what fraction of the database comes from curated RefSeq entries (`NP_`) versus computationally predicted models (`XP_`), and how many peptides involve predicted models |

Results included in `results/`:

- `k_sensitivity_summary.txt`: database size, coverage, feature calls and Step 3 failure rate for k = 6–12.
- `npxp_results.tsv`: curated vs predicted accessions in the database and their involvement in the example peptides.

---

## Limitations

Please keep these in mind when interpreting PEXMap results:

- **Annotated features only.** PEXMap can only report exons and junctions that exist in its reference. It does not discover novel exons, junctions or variant peptides; a peptide from an unannotated junction will be unmapped or mapped to the closest annotated feature.
- **Predicted models in the reference.** The reference includes computationally predicted RefSeq models (`XP_`) as well as curated ones (`NP_`). A peptide specific to one curated isoform can appear shared because a predicted model contains the same sequence (see `npxp_split.py`).
- **Isoleucine and leucine are treated as different amino acids**, although mass spectrometry generally cannot tell them apart. Claims that a peptide is unique to one exon or isoform should take this into account.
- **No modifications.** Post-translational modifications and N-terminal methionine removal are not modelled; input peptides must be plain sequences.
- **Very short exons.** Junction k-mers are built only between consecutive exons. Where an internal exon is shorter than k − 1 residues (≤ 6 aa for k = 8), k-mers spanning three exons are not built.
- **Use Step 3.** k-mer matching alone produces a small fraction (~0.9% on the benchmark) of assignments where the full peptide is not actually present in the assigned isoform.

---

## Runtime and resources

| Task | Approximate runtime |
| ---- | ------------------- |
| Example dataset (1,000 peptides), Steps 1–3 | ~2–3 min, mostly spent loading the database |
| Complete human PeptideAtlas build (Step 2) | ~30 min |

Runtime is dominated by loading the ~1 GB database into memory; the lookup itself scales linearly with the number of peptides. Full-scale runs need enough RAM to hold the database (several times its size on disk). Timings depend on hardware.

---

## Troubleshooting

| Problem | Cause / solution |
| ------- | ---------------- |
| `FileNotFoundError: data/octamerDB.pkl` | The database is not part of the repository. Download it (see [Reference database](#reference-database)) into `data/`, or give its full path to `--database`. |
| `MemoryError`, or the process is killed while loading the database | Not enough RAM. Use a machine with more memory, close other programs, or build a smaller organism-specific database with `build_database.py`. |
| `_pickle.UnpicklingError` / `EOFError` when loading | The download is incomplete or corrupted (Google Drive sometimes saves an HTML page instead). Download again and check the size (~1 GB). |
| `Total peptides processed: 0` | The k-mer file is empty or not in `peptide<TAB>k-mer` format. Check that Step 1 ran on a file with one peptide per line. |
| Many peptides excluded as shorter than 8 aa | Expected: they cannot form an 8-mer. Use `--excluded_output` to see them. |
| Many unmapped peptides | Check that the database matches your organism and `--organism`, that sequences are uppercase without modification symbols, and that `--kmer_length` matches the database. |
| Many low-coverage mappings | Often contaminants or sequence variants; remove contaminants beforehand and/or set `--min_coverage`. |
| Multi-line commands fail on Windows | The `\` line breaks are for Linux/macOS. On Windows, write the command on one line, or end lines with `^` (Command Prompt) or `` ` `` (PowerShell). |
| `python: command not found` | Use `python3`, or activate the environment created during installation. |

---

## Repository structure

```
PEXMap
├── scripts
│   ├── generate_kmers.py        # Step 1: peptides -> k-mers
│   ├── annotate_peptides.py     # Step 2: k-mer annotation and QC summary
│   ├── extract_sequences.py     # builds transcript sequences for Step 3
│   ├── verify_mapping.py        # Step 3: exact-substring verification
│   ├── build_database.py        # build your own k-mer database
│   ├── substring_search.py      # benchmark: direct exact substring search
│   ├── simulate_peptides.py     # benchmark: synthetic peptides of known origin
│   └── npxp_split.py            # curated (NP_) vs predicted (XP_) composition
├── example
│   ├── example_tandem_MSpeptides.txt
│   ├── example_kmers.txt
│   └── example_output.tsv
├── results
│   ├── k_sensitivity_summary.txt
│   └── npxp_results.tsv
├── data                         # not in the repository; put octamerDB.pkl here
├── README.md
├── requirements.txt
├── LICENSE
└── .gitignore
```

---

## Citation

If you use PEXMap, please cite:

Awasthi D, *et al.* PEXMap: A proteogenomic method for exon and isoform level mapping of mass spectrometry derived peptides. *bioRxiv* (2026). https://doi.org/10.64898/2026.04.29.721330

PEXMap uses the ENACT exon nomenclature: Verma P, *et al.* Exon Nomenclature And Classification of Transcripts (ENACT) provides a systematic framework to annotate exon attributes. *Genome Research* 35(6):1440–1455 (2025). https://doi.org/10.1101/gr.279878.124

## License

PEXMap is released under the MIT License; see [LICENSE](LICENSE).

## Author

Deepanshi Awasthi, Computational Biology Group, Indian Institute of Science Education and Research (IISER) Mohali, India
