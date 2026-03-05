# Peptide Annotation Tool

This repository provides a pipeline to annotate experimental peptides
against a reference exon / exon-junction peptide database.

## Workflow

1. Generate 8-mer peptides
2. Search against reference database
3. Retrieve gene and transcript annotations

## Installation

git clone https://github.com/username/peptide-annotation-tool
cd peptide-annotation-tool

pip install -r requirements.txt

## Usage

Generate 8-mers:

python scripts/generate_kmers.py input_peptides.txt kmers.txt

Annotate peptides:

python scripts/annotate_peptides.py kmers.txt data/peptide_dataset.pkl annotations.tsv
