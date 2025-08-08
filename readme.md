[![CBC Analysis CI](https://github.com/alikoraykoc/ITS2-CBC-Analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/alikoraykoc/ITS2-CBC-Analysis/actions/workflows/ci.yml) [![DOI](https://zenodo.org/badge/1034604068.svg)](https://doi.org/10.5281/zenodo.16782240)


# CBC Analysis for ITS2

This Python script identifies **Compensatory Base Changes (CBCs)** and **hemi-CBCs (hCBCs)** in the ITS2 secondary structure between two or more taxa, based on a **PAUP\*** "Character change list" output and an aligned XFasta file containing both sequence and secondary structure (dot-bracket notation).

## Features
- Supports **pairwise comparison** of two taxa or **batch mode** for multiple taxon pairs.
- Reports **aligned** and **ungapped** positions for each site and its paired base.
- Distinguishes between CBC, hCBC, unpaired, and gap positions.
- Outputs a tab-delimited table (`.tsv`) for easy downstream analysis.

## Input Requirements
1. **PAUP\*** change list file (e.g., `changes_output.txt`).
2. **Aligned XFasta file** with:
   - Line 1: `>Taxon_Name`
   - Line 2: Aligned nucleotide sequence (with `-` for gaps)
   - Line 3: Aligned secondary structure in dot-bracket notation
3. Names of the taxa to compare (`--seq1` and `--seq2`) or a file with all pairs (`--pairs all`).

## Usage

### Compare all taxon pairs:
```bash
python cbc_analysis.py \
  --changes changes_output.txt \
  --xfasta tetrigidae_aligned.xfasta \
  --pairs all \
  --out cbc_results.tsv

Compare a specific pair:

python cbc_analysis.py \
  --changes changes_output.txt \
  --xfasta tetrigidae_aligned.xfasta \
  --seq1 Tetrix_japonica \
  --seq2 tetrix_bolivari_ITS2 \
  --out cbc_results.tsv

Output

The TSV file includes:
	•	Aligned and ungapped positions for each base and its partner.
	•	Nucleotide identity and structure symbols.
	•	CBC/hCBC classification.
	•	PAUP change events and associated nodes.

Example output columns:

aligned_pos	seq1_base	seq1_struct	seq1_ungapped_pos	seq1_pair_aligned_pos	seq1_pair_ungapped_pos	seq1_pair_base	seq2_base	seq2_struct	seq2_ungapped_pos	seq2_pair_aligned_pos	seq2_pair_ungapped_pos	seq2_pair_base	type	change	nodes	events
