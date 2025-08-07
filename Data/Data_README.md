# Mutation Data: Raw vs. Processed

This folder contains two key files related to mutation data for the gene **ENSG00000198947** from the GnomAD database, along with a processed version that we generated to facilitate downstream analysis.

---

## Files

### 1. Raw GnomAD Data  
**Filename:** `gnomAD_v4.1.0_ENSG00000198947_2025_04_23_11_02_58.csv`

- This is the original mutation dataset downloaded from GnomAD version 4.1.0.
- It contains comprehensive variant information for the specified gene.
- Includes many columns with detailed annotations, allele frequencies, quality metrics, and other metadata.
- The raw format is large and contains much information that may not be necessary for specific analysis pipelines.

### 2. Processed Mutation Data  
**Filename:** `df_final_all (1).csv`

- This file is a curated subset derived from the raw GnomAD data.
- It contains only the essential columns needed for our analyses:
  - `gene`
  - `Chromosome`
  - `Position`
  - `ref` (reference allele)
  - `alt` (alternate allele)
- The preprocessing involved filtering and selecting relevant columns to reduce complexity and file size.
- This streamlined dataset is easier to integrate with downstream bioinformatics pipelines and reduces memory overhead.

---
