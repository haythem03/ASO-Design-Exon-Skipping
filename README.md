# 🔬 Exon Skipping ASO Design Pipeline

This repository contains a Python pipeline to **design exon-skipping antisense oligonucleotides (ASOs)** for therapeutic purposes. It supports **mutation-aware frame disruption simulation**, **NMD prediction**, and **protein domain impact analysis**, based on **Ensembl GTF/FASTA** and **UniProt annotations**.

---

## 🧭 Overview

This tool is built to streamline the discovery of ASOs that promote exon skipping in disease-relevant transcripts — with the goal of restoring open reading frame (ORF), reducing harmful protein expression, or bypassing pathogenic mutations.

It performs:

- CDS extraction from reference GTF and FASTA
- Simulation of exon skipping
- Detection of frameshift and NMD
- Domain analysis using UniProt XML
- Design of internal (ASO-E) and junction-spanning (ASO-J) ASOs
- Integration of mutation data (VCF or TSV)
- Multi-exon skipping to restore reading frame
- Generation of exportable summary tables

## 📌 Features

- ✅ Parses Ensembl GTF and FASTA files
- ✅ Identifies canonical transcript, coding exons, and CDS structure
- ✅ Simulates exon skipping and evaluates:
  - Frame preservation/disruption
  - NMD potential
  - Protein domain loss
- ✅ Designs both internal (ASO-E) and junction-spanning (ASO-J) ASOs
- ✅ Supports mutation injection before simulation
- ✅ Prioritizes frame-rescuing exon combinations
- ✅ Detects loss of protein function using UniProt domain annotations
- ✅ Outputs comprehensive and structured results for downstream analysis

---


---

## 🛠️ Installation
Clone this repository and install dependencies:  

```bash
# Clone repository
git clone https://github.com/haythem03/ASO-Pipeline.git
cd ASO-Pipeline

# Install dependencies
pip install -r Requirements/requirements.txt

