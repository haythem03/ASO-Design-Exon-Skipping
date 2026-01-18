<div align="center">

# 🔬 Exon Skipping ASO Design Pipeline

**A comprehensive framework for designing therapeutic antisense oligonucleotides.**

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)](LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg?style=for-the-badge)](https://github.com/haythem03/ASO-Pipeline/graphs/commit-activity)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-orange.svg?style=for-the-badge)](http://makeapullrequest.com)

[Overview](#-overview) •
[Features](#-key-features) •
[Installation](#-installation) •
[Methodology](#-methodology) •
[Output](#-outputs)

</div>

---

## 🧭 Overview

**ASO Design Pipeline** is a specialized bioinformatics tool built to streamline the discovery of Antisense Oligonucleotides (ASOs) for therapeutic applications. 

Inspired by the need for precision medicine, this tool automates the complex process of designing ASOs that promote **exon skipping** to restore Open Reading Frames (ORFs), bypass pathogenic mutations, or reduce harmful protein expression. By integrating genomic data (Ensembl) with proteomic insights (UniProt), it ensures that designs are not just genetically valid, but proteomically viable.

### Why use this pipeline?
* **Precision:** Designs ASOs specifically for internal exons (ASO-E) or splice junctions (ASO-J).
* **Insight:** Predicts the downstream impact on protein domains and NMD (Nonsense-Mediated Decay).
* **Flexibility:** Handles custom mutation data (VCF/TSV) to simulate patient-specific scenarios.

---

## 🚀 Key Features

This pipeline offers a robust suite of tools for ASO design and simulation:

| Feature | Description |
| :--- | :--- |
| **🧬 Advanced Parsing** | Seamlessly parses **Ensembl GTF** and **FASTA** files to identify canonical transcripts, coding exons, and CDS structures. |
| **⚡ Smart Simulation** | Simulates exon skipping events to evaluate **frame preservation**, **NMD potential**, and **protein stability**. |
| **🎯 Dual Design Mode** | Generates sequences for both **Internal (ASO-E)** and **Junction-Spanning (ASO-J)** oligonucleotides. |
| **🩺 Mutation-Aware** | Injects custom mutations (via VCF/TSV) before simulation to assess therapeutic rescue of specific variants. |
| **🛡️ Domain Analysis** | Integrates **UniProt XML** data to detect potential loss of critical protein functional domains. |
| **🔄 Multi-Exon Skip** | Capable of simulating multi-exon skipping to find the optimal combination for frame restoration. |

---

## 🛠️ Installation

Get the pipeline running on your local machine in minutes.

### Prerequisites
Ensure you have **Python 3.8+** installed.

### Step-by-Step Guide

1.  **Clone the repository**
    ```bash
    git clone [https://github.com/haythem03/ASO-Pipeline.git](https://github.com/haythem03/ASO-Pipeline.git)
    cd ASO-Pipeline
    ```

2.  **Create a Virtual Environment (Optional but Recommended)**
    ```bash
    python -m venv venv
    # Windows
    .\venv\Scripts\activate
    # Linux/Mac
    source venv/bin/activate
    ```

3.  **Install Dependencies**
    ```bash
    pip install -r Requirements/requirements.txt
    ```

---

## 🔬 Methodology

The pipeline follows a rigorous biological logic flow to ensure high-quality candidates:

1.  **Data Ingestion:** Reads reference genome (GTF/FASTA) and protein annotations (UniProt).
2.  **Variant Injection:** (Optional) Applies user-defined mutations to the sequence.
3.  **Skipping Simulation:** Systematically skips exons and recalculates the reading frame.
4.  **Impact Assessment:** * Checks for Premature Termination Codons (PTC).
    * Predicts NMD susceptibility (50-55nt rule).
    * Maps missing regions to UniProt domains.
5.  **Candidate Generation:** Outputs optimized ASO sequences for the valid targets.

---

## 📊 Outputs

The tool generates comprehensive exportable summary tables (CSV/Excel) containing:

* **Transcript Metadata:** Gene ID, Transcript ID, Exon counts.
* **Skipping Metrics:** Frame check (Preserved/Disrupted), NMD prediction.
* **Protein Impact:** Percent identity retained, specific domains lost/retained.
* **ASO Sequences:** * *ASO-E:* Target sequences within the exon.
    * *ASO-J:* Target sequences spanning exon-exon junctions.

---

## 🤝 Contributing

Contributions make the open-source community an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1.  Fork the Project
2.  Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3.  Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4.  Push to the Branch (`git push origin feature/AmazingFeature`)
5.  Open a Pull Request

---

<div align="center">

**[Report Bug](https://github.com/haythem03/ASO-Pipeline/issues) • [Request Feature](https://github.com/haythem03/ASO-Pipeline/issues)**

</div>
