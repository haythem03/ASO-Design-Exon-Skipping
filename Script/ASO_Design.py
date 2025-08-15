#!/usr/bin/env python3
"""
Pipeline for analyzing exon skipping, ASO design, and protein domain effects in the DMD gene.
Downloads GTF and FASTA files, processes mutation data, simulates exon skipping, predicts NMD,
designs ASOs, and analyzes BLAST results.
"""

import pandas as pd
import numpy as np
from pyfaidx import Fasta
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
import re
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import subprocess
import os

# ---- Helper functions ----

def extract_transcript_id(attr):
    match = re.search(r'transcript_id "([^"]+)"', attr)
    return match.group(1) if match else None

def extract_exon_number(attr):
    match = re.search(r'exon_number "(\d+)"', attr)
    return int(match.group(1)) if match else None

def get_cds_sequence(gtf_cds_blocks, genome):
    seq_parts = []
    strand = gtf_cds_blocks.iloc[0]["strand"]
    chrom = gtf_cds_blocks.iloc[0]["seqname"].replace("chr", "")
    for _, row in gtf_cds_blocks.iterrows():
        start, end = row["start"], row["end"]
        seq = genome[chrom][start - 1:end].seq.upper()
        seq_parts.append(str(seq))
    full_seq = ''.join(seq_parts)
    if strand == '-':
        full_seq = str(Seq(full_seq).reverse_complement())
    return full_seq

def replace_stop_codons_in_cds(cds_seq):
    codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq), 3)]
    modified_codons = ["***" if codon in {"TAA", "TAG", "TGA"} else codon for codon in codons]
    return "".join(modified_codons)

def translate_cds_with_star(cds_seq):
    cds_seq_obj = Seq(cds_seq)
    protein_seq = cds_seq_obj.translate(to_stop=False)
    return str(protein_seq)

def detect_frame_status(original_seq, skipped_seq, original_protein, skipped_protein):
    is_multiple_of_three = len(skipped_seq) % 3 == 0
    if skipped_protein == original_protein:
        return "in-frame (no change)"
    if "*" in skipped_protein:
        if skipped_protein.count("*") == 1 and skipped_protein[-1] == "*":
            if len(skipped_protein) == len(original_protein):
                return "in-frame (same length, altered amino acids)"
            elif len(skipped_protein) < len(original_protein):
                return "in-frame (truncated, ends with stop)"
            else:
                return "in-frame (elongated, ends with stop)"
        else:
            return "frameshift (premature stop)"
    if len(skipped_protein) != len(original_protein):
        if is_multiple_of_three:
            return "in-frame (altered length, no stop)"
        else:
            return "frameshift (length mismatch, no stop)"
    if skipped_protein != original_protein and len(skipped_protein) == len(original_protein):
        if is_multiple_of_three:
            return "in-frame (altered amino acids)"
        else:
            return "frameshift (same length, altered amino acids)"
    if not is_multiple_of_three:
        return "frameshift (likely)"
    return "uncertain"

def simulate_exon_skip_and_translate(cds_sequence: str, exon_cds_start: int, exon_cds_end: int):
    full_seq = cds_sequence.upper()
    skipped_seq = full_seq[:exon_cds_start] + full_seq[exon_cds_end + 1:]
    valid_bases = set('ATCG')
    if not set(full_seq).issubset(valid_bases):
        print("❌ Invalid bases in full sequence")
        return None
    if not set(skipped_seq).issubset(valid_bases):
        print("❌ Invalid bases in skipped sequence")
        return None
    try:
        full_protein = str(Seq(full_seq).translate(to_stop=False))
        skipped_protein = str(Seq(skipped_seq).translate(to_stop=False))
    except Exception as e:
        print(f"❌ Translation error: {e}")
        return None
    if skipped_protein == full_protein:
        disruption_type = "in-frame (no change)"
    elif "*" in skipped_protein and skipped_protein[-1] == "*" and skipped_protein.count("*") == 1:
        disruption_type = "in-frame (ends with stop)"
    elif "*" in skipped_protein:
         disruption_type = "frameshift (premature stop)"
    elif len(skipped_seq) % 3 == 0:
        disruption_type = "in-frame (altered protein)"
    else:
        disruption_type = "frameshift (likely)"
    return {
        "original_len": len(full_seq),
        "skipped_len": len(skipped_seq),
        "original_seq": full_seq,
        "skipped_seq": skipped_seq,
        "original_protein": full_protein,
        "skipped_protein": skipped_protein,
        "disruption_type": disruption_type
    }

def predict_nmd(frame_status, exon_start, exon_end, cds_sequence, genomic_to_cds_index, chrom, strand, genome):
    if frame_status not in {"frameshift (premature stop)", "frameshift (likely)"}:
        return "NA"
    try:
        protein = Seq(cds_sequence.upper()).translate(to_stop=False)
    except:
        return "Translation error"
    stop_index = protein.find("*")
    if stop_index == -1:
        return "No premature stop"
    ptc_cds_pos = stop_index * 3
    reverse_map = {v: k for k, v in genomic_to_cds_index.items()}
    if ptc_cds_pos in reverse_map:
        ptc_genomic_pos = reverse_map[ptc_cds_pos]
    else:
        return "Could not map PTC to genome"
    if ptc_genomic_pos < exon_end - 55:
        return "NMD likely"
    elif ptc_genomic_pos >= exon_end - 55:
        return "NMD escape (50–55 rule)"
    elif ptc_cds_pos < 150:
        return "NMD escape (start-proximal)"
    else:
        return "Uncertain"

def sliding_window(s, k=25):
    return [s[i:i + k] for i in range(len(s) - k + 1)]

def filter_aso(s):
    if len(s) != 25:
        return False
    if "GGGG" in s:
        return False
    gc = gc_fraction(s) * 100
    g_pct = s.count("G") / 25 * 100
    return 40 <= gc <= 60 and g_pct <= 36

def custom_tm_formula(seq):
    seq = seq.upper()
    wA = seq.count('A')
    xT = seq.count('U')
    yG = seq.count('G')
    zC = seq.count('C')
    numerator = (yG + zC) - 16.4
    denominator = (wA + xT + yG + zC)
    if denominator == 0:
        return None
    tm = 64.9 + 41 * (numerator / denominator)
    return tm


# ---- Main execution ----

def main():
    # Load files (make sure to have downloaded the files beforehand)
    gtf_path = "Homo_sapiens.GRCh38.104.gtf"
    genome_fasta_path = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    mutation_file = "df_final_all.csv"

    gtf_cols = [
        "seqname", "source", "feature", "start", "end", "score",
        "strand", "frame", "attribute"
    ]

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        names=gtf_cols,
    )

    genome = Fasta(genome_fasta_path)
    mutation_df = pd.read_csv(mutation_file)
    mutation_df["Chromosome"] = mutation_df["Chromosome"].astype(str).str.replace("chr", "").str.strip()

    canonical_tx = "ENST00000357033"

    gtf_cds = gtf[
        (gtf["feature"].isin(["CDS", "stop_codon"])) &
        gtf["attribute"].astype(str).str.contains(f'transcript_id "{canonical_tx}"') &
        gtf["attribute"].astype(str).str.contains('transcript_biotype "protein_coding"')
    ].copy()
    gtf_cds["transcript_id"] = gtf_cds["attribute"].apply(extract_transcript_id)

    gtf_exons = gtf[
        (gtf["feature"] == "exon") &
        gtf["attribute"].astype(str).str.contains(f'transcript_id "{canonical_tx}"') &
        gtf["attribute"].astype(str).str.contains('transcript_biotype "protein_coding"')
    ].copy()
    gtf_exons["transcript_id"] = gtf_exons["attribute"].apply(extract_transcript_id)
    gtf_exons["exon_number"] = gtf_exons["attribute"].apply(extract_exon_number)

    tx_exons = gtf_exons[gtf_exons["transcript_id"] == canonical_tx].sort_values("start")
    tx_cds = gtf_cds[gtf_cds["transcript_id"] == canonical_tx].sort_values("start")

    if len(tx_exons) <= 2 or len(tx_cds) == 0:
        print(f"❌ Not enough exons or CDS entries for {canonical_tx}")
        return

    strand = tx_cds.iloc[0]["strand"]
    chrom = tx_cds.iloc[0]["seqname"].replace("chr", "")

    if strand == "+":
        cds_blocks_sorted = tx_cds.sort_values("start")
    else:
        cds_blocks_sorted = tx_cds.sort_values("start", ascending=False)

    cds_seq_parts = []
    cds_coords = []

    for _, row in cds_blocks_sorted.iterrows():
        start, end = row["start"], row["end"]
        seq = genome[chrom][start - 1:end].seq.upper()
        cds_seq_parts.append(str(seq))
        cds_coords.append((start, end))

    full_cds_sequence = ''.join(cds_seq_parts)
    if strand == '-':
        full_cds_sequence = str(Seq(full_cds_sequence).reverse_complement())

    print(f"✅ CDS sequence start (first 30 bases): {full_cds_sequence[:30]}")
    print(f"✅ Translated protein (first 20 aa): {str(Seq(full_cds_sequence).translate(to_stop=False)[:20])}")

    genomic_to_cds_index = {}
    cds_index = 0
    for start, end in cds_coords:
        rng = range(start, end + 1) if strand == "+" else range(end, start - 1, -1)
        for pos in rng:
            genomic_to_cds_index[pos] = cds_index
            cds_index += 1
coding_exons = []
    for _, exon in tx_exons.iterrows():
        exon_start, exon_end = exon["start"], exon["end"]
        overlap = tx_cds[
            (tx_cds["start"] <= exon_end) & (tx_cds["end"] >= exon_start)
        ]
        if not overlap.empty:
            coding_exons.append(exon)

    results = []

    for idx, exon_row in enumerate(coding_exons[1:-1], start=1):
        exon_start = exon_row["start"]
        exon_end = exon_row["end"]
        original_exon_number = exon_row.get("exon_number", None)

        exon_cds_positions = [
            genomic_to_cds_index[pos] for pos in range(exon_start, exon_end + 1)
            if pos in genomic_to_cds_index
        ]
        if not exon_cds_positions:
            print(f"❌ Exon {exon_start}-{exon_end} has no CDS positions, skipping")
            continue

        exon_cds_start = min(exon_cds_positions)
        exon_cds_end = max(exon_cds_positions)

        sim_result = simulate_exon_skip_and_translate(
            cds_sequence=full_cds_sequence,
            exon_cds_start=exon_cds_start,
            exon_cds_end=exon_cds_end
        )

        if sim_result is None:
            print(f"❌ Simulation failed for exon {exon_start}-{exon_end}")
            continue

        tx_id = canonical_tx
        frame = sim_result["disruption_type"]

        print(f"Transcript {tx_id} | Exon {exon_start}-{exon_end} | Frame: {frame} | Original length: {sim_result['original_len']} | Skipped length: {sim_result['skipped_len']}")

        nmd = predict_nmd(
            frame_status=frame,
            exon_start=exon_start,
            exon_end=exon_end,
            cds_sequence=sim_result["skipped_seq"],
            genomic_to_cds_index=genomic_to_cds_index,
            chrom=chrom,
            strand=strand,
            genome=genome
        )
muts = mutation_df[
            (mutation_df["gene"] == tx_id) &
            (mutation_df["Chromosome"].astype(str) == chrom) &
            (mutation_df["Position"] >= exon_start) &
            (mutation_df["Position"] <= exon_end)
        ]
        mut_count = len(muts)

        try:
            seq = genome[chrom][exon_start - 1:exon_end].seq.upper()
            if strand == "-":
                seq = str(Seq(seq).reverse_complement())
        except:
            continue

        raw_e = sliding_window(seq, 25)
        valid_e = [str(Seq(a).reverse_complement()).replace("T", "U") for a in raw_e if filter_aso(a)]

        junction_seqs = []
        try:
            if strand == "+":
                intron_5p = genome[chrom][exon_end:exon_end + 24].seq.upper()
                intron_3p = genome[chrom][exon_start - 25:exon_start].seq.upper()
                exon_seq_start = seq[:7]
                exon_seq_end = seq[-7:]
            else:
                intron_5p = str(Seq(genome[chrom][exon_start - 25:exon_start].seq.upper()).reverse_complement())
                intron_3p = str(Seq(genome[chrom][exon_end:exon_end + 24].seq.upper()).reverse_complement())
                exon_seq_start = seq[-7:]
                exon_seq_end = seq[:7]

            for i in range(7):
                a_seq = intron_3p + exon_seq_start
                d_seq = exon_seq_end + intron_5p
                if len(a_seq[i:i + 25]) == 25:
                    junction_seqs.append(a_seq[i:i + 25])
                if len(d_seq[i:i + 25]) == 25:
                    junction_seqs.append(d_seq[i:i + 25])
        except Exception as e:
            print(f"❌ Junction ASO error at Exon {exon_start}-{exon_end}: {e}")
            junction_seqs = []

        valid_j = [str(Seq(j).reverse_complement()).replace("T", "U") for j in junction_seqs if filter_aso(j)]

results.append({
            "Transcript": tx_id,
            "Chromosome": chrom,
            "Exon Number": original_exon_number,
            "Exon Seq": str(seq),
            "Exon length": len(seq),
            "Exon Start": exon_start,
            "Exon End": exon_end,
            "Mutation Count": mut_count,
            "Frame": frame,
            "NMD": nmd,
            "ASO-E Count": len(valid_e),
            "ASO-E seq": valid_e,
            "ASO-J seq": valid_j,
            "ASO-J Count": len(valid_j),
            "Original Protein": sim_result["original_protein"],
            "Skipped Protein": sim_result["skipped_protein"]
        })

    df = pd.DataFrame(results)
    ranked_df = df.sort_values(by="Mutation Count", ascending=False).reset_index(drop=True)
    ranked_df.to_csv("ASO_Count.csv", index=False)
    print(ranked_df.head())

# ======================
# --- Build exon map ---
# ======================
exon_blocks = []
for _, exon_row in enumerate(coding_exons):
    exon_number = exon_row["exon_number"]
    exon_start = exon_row["start"]
    exon_end = exon_row["end"]

    cds_positions = [
        genomic_to_cds_index[pos]
        for pos in range(exon_start, exon_end + 1)
        if pos in genomic_to_cds_index
    ]
    if not cds_positions:
        continue

    exon_blocks.append((exon_number, min(cds_positions), max(cds_positions)))

# ========================================
# --- Define the exon-skipping function ---
# ========================================
def skip_exons_and_translate(cds_sequence, exon_blocks, skip_exon_numbers):
    """Skip selected exons and translate the modified CDS sequence."""
    exon_blocks_sorted = sorted(exon_blocks, key=lambda x: x[1])  # sort by CDS start
    new_cds = "".join(
        cds_sequence[start:end+1]
        for exon_number, start, end in exon_blocks_sorted
        if exon_number not in skip_exon_numbers
    )

    full_protein = str(Seq(cds_sequence).translate(to_stop=False))
    skipped_protein = str(Seq(new_cds).translate(to_stop=False))

    if skipped_protein == full_protein:
        frame_status = "in-frame (no change)"
    elif "*" in skipped_protein and skipped_protein[-1] == "*" and skipped_protein.count("*") == 1:
        frame_status = "in-frame (ends with stop)"
    elif "*" in skipped_protein:
        frame_status = "frameshift (premature stop)"
    elif len(new_cds) % 3 == 0:
        frame_status = "in-frame (altered protein)"
    else:
        frame_status = "frameshift (likely)"

    return {
        "skipped_cds": new_cds,
        "original_protein": full_protein,
        "skipped_protein": skipped_protein,
        "frame_status": frame_status,
        "length_diff": len(cds_sequence) - len(new_cds)
    }

# ==================================
# --- Simulate skipping examples ---
# ==================================
result_50_52 = skip_exons_and_translate(full_cds_sequence, exon_blocks, [48, 49, 50])
print("🧬 Skipping Exons 50–52")
print("Frame effect:", result_50_52["frame_status"])
print("Length diff:", result_50_52["length_diff"])
print("Protein (first 60 aa):", result_50_52["skipped_protein"])
print()

result_50_53 = skip_exons_and_translate(full_cds_sequence, exon_blocks, [48, 49, 50, 51])
print("🧬 Skipping Exons 50–53")
print("Frame effect:", result_50_53["frame_status"])
print("Length diff:", result_50_53["length_diff"])
print("Protein (first 60 aa):", result_50_53["skipped_protein"])
