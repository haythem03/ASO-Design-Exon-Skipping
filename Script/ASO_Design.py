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
