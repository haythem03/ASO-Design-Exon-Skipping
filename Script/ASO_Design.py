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
