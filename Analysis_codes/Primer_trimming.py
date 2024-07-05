#!/usr/bin/env python

# Manuscript: Assembly and maturation of calf gut microbiome from neonate to post-puberty
# Code written by Jae-Yun Lee
# 2024-07-02
# Job description: Trimming 338F and 805R PCR primer sequences from the raw reads.

#### Program information ####
# Cutadapt v4.0 (Marcel Martin, EMBnet.journal, 2011. doi:10.14806/ej.17.1.200 (https://cutadapt.readthedocs.io/en/stable/index.html))

import argparse
import os
import re
import subprocess
from pathlib import Path

def run_cutadapt(forward_read, reverse_read, outdir, f_primer_seq, r_primer_seq):

    acc_num = re.sub("_1\.fastq\.gz$", "", Path(forward_read).name)
    print(f'Accession number: {acc_num}')
    print(f'\tForward Primer sequence: {f_primer_seq}')
    print(f'\tReverse Primer sequence: {r_primer_seq}')

    cmd = [
        'cutadapt',
        '--cores', '0',
        '-o', os.path.join(outdir, Path(forward_read).name),
        "-p", os.path.join(outdir, Path(reverse_read).name),
        '-g', f_primer_seq,
        '-G', r_primer_seq,
        "--discard-untrimmed",
        forward_read,
        reverse_read
    ]

    print('Command:', ' '.join(cmd), '\n')

    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error occurred: {e.stderr}")
        
    log_name = os.path.join(outdir, f"{acc_num}.log")
    with open(log_name, "wt") as f:
        f.write(out.stdout)

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description="Run cutadapt for primer trimming")
    parser.add_argument("-1", "--forward_fastq", type=str, required=True, help="Forward fastq file path")
    parser.add_argument("-2", "--reverse_fastq", type=str, required=True, help="Reverse fastq file path")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory path")
    parser.add_argument("-f", "--forward_primer", type=str, required=True, default="CCTACGGGNGGCWGCAG", help="Forward primer sequence")
    parser.add_argument("-r", "--reverse_primer", type=str, required=True, default="GACTACHVGGGTATCTAATCC", help="Reverse primer sequence")
    args = parser.parse_args()

    # Run program
    run_cutadapt(args.forward_fastq, args.reverse_fastq, args.outdir, args.forward_primer, args.reverse_primer)

    # Output files had been imported into the Qiime2 platform (v2023.5) (https://qiime2.org/) and underwent further analysis as described in Methods section of our manuscript.
