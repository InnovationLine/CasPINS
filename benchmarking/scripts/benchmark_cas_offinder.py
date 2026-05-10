#!/usr/bin/env python
"""
Cas-OFFinder Off-Target Benchmark

Addresses reviewer comment R3-Mand-3:
  "While the authors mention off-target scoring for designed gRNAs, very little
  detail is provided how this is calculated. The manuscript also doesn't mention
  any existing tools specifically for gRNA off-target estimation (e.g.
  Cas-OFFinder). Evaluating off-target effects would be another important
  comparison for the benchmarking section."

WORKFLOW (two-stage — requires manual Cas-OFFinder execution):
================================================================

  Stage 1 — Run this script (generates input files):
      python benchmark_cas_offinder.py --stage prepare

  Stage 2 — Run Cas-OFFinder manually (see instructions below):
      cas-offinder CasOFFinder_input.txt C benchmarking/data/cas_offinder/output.txt

  Stage 3 — Compare results:
      python benchmark_cas_offinder.py --stage compare

  OR run both automatically if Cas-OFFinder binary is on your PATH:
      python benchmark_cas_offinder.py --stage auto --cas-offinder-path cas-offinder

Cas-OFFinder download:
    https://github.com/snugel/cas-offinder
    Web version: https://www.rgenome.net/cas-offinder/

WHAT IS COMPARED:
  For the top-5 gRNAs per benchmark gene (TP53, ATE1, VEGFA, DBH, EMX1):
  - CasPINS off-target count at MM0 / MM1 / MM2 / MM3
  - Cas-OFFinder off-target count at MM0 / MM1 / MM2 / MM3
  Agreement between the two is reported as Table S6.

NOTE on CasPINS off-target methodology:
  CasPINS implements the Cutting Frequency Determination (CFD) score from
  Doench et al. 2016, combined with a sequence-alignment scan (MM0-3) against
  a locally cached reference. Because CasPINS targets locally-stored genome
  sequences (not the full genome index), its off-target counts are approximate
  for gRNAs whose sequences match elsewhere in the genome. Cas-OFFinder performs
  an exhaustive genome-wide search and serves as the ground truth here.
"""

import sys
import os
import json
import csv
import argparse
import subprocess
import shutil
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

BENCHMARK_DIR = os.path.join(os.path.dirname(__file__), '..')
RESULTS_DIR = os.path.join(BENCHMARK_DIR, 'results', 'grna_design')
DATA_DIR = os.path.join(BENCHMARK_DIR, 'data')
OFFINDER_DIR = os.path.join(DATA_DIR, 'cas_offinder')
OUTPUT_DIR = os.path.join(RESULTS_DIR, 'comparison')

# Top N gRNAs per gene to submit to Cas-OFFinder
TOP_N = 5

# Genes and SpCas9 NGG PAM
BENCHMARK_GENES = ['TP53', 'ATE1', 'VEGFA', 'DBH', 'EMX1']
PAM_PATTERN = 'NGG'
GENOME = 'hg38'


# ── helpers ──────────────────────────────────────────────────────────────────

def load_caspins_top_grnas(gene: str, top_n: int = TOP_N) -> list:
    """
    Load top-N gRNAs for a gene from CasPINS benchmark results.
    Returns list of dicts: [{sequence, composite_score, off_target_count, cfd_score}, ...]
    """
    json_path = os.path.join(RESULTS_DIR, 'caspins_grna_benchmark_results.json')
    if not os.path.exists(json_path):
        print(f"ERROR: {json_path} not found. Run benchmark_grna_design.py first.")
        sys.exit(1)

    with open(json_path) as f:
        data = json.load(f)

    for entry in data:
        if entry.get('gene', '').upper() == gene.upper():
            grnas = entry.get('grnas', [])
            result = []
            for g in grnas[:top_n]:
                scores = g.get('scores', {})
                result.append({
                    'gene': gene,
                    'sequence': g['sequence'].upper(),
                    'pam': g.get('pam', 'NGG'),
                    'composite_score': g.get('composite_score', scores.get('composite', 0)),
                    'doench_2016': scores.get('doench_2016', 0),
                    'cfd_score': scores.get('cfd', None),
                    'caspins_ot_count': g.get('off_targets_count', None),
                })
            return result

    print(f"WARNING: No CasPINS results found for gene {gene}")
    return []


# ── Stage 1: Prepare ─────────────────────────────────────────────────────────

def prepare_cas_offinder_input() -> str:
    """
    Generate the Cas-OFFinder input file for the top-N gRNAs of each benchmark gene.

    Cas-OFFinder input format (tab-separated):
      Line 1: genome_db_dir
      Line 2: PAM + N×20 for the gRNA pattern  (e.g.  NNNNNNNNNNNNNNNNNNNNNGG)
      Lines 3+: grna_sequence + PAM_bases  mismatches_allowed  [name]

    We request MM0 through MM3 by writing 4 input lines per gRNA,
    each with a different mismatch count.
    """
    os.makedirs(OFFINDER_DIR, exist_ok=True)

    # Cas-OFFinder needs a local genome DB directory.
    # For hg38, download from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
    # This script generates the input; the user must set GENOME_DB_DIR.
    genome_db_env = os.environ.get('CAS_OFFINDER_GENOME_DB', '')
    if not genome_db_env:
        genome_db_env = '/path/to/hg38_genome_db'  # placeholder
        print("NOTE: Set CAS_OFFINDER_GENOME_DB env var to the directory containing")
        print("      hg38 chromosome FASTA files for Cas-OFFinder.")
        print(f"      Placeholder written: {genome_db_env}")

    # Collect all gRNAs
    all_grnas = []
    for gene in BENCHMARK_GENES:
        grnas = load_caspins_top_grnas(gene, TOP_N)
        all_grnas.extend(grnas)
        print(f"  {gene}: {len(grnas)} gRNAs selected")

    if not all_grnas:
        print("ERROR: No gRNAs found. Run benchmark_grna_design.py first.")
        sys.exit(1)

    # Write Cas-OFFinder input file
    input_path = os.path.join(OFFINDER_DIR, 'CasOFFinder_input.txt')
    with open(input_path, 'w') as f:
        # Line 1: genome DB path
        f.write(genome_db_env + '\n')
        # Line 2: PAM pattern (Cas-OFFinder: 20 N's for gRNA + NGG)
        f.write('N' * 20 + PAM_PATTERN + '\n')
        # One line per gRNA × one for each mismatch level 0-3
        for entry in all_grnas:
            seq = entry['sequence']
            grna_id = f"{entry['gene']}_r{all_grnas.index(entry) + 1}"
            # Cas-OFFinder expects: gRNA_seq (20nt) + PAM_sequence (3nt) mismatch_count id
            # For NGG: append NNN as wildcard PAM (Cas-OFFinder matches NGG internally)
            for mm in range(4):
                f.write(f"{seq}NNN\t{mm}\t{grna_id}_mm{mm}\n")

    print(f"\nCas-OFFinder input written to: {input_path}")

    # Write a metadata mapping for the compare stage
    meta_path = os.path.join(OFFINDER_DIR, 'grna_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump(all_grnas, f, indent=2)

    # Write manual execution instructions
    instr_path = os.path.join(OFFINDER_DIR, 'INSTRUCTIONS.md')
    with open(instr_path, 'w') as f:
        f.write("# Cas-OFFinder Execution Instructions\n\n")
        f.write("## Download Cas-OFFinder\n")
        f.write("  https://github.com/snugel/cas-offinder/releases\n")
        f.write("  OR use the web version: https://www.rgenome.net/cas-offinder/\n\n")
        f.write("## Download hg38 genome\n")
        f.write("  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\n")
        f.write("  Unzip chromosomes into a directory, e.g.: ~/hg38_db/\n\n")
        f.write("## Set environment variable\n")
        f.write("  export CAS_OFFINDER_GENOME_DB=~/hg38_db\n\n")
        f.write("## Re-run the prepare step\n")
        f.write("  python benchmarking/scripts/benchmark_cas_offinder.py --stage prepare\n")
        f.write("  (This regenerates CasOFFinder_input.txt with the correct genome path)\n\n")
        f.write("## Run Cas-OFFinder (CPU mode)\n")
        f.write("  cas-offinder CasOFFinder_input.txt C CasOFFinder_output.txt\n\n")
        f.write("## Run comparison\n")
        f.write("  python benchmarking/scripts/benchmark_cas_offinder.py --stage compare\n\n")
        f.write("## Expected output\n")
        f.write(f"  gRNAs submitted: {len(all_grnas)} ({TOP_N} per gene × {len(BENCHMARK_GENES)} genes)\n")
        f.write("  Mismatch levels: 0, 1, 2, 3\n")
        f.write("  Results: benchmarking/results/grna_design/comparison/offtarget_comparison.csv\n")

    print(f"Instructions written to: {instr_path}")
    return input_path


# ── Stage 2: Parse Cas-OFFinder output ───────────────────────────────────────

def parse_cas_offinder_output(output_path: str) -> dict:
    """
    Parse Cas-OFFinder output file.

    Cas-OFFinder output format (tab-separated):
      sequence  chromosome  position  direction  mismatches  [id]

    Returns dict: {grna_id: {mm0: count, mm1: count, mm2: count, mm3: count}}
    """
    if not os.path.exists(output_path):
        print(f"ERROR: Cas-OFFinder output not found at {output_path}")
        print("Run Cas-OFFinder first (see INSTRUCTIONS.md)")
        return {}

    counts = {}   # grna_id_base -> {0: int, 1: int, 2: int, 3: int}

    with open(output_path, newline='') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue

            # Cas-OFFinder output columns:
            # 0: guide+PAM sequence submitted
            # 1: chromosome
            # 2: position
            # 3: direction
            # 4: mismatch count
            # 5 (optional): ID we appended in the input
            mm = int(parts[4])
            grna_id_full = parts[5].strip() if len(parts) > 5 else ''
            # Strip the _mm{N} suffix to get base ID
            grna_base = '_'.join(grna_id_full.split('_')[:-1]) if '_mm' in grna_id_full else grna_id_full

            if grna_base not in counts:
                counts[grna_base] = {0: 0, 1: 0, 2: 0, 3: 0}
            if mm <= 3:
                counts[grna_base][mm] += 1

    # MM0 count should exclude the on-target itself (1 perfect match = the guide itself)
    for gid in counts:
        counts[gid][0] = max(0, counts[gid][0] - 1)

    return counts


# ── Stage 3: Compare ─────────────────────────────────────────────────────────

def compare_offtargets(cas_offinder_counts: dict) -> list:
    """
    Compare CasPINS off-target estimates with Cas-OFFinder ground truth.
    Returns list of dicts for Supplementary Table S6.
    """
    meta_path = os.path.join(OFFINDER_DIR, 'grna_metadata.json')
    if not os.path.exists(meta_path):
        print("ERROR: grna_metadata.json not found. Run --stage prepare first.")
        sys.exit(1)

    with open(meta_path) as f:
        all_grnas = json.load(f)

    records = []
    for i, entry in enumerate(all_grnas, 1):
        gene = entry['gene']
        seq = entry['sequence']
        grna_id = f"{gene}_r{i}"

        co_counts = cas_offinder_counts.get(grna_id, {})

        record = {
            'gene': gene,
            'grna_sequence': seq,
            'composite_score': round(entry.get('composite_score', 0), 4),
            'doench_2016': round(entry.get('doench_2016', 0), 4),
            # CasPINS off-target (approximate, sequence-alignment based)
            'caspins_ot_count': entry.get('caspins_ot_count', 'N/A'),
            'caspins_cfd_score': round(entry.get('cfd_score', 0), 4) if entry.get('cfd_score') else 'N/A',
            # Cas-OFFinder (genome-wide, ground truth)
            'casoffinder_mm0': co_counts.get(0, 'N/A'),
            'casoffinder_mm1': co_counts.get(1, 'N/A'),
            'casoffinder_mm2': co_counts.get(2, 'N/A'),
            'casoffinder_mm3': co_counts.get(3, 'N/A'),
        }

        # Compute agreement if both tools have data
        if (record['caspins_ot_count'] not in (None, 'N/A') and
                record['casoffinder_mm3'] not in (None, 'N/A')):
            # Agreement: CasPINS count vs Cas-OFFinder MM3 total
            diff = abs(int(record['caspins_ot_count']) - int(record['casoffinder_mm3']))
            record['count_difference'] = diff
            record['agreement'] = 'HIGH' if diff <= 2 else ('MODERATE' if diff <= 10 else 'LOW')
        else:
            record['count_difference'] = 'N/A'
            record['agreement'] = 'N/A'

        records.append(record)

    return records


def write_comparison_table(records: list):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    csv_path = os.path.join(OUTPUT_DIR, 'offtarget_comparison.csv')
    fieldnames = [
        'gene', 'grna_sequence', 'composite_score', 'doench_2016',
        'caspins_ot_count', 'caspins_cfd_score',
        'casoffinder_mm0', 'casoffinder_mm1', 'casoffinder_mm2', 'casoffinder_mm3',
        'count_difference', 'agreement',
    ]
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)

    print(f"\nSupplementary Table S6 saved to: {csv_path}")

    # Summary statistics
    n = len(records)
    high = sum(1 for r in records if r.get('agreement') == 'HIGH')
    moderate = sum(1 for r in records if r.get('agreement') == 'MODERATE')
    low = sum(1 for r in records if r.get('agreement') == 'LOW')
    na = sum(1 for r in records if r.get('agreement') == 'N/A')

    print(f"\nSummary ({n} gRNAs across {len(BENCHMARK_GENES)} genes):")
    print(f"  High agreement (diff ≤2):     {high} ({high/n*100:.0f}%)")
    print(f"  Moderate agreement (diff ≤10): {moderate} ({moderate/n*100:.0f}%)")
    print(f"  Low agreement:                {low} ({low/n*100:.0f}%)")
    if na:
        print(f"  N/A (no Cas-OFFinder data):   {na}")


# ── Auto mode: run Cas-OFFinder if binary is available ───────────────────────

def run_cas_offinder_auto(input_path: str, cas_offinder_bin: str) -> str:
    """Run Cas-OFFinder binary and return output path."""
    output_path = os.path.join(OFFINDER_DIR, 'CasOFFinder_output.txt')
    cmd = [cas_offinder_bin, input_path, 'C', output_path]
    print(f"\nRunning: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            print(f"Cas-OFFinder error:\n{result.stderr}")
            sys.exit(1)
        print(f"Cas-OFFinder completed. Output: {output_path}")
        return output_path
    except FileNotFoundError:
        print(f"ERROR: Cas-OFFinder binary not found at '{cas_offinder_bin}'")
        print("Download from https://github.com/snugel/cas-offinder/releases")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: Cas-OFFinder timed out after 1 hour.")
        sys.exit(1)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Cas-OFFinder off-target benchmark for CasPINS R1 revision'
    )
    parser.add_argument(
        '--stage',
        choices=['prepare', 'compare', 'auto'],
        default='prepare',
        help='prepare: generate input files; compare: analyse Cas-OFFinder output; '
             'auto: prepare + run Cas-OFFinder + compare (requires binary on PATH)',
    )
    parser.add_argument(
        '--cas-offinder-path',
        default='cas-offinder',
        help='Path to Cas-OFFinder binary (for --stage auto)',
    )
    parser.add_argument(
        '--output-file',
        default=None,
        help='Cas-OFFinder output file to parse (for --stage compare)',
    )
    args = parser.parse_args()

    print('=' * 60)
    print('Cas-OFFinder Off-Target Benchmark')
    print(f'Stage: {args.stage}')
    print('=' * 60)

    if args.stage in ('prepare', 'auto'):
        input_path = prepare_cas_offinder_input()

    if args.stage == 'auto':
        output_path = run_cas_offinder_auto(input_path, args.cas_offinder_path)
    elif args.stage == 'compare':
        output_path = args.output_file or os.path.join(OFFINDER_DIR, 'CasOFFinder_output.txt')
    else:
        print(f"\nNext step: Run Cas-OFFinder manually, then:")
        print(f"  python benchmark_cas_offinder.py --stage compare")
        return

    if args.stage in ('compare', 'auto'):
        print(f"\nParsing Cas-OFFinder output: {output_path}")
        counts = parse_cas_offinder_output(output_path)
        if not counts:
            print("No Cas-OFFinder counts loaded. Check output file format.")
            return
        records = compare_offtargets(counts)
        write_comparison_table(records)

    print('\n' + '=' * 60)
    print('DONE')
    print('=' * 60)


if __name__ == '__main__':
    main()
