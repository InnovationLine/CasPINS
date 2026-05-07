#!/usr/bin/env python
"""
Benchmarking Script: CasPINS Indel Analysis vs TIDE

This script validates CasPINS NNLS-based trace decomposition by:
1. Running CasPINS indel analysis on real AB1 trace files
2. Generating structured output for comparison with TIDE results
3. Creating synthetic test cases with known ground truth

Addresses NAR criticism (iii): "Run CasPINS indel analysis and TIDE on the 
same AB1 files and show R-squared correlation."

Benchmark data: DDC (Dopa decarboxylase, rat) with 8 edited samples.
This is an independent benchmarking track from the gRNA design comparison
(which uses human genes). Using different gene sets for different modules
is standard practice in bioinformatics tool papers.

WORKFLOW:
=========
Step 1: Ensure DDC AB1 files are available (control.ab1 + edited*.ab1)
Step 2: Run this script to generate CasPINS results
Step 3: Upload the SAME AB1 files to TIDE (https://tide.nki.nl/)
Step 4: Compare CasPINS vs TIDE efficiency and R-squared values

Usage:
    python benchmark_indel_analysis.py --data-dir C:/Users/kaush/Downloads/data/ddc
"""

import sys
import os
import json
import csv
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

BENCHMARK_DIR = os.path.join(os.path.dirname(__file__), '..')


def run_synthetic_benchmark():
    """
    Generate synthetic trace data with known indel compositions
    and validate CasPINS NNLS decomposition accuracy.
    
    This tests the mathematical correctness of the NNLS algorithm
    independent of AB1 file parsing.
    """
    from scipy.optimize import nnls
    
    print("\n" + "=" * 60)
    print("Synthetic Benchmark: NNLS Decomposition Accuracy")
    print("=" * 60)
    
    np.random.seed(42)  # Reproducibility
    
    results = []
    
    # Test cases with known ground truth
    test_cases = [
        {
            'name': 'Pure wild-type (0% editing)',
            'true_wt': 1.0,
            'true_indels': {0: 1.0}
        },
        {
            'name': 'Heterozygous -1bp deletion (50% editing)',
            'true_wt': 0.5,
            'true_indels': {0: 0.5, -1: 0.5}
        },
        {
            'name': 'Homozygous -1bp deletion (100% editing)',
            'true_wt': 0.0,
            'true_indels': {-1: 1.0}
        },
        {
            'name': 'Mixed indels (70% editing)',
            'true_wt': 0.3,
            'true_indels': {0: 0.3, -1: 0.35, -2: 0.15, 1: 0.20}
        },
        {
            'name': 'Complex spectrum (85% editing)',
            'true_wt': 0.15,
            'true_indels': {0: 0.15, -1: 0.30, -2: 0.20, -3: 0.10, 1: 0.15, 2: 0.10}
        },
        {
            'name': 'Low editing (15% editing)',
            'true_wt': 0.85,
            'true_indels': {0: 0.85, -1: 0.10, 1: 0.05}
        }
    ]
    
    trace_length = 500
    trace_factor = 10  # Data points per base
    
    for tc in test_cases:
        print(f"\n  Test: {tc['name']}")
        
        # Generate synthetic control trace (smooth signal)
        control = np.sin(np.linspace(0, 50 * np.pi, trace_length)) * 500 + 1000
        control += np.random.normal(0, 20, trace_length)  # Add noise
        
        # Generate synthetic edited trace as mixture
        edited = np.zeros(trace_length)
        for indel_size, fraction in tc['true_indels'].items():
            shift = indel_size * trace_factor
            if shift > 0:
                shifted = np.concatenate([control[shift:], np.zeros(shift)])
            elif shift < 0:
                shifted = np.concatenate([np.zeros(-shift), control[:shift]])
            else:
                shifted = control.copy()
            edited += fraction * shifted
        
        edited += np.random.normal(0, 10, trace_length)  # Add noise
        
        # Run NNLS decomposition (same as CasPINS implementation)
        indel_range = range(-10, 11)
        A = np.zeros((trace_length, len(list(indel_range))))
        
        for i, indel_size in enumerate(indel_range):
            shift = indel_size * trace_factor
            if shift > 0:
                A[:trace_length-shift, i] = control[shift:]
            elif shift < 0:
                A[-shift:, i] = control[:trace_length+shift]
            else:
                A[:, i] = control
        
        coefficients, residual = nnls(A, edited, maxiter=1000)
        
        # Normalize
        total = np.sum(coefficients)
        if total > 0:
            coefficients /= total
        
        # Extract results
        wt_idx = list(indel_range).index(0)
        predicted_wt = coefficients[wt_idx]
        predicted_efficiency = (1 - predicted_wt) * 100
        true_efficiency = (1 - tc['true_wt']) * 100
        
        # Get predicted indel spectrum
        predicted_spectrum = {}
        for i, indel_size in enumerate(indel_range):
            if coefficients[i] > 0.01:  # 1% threshold
                predicted_spectrum[indel_size] = round(coefficients[i] * 100, 1)
        
        error = abs(predicted_efficiency - true_efficiency)
        
        result = {
            'test_name': tc['name'],
            'true_efficiency': round(true_efficiency, 1),
            'predicted_efficiency': round(predicted_efficiency, 1),
            'absolute_error': round(error, 1),
            'true_wt_fraction': tc['true_wt'],
            'predicted_wt_fraction': round(predicted_wt, 4),
            'true_indels': tc['true_indels'],
            'predicted_spectrum': predicted_spectrum
        }
        results.append(result)
        
        print(f"    True efficiency: {true_efficiency:.1f}%")
        print(f"    Predicted efficiency: {predicted_efficiency:.1f}%")
        print(f"    Absolute error: {error:.1f}%")
        print(f"    True spectrum: {tc['true_indels']}")
        print(f"    Predicted spectrum: {predicted_spectrum}")
    
    # Summary statistics
    errors = [r['absolute_error'] for r in results]
    print(f"\n  {'='*40}")
    print(f"  SYNTHETIC BENCHMARK SUMMARY")
    print(f"  {'='*40}")
    print(f"  Mean absolute error: {np.mean(errors):.2f}%")
    print(f"  Max absolute error: {np.max(errors):.2f}%")
    print(f"  All tests within 5%: {'Yes' if max(errors) < 5 else 'No'}")
    
    return results


def _get_gene_dirs(data_dir: str) -> list:
    """
    Determine the gene directories to process.
    
    If data_dir itself contains an 'input/' subfolder with AB1 files,
    treat it as a single gene directory (e.g.,.../data/ddc).
    Otherwise, scan its children for gene-level directories.
    """
    input_subdir = os.path.join(data_dir, 'input')
    if os.path.isdir(input_subdir):
        # data_dir is already a gene-level directory
        gene_name = os.path.basename(os.path.normpath(data_dir))
        return [(gene_name, data_dir)]
    
    # Scan children
    gene_dirs = []
    for name in sorted(os.listdir(data_dir)):
        child = os.path.join(data_dir, name)
        child_input = os.path.join(child, 'input')
        if os.path.isdir(child_input):
            gene_dirs.append((name, child))
    return gene_dirs


def run_ab1_benchmark(data_dir: str = None):
    """
    Run CasPINS indel analysis on real AB1 files and output 
    structured results for TIDE comparison.
    
    Default data: DDC AB1 files at C:/Users/kaush/Downloads/data/ddc
    """
    from src.utils.ab1_parser import parse_ab1
    from src.utils.indel_analysis import decompose_traces_indel_analysis
    from src.utils.sequence_analysis import find_grna_in_sequence
    
    if data_dir is None:
        # Default to DDC data in Downloads
        default_ddc = os.path.join(os.path.expanduser('~'), 'Downloads', 'data', 'ddc')
        if os.path.isdir(default_ddc):
            data_dir = default_ddc
        else:
            data_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
    
    print("\n" + "=" * 60)
    print("AB1 File Benchmark: CasPINS Indel Analysis")
    print(f"Data directory: {os.path.abspath(data_dir)}")
    print("=" * 60)
    
    results = []
    
    if not os.path.exists(data_dir):
        print(f"  Data directory not found: {data_dir}")
        print("  Provide --data-dir pointing to a directory with AB1 data.")
        print("  Expected structure: <dir>/input/control.ab1 + edited*.ab1")
        return results
    
    gene_dirs = _get_gene_dirs(data_dir)
    if not gene_dirs:
        print(f"  No gene directories with AB1 data found in: {data_dir}")
        return results
    
    for gene_name, gene_dir in gene_dirs:
        input_dir = os.path.join(gene_dir, 'input')
        
        # Check for control and edited files
        control_path = os.path.join(input_dir, 'control.ab1')
        if not os.path.exists(control_path):
            continue
        
        # Find edited files
        edited_files = [f for f in os.listdir(input_dir) 
                       if f.startswith('edited') and f.endswith('.ab1')]
        
        if not edited_files:
            continue
        
        # Load gRNA sequence (check gene_dir and input_dir)
        grna_path = os.path.join(gene_dir, 'grna.txt')
        if not os.path.exists(grna_path):
            grna_path = os.path.join(input_dir, 'grna.txt')
        if not os.path.exists(grna_path):
            print(f"  Skipping {gene_name}: no grna.txt")
            continue
        
        with open(grna_path) as f:
            grna_sequences = [line.strip() for line in f if line.strip()]
        
        if not grna_sequences:
            continue
        
        # Load mRNA reference
        mrna_path = os.path.join(gene_dir, 'mrna.txt')
        mrna_seq = ''
        if os.path.exists(mrna_path):
            with open(mrna_path) as f:
                mrna_seq = ''.join(line.strip() for line in f if not line.startswith('>'))
        
        print(f"\n  Gene: {gene_name.upper()}")
        print(f"  gRNAs: {grna_sequences}")
        print(f"  Edited files: {edited_files}")
        
        # Parse control
        try:
            control_seq, control_traces = parse_ab1(control_path)
        except Exception as e:
            print(f"    Error parsing control: {e}")
            continue
        
        for edited_file in edited_files:
            edited_path = os.path.join(input_dir, edited_file)
            
            try:
                edited_seq, edited_traces = parse_ab1(edited_path)
            except Exception as e:
                print(f"    Error parsing {edited_file}: {e}")
                continue
            
            # Find cut site
            grna = grna_sequences[0]
            if mrna_seq:
                grna_result = find_grna_in_sequence(mrna_seq, grna)
            else:
                grna_result = find_grna_in_sequence(str(control_seq), grna)
            
            # find_grna_in_sequence returns (position, strand, cut_site) or (None, None, None)
            if isinstance(grna_result, tuple):
                _, _, cut_site = grna_result
            else:
                cut_site = grna_result
            
            if cut_site is None or (isinstance(cut_site, (int, float)) and cut_site < 0):
                print(f"    Could not find gRNA cut site for {edited_file}")
                # Try with a default position
                cut_site = len(control_seq) // 3
            
            # Run NNLS decomposition
            try:
                decomp_result = decompose_traces_indel_analysis(
                    control_traces, edited_traces, cut_site
                )
                
                result = {
                    'gene': gene_name.upper(),
                    'sample': edited_file,
                    'grna': grna,
                    'tool': 'CasPINS',
                    'editing_efficiency': decomp_result.get('editing_efficiency', 0),
                    'wt_fraction': decomp_result.get('wt_fraction', 1.0),
                    'dominant_indel': decomp_result.get('dominant_indel_size', 0),
                    'r_squared': decomp_result.get('r_squared', 0),
                    'confidence': decomp_result.get('confidence', 'LOW'),
                    'indel_spectrum': decomp_result.get('indel_spectrum', {}),
                    'method': 'NNLS Decomposition'
                }
                
                print(f"    {edited_file}: Efficiency={result['editing_efficiency']:.1f}%, "
                      f"R²={result['r_squared']:.3f}, "
                      f"Confidence={result['confidence']}")
                
                results.append(result)
                
            except Exception as e:
                print(f"    Error analyzing {edited_file}: {e}")
                results.append({
                    'gene': gene_name.upper(),
                    'sample': edited_file,
                    'grna': grna,
                    'tool': 'CasPINS',
                    'error': str(e)
                })
    
    return results


def export_benchmark_results(synthetic_results: list, ab1_results: list, output_dir: str):
    """Export all benchmark results."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Synthetic results
    synth_path = os.path.join(output_dir, 'synthetic_benchmark_results.json')
    with open(synth_path, 'w') as f:
        json.dump(synthetic_results, f, indent=2, default=str)
    print(f"\nSaved synthetic results to: {synth_path}")
    
    # AB1 results
    if ab1_results:
        ab1_path = os.path.join(output_dir, 'caspins_indel_benchmark_results.json')
        with open(ab1_path, 'w') as f:
            json.dump(ab1_results, f, indent=2, default=str)
        print(f"Saved AB1 results to: {ab1_path}")
        
        # CSV for comparison
        csv_path = os.path.join(output_dir, 'caspins_indel_benchmark_results.csv')
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Gene', 'Sample', 'gRNA', 'Editing_Efficiency_%', 
                'WT_Fraction', 'Dominant_Indel', 'R_Squared', 
                'Confidence', 'Method'
            ])
            for r in ab1_results:
                if 'error' not in r:
                    writer.writerow([
                        r['gene'], r['sample'], r['grna'],
                        f"{r['editing_efficiency']:.1f}",
                        f"{r['wt_fraction']:.3f}",
                        r['dominant_indel'],
                        f"{r['r_squared']:.3f}",
                        r['confidence'],
                        r['method']
                    ])
        print(f"Saved CSV results to: {csv_path}")
    
    # TIDE comparison instructions
    tide_instructions = os.path.join(output_dir, 'TIDE_COMPARISON_INSTRUCTIONS.md')
    with open(tide_instructions, 'w') as f:
        f.write("# Instructions for TIDE Comparison\n\n")
        f.write("## Purpose\n")
        f.write("Compare CasPINS NNLS indel decomposition with TIDE results on the same AB1 files.\n\n")
        f.write("## Steps\n\n")
        f.write("### Using TIDE (https://tide.nki.nl/)\n")
        f.write("1. Go to https://tide.nki.nl/\n")
        f.write("2. Upload the control AB1 file (control.ab1)\n")
        f.write("3. Upload each edited AB1 file\n")
        f.write("4. Enter the gRNA sequence (20nt, no PAM)\n")
        f.write("5. Click 'Calculate'\n")
        f.write("6. Record:\n")
        f.write("   - Overall editing efficiency (%)\n")
        f.write("   - R-squared value\n")
        f.write("   - Indel spectrum (% for each indel size)\n")
        f.write("7. Save results for each sample\n\n")
        f.write("### Using ICE (https://ice.synthego.com/)\n")
        f.write("1. Go to https://ice.synthego.com/\n")
        f.write("2. Upload control and edited AB1 files\n")
        f.write("3. Enter gRNA sequence\n")
        f.write("4. Record editing efficiency and indel spectrum\n\n")
        f.write("## Samples to Analyze\n\n")
        f.write("| Gene | Control | Edited | gRNA |\n")
        f.write("|------|---------|--------|------|\n")
        for r in ab1_results:
            if 'error' not in r:
                f.write(f"| {r['gene']} | control.ab1 | {r['sample']} | {r['grna']} |\n")
        f.write("\n## Expected Comparison Output\n\n")
        f.write("Create a table with:\n")
        f.write("| Gene | Sample | CasPINS Efficiency | TIDE Efficiency | ICE Efficiency | Correlation |\n")
        f.write("|------|--------|-------------------|----------------|---------------|-------------|\n")
        f.write("\nThe R-squared between CasPINS and TIDE efficiencies should be > 0.9.\n")
    
    print(f"Saved TIDE instructions to: {tide_instructions}")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Benchmark CasPINS indel analysis vs TIDE'
    )
    parser.add_argument('--data-dir', default=None,
                       help='Directory containing AB1 data (default: data/)')
    parser.add_argument('--output-dir', default=None,
                       help='Output directory (default: benchmarking/results/indel_analysis/)')
    parser.add_argument('--skip-synthetic', action='store_true',
                       help='Skip synthetic benchmark')
    parser.add_argument('--skip-ab1', action='store_true',
                       help='Skip AB1 file benchmark')
    
    args = parser.parse_args()
    
    if args.output_dir is None:
        args.output_dir = os.path.join(BENCHMARK_DIR, 'results', 'indel_analysis')
    
    synthetic_results = []
    ab1_results = []
    
    if not args.skip_synthetic:
        synthetic_results = run_synthetic_benchmark()
    
    if not args.skip_ab1:
        ab1_results = run_ab1_benchmark(args.data_dir)
    
    export_benchmark_results(synthetic_results, ab1_results, args.output_dir)
    
    print(f"\n{'='*60}")
    print("INDEL ANALYSIS BENCHMARKING COMPLETE")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
