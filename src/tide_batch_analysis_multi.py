"""
TIDE Batch Analysis Pipeline for Multiple Edited Samples
Processes multiple edited AB1 files against a single control
"""

import os
import glob
import argparse
from datetime import datetime

# Import utility modules
from utils.ab1_parser import parse_ab1
from utils.ab1_analyzer import analyze_ab1_details
from utils.sequence_analysis import (
    calculate_similarity, find_grna_in_sequence, 
    align_sequences, find_divergence_point
)
from utils.tide_algorithm import (
    decompose_traces_tide, calculate_editing_efficiency_fallback
)
from utils.file_management import (
    archive_output_files, ensure_output_dir, 
    save_results, check_required_files
)
from utils.primer_design import generate_primer_recommendations
from utils.visualization_multi import (
    plot_tide_analysis_multi, create_summary_report
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='CRISPR TIDE Analysis Pipeline for Multiple Samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings
  python tide_batch_analysis_multi.py
  
  # Skip specific genes
  python tide_batch_analysis_multi.py --skip-genes vmat1
  
  # Force analysis even without input folder
  python tide_batch_analysis_multi.py --force
        """
    )
    
    parser.add_argument(
        '--skip-genes',
        nargs='+',
        choices=['vmat1', 'vmat2', 'ddc'],
        help='Genes to skip during analysis'
    )
    
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force analysis even if input folder is missing'
    )
    
    parser.add_argument(
        '--no-archive',
        action='store_true',
        help='Do not archive previous output files'
    )
    
    return parser.parse_args()


def check_input_folder(gene_folder):
    """
    Check if input folder exists and contains required files.
    
    Args:
        gene_folder: Path to gene folder
        
    Returns:
        tuple: (exists, control_file, edited_files)
    """
    input_folder = os.path.join(gene_folder, "input")
    
    if not os.path.exists(input_folder):
        return False, None, []
    
    # Check for control file
    control_file = os.path.join(input_folder, "control.ab1")
    if not os.path.exists(control_file):
        print(f"  WARNING: No control.ab1 found in {input_folder}")
        return False, None, []
    
    # Find all edited files
    edited_pattern = os.path.join(input_folder, "edited*.ab1")
    edited_files = sorted(glob.glob(edited_pattern))
    
    if not edited_files:
        print(f"  WARNING: No edited*.ab1 files found in {input_folder}")
        return False, control_file, []
    
    return True, control_file, edited_files


def validate_sequences(control_seq, edited_seq, gene_name, sample_name):
    """
    Validate sequences before TIDE analysis.
    
    Returns:
        tuple: (is_valid, error_messages)
    """
    errors = []
    
    # Check for failed sequencing (N-only sequences)
    if set(control_seq.upper()) == {'N'}:
        errors.append(f"Control sequence contains only N's (failed sequencing)")
    
    if set(edited_seq.upper()) == {'N'}:
        errors.append(f"Edited sequence ({sample_name}) contains only N's (failed sequencing)")
    
    # Check minimum length
    min_length = 200  # Minimum for reliable TIDE analysis
    if len(control_seq) < min_length:
        errors.append(f"Control sequence too short ({len(control_seq)} bp, need >{min_length} bp)")
    
    if len(edited_seq) < min_length:
        errors.append(f"Edited sequence too short ({len(edited_seq)} bp, need >{min_length} bp)")
    
    return len(errors) == 0, errors


def process_single_edited_sample(control_file, edited_file, gene_name, grna_sequences, 
                                mrna_seq, output_dir, sample_name):
    """
    Process a single edited sample against control.
    
    Args:
        control_file: Path to control AB1
        edited_file: Path to edited AB1
        gene_name: Gene name
        grna_sequences: List of gRNA sequences
        mrna_seq: mRNA sequence
        output_dir: Output directory
        sample_name: Name of the edited sample
        
    Returns:
        dict: Analysis results
    """
    print(f"    Processing {sample_name}...")
    
    # Parse AB1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    # Validate sequences
    is_valid, validation_errors = validate_sequences(control_seq, edited_seq, gene_name, sample_name)
    
    if not is_valid:
        print(f"      ERROR: Sequence validation failed:")
        for error in validation_errors:
            print(f"        - {error}")
        
        # Return error result
        return {
            'editing_efficiency': 0.0,
            'dominant_indel_size': 'N/A',
            'dominant_indel_percent': 0.0,
            'indel_spectrum': {},
            'quality_score': 0.0,
            'method': 'ERROR',
            'sequence_similarity': 0.0,
            'sample_name': sample_name,
            'expected_cut_site': None,
            'alignment_window': (0, 0),
            'decomposition_window': (0, 0),
            'confidence': 'ERROR - Sequence validation failed',
            'error_messages': validation_errors,
            'control_seq_length': len(control_seq),
            'edited_seq_length': len(edited_seq)
        }
    
    # Calculate similarity
    similarity = calculate_similarity(control_seq, edited_seq)
    
    # Find gRNA positions and cut sites
    cut_sites = []
    for grna in grna_sequences:
        pos, strand, cut_site = find_grna_in_sequence(mrna_seq, grna)
        if cut_site is not None:
            cut_sites.append(cut_site)
    
    expected_cut_site = cut_sites[0] if cut_sites else None
    
    # Validate cut site is within sequence bounds
    if expected_cut_site and expected_cut_site > len(control_seq) - 50:
        print(f"      WARNING: Cut site ({expected_cut_site}) is beyond usable sequence length ({len(control_seq)})")
        expected_cut_site = None
    
    # Perform TIDE analysis
    try:
        efficiency_data = decompose_traces_tide(
            control_traces, edited_traces, expected_cut_site
        )
        efficiency_data['method'] = 'TIDE'
    except Exception as e:
        print(f"      WARNING: TIDE decomposition failed, using fallback method")
        efficiency_data = calculate_editing_efficiency_fallback(
            control_seq, control_traces, edited_seq, edited_traces, expected_cut_site
        )
        efficiency_data['method'] = 'Fallback'
    
    # Add additional data
    efficiency_data['sequence_similarity'] = similarity
    efficiency_data['sample_name'] = sample_name
    efficiency_data['expected_cut_site'] = expected_cut_site
    
    # Find alignment window
    if expected_cut_site and expected_cut_site > 100 and expected_cut_site < len(control_seq) - 50:
        alignment_start = max(0, expected_cut_site - 100)
        alignment_end = expected_cut_site - 10
        decomp_start = expected_cut_site + 5
        decomp_end = min(len(control_seq), len(edited_seq), decomp_start + 100)
    else:
        # Default windows if no valid cut site
        seq_len = min(len(control_seq), len(edited_seq))
        alignment_start = max(0, seq_len // 4)
        alignment_end = seq_len // 2
        decomp_start = alignment_end + 10
        decomp_end = min(seq_len, decomp_start + 50)
    
    # Ensure decomposition window is valid
    if decomp_end <= decomp_start:
        decomp_end = min(len(control_seq), len(edited_seq))
        decomp_start = max(0, decomp_end - 50)
    
    efficiency_data['alignment_window'] = (alignment_start, alignment_end)
    efficiency_data['decomposition_window'] = (decomp_start, decomp_end)
    efficiency_data['control_seq_length'] = len(control_seq)
    efficiency_data['edited_seq_length'] = len(edited_seq)
    
    return efficiency_data


def process_gene_folder(gene_folder, gene_name, args):
    """
    Process a single gene folder with multiple edited samples.
    
    Args:
        gene_folder: Path to gene folder
        gene_name: Gene name
        args: Command line arguments
        
    Returns:
        bool: Success status
    """
    print(f"\nProcessing {gene_name}...")
    
    # Check for input folder
    input_exists, control_file, edited_files = check_input_folder(gene_folder)
    
    if not input_exists:
        if args.force:
            print(f"  WARNING: No input folder found for {gene_name}, skipping...")
        else:
            print(f"  INFO: No input folder found for {gene_name}, skipping...")
        return False
    
    print(f"  Found {len(edited_files)} edited samples to process")
    
    # Read gRNA and mRNA sequences
    grna_file = os.path.join(gene_folder, "grna.txt")
    mrna_file = os.path.join(gene_folder, "mrna.txt")
    
    if not os.path.exists(grna_file) or not os.path.exists(mrna_file):
        print(f"  ERROR: Missing gRNA or mRNA file for {gene_name}")
        return False
    
    # Read sequences
    with open(grna_file, 'r') as f:
        grna_sequences = [line.strip().upper() for line in f if line.strip()]
    
    with open(mrna_file, 'r') as f:
        mrna_seq = ''.join(line.strip() for line in f).upper().replace(' ', '')
    
    # Archive existing output
    if not args.no_archive:
        archive_output_files(gene_folder)
    
    # Ensure output directory
    output_dir = ensure_output_dir(gene_folder)
    
    # Process each edited sample
    all_results = []
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    for edited_file in edited_files:
        sample_name = os.path.basename(edited_file).replace('.ab1', '')
        
        results = process_single_edited_sample(
            control_file, edited_file, gene_name, 
            grna_sequences, mrna_seq, output_dir, sample_name
        )
        
        all_results.append(results)
        
        # Generate individual TIDE plot
        plot_path = plot_tide_analysis_multi(
            control_file, edited_file, output_dir, gene_name,
            sample_name, results, grna_sequences, timestamp
        )
        
        print(f"      Plot saved: {plot_path}")
    
    # Generate summary report
    summary_path = create_summary_report(
        all_results, output_dir, gene_name, timestamp
    )
    print(f"  Summary report saved: {summary_path}")
    
    # Save detailed results
    detailed_results = {
        'gene': gene_name,
        'timestamp': timestamp,
        'grna_sequences': grna_sequences,
        'num_samples': len(all_results),
        'samples': all_results
    }
    
    import json
    results_file = os.path.join(output_dir, f"tide_results_{gene_name}_{timestamp}.json")
    with open(results_file, 'w') as f:
        json.dump(detailed_results, f, indent=2)
    
    print(f"  Detailed results saved: {results_file}")
    
    return True


def main():
    """Main function to process gene folders."""
    # Parse command line arguments
    args = parse_arguments()
    
    print("\n" + "="*60)
    print("TIDE ANALYSIS - MULTIPLE SAMPLES PIPELINE")
    print("="*60)
    print(f"Processing clonally expanded CRISPR edited cell lines")
    print(f"Each edited sample will be compared against control")
    
    # Look for gene folders
    data_dir = "data"
    if not os.path.exists(data_dir):
        print(f"\nERROR: {data_dir} directory not found!")
        return
    
    # Determine which genes to process
    all_genes = ['vmat1', 'vmat2', 'ddc']
    genes_to_process = [g for g in all_genes if g not in (args.skip_genes or [])]
    
    gene_folders = []
    for gene in genes_to_process:
        gene_path = os.path.join(data_dir, gene)
        if os.path.isdir(gene_path):
            gene_folders.append((gene_path, gene))
    
    if not gene_folders:
        print(f"\nNo gene folders found to process!")
        return
    
    print(f"\nFound {len(gene_folders)} gene folders to process")
    
    # Process each gene folder
    successful = 0
    for folder_path, gene_name in gene_folders:
        if process_gene_folder(folder_path, gene_name, args):
            successful += 1
    
    print("\n" + "="*60)
    print(f"Pipeline complete! Successfully processed {successful}/{len(gene_folders)} genes")
    print("="*60)


if __name__ == '__main__':
    main() 