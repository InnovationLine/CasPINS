"""
Indel Analysis Pipeline for Multiple Edited Samples
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
from utils.indel_analysis import (
    decompose_traces_indel_analysis, calculate_editing_efficiency_fallback
)
from utils.file_management import (
    archive_output_files, ensure_output_dir, 
    save_results, check_required_files
)
from utils.primer_design import generate_primer_recommendations
from utils.visualization_multi import (
    plot_indel_analysis_multi, create_summary_report
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='CRISPR Indel Analysis Pipeline for Multiple Samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings (uses configured data directory)
  python indel_analysis_multi.py
  
  # Specify data directory
  python indel_analysis_multi.py --data-dir /path/to/data
  
  # Skip specific genes
  python indel_analysis_multi.py --skip-genes gene1 gene2
  
  # Force analysis even without input folder
  python indel_analysis_multi.py --force
        """
    )
    
    parser.add_argument(
        '--data-dir',
        help='Path to data directory containing gene folders (overrides config)'
    )
    
    parser.add_argument(
        '--skip-genes',
        nargs='+',
        help='Gene folder names to skip during analysis (e.g., gene1 gene2)'
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





def process_single_edited_sample(control_file, edited_file, gene_name, grna_sequences, 
                                mrna_seq, output_dir, sample_name, r_squared_correction=True):
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
    
    # Note any sequence issues but don't block analysis
    if set(control_seq.upper()) == {'N'} or len(control_seq) < 200:
        print(f"      WARNING: Control sequence quality issue (length: {len(control_seq)})")
    if set(edited_seq.upper()) == {'N'} or len(edited_seq) < 200:
        print(f"      WARNING: Edited sequence quality issue (length: {len(edited_seq)})")
    
    # Calculate similarity
    similarity = calculate_similarity(control_seq, edited_seq)
    
    # Find gRNA positions and cut sites
    cut_sites = []
    for grna in grna_sequences:
        pos, strand, cut_site = find_grna_in_sequence(mrna_seq, grna)
        if cut_site is not None:
            cut_sites.append(cut_site)
    
    expected_cut_site = cut_sites[0] if cut_sites else None
    
    # Note if cut site is beyond sequence bounds but don't block it
    if expected_cut_site and expected_cut_site > len(control_seq):
        print(f"      WARNING: Cut site ({expected_cut_site}) is beyond control sequence length ({len(control_seq)})")
    
    # Perform indel analysis
    try:
        efficiency_data = decompose_traces_indel_analysis(
            control_traces, edited_traces, expected_cut_site, r_squared_correction=r_squared_correction
        )
        efficiency_data['method'] = 'Trace Decomposition'
    except Exception as e:
        print(f"      WARNING: Trace decomposition failed, using fallback method")
        efficiency_data = calculate_editing_efficiency_fallback(
            control_seq, control_traces, edited_seq, edited_traces, expected_cut_site
        )
        efficiency_data['method'] = 'Fallback'
    
    # Add additional data
    efficiency_data['sequence_similarity'] = similarity
    efficiency_data['sample_name'] = sample_name
    efficiency_data['expected_cut_site'] = expected_cut_site
    
    # Find alignment window - similar to archived results
    if expected_cut_site and expected_cut_site > 100:
        alignment_start = max(0, expected_cut_site - 100)
        alignment_end = expected_cut_site - 10
        decomp_start = expected_cut_site + 5
        decomp_end = min(len(control_seq), len(edited_seq), decomp_start + 100)
    else:
        # Default windows if no cut site found
        alignment_start = 50
        alignment_end = 150
        decomp_start = 160
        decomp_end = min(len(control_seq), len(edited_seq), 260)
    
    efficiency_data['alignment_window'] = (alignment_start, alignment_end)
    efficiency_data['decomposition_window'] = (decomp_start, decomp_end) if decomp_start < decomp_end else (decomp_end, decomp_start)
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
    print(f"\n{'='*60}")
    print(f"Processing {gene_name}...")
    print(f"  Folder: {gene_folder}")

    # Check for input folder
    input_exists, control_file, edited_files = check_input_folder(gene_folder)
    print(f"  Input check: exists={input_exists}, control={control_file is not None}, edited_count={len(edited_files)}")

    if not input_exists:
        if args.force:
            print(f"  WARNING: No input folder found for {gene_name}, skipping...")
        else:
            print(f"  INFO: No input folder found for {gene_name}, skipping...")
        return False

    print(f"  Found {len(edited_files)} edited samples to process")
    for ef in edited_files[:5]:  # Show first 5
        print(f"    - {os.path.basename(ef)}")
    if len(edited_files) > 5:
        print(f"   ... and {len(edited_files) - 5} more")

    # Read gRNA and mRNA sequences
    grna_file = os.path.join(gene_folder, "grna.txt")
    mrna_file = os.path.join(gene_folder, "mrna.txt")
    
    print(f"  Checking grna.txt: {os.path.exists(grna_file)}")
    print(f"  Checking mrna.txt: {os.path.exists(mrna_file)}")

    if not os.path.exists(grna_file) or not os.path.exists(mrna_file):
        print(f"  ERROR: Missing gRNA or mRNA file for {gene_name}")
        print(f"    grna_file exists: {os.path.exists(grna_file)}")
        print(f"    mrna_file exists: {os.path.exists(mrna_file)}")
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
        
        r_squared_correction = getattr(args, 'r_squared_correction', True)
        results = process_single_edited_sample(
            control_file, edited_file, gene_name, 
            grna_sequences, mrna_seq, output_dir, sample_name,
            r_squared_correction=r_squared_correction
        )
        
        all_results.append(results)
        
        # Generate individual indel analysis plot
        plot_path = plot_indel_analysis_multi(
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
    results_file = os.path.join(output_dir, f"indel_analysis_{gene_name}_{timestamp}.json")
    with open(results_file, 'w') as f:
        json.dump(detailed_results, f, indent=2)
    
    print(f"  Detailed results saved: {results_file}")
    
    return True


def get_data_directory():
    """Get data directory from config or environment."""
    # Try to load from config
    try:
        from config.settings import get_data_directory as config_get_data_dir
        data_dir = config_get_data_dir()
        if data_dir:
            return data_dir
    except ImportError:
        pass
    
    # Check environment variable
    env_data_dir = os.environ.get('CRISPR_DATA_DIR', '')
    if env_data_dir and os.path.isdir(env_data_dir):
        return env_data_dir
    
    # Check for local data folder as fallback
    if os.path.isdir('data'):
        return 'data'
    
    return None


def main():
    """Main function to process gene folders."""
    # Parse command line arguments
    args = parse_arguments()
    
    print("\n" + "="*60)
    print("INDEL ANALYSIS - MULTIPLE SAMPLES PIPELINE")
    print("="*60)
    print(f"Processing clonally expanded CRISPR edited cell lines")
    print(f"Each edited sample will be compared against control")
    
    # Get data directory - priority: command line > config > env > local
    data_dir = args.data_dir if hasattr(args, 'data_dir') and args.data_dir else get_data_directory()
    
    if not data_dir:
        print(f"\nERROR: Data directory not configured!")
        print("Please set data directory using one of:")
        print("  1. GUI Settings panel")
        print("  2. Environment variable CRISPR_DATA_DIR")
        print("  3. Create a 'data' folder in current directory")
        return
    
    if not os.path.exists(data_dir):
        print(f"\nERROR: Data directory '{data_dir}' not found!")
        return
    
    print(f"\nUsing data directory: {data_dir}")
    
    # Dynamically find all gene folders (not hardcoded list)
    gene_folders = []
    skip_genes = [g.lower() for g in (args.skip_genes or [])]
    
    for item in os.listdir(data_dir):
        item_path = os.path.join(data_dir, item)
        if os.path.isdir(item_path) and not item.startswith('.'):
            # Skip if in skip list
            if item.lower() in skip_genes:
                continue
            # Check if it looks like a gene folder (has input folder)
            if os.path.isdir(os.path.join(item_path, 'input')):
                gene_folders.append((item_path, item))
    
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