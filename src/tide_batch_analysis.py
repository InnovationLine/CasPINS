"""
TIDE Batch Analysis Pipeline - Main Module
Coordinates CRISPR TIDE analysis with modular architecture
"""

import os
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
from utils.visualization import (
    validate_for_plotting, plot_tide_analysis, 
    plot_error_figure, ValidationIssue
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='CRISPR TIDE Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default settings (70% similarity threshold)
  python tide_batch_analysis.py
  
  # Force plotting even with validation issues
  python tide_batch_analysis.py --force-plot
  
  # Use custom similarity threshold
  python tide_batch_analysis.py --similarity-threshold 50
  
  # Interactive mode - ask before each analysis
  python tide_batch_analysis.py --interactive
  
  # Skip all validation checks (not recommended)
  python tide_batch_analysis.py --skip-validation
        """
    )
    
    parser.add_argument(
        '--force-plot', 
        action='store_true',
        help='Generate plots even when validation fails (shows warnings)'
    )
    
    parser.add_argument(
        '--similarity-threshold', 
        type=float, 
        default=70.0,
        help='Minimum sequence similarity percentage for plotting (default: 70%%)'
    )
    
    parser.add_argument(
        '--skip-validation',
        action='store_true',
        help='Skip all validation checks (not recommended)'
    )
    
    parser.add_argument(
        '--interactive',
        action='store_true',
        help='Ask for confirmation before processing each gene'
    )
    
    parser.add_argument(
        '--no-archive',
        action='store_true',
        help='Do not archive previous output files'
    )
    
    parser.add_argument(
        '--genes',
        nargs='+',
        choices=['vmat1', 'vmat2', 'ddc'],
        help='Specific genes to analyze (default: all)'
    )
    
    return parser.parse_args()


def analyze_grna_in_mrna(mrna_seq, grna_sequences, grna_positions, gene_name, efficiency_data=None):
    """
    Analyze gRNA sequences in mRNA and generate recommendations.
    
    Args:
        mrna_seq: mRNA sequence
        grna_sequences: List of gRNA sequences
        grna_positions: List of (position, strand, cut_site) tuples
        gene_name: Gene name
        efficiency_data: Optional TIDE analysis results
        
    Returns:
        str: Formatted recommendations
    """
    recommendations = []
    
    # Header
    recommendations.append("\n" + "="*60)
    recommendations.append(f"CRISPR ANALYSIS REPORT - {gene_name.upper()}")
    recommendations.append("="*60)
    recommendations.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # gRNA Summary
    recommendations.append("\n" + "-"*40)
    recommendations.append("gRNA SUMMARY:")
    recommendations.append("-"*40)
    
    for i, (grna, pos_info) in enumerate(zip(grna_sequences, grna_positions)):
        recommendations.append(f"\ngRNA {i+1}: {grna}")
        if pos_info[0] is not None:
            recommendations.append(f"  Position in mRNA: {pos_info[0]}")
            recommendations.append(f"  Strand: {pos_info[1]}")
            recommendations.append(f"  Cut site: {pos_info[2]}")
        else:
            recommendations.append("  *** NOT FOUND IN mRNA ***")
    
    # TIDE Analysis Results
    if efficiency_data:
        recommendations.append("\n" + "-"*40)
        recommendations.append("TIDE ANALYSIS RESULTS:")
        recommendations.append("-"*40)
        recommendations.append(f"Editing Efficiency: {efficiency_data['editing_efficiency']}%")
        
        if efficiency_data.get('dominant_indel_size') != 'Unknown':
            recommendations.append(f"Dominant Indel: {efficiency_data['dominant_indel_size']}bp ({efficiency_data['dominant_indel_percent']}%)")
        
        recommendations.append(f"Quality Score: {efficiency_data['quality_score']}%")
        recommendations.append(f"Confidence: {efficiency_data.get('confidence', 'Unknown')}")
        
        if efficiency_data.get('error'):
            recommendations.append("\n*** ANALYSIS ERRORS: ***")
            for error in efficiency_data.get('error_messages', []):
                recommendations.append(f"  - {error}")
    
    # Primer Recommendations
    primer_recs = generate_primer_recommendations(grna_positions, mrna_seq, gene_name)
    recommendations.append(primer_recs)
    
    return "\n".join(recommendations)


def perform_tide_analysis(control_file, edited_file, gene_name, expected_cut_site, args):
    """
    Perform TIDE analysis with validation.
    
    Args:
        control_file: Path to control AB1
        edited_file: Path to edited AB1
        gene_name: Gene name
        expected_cut_site: Expected cut site position
        args: Command line arguments
        
    Returns:
        tuple: (efficiency_data, plot_path, timestamp)
    """
    print(f"\nPerforming TIDE Analysis for {gene_name}...")
    
    # Parse AB1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    print(f"  Control sequence length: {len(control_seq)} bp")
    print(f"  Edited sequence length: {len(edited_seq)} bp")
    
    # Calculate similarity
    similarity = calculate_similarity(control_seq, edited_seq)
    print(f"  Sequence similarity: {similarity:.1f}%")
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Validate sequences
    should_plot, validation_issues = validate_for_plotting(
        control_seq, edited_seq, expected_cut_site, similarity, args
    )
    
    if args.skip_validation:
        should_plot = True
        validation_issues = []
    
    # Calculate efficiency regardless of validation
    try:
        # Try proper TIDE decomposition
        efficiency_data = decompose_traces_tide(
            control_traces, edited_traces, expected_cut_site
        )
        
        # Add sequence-based validation
        alignment_pos, aligned_ctrl, aligned_edit = align_sequences(
            control_seq[max(0, expected_cut_site-50):expected_cut_site+50],
            edited_seq[max(0, expected_cut_site-50):expected_cut_site+50]
        )
        
        efficiency_data['sequence_similarity'] = similarity
        
        # Adjust confidence based on issues
        if validation_issues:
            if any(i.severity == 'error' for i in validation_issues):
                efficiency_data['confidence'] = 'LOW - Validation errors'
            else:
                efficiency_data['confidence'] = 'MEDIUM - Validation warnings'
                
    except Exception as e:
        print(f"  WARNING: TIDE decomposition failed, using fallback method")
        efficiency_data = calculate_editing_efficiency_fallback(
            control_seq, control_traces, edited_seq, edited_traces, expected_cut_site
        )
    
    # Generate plot or error figure
    output_dir = ensure_output_dir(os.path.dirname(control_file))
    
    if should_plot:
        plot_path = plot_tide_analysis(
            control_file, edited_file, output_dir, gene_name,
            expected_cut_site, efficiency_data, validation_issues, args
        )
        print(f"  [SUCCESS] Analysis plot saved: {plot_path}")
    else:
        # Create error plot
        error_messages = [issue.message for issue in validation_issues if issue.severity == 'error']
        plot_path = plot_error_figure(output_dir, gene_name, error_messages)
        print(f"  [ERROR] Error plot saved: {plot_path}")
        
        # Add error info to efficiency data
        efficiency_data['error'] = True
        efficiency_data['error_messages'] = error_messages
    
    return efficiency_data, plot_path, timestamp


def process_gene_folder(gene_folder, gene_name, args):
    """
    Process a single gene folder.
    
    Args:
        gene_folder: Path to gene folder
        gene_name: Gene name
        args: Command line arguments
        
    Returns:
        bool: Success status
    """
    print(f"\nProcessing {gene_name}...")
    
    # Interactive mode check
    if args.interactive:
        response = input(f"  Process {gene_name}? (y/n): ").lower()
        if response != 'y':
            print(f"  Skipping {gene_name}")
            return True
    
    # Archive existing output files
    if not args.no_archive:
        archive_output_files(gene_folder)
    
    # Check required files
    files_exist, missing_files = check_required_files(gene_folder)
    if not files_exist:
        print(f"  ERROR: Missing required files in {gene_folder}:")
        for file in missing_files:
            print(f"    - {file}")
        return False
    
    # Read input files
    control_file = os.path.join(gene_folder, "control.ab1")
    edited_file = os.path.join(gene_folder, "edited.ab1")
    grna_file = os.path.join(gene_folder, "grna.txt")
    mrna_file = os.path.join(gene_folder, "mrna.txt")
    
    # Read gRNA sequences
    with open(grna_file, 'r') as f:
        grna_sequences = [line.strip().upper() for line in f if line.strip()]
    
    print(f"  Found {len(grna_sequences)} gRNA sequences")
    
    # Read mRNA sequence
    with open(mrna_file, 'r') as f:
        mrna_seq = ''.join(line.strip() for line in f).upper().replace(' ', '')
    
    # Ensure output directory exists
    output_dir = ensure_output_dir(gene_folder)
    
    # Perform detailed AB1 analysis
    ab1_analysis = analyze_ab1_details(control_file, edited_file, grna_sequences, mrna_seq, gene_name)
    ab1_file = save_results(output_dir, f"ab1_analysis_{gene_name}.txt", ab1_analysis)
    print(f"  [SUCCESS] AB1 analysis saved: {ab1_file}")
    
    # Find gRNAs in mRNA
    grna_positions = []
    expected_cut_sites = []
    
    for i, grna in enumerate(grna_sequences, 1):
        pos, strand, cut_site = find_grna_in_sequence(mrna_seq, grna)
        
        if pos is not None:
            grna_positions.append((pos, strand, cut_site))
            expected_cut_sites.append(cut_site)
            print(f"  gRNA {i} found on {strand} strand at position {pos}, cut site: {cut_site}")
        else:
            grna_positions.append((None, None, None))
            print(f"  ERROR: gRNA {i} ({grna}) not found in mRNA sequence!")
    
    # Process based on whether gRNAs were found
    if not any(pos[0] is not None for pos in grna_positions):
        print("  ERROR: No gRNAs found in the mRNA sequence!")
        
        # Generate error recommendations
        recommendations = analyze_grna_in_mrna(mrna_seq, grna_sequences, grna_positions, gene_name)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        rec_file = save_results(output_dir, f"recommendations_{gene_name}_{timestamp}.txt", recommendations)
        print(f"  [ERROR] Analysis failed. Error report saved: {rec_file}")
        
        return False
    
    # Perform TIDE analysis
    try:
        expected_cut_site = expected_cut_sites[0] if expected_cut_sites else None
        efficiency_data, plot_path, timestamp = perform_tide_analysis(
            control_file, edited_file, gene_name, expected_cut_site, args
        )
        
        # Generate recommendations
        recommendations = analyze_grna_in_mrna(
            mrna_seq, grna_sequences, grna_positions, gene_name, efficiency_data
        )
        
        # Save recommendations
        rec_file = save_results(
            output_dir, 
            f"recommendations_{gene_name}_{timestamp}.txt", 
            recommendations
        )
        print(f"  [SUCCESS] Recommendations saved: {rec_file}")
        
        return True
        
    except Exception as e:
        print(f"  ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main function to process gene folders."""
    # Parse command line arguments
    args = parse_arguments()
    
    print("\n" + "="*60)
    print("TIDE ANALYSIS BATCH PIPELINE")
    print("="*60)
    print(f"Settings:")
    print(f"  Similarity threshold: {args.similarity_threshold}%")
    print(f"  Force plot: {args.force_plot}")
    print(f"  Skip validation: {args.skip_validation}")
    print(f"  Interactive mode: {args.interactive}")
    print(f"  Archive previous files: {not args.no_archive}")
    
    # Look for gene folders
    data_dir = "data"
    if not os.path.exists(data_dir):
        print(f"\nERROR: {data_dir} directory not found!")
        return
    
    # Determine which genes to process
    if args.genes:
        genes_to_process = args.genes
    else:
        genes_to_process = ['vmat1', 'vmat2', 'ddc']
    
    gene_folders = []
    for gene in genes_to_process:
        gene_path = os.path.join(data_dir, gene)
        if os.path.isdir(gene_path):
            gene_folders.append((gene_path, gene))
        else:
            print(f"\nWARNING: Gene folder '{gene}' not found in {data_dir}")
    
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