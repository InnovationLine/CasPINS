import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

def parse_ab1(file_path):
    """Parse an .ab1 file to extract sequence and chromatogram traces."""
    record = SeqIO.read(file_path, "abi")
    sequence = str(record.seq)
    traces = {
        'A': record.annotations['abif_raw']['DATA9'],
        'C': record.annotations['abif_raw']['DATA10'],
        'G': record.annotations['abif_raw']['DATA11'],
        'T': record.annotations['abif_raw']['DATA12'],
    }
    return sequence, traces

def find_divergence_point(seq1, seq2):
    """Find the point where two sequences start to diverge significantly."""
    min_len = min(len(seq1), len(seq2))
    window_size = 10
    
    for i in range(min_len - window_size):
        window1 = seq1[i:i+window_size]
        window2 = seq2[i:i+window_size]
        mismatches = sum(1 for a, b in zip(window1, window2) if a != b)
        
        if mismatches >= 3:
            return i
    
    return min_len // 2

def plot_tide_analysis(control_file, edited_file, output_dir, gene_name):
    """Create a TIDE-style analysis plot comparing control and edited samples."""
    print(f"\nPerforming TIDE Analysis for {gene_name}...")
    
    # Parse both .ab1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    print(f"  Control sequence length: {len(control_seq)} bp")
    print(f"  Edited sequence length: {len(edited_seq)} bp")
    
    # Find where sequences diverge
    divergence_point = find_divergence_point(control_seq, edited_seq)
    print(f"  Divergence point: position {divergence_point}")
    
    # Define viewing window
    window_start = max(0, divergence_point - 50)
    window_end = min(len(control_seq), divergence_point + 100)
    
    # Create the plot
    fig = plt.figure(figsize=(16, 10))
    
    # Colors for each base
    colors = {'A': 'green', 'C': 'blue', 'G': 'black', 'T': 'red'}
    
    # Plot 1: Control chromatogram
    ax1 = plt.subplot(3, 1, 1)
    ax1.set_title(f'{gene_name.upper()} - Control Sample (WT) Chromatogram', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Signal Intensity', fontsize=12)
    
    trace_factor = 12
    trace_start = window_start * trace_factor
    trace_end = window_end * trace_factor
    
    max_signal = 0
    for base, color in colors.items():
        if trace_end <= len(control_traces[base]):
            trace_data = control_traces[base][trace_start:trace_end]
            x_vals = np.arange(len(trace_data))
            ax1.plot(x_vals, trace_data, color=color, label=base, alpha=0.8, linewidth=1)
            if len(trace_data) > 0:
                max_signal = max(max_signal, np.max(trace_data))
    
    cut_trace_pos = (divergence_point - window_start) * trace_factor
    ax1.axvline(x=cut_trace_pos, color='red', linestyle='--', alpha=0.5, label='Divergence point')
    ax1.legend(loc='upper right')
    ax1.set_ylim(0, max_signal * 1.1)
    
    # Plot 2: Edited chromatogram
    ax2 = plt.subplot(3, 1, 2)
    ax2.set_title(f'{gene_name.upper()} - Edited Sample Chromatogram', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Signal Intensity', fontsize=12)
    
    for base, color in colors.items():
        if trace_end <= len(edited_traces[base]):
            trace_data = edited_traces[base][trace_start:trace_end]
            x_vals = np.arange(len(trace_data))
            ax2.plot(x_vals, trace_data, color=color, label=base, alpha=0.8, linewidth=1)
    
    ax2.axvline(x=cut_trace_pos, color='red', linestyle='--', alpha=0.5, label='Divergence point')
    ax2.legend(loc='upper right')
    ax2.set_ylim(0, max_signal * 1.1)
    
    # Plot 3: Overlay comparison
    ax3 = plt.subplot(3, 1, 3)
    ax3.set_title(f'{gene_name.upper()} - Overlay: Control (solid) vs Edited (dashed)', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Signal Intensity', fontsize=12)
    ax3.set_xlabel('Trace Position', fontsize=12)
    
    # Focus on region around divergence
    focus_start = max(0, cut_trace_pos - 200)
    focus_end = min(cut_trace_pos + 400, (window_end - window_start) * trace_factor)
    
    for base, color in colors.items():
        if trace_end <= len(control_traces[base]) and trace_end <= len(edited_traces[base]):
            try:
                control_data = control_traces[base][trace_start:trace_end][focus_start:focus_end]
                edited_data = edited_traces[base][trace_start:trace_end][focus_start:focus_end]
                
                if len(control_data) > 0 and len(edited_data) > 0:
                    x_vals = np.arange(len(control_data))
                    ax3.plot(x_vals, control_data, color=color, label=f'{base} (WT)', 
                            alpha=0.7, linewidth=1.5)
                    ax3.plot(x_vals, edited_data, color=color, linestyle='--', 
                            alpha=0.7, linewidth=1.5)
            except:
                continue
    
    ax3.axvline(x=min(200, focus_end-focus_start-10), color='red', linestyle='--', alpha=0.5, label='Divergence point')
    ax3.legend(loc='upper right', ncol=2)
    ax3.set_ylim(0, max_signal * 1.1)
    
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    output_subdir = os.path.join(output_dir, "output")
    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir)
    
    # Save the plot with timestamp in output subfolder
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_subdir, f"tide_analysis_{gene_name}_{timestamp}.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  [SUCCESS] Plot saved: {output_file}")
    
    return divergence_point, timestamp

def analyze_grna_in_mrna(mrna_seq, grna_sequences, gene_name):
    """Analyze gRNA positions in mRNA and recommend primers."""
    recommendations = []
    recommendations.append(f"PRIMER RECOMMENDATIONS FOR {gene_name.upper()}")
    recommendations.append("=" * 60)
    recommendations.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    recommendations.append(f"mRNA length: {len(mrna_seq)} bp")
    recommendations.append("")
    
    grna_positions = []
    
    # Find each gRNA in the mRNA
    for i, grna in enumerate(grna_sequences, 1):
        recommendations.append(f"\nGuide RNA {i}: {grna}")
        grna_rc = str(Seq(grna).reverse_complement())
        
        if grna in mrna_seq:
            pos = mrna_seq.find(grna)
            grna_positions.append((pos, 'forward'))
            recommendations.append(f"  [FOUND] on forward strand at position {pos}")
            recommendations.append(f"  Cut site (3bp before PAM): position {pos + 17}")
        elif grna_rc in mrna_seq:
            pos = mrna_seq.find(grna_rc)
            grna_positions.append((pos, 'reverse'))
            recommendations.append(f"  [FOUND] on reverse strand at position {pos}")
            recommendations.append(f"  Cut site: position {pos + 3}")
        else:
            recommendations.append(f"  [NOT FOUND] in mRNA sequence!")
    
    if grna_positions:
        # Sort positions to find the range
        positions = [pos for pos, strand in grna_positions]
        min_pos = min(positions)
        max_pos = max(positions)
        
        recommendations.append("\n" + "=" * 60)
        recommendations.append("RECOMMENDED SEQUENCING PRIMERS:")
        recommendations.append("=" * 60)
        
        # Design primers ~200bp away from the gRNA region
        # Forward primer
        primer_f_start = max(0, min_pos - 200)
        primer_f_seq = mrna_seq[primer_f_start:primer_f_start+20]
        
        recommendations.append(f"\nPrimary Primer Set (covers all gRNAs):")
        recommendations.append(f"Forward Primer: 5'-{primer_f_seq}-3'")
        recommendations.append(f"  Position: {primer_f_start}-{primer_f_start+20}")
        recommendations.append(f"  Tm: ~{calculate_tm(primer_f_seq)}°C")
        
        # Reverse primer
        primer_r_start = min(len(mrna_seq)-20, max_pos + 23 + 180)
        primer_r_seq = str(Seq(mrna_seq[primer_r_start:primer_r_start+20]).reverse_complement())
        
        recommendations.append(f"\nReverse Primer: 5'-{primer_r_seq}-3'")
        recommendations.append(f"  Position: {primer_r_start}-{primer_r_start+20}")
        recommendations.append(f"  Tm: ~{calculate_tm(primer_r_seq)}°C")
        
        recommendations.append(f"\nExpected amplicon size: ~{primer_r_start - primer_f_start + 20} bp")
        recommendations.append(f"gRNA region spans: {min_pos}-{max_pos + 23} ({max_pos - min_pos + 23} bp)")
        
        # Alternative closer primers
        alt_f_start = max(0, min_pos - 100)
        alt_f_seq = mrna_seq[alt_f_start:alt_f_start+20]
        alt_r_start = min(len(mrna_seq)-20, max_pos + 23 + 100)
        alt_r_seq = str(Seq(mrna_seq[alt_r_start:alt_r_start+20]).reverse_complement())
        
        recommendations.append("\n" + "-" * 40)
        recommendations.append("Alternative Primer Set (closer to gRNAs):")
        recommendations.append(f"Forward: 5'-{alt_f_seq}-3' (pos {alt_f_start})")
        recommendations.append(f"Reverse: 5'-{alt_r_seq}-3' (pos {alt_r_start})")
        recommendations.append(f"Amplicon: ~{alt_r_start - alt_f_start + 20} bp")
        
    else:
        recommendations.append("\n*** WARNING: No gRNAs found in the provided mRNA sequence! ***")
        recommendations.append("Please verify:")
        recommendations.append("  1. The mRNA sequence is correct")
        recommendations.append("  2. The gRNA sequences are correct")
        recommendations.append("  3. The gRNAs target this specific gene")
    
    return "\n".join(recommendations)

def calculate_tm(sequence):
    """Simple Tm calculation (Wallace rule)."""
    return (sequence.count('G') + sequence.count('C')) * 4 + (sequence.count('A') + sequence.count('T')) * 2

def process_gene_folder(gene_folder, gene_name):
    """Process a single gene folder."""
    print(f"\nProcessing {gene_name}...")
    
    # Check required files
    control_file = os.path.join(gene_folder, "control.ab1")
    edited_file = os.path.join(gene_folder, "edited.ab1")
    grna_file = os.path.join(gene_folder, "grna.txt")
    mrna_file = os.path.join(gene_folder, "mrna.txt")
    
    # Check if all required files exist (mrna.txt is now mandatory)
    missing_files = []
    for file_path, file_name in [(control_file, "control.ab1"), 
                                  (edited_file, "edited.ab1"), 
                                  (grna_file, "grna.txt"),
                                  (mrna_file, "mrna.txt")]:
        if not os.path.exists(file_path):
            missing_files.append(file_name)
    
    if missing_files:
        print(f"  ERROR: Missing required files in {gene_folder}:")
        for file in missing_files:
            print(f"    - {file}")
        print("  All files (control.ab1, edited.ab1, grna.txt, mrna.txt) are required.")
        return False
    
    # Read gRNA sequences
    with open(grna_file, 'r') as f:
        grna_sequences = [line.strip().upper() for line in f if line.strip()]
    
    print(f"  Found {len(grna_sequences)} gRNA sequences")
    
    # Perform TIDE analysis
    try:
        divergence_point, timestamp = plot_tide_analysis(control_file, edited_file, gene_folder, gene_name)
        
        # Read mRNA sequence
        with open(mrna_file, 'r') as f:
            mrna_seq = ''.join(line.strip() for line in f).upper().replace(' ', '')
        
        # Generate recommendations
        recommendations = analyze_grna_in_mrna(mrna_seq, grna_sequences, gene_name)
        
        # Create output directory if it doesn't exist
        output_subdir = os.path.join(gene_folder, "output")
        if not os.path.exists(output_subdir):
            os.makedirs(output_subdir)
        
        # Save recommendations with timestamp in output subfolder
        rec_file = os.path.join(output_subdir, f"recommendations_{gene_name}_{timestamp}.txt")
        with open(rec_file, 'w', encoding='utf-8') as f:
            f.write(recommendations)
        
        print(f"  [SUCCESS] Recommendations saved: {rec_file}")
        
        return True
        
    except Exception as e:
        print(f"  ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main function to process all gene folders."""
    print("\n" + "="*60)
    print("TIDE ANALYSIS BATCH PIPELINE")
    print("="*60)
    
    # Look for gene folders directly in data/
    data_dir = "data"
    if not os.path.exists(data_dir):
        print(f"ERROR: {data_dir} directory not found!")
        return
    
    gene_folders = []
    
    # Check for standard gene folders
    for gene in ['vmat1', 'vmat2', 'ddc']:
        gene_path = os.path.join(data_dir, gene)
        if os.path.isdir(gene_path):
            gene_folders.append((gene_path, gene))
    
    if not gene_folders:
        print(f"No gene folders found in {data_dir}!")
        print("Please create folders named 'vmat1', 'vmat2', or 'ddc' in the data/ directory")
        print("\nExpected structure:")
        print("  data/")
        print("    vmat1/")
        print("      control.ab1")
        print("      edited.ab1")
        print("      grna.txt")
        print("      mrna.txt")
        return
    
    print(f"\nFound {len(gene_folders)} gene folders to process")
    print("Note: mrna.txt is required for primer recommendations")
    
    # Process each gene folder
    successful = 0
    for folder_path, gene_name in gene_folders:
        if process_gene_folder(folder_path, gene_name):
            successful += 1
    
    print("\n" + "="*60)
    print(f"Pipeline complete! Successfully processed {successful}/{len(gene_folders)} genes")
    print("="*60)

if __name__ == '__main__':
    main() 