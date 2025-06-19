import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

def parse_ab1(file_path):
    """
    Parses an .ab1 file to extract sequence and chromatogram traces.
    """
    record = SeqIO.read(file_path, "abi")
    sequence = str(record.seq)
    
    # Extract trace data
    traces = {
        'A': record.annotations['abif_raw']['DATA9'],
        'C': record.annotations['abif_raw']['DATA10'],
        'G': record.annotations['abif_raw']['DATA11'],
        'T': record.annotations['abif_raw']['DATA12'],
    }
    return sequence, traces

def find_divergence_point(seq1, seq2):
    """
    Find the point where two sequences start to diverge significantly.
    """
    min_len = min(len(seq1), len(seq2))
    window_size = 10
    
    for i in range(min_len - window_size):
        window1 = seq1[i:i+window_size]
        window2 = seq2[i:i+window_size]
        mismatches = sum(1 for a, b in zip(window1, window2) if a != b)
        
        if mismatches >= 3:  # If 3 or more mismatches in 10bp window
            return i
    
    return min_len // 2  # Default to middle if no clear divergence

def plot_tide_analysis(control_file, edited_file):
    """
    Create a TIDE-style analysis plot comparing control and edited samples.
    """
    print("\nPerforming TIDE Analysis...")
    
    # Parse both .ab1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    print(f"Control sequence length: {len(control_seq)}")
    print(f"Edited sequence length: {len(edited_seq)}")
    
    # Find where sequences diverge
    divergence_point = find_divergence_point(control_seq, edited_seq)
    print(f"Sequences appear to diverge around position: {divergence_point}")
    
    # Define viewing window
    window_start = max(0, divergence_point - 50)
    window_end = min(len(control_seq), divergence_point + 100)
    
    # Create the plot
    fig = plt.figure(figsize=(16, 10))
    
    # Colors for each base
    colors = {'A': 'green', 'C': 'blue', 'G': 'black', 'T': 'red'}
    
    # Plot 1: Control chromatogram
    ax1 = plt.subplot(3, 1, 1)
    ax1.set_title('Control Sample (WT) Chromatogram', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Signal Intensity', fontsize=12)
    
    # Approximate trace positions (usually ~10-15 trace points per base)
    trace_factor = 12
    trace_start = window_start * trace_factor
    trace_end = window_end * trace_factor
    
    max_signal = 0
    for base, color in colors.items():
        if trace_end <= len(control_traces[base]):
            trace_data = control_traces[base][trace_start:trace_end]
            x_vals = np.arange(len(trace_data))
            ax1.plot(x_vals, trace_data, color=color, label=base, alpha=0.8, linewidth=1)
            max_signal = max(max_signal, np.max(trace_data))
    
    # Mark potential cut site
    cut_trace_pos = (divergence_point - window_start) * trace_factor
    ax1.axvline(x=cut_trace_pos, color='red', linestyle='--', alpha=0.5, label='Divergence point')
    ax1.legend(loc='upper right')
    ax1.set_ylim(0, max_signal * 1.1)
    
    # Plot 2: Edited chromatogram
    ax2 = plt.subplot(3, 1, 2)
    ax2.set_title('Edited Sample (A2) Chromatogram', fontsize=14, fontweight='bold')
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
    ax3.set_title('Overlay: Control (solid) vs Edited (dashed)', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Signal Intensity', fontsize=12)
    ax3.set_xlabel('Trace Position', fontsize=12)
    
    # Focus on region around divergence
    focus_start = max(0, cut_trace_pos - 200)
    focus_end = min(len(control_traces['A']), cut_trace_pos + 400)
    
    for base, color in colors.items():
        # Control traces (solid)
        if focus_end <= len(control_traces[base]):
            control_data = control_traces[base][trace_start:trace_end][focus_start:focus_end]
            x_vals = np.arange(len(control_data))
            ax3.plot(x_vals, control_data, color=color, label=f'{base} (WT)', 
                    alpha=0.7, linewidth=1.5)
        
        # Edited traces (dashed)
        if focus_end <= len(edited_traces[base]):
            edited_data = edited_traces[base][trace_start:trace_end][focus_start:focus_end]
            ax3.plot(x_vals, edited_data, color=color, linestyle='--', 
                    alpha=0.7, linewidth=1.5)
    
    ax3.axvline(x=200, color='red', linestyle='--', alpha=0.5, label='Divergence point')
    ax3.legend(loc='upper right', ncol=2)
    ax3.set_ylim(0, max_signal * 1.1)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = "tide_analysis_demo.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nTIDE analysis plot saved to: {output_file}")
    
    # Show the plot
    plt.show()
    
    # Print sequence comparison around divergence point
    print(f"\nSequence comparison around position {divergence_point}:")
    comp_start = max(0, divergence_point - 20)
    comp_end = min(len(control_seq), len(edited_seq), divergence_point + 30)
    
    control_substr = control_seq[comp_start:comp_end]
    edited_substr = edited_seq[comp_start:comp_end]
    
    print(f"Control: {control_substr}")
    print(f"Edited:  {edited_substr}")
    print(f"         {''.join([' ' if c==e else '^' for c,e in zip(control_substr, edited_substr)])}")

if __name__ == '__main__':
    CONTROL_FILE = "workdir/control.ab1"
    EDITED_FILE = "workdir/edited.ab1"
    
    if os.path.exists(CONTROL_FILE) and os.path.exists(EDITED_FILE):
        plot_tide_analysis(CONTROL_FILE, EDITED_FILE)
    else:
        print("Error: Could not find the required .ab1 files in workdir/") 