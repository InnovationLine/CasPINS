"""
Visualization Module for Multiple Sample TIDE Analysis
Generates TIDE-style plots similar to the web tool
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
from datetime import datetime


def plot_tide_analysis_multi(control_file, edited_file, output_dir, gene_name,
                           sample_name, analysis_results, grna_sequences, timestamp):
    """
    Create TIDE-style analysis plot for a single edited sample.
    
    Args:
        control_file: Path to control AB1
        edited_file: Path to edited AB1
        output_dir: Output directory
        gene_name: Gene name
        sample_name: Sample name (e.g., editedA2)
        analysis_results: TIDE analysis results
        grna_sequences: List of gRNA sequences
        timestamp: Timestamp for file naming
        
    Returns:
        str: Path to saved plot
    """
    from .ab1_parser import parse_ab1
    
    # Parse AB1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    # Create figure with TIDE-style layout
    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.5, 1.5, 1], hspace=0.3, wspace=0.3)
    
    # Title
    fig.suptitle(f"TIDE Analysis - {gene_name.upper()} - {sample_name}", fontsize=18, fontweight='bold')
    
    # Get analysis data
    expected_cut_site = analysis_results.get('expected_cut_site', None)
    alignment_window = analysis_results.get('alignment_window', (0, 200))
    decomp_window = analysis_results.get('decomposition_window', (210, 310))
    
    # Top panel: Input visualization (alignment + decomposition windows)
    ax_input = fig.add_subplot(gs[0, :])
    plot_input_windows(ax_input, control_seq, edited_seq, control_traces, edited_traces,
                      expected_cut_site, alignment_window, decomp_window, grna_sequences[0])
    
    # Middle left: Alignment window detail
    ax_align = fig.add_subplot(gs[1, 0])
    plot_alignment_window(ax_align, control_seq, edited_seq, control_traces, edited_traces,
                         expected_cut_site, alignment_window)
    
    # Middle right: Decomposition window detail
    ax_decomp = fig.add_subplot(gs[1, 1])
    plot_decomposition_window(ax_decomp, control_seq, edited_seq, control_traces, edited_traces,
                            expected_cut_site, decomp_window)
    
    # Bottom left: Indel spectrum
    ax_indel = fig.add_subplot(gs[2, 0])
    plot_indel_spectrum_tide(ax_indel, analysis_results)
    
    # Bottom right: Efficiency pie chart
    ax_pie = fig.add_subplot(gs[2, 1])
    plot_efficiency_pie_tide(ax_pie, analysis_results)
    
    # Add analysis info text
    info_text = f"Editing Efficiency: {analysis_results['editing_efficiency']}%\n"
    info_text += f"Wild-type: {analysis_results.get('wt_fraction', 100 - analysis_results['editing_efficiency']):.1f}%\n"
    info_text += f"Quality Score: {analysis_results.get('quality_score', 0):.1f}%\n"
    info_text += f"Confidence: {analysis_results.get('confidence', 'Unknown')}\n"
    info_text += f"Analysis Method: {analysis_results.get('method', 'TIDE')}"
    
    fig.text(0.02, 0.02, info_text, fontsize=10, bbox=dict(boxstyle="round,pad=0.3", 
                                                           facecolor="lightgray", alpha=0.5))
    
    # Save plot
    plot_filename = f"tide_{gene_name}_{sample_name}_{timestamp}.png"
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')  # Reduced DPI for faster generation
    plt.close(fig)  # Explicitly close the figure
    plt.clf()  # Clear the current figure
    
    # Force garbage collection
    import gc
    gc.collect()
    
    return plot_path


def plot_input_windows(ax, control_seq, edited_seq, control_traces, edited_traces,
                      cut_site, alignment_window, decomp_window, grna_seq):
    """Plot the input visualization showing alignment and decomposition windows."""
    
    # Plot parameters
    trace_factor = 10  # Scaling factor for trace to sequence position
    window_buffer = 50
    
    # Handle case where cut_site is None or invalid
    if not cut_site or cut_site <= 0:
        # Use middle of alignment window as reference
        cut_site = (alignment_window[0] + alignment_window[1]) // 2
    
    if True:  # Always plot, even without valid cut site
        # Calculate display range
        display_start = max(0, alignment_window[0] - window_buffer)
        display_end = min(len(control_seq), decomp_window[1] + window_buffer)
        
        # Convert to trace positions
        trace_start = display_start * trace_factor
        trace_end = display_end * trace_factor
        
        # Ensure within bounds
        trace_end = min(trace_end, len(control_traces['A']))
        
        x_vals = np.arange(trace_start, trace_end)
        x_vals_bp = x_vals / trace_factor
        
        # Plot control trace (combined signal)
        if len(x_vals) > 0:
            control_signal = np.zeros(len(x_vals))
            for base in ['A', 'C', 'G', 'T']:
                if trace_end <= len(control_traces[base]):
                    control_signal += control_traces[base][trace_start:trace_end]
            
            # Normalize
            if len(control_signal) > 0 and np.max(control_signal) > 0:
                control_signal = control_signal / np.max(control_signal)
                ax.plot(x_vals_bp, control_signal, 'b-', alpha=0.7, linewidth=1, label='Control')
        else:
            # Handle case when x_vals is empty
            ax.text(0.5, 0.5, 'No data available in this range', 
                   ha='center', va='center', transform=ax.transAxes)
        
        # Mark alignment window
        ax.axvspan(alignment_window[0], alignment_window[1], alpha=0.2, color='green', label='Alignment Window')
        
        # Mark decomposition window
        ax.axvspan(decomp_window[0], decomp_window[1], alpha=0.2, color='orange', label='Decomposition Window')
        
        # Mark cut site
        ax.axvline(x=cut_site, color='red', linestyle='--', linewidth=2, label='Expected Cut Site')
        
        # Add gRNA annotation
        grna_start = cut_site - 20  # gRNA typically 20bp upstream of cut
        grna_end = cut_site - 3  # PAM is 3bp
        ax.axhspan(0.9, 1.0, xmin=(grna_start-display_start)/(display_end-display_start),
                  xmax=(grna_end-display_start)/(display_end-display_start), 
                  color='purple', alpha=0.3)
        ax.text((grna_start + grna_end) / 2, 0.95, f"gRNA: {grna_seq}", 
               ha='center', va='center', fontsize=8, color='purple')
    
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Signal Intensity', fontsize=12)
    ax.set_title('Input: Alignment and Decomposition Windows', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.1)


def plot_alignment_window(ax, control_seq, edited_seq, control_traces, edited_traces,
                         cut_site, alignment_window):
    """Plot detailed view of alignment window."""
    
    trace_factor = 10
    start_bp, end_bp = alignment_window
    
    # Convert to trace positions
    start_trace = start_bp * trace_factor
    end_trace = end_bp * trace_factor
    
    # Ensure within bounds
    end_trace = min(end_trace, len(control_traces['A']), len(edited_traces['A']))
    
    if end_trace > start_trace:
        x_vals = np.arange(end_trace - start_trace)
        x_vals_bp = start_bp + (x_vals / trace_factor)
        
        # Plot control and edited traces
        for traces, color, label in [(control_traces, 'blue', 'Control'), 
                                     (edited_traces, 'red', 'Edited')]:
            signal = np.zeros(len(x_vals))
            for base in ['A', 'C', 'G', 'T']:
                signal += traces[base][start_trace:end_trace]
            
            if np.max(signal) > 0:
                signal = signal / np.max(signal)
            
            ax.plot(x_vals_bp, signal, color=color, alpha=0.7, linewidth=1.5, label=label)
        
        # Add sequence alignment indicator at bottom
        seq_start = max(0, start_bp)
        seq_end = min(len(control_seq), len(edited_seq), end_bp)
        
        matches = []
        for i in range(seq_start, seq_end):
            if control_seq[i] == edited_seq[i]:
                matches.append(1)
            else:
                matches.append(0)
        
        if matches:
            match_x = np.linspace(start_bp, end_bp, len(matches))
            ax2 = ax.twinx()
            ax2.bar(match_x, matches, width=(end_bp-start_bp)/len(matches), 
                   alpha=0.3, color='green', label='Sequence Match')
            ax2.set_ylim(0, 2)
            ax2.set_ylabel('Sequence Match', fontsize=10)
            ax2.set_yticks([0, 1])
            ax2.set_yticklabels(['Mismatch', 'Match'])
    
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Signal Intensity', fontsize=12)
    ax.set_title('Alignment Window (Before Cut Site)', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.1)


def plot_decomposition_window(ax, control_seq, edited_seq, control_traces, edited_traces,
                            cut_site, decomp_window):
    """Plot detailed view of decomposition window."""
    
    trace_factor = 10
    start_bp, end_bp = decomp_window
    
    # Convert to trace positions
    start_trace = start_bp * trace_factor
    end_trace = end_bp * trace_factor
    
    # Ensure within bounds
    end_trace = min(end_trace, len(control_traces['A']), len(edited_traces['A']))
    
    if end_trace > start_trace:
        x_vals = np.arange(end_trace - start_trace)
        x_vals_bp = start_bp + (x_vals / trace_factor)
        
        # Plot control and edited traces
        control_signal = np.zeros(len(x_vals))
        edited_signal = np.zeros(len(x_vals))
        
        for base in ['A', 'C', 'G', 'T']:
            control_signal += control_traces[base][start_trace:end_trace]
            edited_signal += edited_traces[base][start_trace:end_trace]
        
        # Normalize
        if np.max(control_signal) > 0:
            control_signal = control_signal / np.max(control_signal)
        if np.max(edited_signal) > 0:
            edited_signal = edited_signal / np.max(edited_signal)
        
        ax.plot(x_vals_bp, control_signal, 'b-', alpha=0.7, linewidth=1.5, label='Control')
        ax.plot(x_vals_bp, edited_signal, 'r-', alpha=0.7, linewidth=1.5, label='Edited')
        
        # Highlight differences
        diff_signal = np.abs(control_signal - edited_signal)
        ax.fill_between(x_vals_bp, 0, diff_signal, alpha=0.3, color='yellow', label='Difference')
    
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_ylabel('Signal Intensity', fontsize=12)
    ax.set_title('Decomposition Window (After Cut Site)', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.1)


def plot_indel_spectrum_tide(ax, analysis_results):
    """Plot indel spectrum similar to TIDE web tool."""
    
    indel_spectrum = analysis_results.get('indel_spectrum', {})
    
    # Check for error conditions
    if analysis_results.get('confidence', '').startswith('ERROR'):
        ax.text(0.5, 0.5, 'No editing detected\n\nPossible reasons:\n• No indels present\n• Poor sequence quality\n• Wrong cut site location', 
               ha='center', va='center', transform=ax.transAxes, fontsize=12,
               bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcoral", alpha=0.5))
        ax.set_title('Indel Spectrum - No Signal', fontsize=14, fontweight='bold')
        ax.set_xlim(-11, 11)
        ax.set_ylim(0, 1)
        return
    
    if indel_spectrum:
        # Sort by indel size
        sizes = sorted(indel_spectrum.keys())
        frequencies = [indel_spectrum[s] for s in sizes]
        
        # Create color map (red for insertions, blue for deletions)
        colors = ['red' if s < 0 else 'blue' for s in sizes]
        
        # Create bar plot
        bars = ax.bar(sizes, frequencies, color=colors, alpha=0.7, edgecolor='black', linewidth=1)
        
        # Highlight dominant indel
        dominant_size = analysis_results.get('dominant_indel_size', None)
        if dominant_size in sizes:
            idx = sizes.index(dominant_size)
            bars[idx].set_edgecolor('gold')
            bars[idx].set_linewidth(3)
        
        ax.set_xlabel('Indel Size (bp)', fontsize=12)
        ax.set_ylabel('Frequency (%)', fontsize=12)
        ax.set_title('Indel Spectrum', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add legend
        red_patch = mpatches.Patch(color='red', alpha=0.7, label='Insertions')
        blue_patch = mpatches.Patch(color='blue', alpha=0.7, label='Deletions')
        ax.legend(handles=[red_patch, blue_patch], loc='upper right')
        
        # Set x-axis limits
        ax.set_xlim(-11, 11)
        ax.set_xticks(range(-10, 11, 2))
    else:
        # No indel spectrum available
        ax.text(0.5, 0.5, 'Indel spectrum not available\n(Fallback method used)', 
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
    
    ax.set_title('Indel Spectrum', fontsize=14, fontweight='bold')


def plot_efficiency_pie_tide(ax, analysis_results):
    """Plot editing efficiency pie chart similar to TIDE web tool."""
    
    efficiency = analysis_results['editing_efficiency']
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # For clonal lines, we might have:
    # - Wild-type fraction
    # - Various indel fractions
    # - Unaccounted signal (noise/complex indels)
    
    # Create pie chart with TIDE-style colors
    colors = ['#ff4444', '#4444ff']  # Red for edited, blue for wild-type
    wedges, texts, autotexts = ax.pie([efficiency, wt_fraction], 
                                      labels=['Edited', f'Wild-type'],
                                      colors=colors,
                                      autopct='%1.1f%%',
                                      startangle=90,
                                      explode=(0.05, 0))  # Slightly separate edited slice
    
    # Enhance text
    for text in texts:
        text.set_fontsize(12)
        text.set_fontweight('bold')
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(11)
        autotext.set_fontweight('bold')
    
    # Add title with efficiency
    title = f"Editing Efficiency: {efficiency}%"
    if analysis_results.get('dominant_indel_size') not in [None, 'Unknown']:
        title += f"\nDominant: {analysis_results['dominant_indel_size']}bp ({analysis_results['dominant_indel_percent']}%)"
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)


def create_summary_report(all_results, output_dir, gene_name, timestamp):
    """
    Create a summary report for all analyzed samples.
    
    Args:
        all_results: List of analysis results for all samples
        output_dir: Output directory
        gene_name: Gene name
        timestamp: Timestamp
        
    Returns:
        str: Path to summary report
    """
    # Create summary figure
    n_samples = len(all_results)
    n_cols = min(4, n_samples)
    n_rows = (n_samples + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    if n_samples == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    fig.suptitle(f"TIDE Analysis Summary - {gene_name.upper()}\n{n_samples} Clonal Cell Lines", 
                fontsize=16, fontweight='bold')
    
    # Plot efficiency pie charts for each sample
    for idx, (result, ax) in enumerate(zip(all_results, axes.flat)):
        plot_efficiency_pie_tide(ax, result)
        ax.set_title(f"{result['sample_name']}\nEfficiency: {result['editing_efficiency']}%", 
                    fontsize=12, fontweight='bold')
    
    # Hide unused subplots
    for idx in range(n_samples, n_rows * n_cols):
        axes.flat[idx].axis('off')
    
    # Add summary statistics
    avg_efficiency = np.mean([r['editing_efficiency'] for r in all_results])
    std_efficiency = np.std([r['editing_efficiency'] for r in all_results])
    
    # Count different efficiency categories
    homozygous = sum(1 for r in all_results if r['editing_efficiency'] > 80)
    heterozygous = sum(1 for r in all_results if 40 <= r['editing_efficiency'] <= 60)
    low_edit = sum(1 for r in all_results if 0 < r['editing_efficiency'] < 40)
    no_edit = sum(1 for r in all_results if r['editing_efficiency'] == 0)
    
    summary_text = f"Average Efficiency: {avg_efficiency:.1f}% ± {std_efficiency:.1f}%\n"
    summary_text += f"Range: {min(r['editing_efficiency'] for r in all_results):.1f}% - "
    summary_text += f"{max(r['editing_efficiency'] for r in all_results):.1f}%\n"
    summary_text += f"Likely homozygous (>80%): {homozygous} | Heterozygous (40-60%): {heterozygous} | "
    summary_text += f"Low editing (<40%): {low_edit} | No editing: {no_edit}"
    
    fig.text(0.5, 0.02, summary_text, ha='center', fontsize=12, 
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.5))
    
    plt.tight_layout()
    
    # Save summary
    summary_path = os.path.join(output_dir, f"tide_summary_{gene_name}_{timestamp}.png")
    plt.savefig(summary_path, dpi=150, bbox_inches='tight')  # Reduced DPI
    plt.close('all')  # Close all figures
    
    # Force cleanup
    import gc
    gc.collect()
    
    return summary_path 