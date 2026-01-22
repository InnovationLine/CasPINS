"""
Visualization Module for Multiple Sample Indel Analysis
Generates analysis plots for trace decomposition results
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
from datetime import datetime


def ensure_numpy_traces(traces):
    """
    Ensure trace data is in numpy array format.
    
    AB1 files can return traces as tuples or other iterables.
    This function converts them to numpy arrays for consistent processing.
    
    Args:
        traces: Dictionary with keys A, C, G, T containing trace data
        
    Returns:
        Dictionary with numpy arrays for each base
    """
    converted = {}
    for base in ['A', 'C', 'G', 'T']:
        if base in traces:
            data = traces[base]
            if isinstance(data, np.ndarray):
                converted[base] = data
            elif isinstance(data, (tuple, list)):
                converted[base] = np.array(data, dtype=np.float64)
            else:
                # Try to convert whatever it is
                try:
                    converted[base] = np.array(data, dtype=np.float64)
                except Exception:
                    # Create empty array as fallback
                    converted[base] = np.array([], dtype=np.float64)
        else:
            converted[base] = np.array([], dtype=np.float64)
    return converted


def plot_indel_analysis_multi(control_file, edited_file, output_dir, gene_name,
                             sample_name, analysis_results, grna_sequences, timestamp):
    """
    Create indel analysis plot for a single edited sample.
    
    Args:
        control_file: Path to control AB1 file
        edited_file: Path to edited AB1 file
        output_dir: Output directory for plots
        gene_name: Gene name
        sample_name: Sample name (e.g., editedA1)
        analysis_results: Indel analysis results
        grna_sequences: List of gRNA sequences
        timestamp: Timestamp string
        
    Returns:
        str: Path to saved plot
    """
    from .ab1_parser import parse_ab1
    
    # Parse AB1 files
    control_seq, control_traces_raw = parse_ab1(control_file)
    edited_seq, edited_traces_raw = parse_ab1(edited_file)
    
    # Ensure traces are numpy arrays (AB1 files can return tuples)
    control_traces = ensure_numpy_traces(control_traces_raw)
    edited_traces = ensure_numpy_traces(edited_traces_raw)
    
    # Create figure with analysis layout
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.5, 1.5, 1], hspace=0.3, wspace=0.3)
    
    # Add main title
    fig.suptitle(f"Indel Analysis - {gene_name.upper()} - {sample_name}", fontsize=18, fontweight='bold')
    
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
    plot_indel_spectrum(ax_indel, analysis_results)
    
    # Bottom right: Efficiency pie chart
    ax_pie = fig.add_subplot(gs[2, 1])
    plot_efficiency_pie(ax_pie, analysis_results)
    
    # Add analysis info text
    info_text = f"Editing Efficiency: {analysis_results['editing_efficiency']}%\n"
    info_text += f"Wild-type: {analysis_results.get('wt_fraction', 100 - analysis_results['editing_efficiency']):.1f}%\n"
    info_text += f"Quality Score: {analysis_results.get('quality_score', 0):.1f}%\n"
    info_text += f"Confidence: {analysis_results.get('confidence', 'Unknown')}\n"
    info_text += f"Analysis Method: {analysis_results.get('method', 'Trace Decomposition')}"
    
    fig.text(0.02, 0.02, info_text, fontsize=10, bbox=dict(boxstyle="round,pad=0.3", 
                                                           facecolor="lightgray", alpha=0.5))
    
    # Save plot
    plot_filename = f"indel_analysis_{gene_name}_{sample_name}_{timestamp}.png"
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
    """Plot the input visualization showing alignment and decomposition windows with ATCG colors."""
    
    # Ensure traces are numpy arrays (in case they were passed as tuples)
    control_traces = ensure_numpy_traces(control_traces)
    edited_traces = ensure_numpy_traces(edited_traces)
    
    # ATCG color scheme
    BASE_COLORS = {'A': '#00AA00', 'T': '#FF0000', 'G': '#333333', 'C': '#0000FF'}
    
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
        
        # Ensure within bounds - check if trace data exists
        if len(control_traces['A']) > 0:
            trace_end = min(trace_end, len(control_traces['A']))
        else:
            # No trace data available
            ax.text(0.5, 0.5, 'No trace data available', ha='center', va='center', transform=ax.transAxes)
            return
        
        x_vals = np.arange(trace_start, trace_end)
        x_vals_bp = x_vals / trace_factor
        
        # Plot individual base traces with ATCG colors
        if len(x_vals) > 0:
            max_signal = 1
            for base in ['A', 'T', 'G', 'C']:
                if len(control_traces[base]) > 0 and trace_end <= len(control_traces[base]):
                    signal = control_traces[base][trace_start:trace_end].astype(float)
                    if len(signal) > 0:
                        max_signal = max(max_signal, np.max(signal))
            
            # Plot each base channel with its color
            for base in ['A', 'T', 'G', 'C']:
                if len(control_traces[base]) > 0 and trace_end <= len(control_traces[base]):
                    signal = control_traces[base][trace_start:trace_end].astype(float)
                    if len(signal) > 0:
                        signal = signal / max_signal if max_signal > 0 else signal
                        ax.plot(x_vals_bp, signal, color=BASE_COLORS[base], alpha=0.6, 
                               linewidth=0.6, label=base)
            
            # Add sequence text at top with colors
            seq_start = max(0, display_start)
            seq_end = min(len(control_seq), display_end)
            step = max(1, (seq_end - seq_start) // 60)  # Show ~60 bases max
            
            for i, pos in enumerate(range(seq_start, min(seq_end, seq_start + 60), step)):
                if pos < len(control_seq):
                    base = control_seq[pos]
                    x_pos = pos
                    ax.text(x_pos, 1.02, base, fontsize=5, ha='center', va='bottom',
                           color=BASE_COLORS.get(base, '#888888'), fontweight='bold',
                           fontfamily='monospace')
        else:
            # Handle case when x_vals is empty
            ax.text(0.5, 0.5, 'No data available in this range', 
                   ha='center', va='center', transform=ax.transAxes)
        
        # Mark alignment window
        ax.axvspan(alignment_window[0], alignment_window[1], alpha=0.15, color='green', label='Alignment Window')
        
        # Mark decomposition window
        ax.axvspan(decomp_window[0], decomp_window[1], alpha=0.15, color='orange', label='Decomp. Window')
        
        # Mark cut site
        ax.axvline(x=cut_site, color='red', linestyle='--', linewidth=2, label='Cut Site')
        
        # Add gRNA annotation
        grna_start = cut_site - 20  # gRNA typically 20bp upstream of cut
        grna_end = cut_site - 3  # PAM is 3bp
        ax.axhspan(0.92, 1.0, xmin=(grna_start-display_start)/(display_end-display_start),
                  xmax=(grna_end-display_start)/(display_end-display_start), 
                  color='purple', alpha=0.3)
        ax.text((grna_start + grna_end) / 2, 0.96, f"gRNA: {grna_seq[:20]}...", 
               ha='center', va='center', fontsize=7, color='purple')
    
    ax.set_xlabel('Position (bp)', fontsize=10)
    ax.set_ylabel('Signal Intensity', fontsize=10)
    ax.set_title('Input Overview (A=Green, T=Red, G=Black, C=Blue)', fontsize=11, fontweight='bold')
    ax.legend(loc='lower right', fontsize=7, ncol=4)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.1)


def plot_alignment_window(ax, control_seq, edited_seq, control_traces, edited_traces,
                         cut_site, alignment_window):
    """Plot detailed view of alignment window with ATCG color coding."""
    
    # ATCG color scheme
    BASE_COLORS = {'A': '#00AA00', 'T': '#FF0000', 'G': '#333333', 'C': '#0000FF'}
    
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
        
        # Plot individual base traces with ATCG colors for CONTROL
        max_signal = 1
        for base in ['A', 'T', 'G', 'C']:
            signal = control_traces[base][start_trace:end_trace].astype(float)
            max_signal = max(max_signal, np.max(signal))
        
        for base in ['A', 'T', 'G', 'C']:
            signal = control_traces[base][start_trace:end_trace].astype(float)
            signal = signal / max_signal if max_signal > 0 else signal
            ax.plot(x_vals_bp, signal, color=BASE_COLORS[base], alpha=0.7, 
                   linewidth=0.8, label=f'{base}')
        
        # Add sequence text at top with colors
        seq_start = max(0, start_bp)
        seq_end = min(len(control_seq), end_bp)
        seq_display_start = int(len(x_vals_bp) * 0.1)
        seq_display_end = int(len(x_vals_bp) * 0.9)
        step = max(1, (seq_end - seq_start) // 40)  # Show ~40 bases max
        
        for i, pos in enumerate(range(seq_start, min(seq_end, seq_start + 40), step)):
            if pos < len(control_seq):
                base = control_seq[pos]
                x_pos = start_bp + (pos - seq_start)
                ax.text(x_pos, 1.05, base, fontsize=7, ha='center', va='bottom',
                       color=BASE_COLORS.get(base, '#888888'), fontweight='bold',
                       fontfamily='monospace')
    
    ax.set_xlabel('Position (bp)', fontsize=10)
    ax.set_ylabel('Signal Intensity', fontsize=10)
    ax.set_title('Alignment Window - Control (A=Green, T=Red, G=Black, C=Blue)', 
                fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=7, ncol=4)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.2)


def plot_decomposition_window(ax, control_seq, edited_seq, control_traces, edited_traces,
                            cut_site, decomp_window):
    """Plot detailed view of decomposition window with ATCG colors for edited sample."""
    
    # Ensure traces are numpy arrays
    control_traces = ensure_numpy_traces(control_traces)
    edited_traces = ensure_numpy_traces(edited_traces)

    # ATCG color scheme
    BASE_COLORS = {'A': '#00AA00', 'T': '#FF0000', 'G': '#333333', 'C': '#0000FF'}

    trace_factor = 10
    start_bp, end_bp = decomp_window

    # Convert to trace positions
    start_trace = start_bp * trace_factor
    end_trace = end_bp * trace_factor

    # Ensure within bounds - check if traces have data
    if len(control_traces['A']) == 0 or len(edited_traces['A']) == 0:
        ax.text(0.5, 0.5, 'No trace data available', ha='center', va='center', transform=ax.transAxes)
        return
        
    end_trace = min(end_trace, len(control_traces['A']), len(edited_traces['A']))

    if end_trace > start_trace:
        x_vals = np.arange(end_trace - start_trace)
        x_vals_bp = start_bp + (x_vals / trace_factor)

        # Plot individual base traces with ATCG colors for EDITED sample
        max_signal = 1
        for base in ['A', 'T', 'G', 'C']:
            if len(edited_traces[base]) > 0 and end_trace <= len(edited_traces[base]):
                signal = edited_traces[base][start_trace:end_trace].astype(float)
                if len(signal) > 0:
                    max_signal = max(max_signal, np.max(signal))

        for base in ['A', 'T', 'G', 'C']:
            if len(edited_traces[base]) > 0 and end_trace <= len(edited_traces[base]):
                signal = edited_traces[base][start_trace:end_trace].astype(float)
                if len(signal) > 0:
                    signal = signal / max_signal if max_signal > 0 else signal
                    ax.plot(x_vals_bp, signal, color=BASE_COLORS[base], alpha=0.7,
                           linewidth=0.8, label=f'{base}')
        
        # Add sequence text at top with colors for EDITED sample
        seq_start = max(0, start_bp)
        seq_end = min(len(edited_seq), end_bp)
        step = max(1, (seq_end - seq_start) // 40)
        
        for i, pos in enumerate(range(seq_start, min(seq_end, seq_start + 40), step)):
            if pos < len(edited_seq):
                base = edited_seq[pos]
                x_pos = start_bp + (pos - seq_start)
                ax.text(x_pos, 1.05, base, fontsize=7, ha='center', va='bottom',
                       color=BASE_COLORS.get(base, '#888888'), fontweight='bold',
                       fontfamily='monospace')
    
    ax.set_xlabel('Position (bp)', fontsize=10)
    ax.set_ylabel('Signal Intensity', fontsize=10)
    ax.set_title('Decomposition Window - Edited (A=Green, T=Red, G=Black, C=Blue)', 
                fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=7, ncol=4)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.2)


def plot_indel_spectrum(ax, analysis_results):
    """Plot indel spectrum from trace decomposition with TIDE-style colors."""
    
    indel_spectrum = analysis_results.get('indel_spectrum', {})
    efficiency = analysis_results.get('editing_efficiency', 0)
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # TIDE-style colors
    DELETION_COLOR = '#DC143C'   # Crimson red
    WILDTYPE_COLOR = '#228B22'   # Forest green
    INSERTION_COLOR = '#4169E1'  # Royal blue
    
    # Check for error conditions
    if analysis_results.get('confidence', '').startswith('ERROR'):
        ax.text(0.5, 0.5, 'No editing detected\n\nPossible reasons:\n• No indels present\n• Poor sequence quality\n• Wrong cut site location', 
               ha='center', va='center', transform=ax.transAxes, fontsize=12,
               bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcoral", alpha=0.5))
        ax.set_title('Indel Spectrum - No Signal', fontsize=14, fontweight='bold')
        ax.set_xlim(-11, 11)
        ax.set_ylim(0, 1)
        return
    
    # Always show full spectrum from -10 to +10
    all_sizes = list(range(-10, 11))
    frequencies = []
    colors = []
    
    for size in all_sizes:
        if size == 0:
            # Wild-type frequency
            frequencies.append(wt_fraction)
            colors.append(WILDTYPE_COLOR)
        else:
            # Get frequency from spectrum, handle both int and str keys
            freq = indel_spectrum.get(size, indel_spectrum.get(str(size), 0))
            frequencies.append(freq)
            if size < 0:
                colors.append(DELETION_COLOR)  # Red for deletions
            else:
                colors.append(INSERTION_COLOR)  # Blue for insertions
    
    # Create bar plot
    bars = ax.bar(all_sizes, frequencies, color=colors, alpha=0.85, edgecolor='black', linewidth=0.5)
    
    # Add value labels for significant bars
    for i, (size, freq) in enumerate(zip(all_sizes, frequencies)):
        if freq > 3:
            ax.text(size, freq + 1.5, f'{freq:.0f}%', ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    # Highlight dominant indel
    dominant_size = analysis_results.get('dominant_indel_size', None)
    if dominant_size is not None and dominant_size in all_sizes:
        idx = all_sizes.index(dominant_size)
        bars[idx].set_edgecolor('gold')
        bars[idx].set_linewidth(2)
    
    ax.set_xlabel('Indel Size (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency (%)', fontsize=12, fontweight='bold')
    ax.set_title('Indel Spectrum', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add legend with TIDE colors
    del_patch = mpatches.Patch(color=DELETION_COLOR, alpha=0.85, label='Deletions')
    wt_patch = mpatches.Patch(color=WILDTYPE_COLOR, alpha=0.85, label='Wild-type (0)')
    ins_patch = mpatches.Patch(color=INSERTION_COLOR, alpha=0.85, label='Insertions')
    ax.legend(handles=[del_patch, wt_patch, ins_patch], loc='upper right', fontsize=9)
    
    # Set x-axis limits
    ax.set_xlim(-11, 11)
    ax.set_xticks(range(-10, 11, 2))
    
    # Set y-axis dynamically
    max_freq = max(frequencies) if frequencies else 100
    ax.set_ylim(0, min(100, max_freq * 1.2 + 5))


def plot_efficiency_pie(ax, analysis_results):
    """Plot editing efficiency pie chart with TIDE-style colors."""
    
    efficiency = analysis_results['editing_efficiency']
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # Ensure values are valid
    if efficiency < 0:
        efficiency = 0
    if wt_fraction < 0:
        wt_fraction = 100 - efficiency
    
    # Normalize to 100%
    total = efficiency + wt_fraction
    if total > 0:
        efficiency = (efficiency / total) * 100
        wt_fraction = (wt_fraction / total) * 100
    
    # TIDE-style colors: Red for edited, Green for wild-type
    colors = ['#DC143C', '#228B22']  # Crimson red, Forest green
    
    wedges, texts, autotexts = ax.pie([efficiency, wt_fraction], 
                                      labels=['Edited', 'Wild-type'],
                                      colors=colors,
                                      autopct='%1.1f%%',
                                      startangle=90,
                                      explode=(0.02, 0),
                                      wedgeprops={'edgecolor': 'white', 'linewidth': 1})
    
    # Enhance text
    for text in texts:
        text.set_fontsize(11)
        text.set_fontweight('bold')
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(10)
        autotext.set_fontweight('bold')
    
    # Add title with efficiency
    title = f"Editing Efficiency: {efficiency:.1f}%"
    dominant_size = analysis_results.get('dominant_indel_size')
    dominant_pct = analysis_results.get('dominant_indel_percent', 0)
    if dominant_size not in [None, 'Unknown', 0, 'N/A']:
        title += f"\nDominant: {dominant_size}bp ({dominant_pct}%)"
    
    ax.set_title(title, fontsize=13, fontweight='bold', pad=15)


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
    
    # Main title
    fig.suptitle(f"Indel Analysis Summary - {gene_name.upper()}\n{n_samples} Clonal Cell Lines", 
                fontsize=16, fontweight='bold')
    
    for i, result in enumerate(all_results):
        ax = plt.subplot(n_rows, n_cols, i+1)
        plot_efficiency_pie(ax, result)
        sample_name = result.get('sample_name', f'Sample {i+1}')
        ax.set_title(sample_name, fontsize=12)
        
        # Add efficiency text
        eff = result.get('editing_efficiency', 0)
        ax.text(0.5, -0.2, f"{eff}% edited", 
                transform=ax.transAxes, ha='center', fontsize=10)
    
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
    
    # Save summary plot
    summary_path = os.path.join(output_dir, f"indel_analysis_summary_{gene_name}_{timestamp}.png")
    plt.savefig(summary_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Force cleanup
    import gc
    gc.collect()
    
    return summary_path 