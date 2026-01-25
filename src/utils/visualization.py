"""
Visualization Module for Indel Analysis
Handles plotting trace decomposition analysis graphs with validation
"""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime


class ValidationIssue:
    """Class to track validation issues"""
    def __init__(self, severity, message):
        self.severity = severity  # 'error', 'warning', 'info'
        self.message = message


def validate_for_plotting(control_seq, edited_seq, expected_cut_site, similarity, args):
    """
    Validate sequences for plotting based on criteria and user arguments.
    
    Args:
        control_seq: Control sequence
        edited_seq: Edited sequence
        expected_cut_site: Expected cut site position
        similarity: Calculated similarity percentage
        args: Command line arguments
        
    Returns:
        tuple: (should_plot, validation_issues)
    """
    validation_issues = []
    
    # Check cut site within bounds
    min_seq_len = min(len(control_seq), len(edited_seq))
    if expected_cut_site and expected_cut_site >= min_seq_len:
        validation_issues.append(ValidationIssue(
            'error',
            f"Expected cut site ({expected_cut_site}) is beyond sequence length ({min_seq_len})"
        ))
    
    # Check similarity
    if similarity < args.similarity_threshold:
        validation_issues.append(ValidationIssue(
            'error' if similarity < 50 else 'warning',
            f"Low similarity between sequences ({similarity:.1f}% < {args.similarity_threshold}% threshold)"
        ))
    
    # Determine if we should plot
    has_errors = any(issue.severity == 'error' for issue in validation_issues)
    should_plot = args.force_plot or not has_errors
    
    return should_plot, validation_issues


def plot_indel_analysis(control_file, edited_file, output_dir, gene_name,
                       grna_info, efficiency_data, timestamp=None,
                       validation_issues=None, expected_cut_site=None):
    """
    Create indel analysis plot with validation warnings.
    
    Args:
        control_file: Path to control AB1
        edited_file: Path to edited AB1
        output_dir: Output directory
        gene_name: Gene name
        grna_info: gRNA information dict
        efficiency_data: Indel analysis results
        timestamp: Optional timestamp
        validation_issues: Optional list of ValidationIssue objects
        expected_cut_site: Optional expected cut site position
    """
    from .ab1_parser import parse_ab1
    from .sequence_analysis import find_divergence_point
    
    # Initialize validation_issues if not provided
    if validation_issues is None:
        validation_issues = []
    
    # Get expected_cut_site from grna_info or efficiency_data if not provided
    if expected_cut_site is None:
        expected_cut_site = grna_info.get('cut_site') if grna_info else None
        if expected_cut_site is None:
            expected_cut_site = efficiency_data.get('cut_site')
    
    # Parse AB1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Title with warnings if present
    title = f"Indel Analysis - {gene_name.upper()}"
    if validation_issues:
        title += " [ISSUES DETECTED]"
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    # Add validation warnings box if present
    if validation_issues:
        warnings_text = "VALIDATION ISSUES:\n"
        for issue in validation_issues:
            symbol = "⚠️" if issue.severity == 'warning' else "❌"
            warnings_text += f"{symbol} {issue.message}\n"
        
        # Add text box with warnings
        props = dict(boxstyle='round', facecolor='yellow' if any(i.severity == 'warning' for i in validation_issues) else 'lightcoral', alpha=0.8)
        fig.text(0.5, 0.94, warnings_text.strip(), transform=fig.transFigure, 
                fontsize=10, verticalalignment='top', horizontalalignment='center', bbox=props)
    
    # Subplot 1: Chromatogram overlay
    plot_chromatogram_overlay(ax1, control_traces, edited_traces, expected_cut_site, "Chromatogram Overlay")
    
    # Subplot 2: Sequence alignment visualization
    divergence_point = find_divergence_point(control_seq, edited_seq)
    plot_sequence_alignment(ax2, control_seq, edited_seq, divergence_point, expected_cut_site, "Sequence Alignment")
    
    # Subplot 3: Trace decomposition or signal decay
    if efficiency_data.get('signal_decay_ratio') is not None:
        ax3.set_title('Signal Decay Analysis (Fallback Method)', fontsize=12)
    else:
        plot_signal_analysis(ax3, efficiency_data, "Signal Analysis")
    
    # Subplot 4: Editing efficiency pie chart
    plot_efficiency_pie(ax4, efficiency_data['editing_efficiency'], efficiency_data)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plot_filename = f"indel_analysis_{gene_name}_{timestamp}.png"
    if validation_issues and any(i.severity == 'error' for i in validation_issues):
        plot_filename = plot_filename.replace('.png', '_WARNING.png')
    
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_path


def plot_error_figure(output_dir, gene_name, validation_errors):
    """
    Create an error plot when analysis cannot proceed.
    
    Args:
        output_dir: Output directory
        gene_name: Gene name
        validation_errors: List of error messages
        
    Returns:
        str: Path to saved error plot
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    
    error_text = f"Indel Analysis Failed for {gene_name.upper()}\n\n"
    error_text += "\n".join(f"• {error}" for error in validation_errors)
    error_text += "\n\nPlease ensure:\n"
    error_text += "• Sequencing covers the gRNA target region\n"
    error_text += "• Control and edited samples are from the same amplicon\n"
    error_text += "• Sequence quality is sufficient\n"
    error_text += "\nUse --force-plot to generate plots anyway (not recommended)"
    
    ax.text(0.5, 0.5, error_text,
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=14,
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcoral", alpha=0.5))
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    # Save error plot
    error_plot_path = os.path.join(output_dir, f"indel_analysis_{gene_name}_{timestamp}_ERROR.png")
    plt.savefig(error_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return error_plot_path


def plot_chromatogram_overlay(ax, control_traces, edited_traces, cut_site, title):
    """Plot chromatogram overlay around cut site."""
    trace_factor = 12
    window = 100
    
    if cut_site and control_traces and 'A' in control_traces:
        # Check if cut site is within trace bounds
        max_trace_pos = len(control_traces['A']) // trace_factor
        if cut_site > max_trace_pos:
            # Cut site is beyond trace length
            ax.text(0.5, 0.5, f"Cut site ({cut_site}) is beyond trace length ({max_trace_pos} bp)",
                   ha='center', va='center', transform=ax.transAxes, fontsize=12, color='red')
            ax.set_title(title + " - OUT OF BOUNDS")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            return
        
        start = max(0, (cut_site - window // 2) * trace_factor)
        end = min(len(control_traces['A']), (cut_site + window // 2) * trace_factor)
        
        x_vals = np.arange(start, end)
        x_vals_norm = (x_vals - start) / trace_factor
        
        # Plot control in blue
        control_sum = np.zeros(end - start)
        for base in ['A', 'C', 'G', 'T']:
            control_sum += control_traces[base][start:end]
        ax.plot(x_vals_norm, control_sum, 'b-', alpha=0.7, label='Control')
        
        # Plot edited in red
        if end <= len(edited_traces['A']):
            edited_sum = np.zeros(end - start)
            for base in ['A', 'C', 'G', 'T']:
                edited_sum += edited_traces[base][start:end]
            ax.plot(x_vals_norm, edited_sum, 'r-', alpha=0.7, label='Edited')
        
        # Mark cut site
        cut_pos_norm = window // 2
        ax.axvline(x=cut_pos_norm, color='green', linestyle='--', label='Expected Cut Site')
    
    ax.set_title(title)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Signal Intensity')
    ax.legend()
    ax.grid(True, alpha=0.3)


def plot_sequence_alignment(ax, control_seq, edited_seq, divergence_point, cut_site, title):
    """Plot sequence alignment visualization."""
    window_size = 50
    
    if divergence_point >= 0:
        start = max(0, divergence_point - window_size)
        end = min(len(control_seq), divergence_point + window_size)
        
        # Create color map for matches/mismatches
        colors = []
        for i in range(start, end):
            if i < len(edited_seq) and control_seq[i] == edited_seq[i]:
                colors.append('lightgreen')
            else:
                colors.append('lightcoral')
        
        # Plot as color bars
        ax.bar(range(len(colors)), [1]*len(colors), color=colors, width=1)
        
        # Mark divergence point
        if start <= divergence_point < end:
            ax.axvline(x=divergence_point-start, color='red', linestyle='--', label='Divergence Point')
        
        # Mark cut site if in range
        if cut_site and start <= cut_site < end:
            ax.axvline(x=cut_site-start, color='green', linestyle='--', label='Expected Cut Site')
    
    ax.set_title(title)
    ax.set_xlabel('Position')
    ax.set_ylabel('Match')
    ax.set_ylim(0, 1.5)
    ax.legend()


def plot_indel_spectrum(ax, efficiency_data):
    """Plot indel spectrum from trace decomposition analysis."""
    if efficiency_data.get('indel_spectrum'):
        sizes = list(efficiency_data['indel_spectrum'].keys())
        frequencies = list(efficiency_data['indel_spectrum'].values())
        
        colors = ['red' if s < 0 else 'blue' for s in sizes]
        ax.bar(sizes, frequencies, color=colors, alpha=0.7)
        
        ax.set_title('Indel Spectrum')
        ax.set_xlabel('Indel Size (bp)')
        ax.set_ylabel('Frequency (%)')
        ax.grid(True, alpha=0.3)
        
        # Add legend
        red_patch = mpatches.Patch(color='red', label='Insertions')
        blue_patch = mpatches.Patch(color='blue', label='Deletions')
        ax.legend(handles=[red_patch, blue_patch])


def plot_signal_analysis(ax, efficiency_data, title):
    """Plot signal analysis metrics."""
    metrics = {
        'Editing\nEfficiency': efficiency_data['editing_efficiency'],
        'Quality\nScore': efficiency_data['quality_score'],
        'Sequence\nSimilarity': efficiency_data.get('sequence_similarity', 0)
    }
    
    x_pos = np.arange(len(metrics))
    values = list(metrics.values())
    
    bars = ax.bar(x_pos, values, alpha=0.7)
    
    # Color bars based on value
    for bar, value in zip(bars, values):
        if value < 30:
            bar.set_color('red')
        elif value < 70:
            bar.set_color('orange')
        else:
            bar.set_color('green')
    
    ax.set_title(title)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(metrics.keys())
    ax.set_ylabel('Percentage (%)')
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)


def plot_efficiency_pie(ax, efficiency, efficiency_data):
    """Plot editing efficiency as a pie chart."""
    edited = efficiency
    unedited = 100 - efficiency
    
    # Create pie chart
    wedges, texts, autotexts = ax.pie([edited, unedited], 
                                      labels=['Edited', 'Unedited'],
                                      colors=['#FF6B6B', '#4ECDC4'],
                                      autopct='%1.1f%%',
                                      startangle=90)
    
    # Add title with additional info
    title = f"Editing Efficiency: {efficiency}%"
    if efficiency_data.get('dominant_indel_size') != 'Unknown':
        title += f"\nDominant Indel: {efficiency_data['dominant_indel_size']}bp ({efficiency_data['dominant_indel_percent']}%)"
    
    confidence = efficiency_data.get('confidence', 'Unknown')
    title += f"\nConfidence: {confidence}"
    
    ax.set_title(title, fontsize=12) 