"""
TIDE-Style Visualization Module
Implements visualization similar to TIDE Analysis (https://tide.nki.nl/)

This module creates publication-ready figures with:
1. Sequence alignment view with base calls (A, T, C, G)
2. Chromatogram traces with individual base channels
3. Indel spectrum in TIDE format
4. Decomposition quality metrics
5. Streamlit-compatible display components
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from typing import Dict, List, Tuple, Optional
import io
import base64


# TIDE-style color scheme (matching actual TIDE tool)
TIDE_COLORS = {
    # Base colors for sequence display
    'A': '#00AA00',      # Green (darker for visibility)
    'T': '#FF0000',      # Red
    'G': '#000000',      # Black
    'C': '#0000FF',      # Blue
    # Indel spectrum colors (matching TIDE)
    'deletion': '#DC143C',    # Crimson red for deletions
    'insertion': '#4169E1',   # Royal blue for insertions  
    'wildtype': '#228B22',    # Forest green for wild-type (TIDE uses green)
    # Pie chart colors (matching TIDE)
    'edited': '#DC143C',      # Red for edited
    'unedited': '#228B22',    # Green for wild-type/unedited
    # Background
    'background': '#F8F9FA',
    'grid': '#E0E0E0'
}


def create_tide_style_figure(control_seq: str, edited_seq: str,
                              control_traces: Dict, edited_traces: Dict,
                              cut_site: int, analysis_results: Dict,
                              grna_sequence: str, sample_name: str = "Sample") -> plt.Figure:
    """
    Create a complete TIDE-style analysis figure.
    
    Args:
        control_seq: Control sequence string
        edited_seq: Edited sequence string
        control_traces: Dictionary of control chromatogram traces
        edited_traces: Dictionary of edited chromatogram traces
        cut_site: Expected cut site position
        analysis_results: Indel analysis results dictionary
        grna_sequence: Guide RNA sequence
        sample_name: Name of the sample
        
    Returns:
        matplotlib Figure object
    """
    # Validate inputs
    if not control_seq or len(control_seq) < 10:
        control_seq = "N" * 200  # Placeholder
    if not edited_seq or len(edited_seq) < 10:
        edited_seq = control_seq  # Use control as fallback
    
    # Ensure cut_site is valid
    if not cut_site or cut_site < 0 or cut_site > len(control_seq):
        cut_site = min(100, len(control_seq) // 2)
    
    # Create figure with TIDE-style layout
    fig = plt.figure(figsize=(16, 12), facecolor='white')
    
    # Define grid layout similar to TIDE
    gs = gridspec.GridSpec(4, 2, figure=fig, 
                          height_ratios=[1.0, 1.2, 1.0, 0.8],
                          hspace=0.35, wspace=0.25)
    
    # Title
    fig.suptitle(f'CRISPR Indel Analysis - {sample_name}', 
                fontsize=14, fontweight='bold', y=0.98)
    
    # Panel 1: Sequence Alignment View (Top full width)
    ax_align = fig.add_subplot(gs[0, :])
    _plot_sequence_alignment_matplotlib(ax_align, control_seq, edited_seq, cut_site, grna_sequence)
    
    # Panel 2: Control Chromatogram (Left)
    ax_control = fig.add_subplot(gs[1, 0])
    _plot_chromatogram_matplotlib(ax_control, control_traces, control_seq, cut_site, 
                                 title="Control Sample", window=80)
    
    # Panel 3: Edited Chromatogram (Right)
    ax_edited = fig.add_subplot(gs[1, 1])
    _plot_chromatogram_matplotlib(ax_edited, edited_traces, edited_seq, cut_site,
                                title="Edited Sample", window=80)
    
    # Panel 4: Indel Spectrum (Bottom Left) - TIDE style
    ax_spectrum = fig.add_subplot(gs[2, 0])
    _plot_indel_spectrum_matplotlib(ax_spectrum, analysis_results)
    
    # Panel 5: Efficiency Pie Chart (Bottom Right)
    ax_pie = fig.add_subplot(gs[2, 1])
    _plot_efficiency_pie_matplotlib(ax_pie, analysis_results)
    
    # Panel 6: Summary Statistics (Full width bottom)
    ax_summary = fig.add_subplot(gs[3, :])
    _plot_summary_matplotlib(ax_summary, analysis_results, grna_sequence, sample_name)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    return fig


def _plot_sequence_alignment_matplotlib(ax: plt.Axes, control_seq: str, edited_seq: str,
                                         cut_site: int, grna_sequence: str):
    """
    Plot sequence alignment view using matplotlib (not HTML).
    Shows ATCG bases with colors similar to TIDE.
    """
    ax.set_facecolor('#FFFFFF')
    
    # Calculate display window
    window_before = 25
    window_after = 40
    start = max(0, cut_site - window_before)
    end = min(len(control_seq), cut_site + window_after)
    
    # Prepare sequences
    ctrl_segment = control_seq[start:end] if len(control_seq) > start else "N" * (end - start)
    edit_segment = edited_seq[start:end] if len(edited_seq) > start else ctrl_segment
    
    # Ensure we have content
    if len(ctrl_segment) == 0:
        ctrl_segment = "N" * 50
    if len(edit_segment) == 0:
        edit_segment = ctrl_segment
    
    # Ensure same length
    max_len = max(len(ctrl_segment), len(edit_segment))
    ctrl_segment = ctrl_segment.ljust(max_len, '-')
    edit_segment = edit_segment.ljust(max_len, '-')
    
    # Position spacing
    base_width = 0.8
    y_control = 2.0
    y_edited = 1.0
    y_match = 0.3
    
    # Plot position numbers (every 10 bases)
    for i in range(0, len(ctrl_segment), 10):
        pos = start + i
        x = i * base_width
        ax.text(x, 3.0, str(pos), fontsize=7, ha='center', va='bottom', fontfamily='monospace')
    
    # Plot control sequence with colored bases
    for i, base in enumerate(ctrl_segment):
        x = i * base_width
        color = TIDE_COLORS.get(base.upper(), '#888888')
        ax.text(x, y_control, base, fontsize=9, ha='center', va='center',
               color=color, fontweight='bold', fontfamily='monospace')
    
    # Plot edited sequence with colored bases
    for i, base in enumerate(edit_segment):
        x = i * base_width
        color = TIDE_COLORS.get(base.upper(), '#888888')
        
        # Highlight mismatches with background
        if i < len(ctrl_segment) and base.upper() != ctrl_segment[i].upper():
            rect = plt.Rectangle((x - 0.35, y_edited - 0.3), 0.7, 0.6, 
                                 facecolor='#FFFF99', alpha=0.7, edgecolor='orange', linewidth=0.5)
            ax.add_patch(rect)
        
        ax.text(x, y_edited, base, fontsize=9, ha='center', va='center',
               color=color, fontweight='bold', fontfamily='monospace')
    
    # Plot match indicators
    for i in range(min(len(ctrl_segment), len(edit_segment))):
        x = i * base_width
        if ctrl_segment[i].upper() == edit_segment[i].upper():
            ax.text(x, y_match, '|', fontsize=7, ha='center', va='center', color='#AAAAAA')
        else:
            ax.text(x, y_match, '×', fontsize=9, ha='center', va='center', color='red', fontweight='bold')
    
    # Mark cut site
    cut_x = (cut_site - start) * base_width
    if 0 <= cut_x <= len(ctrl_segment) * base_width:
        ax.axvline(x=cut_x, color='red', linewidth=2, linestyle='--', alpha=0.8)
        ax.text(cut_x, 3.5, '▼ Cut Site', fontsize=8, color='red', ha='center', fontweight='bold')
    
    # Labels
    ax.text(-2, y_control, 'Control:', fontsize=9, ha='right', va='center', fontweight='bold', color='#333333')
    ax.text(-2, y_edited, 'Edited:', fontsize=9, ha='right', va='center', fontweight='bold', color='#333333')
    
    ax.set_xlim(-4, len(ctrl_segment) * base_width + 1)
    ax.set_ylim(-0.2, 4.0)
    ax.axis('off')
    
    # Legend for base colors
    legend_text = 'Base Colors:  A=Green  T=Red  G=Black  C=Blue  |  Yellow highlight = mismatch'
    ax.set_title(f'Sequence Alignment\n{legend_text}', fontsize=10, fontweight='bold', loc='left')


def _plot_chromatogram_matplotlib(ax: plt.Axes, traces: Dict, sequence: str,
                                   cut_site: int, title: str = "Chromatogram",
                                   window: int = 80):
    """
    Plot chromatogram with individual base channels.
    Fixed to handle tuple/list traces properly.
    """
    # Check if traces exist and have data
    if not traces:
        ax.text(0.5, 0.5, 'No trace data available', ha='center', va='center', 
                transform=ax.transAxes, fontsize=11, color='gray')
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_facecolor('#F5F5F5')
        return
    
    # Check for valid trace data
    trace_a = traces.get('A')
    if trace_a is None:
        ax.text(0.5, 0.5, 'No trace data available', ha='center', va='center', 
                transform=ax.transAxes, fontsize=11, color='gray')
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_facecolor('#F5F5F5')
        return
    
    # Convert traces to numpy arrays safely
    try:
        trace_arrays = {}
        for base in ['A', 'T', 'G', 'C']:
            if base in traces and traces[base] is not None:
                trace_data = traces[base]
                # Handle tuple, list, or array
                if isinstance(trace_data, (tuple, list)):
                    trace_arrays[base] = np.array(trace_data, dtype=float)
                elif isinstance(trace_data, np.ndarray):
                    trace_arrays[base] = trace_data.astype(float)
                else:
                    trace_arrays[base] = np.array([trace_data], dtype=float)
        
        if not trace_arrays:
            ax.text(0.5, 0.5, 'Invalid trace data format', ha='center', va='center', 
                    transform=ax.transAxes, fontsize=11, color='gray')
            ax.set_title(title, fontsize=10, fontweight='bold')
            return
            
    except Exception as e:
        ax.text(0.5, 0.5, f'Error processing traces: {str(e)[:30]}', ha='center', va='center', 
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_title(title, fontsize=10, fontweight='bold')
        return
    
    # Get trace length
    trace_len = max(len(arr) for arr in trace_arrays.values())
    if trace_len < 100:
        ax.text(0.5, 0.5, 'Insufficient trace data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=11, color='gray')
        ax.set_title(title, fontsize=10, fontweight='bold')
        return
    
    # Calculate window
    trace_factor = 10
    seq_len = len(sequence) if sequence else trace_len // trace_factor
    
    start_bp = max(0, cut_site - window // 2)
    end_bp = min(seq_len, cut_site + window // 2)
    
    start_trace = start_bp * trace_factor
    end_trace = min(end_bp * trace_factor, trace_len)
    
    if end_trace <= start_trace or end_trace - start_trace < 50:
        # Use full range if window calculation fails
        start_trace = 0
        end_trace = min(1000, trace_len)
        start_bp = 0
    
    x_vals = np.arange(end_trace - start_trace)
    x_bp = start_bp + (x_vals / trace_factor)
    
    # Normalize and plot each channel
    max_signal = 1
    for base, arr in trace_arrays.items():
        if end_trace <= len(arr):
            max_signal = max(max_signal, np.max(arr[start_trace:end_trace]))
    
    # Plot traces with TIDE-style colors
    color_map = {'A': '#00AA00', 'T': '#FF0000', 'G': '#333333', 'C': '#0000FF'}
    
    for base in ['A', 'T', 'G', 'C']:
        if base in trace_arrays:
            arr = trace_arrays[base]
            if end_trace <= len(arr):
                signal = arr[start_trace:end_trace]
                signal = signal / max_signal if max_signal > 0 else signal
                ax.plot(x_bp, signal, color=color_map[base], alpha=0.8, linewidth=0.8, label=base)
    
    # Mark cut site
    if start_bp <= cut_site <= end_bp:
        ax.axvline(x=cut_site, color='red', linewidth=2, linestyle='--', alpha=0.8)
        ax.text(cut_site, 1.05, '▼', fontsize=10, color='red', ha='center')
    
    ax.set_xlabel('Position (bp)', fontsize=9)
    ax.set_ylabel('Signal Intensity', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8, ncol=4, framealpha=0.9)
    ax.set_ylim(-0.05, 1.15)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_facecolor('#FAFAFA')


def _plot_indel_spectrum_matplotlib(ax: plt.Axes, analysis_results: Dict):
    """
    Plot indel spectrum in TIDE style with proper color coding.
    Deletions: Red, Wild-type: Green, Insertions: Blue
    """
    indel_spectrum = analysis_results.get('indel_spectrum', {})
    efficiency = analysis_results.get('editing_efficiency', 0)
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # Prepare data
    all_indels = list(range(-10, 11))
    frequencies = []
    colors = []
    
    for indel in all_indels:
        if indel == 0:
            frequencies.append(wt_fraction)
            colors.append(TIDE_COLORS['wildtype'])  # Green for wild-type
        else:
            freq = indel_spectrum.get(indel, indel_spectrum.get(str(indel), 0))
            frequencies.append(freq)
            if indel < 0:
                colors.append(TIDE_COLORS['deletion'])  # Red for deletions
            else:
                colors.append(TIDE_COLORS['insertion'])  # Blue for insertions
    
    # Create bar plot
    bars = ax.bar(all_indels, frequencies, color=colors, edgecolor='black', linewidth=0.5, alpha=0.85)
    
    # Add value labels for significant bars
    for i, (indel, freq) in enumerate(zip(all_indels, frequencies)):
        if freq > 3:
            ax.text(indel, freq + 1.5, f'{freq:.0f}%', ha='center', va='bottom', fontsize=7, fontweight='bold')
    
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.axvline(x=0, color='gray', linewidth=1, linestyle='-', alpha=0.3)
    
    ax.set_xlabel('Indel Size (bp)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Frequency (%)', fontsize=10, fontweight='bold')
    ax.set_title('Indel Spectrum', fontsize=11, fontweight='bold')
    
    # Legend with TIDE colors
    legend_elements = [
        mpatches.Patch(color=TIDE_COLORS['deletion'], label='Deletions', alpha=0.85),
        mpatches.Patch(color=TIDE_COLORS['wildtype'], label='Wild-type (0)', alpha=0.85),
        mpatches.Patch(color=TIDE_COLORS['insertion'], label='Insertions', alpha=0.85)
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, framealpha=0.9)
    
    ax.set_xlim(-11, 11)
    ax.set_xticks(range(-10, 11, 2))
    ax.grid(True, axis='y', alpha=0.3, linestyle='--')
    ax.set_facecolor('#FAFAFA')
    
    # Set y-axis limit dynamically
    max_freq = max(frequencies) if frequencies else 100
    ax.set_ylim(0, min(100, max_freq * 1.2 + 5))


def _plot_efficiency_pie_matplotlib(ax: plt.Axes, analysis_results: Dict):
    """
    Plot editing efficiency pie chart with TIDE-style colors.
    Edited: Red, Wild-type: Green (matching TIDE)
    """
    efficiency = analysis_results.get('editing_efficiency', 0)
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # Ensure values are valid
    if efficiency < 0:
        efficiency = 0
    if wt_fraction < 0:
        wt_fraction = 100 - efficiency
    
    # Ensure they sum to 100
    total = efficiency + wt_fraction
    if total > 0:
        efficiency = (efficiency / total) * 100
        wt_fraction = (wt_fraction / total) * 100
    
    values = [efficiency, wt_fraction]
    # TIDE uses Red for edited, Green for wild-type
    colors = [TIDE_COLORS['edited'], TIDE_COLORS['unedited']]
    labels = ['Edited', 'Wild-type']
    
    # Create pie chart
    wedges, texts, autotexts = ax.pie(values, labels=labels, colors=colors,
                                      autopct='%1.1f%%', startangle=90,
                                      explode=(0.02, 0), shadow=False,
                                      wedgeprops={'edgecolor': 'white', 'linewidth': 1})
    
    for text in texts:
        text.set_fontsize(10)
        text.set_fontweight('bold')
    for autotext in autotexts:
        autotext.set_fontsize(9)
        autotext.set_fontweight('bold')
        autotext.set_color('white')
    
    # Title with metrics
    r2 = analysis_results.get('quality_score', analysis_results.get('r_squared', 0))
    confidence = analysis_results.get('confidence', 'Unknown')
    
    title = f"Editing Efficiency: {efficiency:.1f}%"
    dominant = analysis_results.get('dominant_indel_size')
    if dominant and dominant not in ['N/A', 'Unknown', None, 0]:
        title += f"\nDominant Indel: {dominant} bp"
    
    ax.set_title(title, fontsize=11, fontweight='bold', pad=15)


def _plot_summary_matplotlib(ax: plt.Axes, analysis_results: Dict, grna_sequence: str, sample_name: str):
    """
    Plot summary statistics panel.
    """
    ax.set_facecolor('#F0F0F0')
    ax.axis('off')
    
    efficiency = analysis_results.get('editing_efficiency', 0)
    quality = analysis_results.get('quality_score', analysis_results.get('r_squared', 0))
    confidence = analysis_results.get('confidence', 'Unknown')
    method = analysis_results.get('method', 'Trace Decomposition')
    wt = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # Determine confidence color
    if 'HIGH' in str(confidence).upper():
        conf_color = '#228B22'  # Green
    elif 'MEDIUM' in str(confidence).upper():
        conf_color = '#FF8C00'  # Orange
    else:
        conf_color = '#DC143C'  # Red
    
    # Build summary text
    summary_lines = [
        f"Sample: {sample_name}",
        f"gRNA: {grna_sequence[:25]}{'...' if len(grna_sequence) > 25 else ''}",
        f"Analysis Method: {method}",
        "",
        f"Editing Efficiency: {efficiency:.1f}%",
        f"Wild-type Fraction: {wt:.1f}%", 
        f"R² Quality Score: {quality:.1f}%",
        f"Confidence: {confidence}"
    ]
    
    summary_text = '\n'.join(summary_lines)
    
    ax.text(0.5, 0.5, summary_text, transform=ax.transAxes, fontsize=10,
           ha='center', va='center', fontfamily='monospace',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='#CCCCCC', alpha=0.95))


def figure_to_base64(fig: plt.Figure) -> str:
    """
    Convert matplotlib figure to base64 string.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode()
    buf.close()
    return img_str


def display_results_streamlit(st, analysis_results: Dict, sample_name: str):
    """
    Display analysis results using native Streamlit components.
    Call this function with Streamlit's st module.
    
    Args:
        st: Streamlit module
        analysis_results: Analysis results dictionary
        sample_name: Name of sample
    """
    efficiency = analysis_results.get('editing_efficiency', 0)
    quality = analysis_results.get('quality_score', analysis_results.get('r_squared', 0))
    confidence = analysis_results.get('confidence', 'Unknown')
    wt = analysis_results.get('wt_fraction', 100 - efficiency)
    dominant = analysis_results.get('dominant_indel_size', 'N/A')
    
    # Display metrics using native Streamlit
    st.subheader(f"📊 {sample_name} - Analysis Results")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Editing Efficiency", f"{efficiency:.1f}%")
    with col2:
        st.metric("Wild-type", f"{wt:.1f}%")
    with col3:
        st.metric("R² Quality", f"{quality:.1f}%")
    with col4:
        if dominant not in ['N/A', 'Unknown', None, 0]:
            st.metric("Dominant Indel", f"{dominant} bp")
        else:
            st.metric("Dominant Indel", "N/A")
    
    # Confidence badge
    conf_str = str(confidence).upper()
    if 'HIGH' in conf_str:
        st.success(f"✅ Confidence: {confidence}")
    elif 'MEDIUM' in conf_str:
        st.warning(f"⚠️ Confidence: {confidence}")
    else:
        st.info(f"ℹ️ Confidence: {confidence}")


def display_indel_spectrum_streamlit(st, analysis_results: Dict):
    """
    Display indel spectrum using matplotlib with TIDE-style colors.
    This creates a properly color-coded bar chart.
    
    Args:
        st: Streamlit module
        analysis_results: Analysis results dictionary
    """
    import matplotlib.pyplot as plt
    
    indel_spectrum = analysis_results.get('indel_spectrum', {})
    efficiency = analysis_results.get('editing_efficiency', 0)
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Prepare data
    all_indels = list(range(-10, 11))
    frequencies = []
    colors = []
    
    for indel in all_indels:
        if indel == 0:
            frequencies.append(wt_fraction)
            colors.append(TIDE_COLORS['wildtype'])  # Green
        else:
            freq = indel_spectrum.get(indel, indel_spectrum.get(str(indel), 0))
            frequencies.append(freq)
            if indel < 0:
                colors.append(TIDE_COLORS['deletion'])  # Red
            else:
                colors.append(TIDE_COLORS['insertion'])  # Blue
    
    # Create bar plot
    bars = ax.bar(all_indels, frequencies, color=colors, edgecolor='black', linewidth=0.5, alpha=0.85)
    
    # Add value labels
    for i, (indel, freq) in enumerate(zip(all_indels, frequencies)):
        if freq > 3:
            ax.text(indel, freq + 1, f'{freq:.0f}%', ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_xlabel('Indel Size (bp)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Frequency (%)', fontsize=10, fontweight='bold')
    ax.set_title('Indel Spectrum', fontsize=11, fontweight='bold')
    
    # Legend
    legend_elements = [
        mpatches.Patch(color=TIDE_COLORS['deletion'], label='Deletions', alpha=0.85),
        mpatches.Patch(color=TIDE_COLORS['wildtype'], label='Wild-type', alpha=0.85),
        mpatches.Patch(color=TIDE_COLORS['insertion'], label='Insertions', alpha=0.85)
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    ax.set_xlim(-11, 11)
    ax.set_xticks(range(-10, 11, 2))
    ax.grid(True, axis='y', alpha=0.3)
    
    max_freq = max(frequencies) if frequencies else 100
    ax.set_ylim(0, min(100, max_freq * 1.2 + 5))
    
    plt.tight_layout()
    
    # Display in Streamlit
    st.subheader("📈 Indel Spectrum")
    st.pyplot(fig)
    plt.close(fig)
    
    # Show data table in expander (only if not already in an expander)
    import pandas as pd
    try:
        with st.expander("📋 View Indel Data"):
            data = []
            for indel, freq, color in zip(all_indels, frequencies, colors):
                if freq > 0.1:
                    if indel == 0:
                        indel_type = 'Wild-type'
                    elif indel < 0:
                        indel_type = 'Deletion'
                    else:
                        indel_type = 'Insertion'
                    data.append({'Indel Size (bp)': indel, 'Frequency (%)': round(freq, 1), 'Type': indel_type})
            
            if data:
                df = pd.DataFrame(data).sort_values('Frequency (%)', ascending=False)
                st.dataframe(df, hide_index=True, use_container_width=True)
    except Exception:
        # If expander fails (e.g., nested), just skip the data table
        pass


def display_indel_spectrum_simple(st, analysis_results: Dict):
    """
    Display indel spectrum chart only (no expanders) - for use in batch results.
    
    Args:
        st: Streamlit module
        analysis_results: Analysis results dictionary
    """
    import matplotlib.pyplot as plt

    indel_spectrum = analysis_results.get('indel_spectrum', {})
    efficiency = analysis_results.get('editing_efficiency', 0)
    wt_fraction = analysis_results.get('wt_fraction', 100 - efficiency)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 3))

    # Prepare data
    all_indels = list(range(-10, 11))
    frequencies = []
    colors = []

    for indel in all_indels:
        if indel == 0:
            frequencies.append(wt_fraction)
            colors.append(TIDE_COLORS['wildtype'])  # Green
        else:
            freq = indel_spectrum.get(indel, indel_spectrum.get(str(indel), 0))
            frequencies.append(freq)
            if indel < 0:
                colors.append(TIDE_COLORS['deletion'])  # Red
            else:
                colors.append(TIDE_COLORS['insertion'])  # Blue

    # Create bar plot
    bars = ax.bar(all_indels, frequencies, color=colors, edgecolor='black', linewidth=0.5, alpha=0.85)

    # Add value labels for significant bars
    for i, (indel, freq) in enumerate(zip(all_indels, frequencies)):
        if freq > 5:
            ax.text(indel, freq + 1, f'{freq:.0f}%', ha='center', va='bottom', fontsize=7, fontweight='bold')

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_xlabel('Indel Size (bp)', fontsize=9)
    ax.set_ylabel('Frequency (%)', fontsize=9)
    ax.set_title('Indel Spectrum (Red=Del, Green=WT, Blue=Ins)', fontsize=9, fontweight='bold')

    ax.set_xlim(-11, 11)
    ax.set_xticks(range(-10, 11, 2))
    ax.grid(True, axis='y', alpha=0.3)

    max_freq = max(frequencies) if frequencies else 100
    ax.set_ylim(0, min(100, max_freq * 1.2 + 5))

    plt.tight_layout()
    st.pyplot(fig)
    plt.close(fig)


# Convenience functions for backward compatibility
def create_streamlit_display_html(analysis_results: Dict, sample_name: str) -> str:
    """Kept for compatibility - returns empty string, use display_results_streamlit instead."""
    return ""


def create_indel_spectrum_html(analysis_results: Dict) -> str:
    """Kept for compatibility - returns empty string, use display_indel_spectrum_streamlit instead."""
    return ""


# Export functions
__all__ = [
    'create_tide_style_figure',
    'figure_to_base64',
    'display_results_streamlit',
    'display_indel_spectrum_streamlit',
    'create_streamlit_display_html',
    'create_indel_spectrum_html',
    'TIDE_COLORS'
]
