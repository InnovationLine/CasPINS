"""
TIDE Algorithm Module
Implements TIDE decomposition and editing efficiency calculations
"""

import numpy as np
from scipy import signal, optimize
from scipy.optimize import nnls


def decompose_traces_tide(control_traces, edited_traces, cut_site, window_size=50):
    """
    Implement proper TIDE decomposition algorithm.
    Decomposes edited trace into sum of deletion/insertion traces.
    
    Args:
        control_traces: Dictionary of control chromatogram traces
        edited_traces: Dictionary of edited chromatogram traces
        cut_site: Expected cut site position
        window_size: Size of window around cut site to analyze
        
    Returns:
        dict: TIDE analysis results
    """
    # Get trace region around cut site - focus on downstream region where indels appear
    trace_factor = 10  # Typical scaling for AB1 files
    
    # Handle case where cut_site might be None or beyond sequence
    if not cut_site or cut_site < 0:
        # Use middle of sequence as fallback
        cut_site = len(control_traces['A']) // (2 * trace_factor)
    
    # Decomposition window: from cut_site + 5bp to cut_site + window_size + 5bp
    start_pos = max(0, (cut_site + 5) * trace_factor)
    end_pos = min(len(control_traces['A']), (cut_site + window_size + 5) * trace_factor)
    
    # Combine all four channels for analysis
    control_signal = np.zeros(end_pos - start_pos)
    edited_signal = np.zeros(end_pos - start_pos)
    
    for base in ['A', 'C', 'G', 'T']:
        if len(control_traces[base]) > end_pos and len(edited_traces[base]) > end_pos:
            control_signal += control_traces[base][start_pos:end_pos]
            edited_signal += edited_traces[base][start_pos:end_pos]
    
    # Skip if no signal
    if len(control_signal) == 0 or np.max(control_signal) == 0:
        return {
            'editing_efficiency': 0.0,
            'dominant_indel_size': 0,
            'dominant_indel_percent': 0.0,
            'indel_spectrum': {},
            'quality_score': 0.0,
            'confidence': 'ERROR - No signal in decomposition window'
        }
    
    # Normalize signals
    control_signal = control_signal / np.max(control_signal)
    edited_signal = edited_signal / np.max(edited_signal) if np.max(edited_signal) > 0 else edited_signal
    
    # Build matrix of all possible indel traces
    # Each column is a shifted version of the control signal
    indel_range = range(-10, 11)  # -10 to +10 bp indels
    num_indels = len(indel_range)
    
    # Create matrix where each column is a shifted control signal
    A = np.zeros((len(control_signal), num_indels))
    
    for idx, indel_size in enumerate(indel_range):
        if indel_size == 0:  # Wild-type
            A[:, idx] = control_signal
        elif indel_size > 0:  # Deletion (shift left)
            shift_amount = indel_size * trace_factor // 10  # Approximate shift
            if shift_amount < len(control_signal):
                A[:-shift_amount, idx] = control_signal[shift_amount:]
        else:  # Insertion (shift right)
            shift_amount = -indel_size * trace_factor // 10
            if shift_amount < len(control_signal):
                A[shift_amount:, idx] = control_signal[:-shift_amount]
    
    # Use Non-Negative Least Squares (NNLS) to find the best fit
    # This ensures all coefficients are >= 0
    coefficients, residual = nnls(A, edited_signal, maxiter=1000)
    
    # Normalize coefficients to sum to 1
    total_signal = np.sum(coefficients)
    if total_signal > 0:
        coefficients = coefficients / total_signal
    
    # Extract indel spectrum
    indel_spectrum = {}
    wt_fraction = 0
    
    for idx, indel_size in enumerate(indel_range):
        fraction = coefficients[idx]
        if fraction > 0.01:  # Only include significant contributions (>1%)
            if indel_size == 0:
                wt_fraction = fraction
            else:
                indel_spectrum[indel_size] = round(fraction * 100, 1)
    
    # Calculate editing efficiency (100% - WT%)
    editing_efficiency = round((1 - wt_fraction) * 100, 1)
    
    # Find dominant indel
    if indel_spectrum:
        dominant_indel = max(indel_spectrum.items(), key=lambda x: x[1])
        dominant_indel_size = dominant_indel[0]
        dominant_indel_percent = dominant_indel[1]
    else:
        dominant_indel_size = 0
        dominant_indel_percent = 0.0
    
    # Calculate quality score based on how well the model fits
    if len(edited_signal) > 0:
        reconstructed = np.dot(A, coefficients)
        mse = np.mean((edited_signal - reconstructed) ** 2)
        # R-squared calculation
        ss_tot = np.sum((edited_signal - np.mean(edited_signal)) ** 2)
        ss_res = np.sum((edited_signal - reconstructed) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        quality_score = round(max(0, r_squared * 100), 1)
    else:
        quality_score = 0.0
    
    # Determine confidence based on quality and efficiency
    if quality_score > 80 and editing_efficiency > 10:
        confidence = 'HIGH'
    elif quality_score > 60:
        confidence = 'MEDIUM'
    else:
        confidence = 'LOW'
    
    return {
        'editing_efficiency': editing_efficiency,
        'dominant_indel_size': dominant_indel_size,
        'dominant_indel_percent': dominant_indel_percent,
        'indel_spectrum': indel_spectrum,
        'quality_score': quality_score,
        'confidence': confidence,
        'wt_fraction': round(wt_fraction * 100, 1)
    }


def calculate_editing_efficiency_fallback(control_seq, control_traces, edited_seq, edited_traces, cut_position):
    """
    Fallback method for calculating editing efficiency using signal decay.
    Used when TIDE decomposition fails.
    
    Args:
        control_seq: Control sequence
        control_traces: Control chromatogram traces
        edited_seq: Edited sequence
        edited_traces: Edited chromatogram traces
        cut_position: Expected cut position
        
    Returns:
        dict: Efficiency calculation results
    """
    # For clonal cell lines, we expect to see clear differences after the cut site
    # Look at signal quality degradation as a proxy for editing
    
    if not cut_position or cut_position < 50:
        cut_position = len(control_seq) // 2
    
    # Define regions for analysis
    upstream_start = max(0, cut_position - 50)
    upstream_end = cut_position - 5
    downstream_start = cut_position + 5
    downstream_end = min(len(control_seq), cut_position + 100)
    
    # Calculate average peak heights in regions
    trace_factor = 10
    
    def calculate_signal_quality(traces, start_pos, end_pos):
        """Calculate signal quality metric."""
        start_trace = start_pos * trace_factor
        end_trace = min(end_pos * trace_factor, len(traces['A']))
        
        if end_trace <= start_trace:
            return 0
        
        # Sum signal across all channels
        total_signal = 0
        for base in ['A', 'C', 'G', 'T']:
            if end_trace <= len(traces[base]):
                segment = traces[base][start_trace:end_trace]
                total_signal += np.sum(segment)
        
        return total_signal / (end_trace - start_trace) if end_trace > start_trace else 0
    
    # Calculate qualities
    control_upstream = calculate_signal_quality(control_traces, upstream_start, upstream_end)
    control_downstream = calculate_signal_quality(control_traces, downstream_start, downstream_end)
    edited_upstream = calculate_signal_quality(edited_traces, upstream_start, upstream_end)
    edited_downstream = calculate_signal_quality(edited_traces, downstream_start, downstream_end)
    
    # Calculate ratios
    control_ratio = control_downstream / control_upstream if control_upstream > 0 else 0
    edited_ratio = edited_downstream / edited_upstream if edited_upstream > 0 else 0
    
    # For clonal lines, if the downstream signal is significantly degraded, it indicates editing
    signal_degradation = max(0, control_ratio - edited_ratio) / control_ratio if control_ratio > 0 else 0
    
    # Estimate efficiency (clonal lines typically show 50-100% efficiency per allele)
    # If signal is degraded by >30%, likely heterozygous edit (~50% efficiency)
    # If signal is degraded by >60%, likely homozygous edit (~100% efficiency)
    if signal_degradation > 0.6:
        editing_efficiency = 95.0  # Likely homozygous
    elif signal_degradation > 0.3:
        editing_efficiency = 50.0  # Likely heterozygous
    elif signal_degradation > 0.1:
        editing_efficiency = 25.0  # Partial editing
    else:
        editing_efficiency = 0.0  # No significant editing
    
    return {
        'editing_efficiency': editing_efficiency,
        'signal_decay_ratio': round(signal_degradation, 3),
        'quality_score': 50.0,  # Fallback method has lower confidence
        'sequence_similarity': 0.0,
        'dominant_indel_size': 'Unknown',
        'dominant_indel_percent': 'Unknown',
        'confidence': 'LOW - Fallback method used',
        'wt_fraction': round(100 - editing_efficiency, 1)
    } 