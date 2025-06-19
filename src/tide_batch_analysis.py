import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from datetime import datetime
import primer3
from scipy import signal, optimize
from scipy.fft import fft, ifft
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

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

def align_sequences(seq1, seq2):
    """Perform global alignment of two sequences to find exact indel locations."""
    # Use Needleman-Wunsch global alignment
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
    
    if alignments:
        best_alignment = alignments[0]
        aligned_seq1, aligned_seq2, score, start, end = best_alignment
        
        # Find first significant difference
        for i, (base1, base2) in enumerate(zip(aligned_seq1, aligned_seq2)):
            if base1 != base2:
                # Count position in original sequence (excluding gaps)
                original_pos = len(aligned_seq1[:i].replace('-', ''))
                return original_pos, aligned_seq1, aligned_seq2
    
    return len(seq1) // 2, seq1, seq2

def decompose_traces_tide(control_traces, edited_traces, cut_site, window_size=50):
    """
    Implement proper TIDE decomposition algorithm.
    Decomposes edited trace into sum of deletion/insertion traces.
    """
    # Get trace region around cut site
    trace_factor = 10  # More accurate scaling based on typical AB1 files
    
    start_pos = max(0, (cut_site - window_size) * trace_factor)
    end_pos = min(len(control_traces['A']), (cut_site + window_size * 2) * trace_factor)
    
    # Combine all four channels for analysis
    control_signal = np.zeros(end_pos - start_pos)
    edited_signal = np.zeros(end_pos - start_pos)
    
    for base in ['A', 'C', 'G', 'T']:
        if len(control_traces[base]) > end_pos and len(edited_traces[base]) > end_pos:
            control_signal += control_traces[base][start_pos:end_pos]
            edited_signal += edited_traces[base][start_pos:end_pos]
    
    # Normalize signals
    control_signal = control_signal / np.max(control_signal) if np.max(control_signal) > 0 else control_signal
    edited_signal = edited_signal / np.max(edited_signal) if np.max(edited_signal) > 0 else edited_signal
    
    # Find optimal shift using cross-correlation
    correlation = signal.correlate(edited_signal, control_signal, mode='same')
    shift = np.argmax(correlation) - len(edited_signal) // 2
    
    # Calculate indel spectrum (-10 to +10 bp)
    indel_spectrum = {}
    best_fit_score = float('inf')
    best_efficiency = 0
    
    for indel_size in range(-10, 11):
        if indel_size == 0:
            continue
            
        # Create shifted version of control
        if indel_size > 0:  # Deletion
            shifted = np.concatenate([control_signal[indel_size:], np.zeros(indel_size)])
        else:  # Insertion
            shifted = np.concatenate([np.zeros(-indel_size), control_signal[:indel_size]])
        
        # Fit linear combination of WT and indel
        def objective(params):
            wt_fraction, indel_fraction = params
            if wt_fraction < 0 or indel_fraction < 0 or wt_fraction + indel_fraction > 1:
                return float('inf')
            reconstructed = wt_fraction * control_signal + indel_fraction * shifted
            return np.sum((edited_signal - reconstructed) ** 2)
        
        result = optimize.minimize(objective, [0.5, 0.5], bounds=[(0, 1), (0, 1)])
        
        if result.fun < best_fit_score:
            best_fit_score = result.fun
            wt_frac, indel_frac = result.x
            indel_spectrum[indel_size] = indel_frac * 100
            best_efficiency = (1 - wt_frac) * 100
    
    # Find dominant indel
    dominant_indel = max(indel_spectrum.items(), key=lambda x: x[1]) if indel_spectrum else (0, 0)
    
    return {
        'editing_efficiency': round(best_efficiency, 1),
        'dominant_indel_size': dominant_indel[0],
        'dominant_indel_percent': round(dominant_indel[1], 1),
        'indel_spectrum': indel_spectrum,
        'quality_score': round((1 - best_fit_score) * 100, 1)
    }

def design_primers_with_primer3(mrna_seq, target_regions, gene_name):
    """
    Design primers using primer3-py with rigorous parameters.
    
    Args:
        mrna_seq: mRNA sequence
        target_regions: List of (start, end) tuples for regions to amplify
        gene_name: Gene name for primer naming
    
    Returns:
        List of primer pair dictionaries
    """
    primer_results = []
    
    for i, (region_start, region_end) in enumerate(target_regions):
        # Define target region with flanking sequences
        target_start = max(0, region_start - 200)
        target_end = min(len(mrna_seq), region_end + 200)
        target_len = region_end - region_start
        
        # Primer3 parameters
        seq_args = {
            'SEQUENCE_ID': f'{gene_name}_region{i+1}',
            'SEQUENCE_TEMPLATE': mrna_seq,
            'SEQUENCE_TARGET': [region_start, target_len],
            'SEQUENCE_INCLUDED_REGION': [target_start, target_end - target_start]
        }
        
        global_args = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 4,
            'PRIMER_MAX_SELF_END': 2,
            'PRIMER_PAIR_MAX_COMPL_ANY': 4,
            'PRIMER_PAIR_MAX_COMPL_END': 2,
            'PRIMER_PRODUCT_SIZE_RANGE': [[300, 800]],
            'PRIMER_NUM_RETURN': 5
        }
        
        try:
            # Design primers
            primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
            
            # Extract primer pairs
            num_primers = primer3_result.get('PRIMER_PAIR_NUM_RETURNED', 0)
            
            for j in range(num_primers):
                left_seq = primer3_result.get(f'PRIMER_LEFT_{j}_SEQUENCE', '')
                right_seq = primer3_result.get(f'PRIMER_RIGHT_{j}_SEQUENCE', '')
                left_tm = primer3_result.get(f'PRIMER_LEFT_{j}_TM', 0)
                right_tm = primer3_result.get(f'PRIMER_RIGHT_{j}_TM', 0)
                left_gc = primer3_result.get(f'PRIMER_LEFT_{j}_GC_PERCENT', 0)
                right_gc = primer3_result.get(f'PRIMER_RIGHT_{j}_GC_PERCENT', 0)
                left_pos = primer3_result.get(f'PRIMER_LEFT_{j}', [0, 0])[0]
                right_pos = primer3_result.get(f'PRIMER_RIGHT_{j}', [0, 0])[0]
                product_size = primer3_result.get(f'PRIMER_PAIR_{j}_PRODUCT_SIZE', 0)
                
                primer_pair = {
                    'rank': j + 1,
                    'forward_seq': left_seq,
                    'reverse_seq': right_seq,
                    'forward_tm': round(left_tm, 1),
                    'reverse_tm': round(right_tm, 1),
                    'forward_gc': round(left_gc, 1),
                    'reverse_gc': round(right_gc, 1),
                    'forward_pos': left_pos,
                    'reverse_pos': right_pos,
                    'product_size': product_size,
                    'target_region': f'{region_start}-{region_end}'
                }
                
                primer_results.append(primer_pair)
                
        except Exception as e:
            print(f"  WARNING: Primer3 failed for region {i+1}: {str(e)}")
            # Fall back to simple primer design
            fallback_primers = design_simple_primers(mrna_seq, region_start, region_end)
            primer_results.extend(fallback_primers)
    
    return primer_results

def design_simple_primers(seq, target_start, target_end):
    """Fallback simple primer design if primer3 fails."""
    primers = []
    
    # Forward primer 150-200bp upstream
    for offset in [150, 175, 200]:
        f_start = max(0, target_start - offset)
        f_seq = seq[f_start:f_start+20]
        
        # Reverse primer 150-200bp downstream
        r_start = min(len(seq)-20, target_end + offset - 20)
        r_seq = str(Seq(seq[r_start:r_start+20]).reverse_complement())
        
        if len(f_seq) == 20 and len(r_seq) == 20:
            primers.append({
                'rank': len(primers) + 1,
                'forward_seq': f_seq,
                'reverse_seq': r_seq,
                'forward_tm': calculate_tm(f_seq),
                'reverse_tm': calculate_tm(r_seq),
                'forward_gc': (f_seq.count('G') + f_seq.count('C')) * 5,
                'reverse_gc': (r_seq.count('G') + r_seq.count('C')) * 5,
                'forward_pos': f_start,
                'reverse_pos': r_start,
                'product_size': r_start + 20 - f_start,
                'target_region': f'{target_start}-{target_end}',
                'note': 'Simple design (Primer3 unavailable)'
            })
    
    return primers

def calculate_tm(sequence):
    """Calculate melting temperature using nearest-neighbor method."""
    # Simplified calculation - primer3 does this better
    gc_count = sequence.count('G') + sequence.count('C')
    at_count = sequence.count('A') + sequence.count('T')
    
    if len(sequence) < 14:
        return (gc_count * 4) + (at_count * 2)
    else:
        return 64.9 + 41 * (gc_count - 16.4) / len(sequence)

def find_divergence_point(seq1, seq2):
    """Find the point where two sequences start to diverge significantly."""
    min_len = min(len(seq1), len(seq2))
    window_size = 10
    
    # First check if sequences are identical
    if seq1[:min_len] == seq2[:min_len]:
        print("    Sequences are identical - no divergence found")
        return -1  # Return -1 to indicate no divergence
    
    # Check overall similarity
    total_matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    similarity = total_matches / min_len * 100
    print(f"    Overall sequence similarity: {similarity:.1f}%")
    
    # If sequences are completely different, return -1
    if similarity < 50:
        print("    WARNING: Sequences have very low similarity - might be from different regions")
        return -1
    
    # Look for divergence point
    for i in range(min_len - window_size):
        window1 = seq1[i:i+window_size]
        window2 = seq2[i:i+window_size]
        mismatches = sum(1 for a, b in zip(window1, window2) if a != b)
        
        if mismatches >= 3:
            # Check if this is a real divergence or just noise
            # Look ahead to see if mismatches continue
            if i + window_size + 10 < min_len:
                next_window1 = seq1[i+5:i+15]
                next_window2 = seq2[i+5:i+15]
                next_mismatches = sum(1 for a, b in zip(next_window1, next_window2) if a != b)
                
                if next_mismatches >= 3:
                    print(f"    Found divergence at position {i} ({mismatches}/10 mismatches)")
                    return i
            else:
                return i
    
    print("    No clear divergence point found")
    return -1

def calculate_editing_efficiency(control_seq, control_traces, edited_seq, edited_traces, cut_position):
    """
    Calculate CRISPR editing efficiency using proper TIDE decomposition.
    
    Returns:
        dict: Contains editing_efficiency (%), dominant indel info, and quality metrics
    """
    # First try proper TIDE decomposition
    try:
        tide_results = decompose_traces_tide(control_traces, edited_traces, cut_position)
        
        # Add sequence-based validation
        alignment_pos, aligned_ctrl, aligned_edit = align_sequences(
            control_seq[max(0, cut_position-50):cut_position+50],
            edited_seq[max(0, cut_position-50):cut_position+50]
        )
        
        # Validate TIDE results with sequence alignment
        seq_similarity = sum(1 for a, b in zip(aligned_ctrl, aligned_edit) if a == b) / len(aligned_ctrl)
        
        tide_results['sequence_similarity'] = round(seq_similarity * 100, 1)
        
        # Adjust confidence based on multiple factors
        if tide_results['quality_score'] < 50:
            tide_results['confidence'] = 'LOW - Poor signal quality'
        elif abs(tide_results['dominant_indel_size']) > 20:
            tide_results['confidence'] = 'MEDIUM - Large indel detected'
        elif tide_results['editing_efficiency'] > 90:
            tide_results['confidence'] = 'HIGH - Clear editing pattern'
        else:
            tide_results['confidence'] = 'MEDIUM'
            
        return tide_results
        
    except Exception as e:
        print(f"  WARNING: TIDE decomposition failed ({str(e)}), using fallback method")
        
        # Fallback to simpler analysis
        # Define regions for analysis
        upstream_start = max(0, cut_position - 50)
        upstream_end = cut_position - 5
        downstream_start = cut_position + 5
        downstream_end = min(len(control_seq), cut_position + 50)
        
        # Calculate average peak heights in upstream and downstream regions
        trace_factor = 12  # Approximate trace to sequence position factor
        
        def calculate_peak_quality(traces, start_pos, end_pos):
            """Calculate average peak height and signal quality in a region."""
            start_trace = start_pos * trace_factor
            end_trace = end_pos * trace_factor
            
            total_signal = 0
            peak_count = 0
            
            for base in ['A', 'C', 'G', 'T']:
                if end_trace <= len(traces[base]):
                    trace_segment = traces[base][start_trace:end_trace]
                    if len(trace_segment) > 0:
                        # Find peaks (local maxima)
                        peaks = []
                        for i in range(1, len(trace_segment)-1):
                            if trace_segment[i] > trace_segment[i-1] and trace_segment[i] > trace_segment[i+1]:
                                peaks.append(trace_segment[i])
                        
                        if peaks:
                            total_signal += sum(peaks)
                            peak_count += len(peaks)
            
            return total_signal / peak_count if peak_count > 0 else 0
        
        # Calculate quality for control sample
        control_upstream_quality = calculate_peak_quality(control_traces, upstream_start, upstream_end)
        control_downstream_quality = calculate_peak_quality(control_traces, downstream_start, downstream_end)
        
        # Calculate quality for edited sample
        edited_upstream_quality = calculate_peak_quality(edited_traces, upstream_start, upstream_end)
        edited_downstream_quality = calculate_peak_quality(edited_traces, downstream_start, downstream_end)
        
        # Calculate signal decay ratio
        control_ratio = control_downstream_quality / control_upstream_quality if control_upstream_quality > 0 else 1
        edited_ratio = edited_downstream_quality / edited_upstream_quality if edited_upstream_quality > 0 else 0
        
        # Estimate editing efficiency based on signal decay
        signal_decay_difference = control_ratio - edited_ratio
        editing_efficiency = min(100, max(0, signal_decay_difference * 100))
        
        # Calculate overall quality score
        quality_score = min(100, (edited_upstream_quality / control_upstream_quality * 100) if control_upstream_quality > 0 else 0)
        
        # Check sequence similarity
        seq_similarity_after_cut = sum(1 for i in range(downstream_start, min(downstream_end, len(edited_seq), len(control_seq))) 
                                      if i < len(control_seq) and i < len(edited_seq) and control_seq[i] == edited_seq[i])
        expected_matches = downstream_end - downstream_start
        similarity_ratio = seq_similarity_after_cut / expected_matches if expected_matches > 0 else 1
        
        # Adjust editing efficiency based on sequence similarity
        if similarity_ratio < 0.7:
            editing_efficiency = max(editing_efficiency, (1 - similarity_ratio) * 100)
        
        return {
            'editing_efficiency': round(editing_efficiency, 1),
            'signal_decay_ratio': round(signal_decay_difference, 3),
            'quality_score': round(quality_score, 1),
            'sequence_similarity': round(similarity_ratio * 100, 1),
            'dominant_indel_size': 'Unknown',
            'dominant_indel_percent': 'Unknown',
            'confidence': 'LOW - Fallback method used'
        }

def plot_tide_analysis_with_validation(control_file, edited_file, output_dir, gene_name, expected_cut_site=None):
    """
    Create a TIDE-style analysis plot with validation against expected cut site.
    
    Args:
        control_file: Path to control AB1 file
        edited_file: Path to edited AB1 file
        output_dir: Output directory
        gene_name: Gene name
        expected_cut_site: Expected cut site position based on gRNA location
    """
    print(f"\nPerforming TIDE Analysis for {gene_name}...")
    
    # Parse both .ab1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    print(f"  Control sequence length: {len(control_seq)} bp")
    print(f"  Edited sequence length: {len(edited_seq)} bp")
    
    # Find where sequences diverge
    divergence_point = find_divergence_point(control_seq, edited_seq)
    
    # Handle cases where no divergence is found
    if divergence_point == -1:
        print("  No divergence detected between control and edited sequences")
        if expected_cut_site:
            print(f"  Using expected cut site from gRNA: position {expected_cut_site}")
            analysis_position = expected_cut_site
            divergence_point = expected_cut_site  # For reporting
        else:
            print("  ERROR: No divergence found and no expected cut site provided!")
            # Use middle of sequence as fallback
            analysis_position = len(control_seq) // 2
            divergence_point = 0  # For reporting
    else:
        print(f"  Detected divergence point: position {divergence_point}")
        
        # Validate against expected cut site
        if expected_cut_site:
            print(f"  Expected cut site from gRNA: position {expected_cut_site}")
            
            # Check if divergence is near expected cut site (within 10bp)
            distance = abs(divergence_point - expected_cut_site)
            if distance > 10:
                print(f"  WARNING: Divergence point is {distance}bp away from expected cut site!")
                print(f"  This may indicate:")
                print(f"    - The wrong gRNA sequence was provided")
                print(f"    - The AB1 sequences don't match the mRNA reference")
                print(f"    - Off-target editing occurred")
                
                # Use expected cut site for analysis if divergence is too far
                if distance > 50:
                    print(f"  Using expected cut site ({expected_cut_site}) for analysis instead of divergence point")
                    analysis_position = expected_cut_site
                else:
                    print(f"  Using detected divergence point for analysis")
                    analysis_position = divergence_point
            else:
                print(f"  ✓ Divergence point matches expected cut site (within {distance}bp)")
                analysis_position = divergence_point
        else:
            # No expected cut site provided, use divergence point
            analysis_position = divergence_point
    
    # Calculate editing efficiency at the validated position
    efficiency_data = calculate_editing_efficiency(
        control_seq, control_traces, 
        edited_seq, edited_traces, 
        analysis_position
    )
    
    # Add validation info to efficiency data
    if expected_cut_site:
        efficiency_data['cut_site_validation'] = {
            'expected': expected_cut_site,
            'detected': divergence_point,
            'distance': abs(divergence_point - expected_cut_site),
            'validated': abs(divergence_point - expected_cut_site) <= 10
        }
    
    print(f"  Editing efficiency: {efficiency_data['editing_efficiency']}%")
    print(f"  Quality score: {efficiency_data['quality_score']}%")
    
    # Check if analysis position is within sequence bounds
    if analysis_position >= len(control_seq):
        print(f"  ERROR: Analysis position ({analysis_position}) is beyond sequence length ({len(control_seq)})")
        print(f"  The AB1 files may be from a different amplicon than expected")
        # Adjust to plot the end of the sequence
        analysis_position = len(control_seq) - 50
        
    # Define viewing window
    window_start = max(0, analysis_position - 50)
    window_end = min(len(control_seq), analysis_position + 100)
    
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
    
    cut_trace_pos = (analysis_position - window_start) * trace_factor
    ax1.axvline(x=cut_trace_pos, color='red', linestyle='--', alpha=0.5, label=f'Cut site (pos {analysis_position})')
    
    # If there's a mismatch, show both positions
    if expected_cut_site and abs(divergence_point - expected_cut_site) > 10:
        div_trace_pos = (divergence_point - window_start) * trace_factor
        ax1.axvline(x=div_trace_pos, color='orange', linestyle=':', alpha=0.5, label=f'Divergence (pos {divergence_point})')
    
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
    
    ax2.axvline(x=cut_trace_pos, color='red', linestyle='--', alpha=0.5, label=f'Cut site (pos {analysis_position})')
    
    if expected_cut_site and abs(divergence_point - expected_cut_site) > 10:
        div_trace_pos = (divergence_point - window_start) * trace_factor
        ax2.axvline(x=div_trace_pos, color='orange', linestyle=':', alpha=0.5, label=f'Divergence (pos {divergence_point})')
    
    ax2.legend(loc='upper right')
    ax2.set_ylim(0, max_signal * 1.1)
    
    # Plot 3: Overlay comparison
    ax3 = plt.subplot(3, 1, 3)
    ax3.set_title(f'{gene_name.upper()} - Overlay: Control (solid) vs Edited (dashed)', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Signal Intensity', fontsize=12)
    ax3.set_xlabel('Trace Position', fontsize=12)
    
    # Focus on region around cut site
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
    
    ax3.axvline(x=min(200, focus_end-focus_start-10), color='red', linestyle='--', alpha=0.5, label='Cut site')
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
    
    return divergence_point, timestamp, efficiency_data


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
    
    # Calculate editing efficiency
    efficiency_data = calculate_editing_efficiency(
        control_seq, control_traces, 
        edited_seq, edited_traces, 
        divergence_point
    )
    
    print(f"  Editing efficiency: {efficiency_data['editing_efficiency']}%")
    print(f"  Quality score: {efficiency_data['quality_score']}%")
    
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
    
    return divergence_point, timestamp, efficiency_data

def analyze_grna_in_mrna(mrna_seq, grna_sequences, gene_name, efficiency_data=None):
    """Analyze gRNA positions in mRNA and recommend primers using primer3."""
    recommendations = []
    recommendations.append(f"PRIMER RECOMMENDATIONS FOR {gene_name.upper()}")
    recommendations.append("=" * 60)
    recommendations.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    recommendations.append(f"mRNA length: {len(mrna_seq)} bp")
    recommendations.append("")
    
    # Add editing efficiency data if available
    if efficiency_data:
        recommendations.append("CRISPR EDITING EFFICIENCY ANALYSIS:")
        recommendations.append("-" * 60)
        
        # Add validation warning if present
        if 'cut_site_validation' in efficiency_data:
            validation = efficiency_data['cut_site_validation']
            if not validation['validated']:
                recommendations.append("⚠️ WARNING: CUT SITE VALIDATION FAILED!")
                recommendations.append(f"  Expected cut site (from gRNA): position {validation['expected']}")
                recommendations.append(f"  Detected divergence point: position {validation['detected']}")
                recommendations.append(f"  Distance: {validation['distance']} bp")
                recommendations.append("")
                recommendations.append("This large discrepancy suggests:")
                recommendations.append("  • The provided gRNA/mRNA sequences may be incorrect")
                recommendations.append("  • The AB1 files may be from a different target region")
                recommendations.append("  • Off-target editing may have occurred")
                recommendations.append("")
        
        recommendations.append(f"Editing Efficiency: {efficiency_data['editing_efficiency']}%")
        
        # Add more detailed TIDE results if available
        if 'dominant_indel_size' in efficiency_data and efficiency_data['dominant_indel_size'] != 'Unknown':
            recommendations.append(f"Dominant Indel: {efficiency_data['dominant_indel_size']} bp ({efficiency_data['dominant_indel_percent']}%)")
        
        recommendations.append(f"Signal Quality Score: {efficiency_data['quality_score']}%")
        recommendations.append(f"Sequence Similarity After Cut: {efficiency_data['sequence_similarity']}%")
        
        if 'confidence' in efficiency_data:
            recommendations.append(f"Analysis Confidence: {efficiency_data['confidence']}")
        
        recommendations.append("")
        
        # Add detailed interpretation and recommendations
        if efficiency_data['editing_efficiency'] >= 70:
            recommendations.append("Interpretation: HIGH editing efficiency detected")
            recommendations.append("Recommendation: This sample shows excellent editing. Proceed with clonal isolation.")
        elif efficiency_data['editing_efficiency'] >= 30:
            recommendations.append("Interpretation: MODERATE editing efficiency detected")
            recommendations.append("Recommendation: Consider enriching edited cells before clonal isolation.")
        else:
            recommendations.append("Interpretation: LOW editing efficiency detected")
            recommendations.append("Recommendation: Consider optimizing transfection conditions or gRNA design.")
            
        # Quality-based recommendations
        if efficiency_data['quality_score'] < 50:
            recommendations.append("\nWARNING: Low signal quality detected!")
            recommendations.append("Recommendation: Re-sequence samples with higher quality DNA or optimize PCR conditions.")
        
        recommendations.append("")

    grna_positions = []
    cut_sites = []
    
    # Find each gRNA in the mRNA
    for i, grna in enumerate(grna_sequences, 1):
        recommendations.append(f"\nGuide RNA {i}: {grna}")
        grna_rc = str(Seq(grna).reverse_complement())
        
        if grna in mrna_seq:
            pos = mrna_seq.find(grna)
            cut_site = pos + 17  # 3bp before PAM
            grna_positions.append((pos, 'forward'))
            cut_sites.append(cut_site)
            recommendations.append(f"  [FOUND] on forward strand at position {pos}")
            recommendations.append(f"  Cut site (3bp before PAM): position {cut_site}")
        elif grna_rc in mrna_seq:
            pos = mrna_seq.find(grna_rc)
            cut_site = pos + 3
            grna_positions.append((pos, 'reverse'))
            cut_sites.append(cut_site)
            recommendations.append(f"  [FOUND] on reverse strand at position {pos}")
            recommendations.append(f"  Cut site: position {cut_site}")
        else:
            recommendations.append(f"  [NOT FOUND] in mRNA sequence!")
            recommendations.append("  WARNING: This gRNA does not match the provided mRNA sequence.")
    
    if grna_positions:
        # Sort positions to find the range
        positions = [pos for pos, strand in grna_positions]
        min_pos = min(positions)
        max_pos = max(positions)
        min_cut = min(cut_sites)
        max_cut = max(cut_sites)
        
        recommendations.append("\n" + "=" * 60)
        recommendations.append("RECOMMENDED SEQUENCING PRIMERS (Designed with Primer3):")
        recommendations.append("=" * 60)
        
        # Define target regions for primer design
        # Target region should span all cut sites with adequate flanking
        target_regions = [(min_cut - 50, max_cut + 50)]
        
        # Use primer3 to design primers
        primer3_results = design_primers_with_primer3(mrna_seq, target_regions, gene_name)
        
        if primer3_results:
            recommendations.append("\nPrimer3-designed primer pairs (ranked by quality):")
            recommendations.append("-" * 60)
            
            for i, primer in enumerate(primer3_results[:3]):  # Show top 3
                recommendations.append(f"\nPrimer Pair {primer['rank']}:")
                recommendations.append(f"  Forward: 5'-{primer['forward_seq']}-3'")
                recommendations.append(f"    Position: {primer['forward_pos']}")
                recommendations.append(f"    Tm: {primer['forward_tm']}°C, GC: {primer['forward_gc']}%")
                recommendations.append(f"  Reverse: 5'-{primer['reverse_seq']}-3'")
                recommendations.append(f"    Position: {primer['reverse_pos']}")
                recommendations.append(f"    Tm: {primer['reverse_tm']}°C, GC: {primer['reverse_gc']}%")
                recommendations.append(f"  Product size: {primer['product_size']} bp")
                
                # Add specific notes for the best primer pair
                if i == 0:
                    recommendations.append(f"  ** RECOMMENDED - Best overall primer pair **")
        else:
            # Fallback to simple primer design if primer3 fails
            recommendations.append("\n[Primer3 unavailable - using simple primer design]")
            
            # Design primers ~200bp away from the gRNA region
            primer_f_start = max(0, min_pos - 200)
            primer_f_seq = mrna_seq[primer_f_start:primer_f_start+20]
            
            recommendations.append(f"\nFallback Primer Set:")
            recommendations.append(f"Forward Primer: 5'-{primer_f_seq}-3'")
            recommendations.append(f"  Position: {primer_f_start}-{primer_f_start+20}")
            recommendations.append(f"  Tm: ~{calculate_tm(primer_f_seq)}°C")
            
            # Reverse primer
            primer_r_start = min(len(mrna_seq)-20, max_pos + 23 + 180)
            primer_r_seq = str(Seq(mrna_seq[primer_r_start:primer_r_start+20]).reverse_complement())
            
            recommendations.append(f"Reverse Primer: 5'-{primer_r_seq}-3'")
            recommendations.append(f"  Position: {primer_r_start}-{primer_r_start+20}")
            recommendations.append(f"  Tm: ~{calculate_tm(primer_r_seq)}°C")
            recommendations.append(f"Expected amplicon size: ~{primer_r_start - primer_f_start + 20} bp")
        
        recommendations.append(f"\ngRNA region spans: {min_pos}-{max_pos + 23} ({max_pos - min_pos + 23} bp)")
        recommendations.append(f"Cut sites span: {min_cut}-{max_cut} ({max_cut - min_cut} bp)")
        
        # Add PCR cycling recommendations
        recommendations.append("\n" + "-" * 40)
        recommendations.append("Recommended PCR conditions:")
        recommendations.append("  - Use high-fidelity polymerase (e.g., Q5, Phusion)")
        recommendations.append("  - Initial denaturation: 98°C for 30s")
        recommendations.append("  - 35 cycles: 98°C 10s, 60°C 20s, 72°C 30s")
        recommendations.append("  - Final extension: 72°C for 2 min")
        
    else:
        recommendations.append("\n*** ERROR: No gRNAs found in the provided mRNA sequence! ***")
        recommendations.append("\nTroubleshooting steps:")
        recommendations.append("  1. Verify the mRNA sequence is correct and complete")
        recommendations.append("  2. Check if gRNAs include PAM sequence (remove if present)")
        recommendations.append("  3. Ensure gRNAs are 20nt long (standard length)")
        recommendations.append("  4. Verify gRNAs target the correct gene")
        recommendations.append("  5. Check if you're using genomic vs cDNA sequence")
        recommendations.append("\nTechnical details:")
        recommendations.append(f"  - Searched for exact matches and reverse complements")
        recommendations.append(f"  - mRNA length: {len(mrna_seq)} bp")
        recommendations.append(f"  - Number of gRNAs: {len(grna_sequences)}")
    
    return "\n".join(recommendations)

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
    
    # Read mRNA sequence FIRST
    with open(mrna_file, 'r') as f:
        mrna_seq = ''.join(line.strip() for line in f).upper().replace(' ', '')
    
    # Verify gRNAs are present in the mRNA
    grna_positions = []
    expected_cut_sites = []
    
    for i, grna in enumerate(grna_sequences, 1):
        original_grna = grna
        grna_length = len(grna)
        
        # Handle non-standard gRNA lengths
        if grna_length > 20:
            print(f"  gRNA {i} is {grna_length}bp (expected 20bp). Checking if it includes PAM...")
            # Try removing potential PAM sequences (last 3 bases if ends with GG)
            if grna.endswith('GG') or grna.endswith('CC'):
                grna_20bp = grna[:20]
                print(f"    Trying without PAM: {grna_20bp}")
            else:
                grna_20bp = grna[:20]
                print(f"    Trimming to 20bp: {grna_20bp}")
        elif grna_length < 20:
            print(f"  WARNING: gRNA {i} is only {grna_length}bp (expected 20bp)")
            grna_20bp = grna
        else:
            grna_20bp = grna
        
        # Search for exact match first
        found = False
        grna_rc = str(Seq(grna_20bp).reverse_complement())
        
        if grna_20bp in mrna_seq:
            pos = mrna_seq.find(grna_20bp)
            cut_site = pos + 17  # 3bp before PAM
            grna_positions.append((pos, 'forward', cut_site))
            expected_cut_sites.append(cut_site)
            print(f"  gRNA {i} found on forward strand at position {pos}, cut site: {cut_site}")
            found = True
        elif grna_rc in mrna_seq:
            pos = mrna_seq.find(grna_rc)
            cut_site = pos + 3
            grna_positions.append((pos, 'reverse', cut_site))
            expected_cut_sites.append(cut_site)
            print(f"  gRNA {i} found on reverse strand at position {pos}, cut site: {cut_site}")
            found = True
        
        # If not found, try partial match (check if part of gRNA is in mRNA)
        if not found:
            # Try searching for the gRNA without first base (common off-by-one error)
            if len(grna_20bp) >= 19:
                partial_grna = grna_20bp[1:]
                if partial_grna in mrna_seq:
                    pos = mrna_seq.find(partial_grna)
                    print(f"  WARNING: Partial match found (missing first base) at position {pos}")
                    print(f"    Expected: {grna_20bp}")
                    print(f"    Found:    {mrna_seq[pos-1:pos+len(partial_grna)]}")
                    found = True
            
            if not found:
                print(f"  ERROR: gRNA {i} ({original_grna}) not found in mRNA sequence!")
                print(f"    Searched for: {grna_20bp} and its reverse complement")
    
    if not grna_positions:
        print("  ERROR: No gRNAs found in the mRNA sequence!")
        print("  Cannot perform TIDE analysis without knowing the expected cut sites.")
        
        # Still generate recommendations but with error message
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        recommendations = analyze_grna_in_mrna(mrna_seq, grna_sequences, gene_name, efficiency_data=None)
        
        # Create output directory if it doesn't exist
        output_subdir = os.path.join(gene_folder, "output")
        if not os.path.exists(output_subdir):
            os.makedirs(output_subdir)
        
        # Save error recommendations
        rec_file = os.path.join(output_subdir, f"recommendations_{gene_name}_{timestamp}.txt")
        with open(rec_file, 'w', encoding='utf-8') as f:
            f.write(recommendations)
        
        print(f"  [ERROR] Analysis failed. Error report saved: {rec_file}")
        return False
    
    # Perform TIDE analysis
    try:
        # Use the expected cut site from gRNA position
        expected_cut_site = expected_cut_sites[0] if expected_cut_sites else None
        divergence_point, timestamp, efficiency_data = plot_tide_analysis_with_validation(
            control_file, edited_file, gene_folder, gene_name, expected_cut_site
        )
        
        # Generate recommendations
        recommendations = analyze_grna_in_mrna(mrna_seq, grna_sequences, gene_name, efficiency_data)
        
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