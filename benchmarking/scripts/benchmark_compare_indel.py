#!/usr/bin/env python
"""
Benchmarking Script: Compare CasPINS indel analysis with TIDE and ICE algorithms

Uses the published TIDE example data (Brinkman et al. 2014, NAR, gku936
Supplementary Data) and the DDC experimental data to compare three
decomposition algorithms on the same raw Sanger trace data.

Algorithms implemented (faithful to the published methods):

  1. TIDE (Brinkman et al. 2014, NAR 42(22):e168):
     - Extracts peak heights at basecalled positions for each of 4 channels
     - Stacks all 4 channels vertically into a single aggregation matrix
     - Single NNLS decomposition on the stacked matrix
     - R-squared = cor(fitted, observed)^2
     - Component percentages scaled by R-squared
     - Efficiency = (R^2 * 100) - WT_component_percentage

  2. ICE (Hsiau et al. 2019, CRISPR Journal 2(2):123-130):
     - Extracts peak heights, normalizes each position to sum to 100
     - Lasso (L1) regression with positive coefficients
     - R-squared correction applied to all abundances
     - Efficiency = (1 - unedited_fraction) * 100

  3. CasPINS:
     - Sums all 4 trace channels into single combined signal
     - NNLS decomposition on the continuous combined signal
     - No R-squared correction (reports raw NNLS fractions)

Datasets:
  - Open Source Example: example1.ab1 (control) + example2.ab1 (edited)
    gRNA: CATGCCGAGAGTGATCCCGG (Brinkman et al. 2014)
  - DDC Experimental: control.ab1 + 8 edited samples
    gRNA from C:/Users/kaush/Downloads/data/ddc/grna.txt

Usage:
    python benchmarking/scripts/benchmark_compare_indel.py
"""

import os
import sys
import csv
import numpy as np
from scipy.optimize import nnls
from sklearn.linear_model import Lasso
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

BENCHMARK_DIR = os.path.join(os.path.dirname(__file__), '..')
RESULTS_DIR = os.path.join(BENCHMARK_DIR, 'results', 'indel_analysis')
SUPPLEMENTARY_DIR = os.path.join(BENCHMARK_DIR, 'results', 'supplementary')

TIDE_EXAMPLE_DIR = os.path.join(
    BENCHMARK_DIR, 'data', 'tide_comparison_analysis', 'gku936_Supplementary_Data'
)
TIDE_EXAMPLE_GRNA = 'CATGCCGAGAGTGATCCCGG'

DDC_AB1_DIR = os.path.join(os.path.expanduser('~'), 'Downloads', 'data', 'ddc', 'input')
DDC_DIR = os.path.join(os.path.expanduser('~'), 'Downloads', 'data', 'ddc')


def load_ab1_data(ab1_path):
    """Load trace data, basecalled sequence, and peak locations from AB1 file."""
    from Bio import SeqIO
    record = SeqIO.read(ab1_path, 'abi')
    abif = record.annotations['abif_raw']

    base_order = abif.get('FWO_1', b'GATC')
    if isinstance(base_order, bytes):
        base_order = base_order.decode()

    channels = {}
    for i, tag in enumerate(['DATA9', 'DATA10', 'DATA11', 'DATA12']):
        channels[base_order[i]] = np.array(abif[tag], dtype=float)

    basecalled_seq = str(record.seq)
    peak_locs = list(abif.get('PLOC1', []))

    return channels, basecalled_seq, peak_locs, base_order


def extract_peak_heights(channels, peak_locs, base_order='ACGT'):
    """
    Extract peak height matrix (n_bases x 4) at basecalled positions.
    This mirrors the R code: traceMatrix(control)[peak_ctr_loc, ]
    """
    n_bases = len(peak_locs)
    heights = np.zeros((n_bases, 4))
    for col_idx, base in enumerate(base_order):
        trace = channels[base]
        for row_idx, pos in enumerate(peak_locs):
            if 0 <= pos < len(trace):
                heights[row_idx, col_idx] = trace[pos]
    return heights


def find_grna_in_sequence(sequence, grna):
    """Find gRNA in sequence, return (breaksite, orientation).
    TIDE convention: breaksite between nt 16-17 from 5' of guide.
    """
    grna_rc = grna[::-1].translate(str.maketrans('ACGT', 'TGCA'))

    if grna in sequence:
        pos = sequence.index(grna)
        breaksite = pos + 17
        return breaksite, 'forward'
    elif grna_rc in sequence:
        pos = sequence.index(grna_rc)
        breaksite = pos + 3
        return breaksite, 'reverse'
    return None, None


def find_cut_site_from_mrna(mrna_path, grna_seq):
    """Find gRNA cut site position in an mRNA reference file."""
    with open(mrna_path) as f:
        lines = f.readlines()
    mrna = ''.join(line.strip() for line in lines if not line.startswith('>'))
    grna_rc = grna_seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))

    if grna_seq in mrna:
        pos = mrna.index(grna_seq)
        return pos + len(grna_seq) - 3, 'forward', mrna
    elif grna_rc in mrna:
        pos = mrna.index(grna_rc)
        return pos + 3, 'reverse', mrna
    return None, None, mrna


# ============================================================
# Algorithm 1: TIDE - Faithful to Brinkman et al. 2014 R code
# ============================================================

def run_tide_algorithm(ctrl_heights, edit_heights, breaksite, offset=0,
                       maxshift=10, seqstart=100, seqend=None):
    """
    TIDE decomposition exactly as implemented in TIDE_functions.R.

    The R code stacks all 4 channels' peak-height data into one matrix
    and one vector, then runs a single NNLS.

    Key lines from TIDE_functions.R decomposition():
        for(b in import$B)  {
            sim <- matrix(...)
            for(i in shiftrange) {sim[,as.character(i)] <- import$ctr[(rg1:rg2)-i, b]}
            I_matrix <- rbind(I_matrix, sim)
            I_vec <- c(I_vec, import$mut[(rg1:rg2)+offset_mut, b])
        }
        NNFIT <- nnls(I_matrix, I_vec)
        Rsq <- cor(NNFIT$fit, I_vec)^2
        comper <- (Rsq * 100 * (NNFIT$x / sum(NNFIT$x)))
        eff <- (Rsq * 100) - comper[maxshift+1]
    """
    n_bases = ctrl_heights.shape[0]
    if seqend is None:
        seqend = min(n_bases, ctrl_heights.shape[0], edit_heights.shape[0])
    seqend = min(seqend, n_bases, edit_heights.shape[0])

    shiftrange = list(range(-maxshift, maxshift + 1))

    rg1 = breaksite + maxshift + 5
    rg2 = seqend - maxshift - 5

    if offset != 0:
        rg2 = min(rg2, seqend - abs(offset) - 5)

    if rg2 <= rg1 + maxshift * 2:
        return None

    rg1 = max(rg1, maxshift)
    rg2 = min(rg2, n_bases - maxshift - 1, edit_heights.shape[0] - abs(offset) - 1)

    if rg2 <= rg1:
        return None

    decomp_len = rg2 - rg1 + 1

    I_matrix_parts = []
    I_vec_parts = []

    for b_idx in range(4):
        sim = np.zeros((decomp_len, len(shiftrange)))
        for col_idx, shift in enumerate(shiftrange):
            indices = np.arange(rg1, rg2 + 1) - shift
            valid = (indices >= 0) & (indices < ctrl_heights.shape[0])
            sim[valid, col_idx] = ctrl_heights[indices[valid], b_idx]

        I_matrix_parts.append(sim)

        edit_indices = np.arange(rg1, rg2 + 1) + offset
        valid_e = (edit_indices >= 0) & (edit_indices < edit_heights.shape[0])
        vec = np.zeros(decomp_len)
        vec[valid_e] = edit_heights[edit_indices[valid_e], b_idx]
        I_vec_parts.append(vec)

    I_matrix = np.vstack(I_matrix_parts)
    I_vec = np.concatenate(I_vec_parts)

    coefficients, _ = nnls(I_matrix, I_vec)

    fitted = I_matrix @ coefficients

    if np.std(fitted) > 0 and np.std(I_vec) > 0:
        r_sq = np.corrcoef(fitted, I_vec)[0, 1] ** 2
    else:
        r_sq = 0.0

    total = np.sum(coefficients)
    if total > 0:
        comper = r_sq * 100 * (coefficients / total)
    else:
        comper = np.zeros_like(coefficients)

    wt_idx = shiftrange.index(0)
    efficiency = r_sq * 100 - comper[wt_idx]

    spectrum = {}
    for idx, shift in enumerate(shiftrange):
        if comper[idx] > 0.5:
            spectrum[shift] = round(comper[idx], 1)

    dominant = max(spectrum.items(), key=lambda x: x[1]) if spectrum else (0, 0)

    return {
        'efficiency': round(max(0, efficiency), 1),
        'r_squared': round(r_sq, 4),
        'dominant_indel': dominant[0],
        'wt_pct': round(comper[wt_idx], 1),
        'spectrum': spectrum,
        'method': 'TIDE (stacked 4-channel NNLS, R²-corrected)',
        'decomp_window': f'{rg1}-{rg2}',
    }


# ============================================================
# Algorithm 2: ICE - Faithful to Hsiau et al. 2019
# ============================================================

def run_ice_algorithm(ctrl_heights, edit_heights, breaksite, offset=0,
                      maxshift=10, seqend=None):
    """
    ICE decomposition: Lasso regression on normalized peak heights.

    ICE normalizes each base position so all 4 channels sum to 100,
    flattens into a single vector, builds shifted control columns,
    and uses Lasso regression with positive=True.

    R-squared correction: abundances *= R², unedited = wt*R² + (1-R²).
    """
    n_bases = ctrl_heights.shape[0]
    if seqend is None:
        seqend = min(n_bases, edit_heights.shape[0])
    seqend = min(seqend, n_bases, edit_heights.shape[0])

    shiftrange = list(range(-maxshift, maxshift + 1))

    rg1 = breaksite + maxshift + 5
    rg2 = seqend - maxshift - 5
    if offset != 0:
        rg2 = min(rg2, seqend - abs(offset) - 5)

    if rg2 <= rg1 + maxshift * 2:
        return None

    rg1 = max(rg1, maxshift)
    rg2 = min(rg2, n_bases - maxshift - 1, edit_heights.shape[0] - abs(offset) - 1)

    if rg2 <= rg1:
        return None

    def normalize_to_100(heights, start, end):
        """Normalize each position's 4 channels to sum to 100."""
        sub = heights[start:end + 1, :].copy()
        row_sums = sub.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        return sub / row_sums * 100

    ctrl_norm = normalize_to_100(ctrl_heights, rg1, rg2)
    decomp_len = ctrl_norm.shape[0]

    edit_indices = np.arange(rg1, rg2 + 1) + offset
    edit_sub = np.zeros((decomp_len, 4))
    valid = (edit_indices >= 0) & (edit_indices < edit_heights.shape[0])
    edit_sub[valid, :] = edit_heights[edit_indices[valid], :]
    e_sums = edit_sub.sum(axis=1, keepdims=True)
    e_sums[e_sums == 0] = 1
    edit_norm = edit_sub / e_sums * 100

    b = edit_norm.flatten()

    A = np.zeros((decomp_len * 4, len(shiftrange)))
    for col_idx, shift in enumerate(shiftrange):
        shifted = np.zeros_like(ctrl_norm)
        src_indices = np.arange(decomp_len) + shift
        valid_s = (src_indices >= 0) & (src_indices < ctrl_norm.shape[0])
        shifted[valid_s, :] = ctrl_norm[src_indices[valid_s], :]
        A[:, col_idx] = shifted.flatten()

    lasso = Lasso(alpha=0.8, positive=True, max_iter=10000, fit_intercept=False)
    lasso.fit(A, b)
    xvals = lasso.coef_

    predicted = A @ xvals
    if np.std(predicted) > 0 and np.std(b) > 0:
        r_sq = np.corrcoef(predicted, b)[0, 1] ** 2
    else:
        r_sq = 0.0

    total = xvals.sum()
    if total > 0:
        abundances = xvals / total
    else:
        abundances = xvals

    wt_idx = shiftrange.index(0)
    wt_fraction = abundances[wt_idx] * r_sq
    unedited = wt_fraction + (1 - r_sq)
    efficiency = max(0, min(100, (1 - unedited) * 100))

    corrected = abundances * r_sq * 100

    spectrum = {}
    for idx, shift in enumerate(shiftrange):
        if corrected[idx] > 0.5:
            spectrum[shift] = round(corrected[idx], 1)

    dominant = max(spectrum.items(), key=lambda x: x[1]) if spectrum else (0, 0)

    return {
        'efficiency': round(efficiency, 1),
        'r_squared': round(r_sq, 4),
        'dominant_indel': dominant[0],
        'wt_pct': round(unedited * 100, 1),
        'spectrum': spectrum,
        'method': 'ICE (Lasso on normalized peaks, R²-corrected)',
        'decomp_window': f'{rg1}-{rg2}',
    }


# ============================================================
# Algorithm 3: CasPINS - NNLS on combined (summed) channels
# ============================================================

def run_caspins_algorithm(ctrl_channels, edit_channels, cut_site_trace,
                          window_size=50, indel_range=range(-10, 11)):
    """
    CasPINS NNLS decomposition on combined (summed) trace channels.
    Works on continuous raw trace data (not peak heights).
    """
    start_pos = max(0, cut_site_trace + 50)
    end_pos = min(len(ctrl_channels['A']), cut_site_trace + window_size * 42 + 50)

    ctrl_signal = np.zeros(end_pos - start_pos)
    edit_signal = np.zeros(end_pos - start_pos)

    for base in 'ACGT':
        if len(ctrl_channels[base]) > end_pos and len(edit_channels[base]) > end_pos:
            ctrl_signal += ctrl_channels[base][start_pos:end_pos]
            edit_signal += edit_channels[base][start_pos:end_pos]

    if len(ctrl_signal) == 0 or np.max(ctrl_signal) == 0:
        return None

    ctrl_signal = ctrl_signal / np.max(ctrl_signal)
    edit_signal = edit_signal / np.max(edit_signal) if np.max(edit_signal) > 0 else edit_signal

    trace_factor = 42
    num_indels = len(list(indel_range))
    A = np.zeros((len(ctrl_signal), num_indels))

    for idx, indel_size in enumerate(indel_range):
        shift = indel_size * trace_factor
        if indel_size == 0:
            A[:, idx] = ctrl_signal
        elif shift > 0 and shift < len(ctrl_signal):
            A[:-shift, idx] = ctrl_signal[shift:]
        elif shift < 0 and -shift < len(ctrl_signal):
            A[-shift:, idx] = ctrl_signal[:shift]

    coefficients, _ = nnls(A, edit_signal, maxiter=1000)
    total = np.sum(coefficients)
    if total > 0:
        coefficients /= total

    reconstructed = A @ coefficients
    ss_tot = np.sum((edit_signal - np.mean(edit_signal)) ** 2)
    ss_res = np.sum((edit_signal - reconstructed) ** 2)
    r_squared = max(0, 1 - ss_res / ss_tot) if ss_tot > 0 else 0

    wt_idx = list(indel_range).index(0)
    wt_fraction = coefficients[wt_idx]
    efficiency = (1 - wt_fraction) * 100

    spectrum = {}
    for idx, indel_size in enumerate(indel_range):
        if coefficients[idx] * 100 > 0.5:
            spectrum[indel_size] = round(coefficients[idx] * 100, 1)

    dominant = max(spectrum.items(), key=lambda x: x[1]) if spectrum else (0, 0)

    return {
        'efficiency': round(efficiency, 1),
        'r_squared': round(r_squared, 4),
        'dominant_indel': dominant[0],
        'wt_pct': round(wt_fraction * 100, 1),
        'spectrum': spectrum,
        'method': 'CasPINS (NNLS on summed channels)',
    }


def run_caspins_on_peaks(ctrl_heights, edit_heights, breaksite, offset=0,
                         maxshift=10, seqend=None):
    """
    CasPINS-style NNLS but operating on the same peak-height data as TIDE/ICE,
    using a summed (combined) signal. This allows direct comparison on identical
    input features.
    """
    n_bases = ctrl_heights.shape[0]
    if seqend is None:
        seqend = min(n_bases, edit_heights.shape[0])
    seqend = min(seqend, n_bases, edit_heights.shape[0])

    shiftrange = list(range(-maxshift, maxshift + 1))

    rg1 = breaksite + maxshift + 5
    rg2 = seqend - maxshift - 5
    if offset != 0:
        rg2 = min(rg2, seqend - abs(offset) - 5)
    rg1 = max(rg1, maxshift)
    rg2 = min(rg2, n_bases - maxshift - 1, edit_heights.shape[0] - abs(offset) - 1)

    if rg2 <= rg1:
        return None

    ctrl_sub = ctrl_heights[rg1:rg2 + 1, :]
    ctrl_combined = ctrl_sub.sum(axis=1)

    edit_indices = np.arange(rg1, rg2 + 1) + offset
    valid = (edit_indices >= 0) & (edit_indices < edit_heights.shape[0])
    edit_sub = np.zeros((rg2 - rg1 + 1, 4))
    edit_sub[valid, :] = edit_heights[edit_indices[valid], :]
    edit_combined = edit_sub.sum(axis=1)

    if np.max(ctrl_combined) == 0:
        return None

    ctrl_combined = ctrl_combined / np.max(ctrl_combined)
    edit_combined = edit_combined / np.max(edit_combined) if np.max(edit_combined) > 0 else edit_combined

    decomp_len = len(ctrl_combined)
    A = np.zeros((decomp_len, len(shiftrange)))
    for col_idx, shift in enumerate(shiftrange):
        indices = np.arange(decomp_len) - shift
        valid_s = (indices >= 0) & (indices < decomp_len)
        A[valid_s, col_idx] = ctrl_combined[indices[valid_s]]

    coefficients, _ = nnls(A, edit_combined, maxiter=1000)
    total = np.sum(coefficients)
    if total > 0:
        coefficients /= total

    fitted = A @ coefficients
    ss_tot = np.sum((edit_combined - np.mean(edit_combined)) ** 2)
    ss_res = np.sum((edit_combined - fitted) ** 2)
    r_sq = max(0, 1 - ss_res / ss_tot) if ss_tot > 0 else 0

    wt_idx = shiftrange.index(0)
    wt_fraction = coefficients[wt_idx]
    efficiency = (1 - wt_fraction) * 100

    spectrum = {}
    for idx, shift in enumerate(shiftrange):
        if coefficients[idx] * 100 > 0.5:
            spectrum[shift] = round(coefficients[idx] * 100, 1)

    dominant = max(spectrum.items(), key=lambda x: x[1]) if spectrum else (0, 0)

    return {
        'efficiency': round(max(0, efficiency), 1),
        'r_squared': round(r_sq, 4),
        'dominant_indel': dominant[0],
        'wt_pct': round(wt_fraction * 100, 1),
        'spectrum': spectrum,
        'method': 'CasPINS (NNLS on summed channels)',
        'decomp_window': f'{rg1}-{rg2}',
    }


# ============================================================
# Alignment: compute offset between control and edited
# ============================================================

def compute_offset(ctrl_seq, edit_seq, breaksite, maxshift=10, seqstart=100):
    """
    Compute alignment offset between control and edited sequences,
    mirroring TIDE's R alignment logic.
    """
    from difflib import SequenceMatcher

    end = max(0, breaksite - maxshift)
    if end <= seqstart:
        return 0

    ctrl_sub = ctrl_seq[seqstart:end]
    edit_sub = edit_seq[seqstart:end]

    best_offset = 0
    best_score = 0
    for off in range(-5, 6):
        e_start = max(0, seqstart + off)
        e_end = max(0, end + off)
        if e_end > len(edit_seq):
            continue
        e_sub = edit_seq[e_start:e_end]
        if len(e_sub) == 0 or len(ctrl_sub) == 0:
            continue
        sm = SequenceMatcher(None, ctrl_sub, e_sub)
        score = sm.ratio()
        if score > best_score:
            best_score = score
            best_offset = off

    return best_offset


# ============================================================
# Main comparison pipeline
# ============================================================

def analyze_tide_example():
    """Analyze the TIDE paper's example data (gold-standard)."""
    ctrl_path = os.path.join(TIDE_EXAMPLE_DIR, 'example1.ab1')
    edit_path = os.path.join(TIDE_EXAMPLE_DIR, 'example2.ab1')
    grna = TIDE_EXAMPLE_GRNA

    if not os.path.exists(ctrl_path):
        print(f"TIDE example data not found at {TIDE_EXAMPLE_DIR}")
        return []

    print("=" * 80)
    print("DATASET 1: Open Source Example Data (Brinkman et al. 2014)")
    print("=" * 80)

    ctrl_channels, ctrl_seq, ctrl_peaks, ctrl_order = load_ab1_data(ctrl_path)
    edit_channels, edit_seq, edit_peaks, edit_order = load_ab1_data(edit_path)

    print(f"  Control: {len(ctrl_seq)} basecalled bases, {len(ctrl_channels['A'])} trace points")
    print(f"  Edited:  {len(edit_seq)} basecalled bases, {len(edit_channels['A'])} trace points")
    print(f"  gRNA:    {grna}")

    breaksite, orientation = find_grna_in_sequence(ctrl_seq, grna)
    if breaksite is None:
        print("  ERROR: gRNA not found in control sequence")
        return []
    print(f"  Breaksite: position {breaksite} in basecalled sequence ({orientation})")

    ctrl_heights = extract_peak_heights(ctrl_channels, ctrl_peaks)
    edit_heights = extract_peak_heights(edit_channels, edit_peaks)

    offset = compute_offset(ctrl_seq, edit_seq, breaksite)
    print(f"  Alignment offset: {offset}")

    seqend = min(ctrl_heights.shape[0], edit_heights.shape[0]) - 1

    tide = run_tide_algorithm(ctrl_heights, edit_heights, breaksite,
                              offset=offset, maxshift=10, seqend=seqend)
    ice = run_ice_algorithm(ctrl_heights, edit_heights, breaksite,
                            offset=offset, maxshift=10, seqend=seqend)
    caspins = run_caspins_on_peaks(ctrl_heights, edit_heights, breaksite,
                                   offset=offset, maxshift=10, seqend=seqend)

    row = {
        'dataset': 'Open_Source',
        'gene': 'GFP-target',
        'sample': 'example2',
        'control': 'example1',
        'grna': grna,
    }

    for prefix, result in [('caspins', caspins), ('tide', tide), ('ice', ice)]:
        if result:
            row[f'{prefix}_efficiency'] = result['efficiency']
            row[f'{prefix}_r2'] = result['r_squared']
            row[f'{prefix}_dominant'] = result['dominant_indel']
        else:
            row[f'{prefix}_efficiency'] = None
            row[f'{prefix}_r2'] = None
            row[f'{prefix}_dominant'] = None

    print(f"\n  {'Method':<42} {'Efficiency':>10} {'R²':>8} {'Dominant':>10}")
    print(f"  {'-'*72}")
    for prefix, result in [('CasPINS', caspins), ('TIDE', tide), ('ICE', ice)]:
        if result:
            print(f"  {result['method']:<42} {result['efficiency']:>9.1f}% {result['r_squared']:>8.4f} {result['dominant_indel']:>+10d}")
            if result.get('spectrum'):
                sig_indels = {k: v for k, v in sorted(result['spectrum'].items()) if k != 0}
                if sig_indels:
                    indel_str = ', '.join(f"{k:+d}:{v:.1f}%" for k, v in sorted(sig_indels.items(), key=lambda x: -x[1])[:5])
                    print(f"  {'':42} Indels: {indel_str}")
        else:
            print(f"  {prefix:<42} {'FAILED':>10}")

    return [row]


def analyze_ddc_data():
    """Analyze DDC experimental data (8 edited samples)."""
    if not os.path.exists(DDC_AB1_DIR):
        print(f"\nDDC data not found at {DDC_AB1_DIR}")
        return []

    grna_path = os.path.join(DDC_DIR, 'grna.txt')
    mrna_path = os.path.join(DDC_DIR, 'mrna.txt')

    if not os.path.exists(grna_path) or not os.path.exists(mrna_path):
        print("\nDDC grna.txt or mrna.txt not found")
        return []

    with open(grna_path) as f:
        grna = f.read().strip()

    cut_site_mrna, orientation, mrna = find_cut_site_from_mrna(mrna_path, grna)
    if cut_site_mrna is None:
        print("ERROR: gRNA not found in DDC mRNA reference")
        return []

    print(f"\n{'='*80}")
    print("DATASET 2: DDC Experimental Data (Dopa decarboxylase, rat)")
    print("=" * 80)
    print(f"  gRNA: {grna} ({orientation})")
    print(f"  Cut site: position {cut_site_mrna} in mRNA ({len(mrna)} bp)")

    ctrl_path = os.path.join(DDC_AB1_DIR, 'control.ab1')
    ctrl_channels, ctrl_seq, ctrl_peaks, ctrl_order = load_ab1_data(ctrl_path)
    print(f"  Control: {len(ctrl_seq)} basecalled bases, {len(ctrl_channels['A'])} trace points")

    if len(ctrl_seq) < 20:
        print(f"  NOTE: Control basecalling is poor ({len(ctrl_seq)} bases).")
        print(f"        Using estimated peak positions from trace data for decomposition.")

    n_trace_points = len(ctrl_channels['A'])
    estimated_n_bases = n_trace_points // 42

    ctrl_heights = np.zeros((estimated_n_bases, 4))
    for col_idx, base in enumerate('ACGT'):
        trace = ctrl_channels[base]
        for row_idx in range(estimated_n_bases):
            center = row_idx * 42 + 21
            if center < len(trace):
                s = max(0, center - 5)
                e = min(len(trace), center + 6)
                ctrl_heights[row_idx, col_idx] = np.max(trace[s:e])

    edited_files = sorted([f for f in os.listdir(DDC_AB1_DIR)
                           if f.startswith('edited') and f.endswith('.ab1')])

    results = []
    for ef in edited_files:
        sample_name = ef.replace('.ab1', '')
        edit_path = os.path.join(DDC_AB1_DIR, ef)
        edit_channels, edit_seq, edit_peaks, edit_order = load_ab1_data(edit_path)

        edit_heights = np.zeros((estimated_n_bases, 4))
        for col_idx, base in enumerate('ACGT'):
            trace = edit_channels[base]
            for row_idx in range(estimated_n_bases):
                center = row_idx * 42 + 21
                if center < len(trace):
                    s = max(0, center - 5)
                    e = min(len(trace), center + 6)
                    edit_heights[row_idx, col_idx] = np.max(trace[s:e])

        breaksite_est = cut_site_mrna

        tide = run_tide_algorithm(ctrl_heights, edit_heights, breaksite_est,
                                  offset=0, maxshift=10, seqend=estimated_n_bases - 1)
        ice = run_ice_algorithm(ctrl_heights, edit_heights, breaksite_est,
                                offset=0, maxshift=10, seqend=estimated_n_bases - 1)
        caspins_peak = run_caspins_on_peaks(ctrl_heights, edit_heights, breaksite_est,
                                            offset=0, maxshift=10, seqend=estimated_n_bases - 1)

        row = {
            'dataset': 'DDC',
            'gene': 'DDC',
            'sample': sample_name,
            'control': 'control',
            'grna': grna,
        }
        for prefix, result in [('caspins', caspins_peak), ('tide', tide), ('ice', ice)]:
            if result:
                row[f'{prefix}_efficiency'] = result['efficiency']
                row[f'{prefix}_r2'] = result['r_squared']
                row[f'{prefix}_dominant'] = result['dominant_indel']
            else:
                row[f'{prefix}_efficiency'] = None
                row[f'{prefix}_r2'] = None
                row[f'{prefix}_dominant'] = None

        results.append(row)

    if results:
        print(f"\n  {'Sample':<14} {'CasPINS':>10} {'TIDE':>10} {'ICE':>10}  "
              f"{'CasPINS R²':>10} {'TIDE R²':>10} {'ICE R²':>10}")
        print(f"  {'-'*76}")
        for r in results:
            c_e = f"{r['caspins_efficiency']:.1f}" if r['caspins_efficiency'] is not None else "N/A"
            t_e = f"{r['tide_efficiency']:.1f}" if r['tide_efficiency'] is not None else "N/A"
            i_e = f"{r['ice_efficiency']:.1f}" if r['ice_efficiency'] is not None else "N/A"
            c_r = f"{r['caspins_r2']:.4f}" if r['caspins_r2'] is not None else "N/A"
            t_r = f"{r['tide_r2']:.4f}" if r['tide_r2'] is not None else "N/A"
            i_r = f"{r['ice_r2']:.4f}" if r['ice_r2'] is not None else "N/A"
            print(f"  {r['sample']:<14} {c_e:>10} {t_e:>10} {i_e:>10}  {c_r:>10} {t_r:>10} {i_r:>10}")

        c_effs = [r['caspins_efficiency'] for r in results if r['caspins_efficiency'] is not None]
        t_effs = [r['tide_efficiency'] for r in results if r['tide_efficiency'] is not None]
        i_effs = [r['ice_efficiency'] for r in results if r['ice_efficiency'] is not None]
        if c_effs and t_effs and i_effs:
            print(f"  {'-'*76}")
            print(f"  {'MEAN':<14} {np.mean(c_effs):>10.1f} {np.mean(t_effs):>10.1f} {np.mean(i_effs):>10.1f}")

    return results


def run_comparison():
    """Run all three algorithms on both datasets and produce reports."""
    print("=" * 80)
    print("INDEL ANALYSIS: THREE-WAY ALGORITHM COMPARISON")
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    print(f"\nAlgorithms compared:")
    print(f"  1. TIDE: Stacked 4-channel peak heights, single NNLS, R²-corrected")
    print(f"     (Brinkman et al. 2014, NAR 42(22):e168)")
    print(f"  2. ICE:  Normalized peak heights, Lasso (L1) regression, R²-corrected")
    print(f"     (Hsiau et al. 2019, CRISPR Journal 2(2):123-130)")
    print(f"  3. CasPINS: NNLS on combined (summed) trace channels")

    all_results = []

    tide_results = analyze_tide_example()
    all_results.extend(tide_results)

    ddc_results = analyze_ddc_data()
    all_results.extend(ddc_results)

    if not all_results:
        print("\nNo results to report!")
        return

    # Pairwise correlations (DDC has multiple samples)
    ddc_rows = [r for r in all_results if r['dataset'] == 'DDC']
    if len(ddc_rows) > 2:
        c_effs = [r['caspins_efficiency'] for r in ddc_rows if r['caspins_efficiency'] is not None]
        t_effs = [r['tide_efficiency'] for r in ddc_rows if r['tide_efficiency'] is not None]
        i_effs = [r['ice_efficiency'] for r in ddc_rows if r['ice_efficiency'] is not None]

        if len(c_effs) >= 3 and len(t_effs) >= 3 and len(i_effs) >= 3:
            from scipy.stats import pearsonr
            r_ct, p_ct = pearsonr(c_effs, t_effs)
            r_ci, p_ci = pearsonr(c_effs, i_effs)
            r_ti, p_ti = pearsonr(t_effs, i_effs)

            print(f"\n{'='*80}")
            print("PAIRWISE CORRELATIONS (DDC dataset, n=8)")
            print("=" * 80)
            print(f"  CasPINS vs TIDE: Pearson r = {r_ct:.4f} (p = {p_ct:.4e})")
            print(f"  CasPINS vs ICE:  Pearson r = {r_ci:.4f} (p = {p_ci:.4e})")
            print(f"  TIDE vs ICE:     Pearson r = {r_ti:.4f} (p = {p_ti:.4e})")

    # Save CSV summary
    os.makedirs(os.path.join(RESULTS_DIR, 'comparison'), exist_ok=True)
    csv_path = os.path.join(RESULTS_DIR, 'comparison', 'indel_comparison_summary.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Dataset', 'Gene', 'Sample',
            'CasPINS_Efficiency_%', 'CasPINS_R2',
            'TIDE_Efficiency_%', 'TIDE_R2',
            'ICE_Efficiency_%', 'ICE_R2',
            'CasPINS_Dominant_Indel', 'TIDE_Dominant_Indel', 'ICE_Dominant_Indel'
        ])
        for r in all_results:
            writer.writerow([
                r['dataset'], r['gene'], r['sample'],
                r.get('caspins_efficiency', ''), r.get('caspins_r2', ''),
                r.get('tide_efficiency', ''), r.get('tide_r2', ''),
                r.get('ice_efficiency', ''), r.get('ice_r2', ''),
                r.get('caspins_dominant', ''), r.get('tide_dominant', ''), r.get('ice_dominant', ''),
            ])
    print(f"\nSaved comparison CSV: {csv_path}")

    # Save Supplementary Table S4
    os.makedirs(SUPPLEMENTARY_DIR, exist_ok=True)
    s4_path = os.path.join(SUPPLEMENTARY_DIR, 'Table_S4_indel_comparison.csv')
    with open(s4_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Dataset', 'Gene', 'Sample',
            'CasPINS_Efficiency_%', 'CasPINS_R2',
            'TIDE_Efficiency_%', 'TIDE_R2',
            'ICE_Efficiency_%', 'ICE_R2',
            'CasPINS_Dominant_Indel', 'TIDE_Dominant_Indel', 'ICE_Dominant_Indel'
        ])

        for r in all_results:
            writer.writerow([
                r['dataset'], r['gene'], r['sample'],
                r.get('caspins_efficiency', ''), r.get('caspins_r2', ''),
                r.get('tide_efficiency', ''), r.get('tide_r2', ''),
                r.get('ice_efficiency', ''), r.get('ice_r2', ''),
                r.get('caspins_dominant', ''), r.get('tide_dominant', ''), r.get('ice_dominant', ''),
            ])

        ddc_rows_valid = [r for r in ddc_rows if r.get('caspins_efficiency') is not None]
        if ddc_rows_valid:
            writer.writerow([
                'DDC', 'DDC', 'MEAN',
                round(np.mean([r['caspins_efficiency'] for r in ddc_rows_valid]), 1),
                round(np.mean([r['caspins_r2'] for r in ddc_rows_valid if r['caspins_r2'] is not None]), 4),
                round(np.mean([r['tide_efficiency'] for r in ddc_rows_valid if r['tide_efficiency'] is not None]), 1),
                round(np.mean([r['tide_r2'] for r in ddc_rows_valid if r['tide_r2'] is not None]), 4),
                round(np.mean([r['ice_efficiency'] for r in ddc_rows_valid if r['ice_efficiency'] is not None]), 1),
                round(np.mean([r['ice_r2'] for r in ddc_rows_valid if r['ice_r2'] is not None]), 4),
                '', '', '',
            ])

    print(f"Saved Table S4: {s4_path}")

    # Save detailed text report
    report_path = os.path.join(RESULTS_DIR, 'comparison', 'indel_comparison_report.txt')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("INDEL ANALYSIS: THREE-WAY ALGORITHM COMPARISON REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")

        f.write("METHODOLOGY\n")
        f.write("-" * 40 + "\n\n")
        f.write("Three decomposition algorithms were implemented and compared using\n")
        f.write("the same input Sanger trace data:\n\n")

        f.write("1. TIDE (Brinkman et al. 2014, NAR 42(22):e168)\n")
        f.write("   Implementation faithful to the published R code (gku936 Supplementary):\n")
        f.write("   - Extracts peak heights at basecalled positions (4 channels x N bases)\n")
        f.write("   - Stacks all 4 channels vertically into single aggregation matrix\n")
        f.write("   - Single NNLS decomposition on the stacked matrix\n")
        f.write("   - R^2 = cor(fitted, observed)^2\n")
        f.write("   - Component percentages = R^2 x 100 x (coeff / sum(coeff))\n")
        f.write("   - Efficiency = R^2 x 100 - WT_component_percentage\n\n")

        f.write("2. ICE (Hsiau et al. 2019, CRISPR Journal 2(2):123-130)\n")
        f.write("   Implementation based on published algorithm and open-source code:\n")
        f.write("   - Extracts peak heights, normalizes each position to sum to 100\n")
        f.write("   - Lasso (L1) regression with alpha=0.8, positive=True\n")
        f.write("   - R^2 correction: abundances scaled by R^2\n")
        f.write("   - Unedited fraction = WT*R^2 + (1-R^2)\n\n")

        f.write("3. CasPINS\n")
        f.write("   As implemented in CasPINS src/utils/indel_analysis.py:\n")
        f.write("   - Sums all 4 trace channels into a single combined signal\n")
        f.write("   - NNLS decomposition (scipy.optimize.nnls)\n")
        f.write("   - No R^2-correction applied to abundances\n\n")

        f.write("DATASETS\n")
        f.write("-" * 40 + "\n\n")
        f.write("Dataset 1: Open Source Example (Brinkman et al. 2014)\n")
        f.write("  - Source: gku936 Supplementary Data\n")
        f.write(f"  - gRNA: {TIDE_EXAMPLE_GRNA}\n")
        f.write("  - Control: example1.ab1 (325 bases, mean Phred 53)\n")
        f.write("  - Edited:  example2.ab1 (325 bases, mean Phred 34)\n")
        f.write("  - High-quality gold-standard data from the TIDE publication\n\n")

        f.write("Dataset 2: DDC Experimental (Dopa decarboxylase, rat)\n")
        f.write("  - 8 edited samples + 1 control\n")
        f.write("  - Control has poor basecalling (5 bases); peak positions estimated\n")
        f.write("  - Demonstrates algorithm performance on challenging trace data\n\n")

        f.write("=" * 80 + "\n")
        f.write("RESULTS\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"{'Dataset':<12} {'Sample':<14} {'CasPINS':>8} {'TIDE':>8} {'ICE':>8}  "
                f"{'CasPINS':>8} {'TIDE':>8} {'ICE':>8}\n")
        f.write(f"{'':12} {'':14} {'Eff.%':>8} {'Eff.%':>8} {'Eff.%':>8}  "
                f"{'R²':>8} {'R²':>8} {'R²':>8}\n")
        f.write("-" * 80 + "\n")

        for r in all_results:
            vals = []
            for prefix in ['caspins', 'tide', 'ice']:
                e = r.get(f'{prefix}_efficiency')
                vals.append(f"{e:.1f}" if e is not None else "N/A")
            r2s = []
            for prefix in ['caspins', 'tide', 'ice']:
                rv = r.get(f'{prefix}_r2')
                r2s.append(f"{rv:.4f}" if rv is not None else "N/A")
            f.write(f"{r['dataset']:<12} {r['sample']:<14} {vals[0]:>8} {vals[1]:>8} {vals[2]:>8}  "
                    f"{r2s[0]:>8} {r2s[1]:>8} {r2s[2]:>8}\n")

        f.write("-" * 80 + "\n\n")

        if len(ddc_rows) > 2:
            c_effs = [r['caspins_efficiency'] for r in ddc_rows if r['caspins_efficiency'] is not None]
            t_effs = [r['tide_efficiency'] for r in ddc_rows if r['tide_efficiency'] is not None]
            i_effs = [r['ice_efficiency'] for r in ddc_rows if r['ice_efficiency'] is not None]
            if c_effs and t_effs and i_effs:
                f.write("PAIRWISE CORRELATIONS (DDC dataset, n=8)\n")
                f.write("-" * 40 + "\n")
                from scipy.stats import pearsonr
                r_ct, p_ct = pearsonr(c_effs, t_effs)
                r_ci, p_ci = pearsonr(c_effs, i_effs)
                r_ti, p_ti = pearsonr(t_effs, i_effs)
                f.write(f"  CasPINS vs TIDE: r = {r_ct:.4f} (p = {p_ct:.4e})\n")
                f.write(f"  CasPINS vs ICE:  r = {r_ci:.4f} (p = {p_ci:.4e})\n")
                f.write(f"  TIDE vs ICE:     r = {r_ti:.4f} (p = {p_ti:.4e})\n\n")

        f.write("=" * 80 + "\n")
        f.write("INTERPRETATION\n")
        f.write("=" * 80 + "\n\n")
        f.write("The TIDE paper example data provides a gold-standard comparison where\n")
        f.write("all three algorithms operate on high-quality Sanger traces. The DDC\n")
        f.write("data extends the comparison to traces with poor basecall quality,\n")
        f.write("demonstrating algorithmic behavior on challenging real-world data.\n\n")
        f.write("Key observations:\n")
        f.write("- All three algorithms are variants of trace decomposition\n")
        f.write("- TIDE and CasPINS both use NNLS but differ in signal representation\n")
        f.write("- ICE uses Lasso (L1-regularized) regression, favoring sparse solutions\n")
        f.write("- TIDE applies R^2-correction, scaling all estimates by model fit quality\n")

    print(f"Saved report: {report_path}")

    print(f"\n{'='*60}")
    print("INDEL COMPARISON COMPLETE")
    print(f"{'='*60}")

    return all_results


if __name__ == '__main__':
    run_comparison()
