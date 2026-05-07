#!/usr/bin/env python
"""
Comparison Script: CasPINS vs CRISPOR vs CHOPCHOP gRNA Rankings

This script compares gRNA design results from CasPINS, CRISPOR, and CHOPCHOP
to generate concordance statistics and visualizations for the manuscript.

Prerequisites:
    1. Run benchmark_grna_design.py first (generates CasPINS results)
    2. Manually query CRISPOR and CHOPCHOP (see CRISPOR_CHOPCHOP_QUERY_INSTRUCTIONS.md)
    3. Place CRISPOR results in: benchmarking/data/CRISPRor/crispror_{gene}.xls
    4. Place CHOPCHOP results in: benchmarking/data/CHOPCHOP/chopchop_{gene}.tsv

Usage:
    python benchmark_compare_grna.py
"""

import sys
import os
import json
import csv
import glob
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

BENCHMARK_DIR = os.path.join(os.path.dirname(__file__), '..')
RESULTS_DIR = os.path.join(BENCHMARK_DIR, 'results', 'grna_design')
DATA_DIR = os.path.join(BENCHMARK_DIR, 'data')
CRISPOR_DIR = os.path.join(DATA_DIR, 'CRISPRor')
CHOPCHOP_DIR = os.path.join(DATA_DIR, 'CHOPCHOP')


def parse_genomic_coord(loc_str: str):
    """Parse a genomic location string like '17:7669567-7669587' and return (chrom, start, end)."""
    if not loc_str or ':' not in loc_str:
        return None, None, None
    parts = loc_str.split(':')
    chrom = parts[0].replace('chr', '')
    range_parts = parts[-1].split('-')
    try:
        start = int(range_parts[0])
        end = int(range_parts[1]) if len(range_parts) > 1 else start + 20
        return chrom, start, end
    except (ValueError, IndexError):
        return None, None, None


def load_caspins_results(results_dir: str) -> dict:
    """Load CasPINS benchmark results with genomic coordinate information."""
    json_path = os.path.join(results_dir, 'caspins_grna_benchmark_results.json')
    if not os.path.exists(json_path):
        print(f"ERROR: CasPINS results not found at {json_path}")
        print("Run benchmark_grna_design.py first.")
        sys.exit(1)
    
    with open(json_path) as f:
        results = json.load(f)
    
    # Convert to gene -> data mapping (including coordinates for region filtering)
    gene_grnas = {}
    for gene_result in results:
        gene = gene_result['gene']
        grnas = gene_result.get('grnas', [])
        
        sequences = []
        coords = []  # parallel list of (chrom, start, end) for each gRNA
        scores = {}
        doench_scores = {}
        
        for g in grnas:
            seq = g['sequence'].upper()
            sequences.append(seq)
            scores[seq] = g.get('composite_score', 0)
            doench_scores[seq] = g.get('doench_2016_score', 0)
            
            loc = g.get('genomic_location', '')
            chrom, start, end = parse_genomic_coord(loc)
            coords.append((chrom, start, end))
        
        gene_grnas[gene] = {
            'sequences': sequences,
            'coords': coords,
            'scores': scores,
            'doench_scores': doench_scores,
            'n_total_found': gene_result.get('n_total_found', len(grnas)),
        }
    
    return gene_grnas


def load_crispor_results(data_dir: str, gene: str) -> dict:
    """
    Load CRISPOR results from.xls Excel file or.tsv file.
    
    Searches in benchmarking/data/CRISPRor/ for files named crispror_{gene}.xls
    Also falls back to benchmarking/data/ for crispor_{gene}.tsv (legacy format).
    
    Expected CRISPOR columns:
    - targetSeq (or guideSeq, Guide Sequence) - the 20nt gRNA sequence
    - doenchScore (or Doench '16) - Doench 2016 score
    """
    # Try.xls files in CRISPRor/ subfolder first (current data layout)
    filepath = None
    file_format = None
    
    for candidate in [
        (os.path.join(CRISPOR_DIR, f'crispror_{gene.lower()}.xls'), 'xls'),
        (os.path.join(CRISPOR_DIR, f'crispror_{gene}.xls'), 'xls'),
        (os.path.join(CRISPOR_DIR, f'crispor_{gene.lower()}.xls'), 'xls'),
        (os.path.join(data_dir, f'crispor_{gene.lower()}.tsv'), 'tsv'),
        (os.path.join(data_dir, f'crispor_{gene}.tsv'), 'tsv'),
    ]:
        if os.path.exists(candidate[0]):
            filepath = candidate[0]
            file_format = candidate[1]
            break
    
    if filepath is None:
        return None
    
    sequences = []
    scores = {}
    doench_scores = {}
    
    seq_columns = ['targetSeq', 'guideSeq', 'Guide Sequence', 'guide_sequence']
    doench_columns = ['doenchScore', "Doench '16", 'Doench2016', 'doench_score']
    
    if file_format == 'xls':
        # Read binary Excel (.xls) using xlrd
        try:
            import xlrd
        except ImportError:
            print(f"  ERROR: xlrd is required to read.xls files. Install with: pip install xlrd")
            return None
        
        try:
            workbook = xlrd.open_workbook(filepath)
            sheet = workbook.sheet_by_index(0)
            
            if sheet.nrows < 2:
                print(f"  WARNING: Empty spreadsheet: {filepath}")
                return None
            
            # CRISPOR.xls files have metadata rows (# Name, # Sequence, etc.)
            # before the actual header row. Find the real header row by looking
            # for a row that contains 'targetSeq' or '#guideId'.
            header_row_idx = None
            for row_idx in range(min(20, sheet.nrows)):
                row_values = [str(sheet.cell_value(row_idx, c)).strip() 
                             for c in range(sheet.ncols)]
                if any(v in seq_columns or v == '#guideId' for v in row_values):
                    header_row_idx = row_idx
                    break
            
            if header_row_idx is None:
                print(f"  WARNING: No header row with gRNA columns found in {filepath}")
                # Show first few rows for debugging
                for r in range(min(5, sheet.nrows)):
                    vals = [str(sheet.cell_value(r, c)).strip() for c in range(min(5, sheet.ncols))]
                    print(f"    Row {r}: {vals}")
                return None
            
            headers = [str(sheet.cell_value(header_row_idx, col)).strip() 
                      for col in range(sheet.ncols)]
            
            # Find relevant column indices
            seq_col_idx = None
            doench_col_idx = None
            
            for i, h in enumerate(headers):
                if seq_col_idx is None and h in seq_columns:
                    seq_col_idx = i
                if doench_col_idx is None and h in doench_columns:
                    doench_col_idx = i
            
            if seq_col_idx is None:
                print(f"  WARNING: No gRNA sequence column found in {filepath}")
                print(f"  Available columns: {headers}")
                return None
            
            # Read data rows (start after header)
            for row_idx in range(header_row_idx + 1, sheet.nrows):
                raw_seq = str(sheet.cell_value(row_idx, seq_col_idx)).upper().replace(' ', '')
                if not raw_seq or len(raw_seq) < 20:
                    continue
                # CRISPOR targetSeq is 23nt (20nt guide + 3nt PAM) - take first 20
                seq = raw_seq[:20]
                sequences.append(seq)
                
                if doench_col_idx is not None:
                    try:
                        doench_scores[seq] = float(sheet.cell_value(row_idx, doench_col_idx))
                    except (ValueError, TypeError):
                        pass
        except Exception as e:
            print(f"  ERROR reading {filepath}: {e}")
            return None
    
    else:
        # Read TSV format (legacy)
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                seq = None
                for col in seq_columns:
                    if col in row:
                        seq = row[col].upper().replace(' ', '')[:20]
                        break
                
                if seq and len(seq) >= 20:
                    sequences.append(seq)
                    
                    for col in doench_columns:
                        if col in row:
                            try:
                                doench_scores[seq] = float(row[col])
                            except (ValueError, TypeError):
                                pass
                            break
    
    if not sequences:
        print(f"  WARNING: No gRNA sequences parsed from {filepath}")
        return None
    
    return {
        'sequences': sequences,
        'scores': scores,
        'doench_scores': doench_scores
    }


def load_chopchop_results(data_dir: str, gene: str) -> dict:
    """
    Load CHOPCHOP results from TSV file.
    
    Searches in benchmarking/data/CHOPCHOP/ subfolder first, then falls back
    to benchmarking/data/ (legacy location). Handles filenames with stray spaces.
    
    Expected CHOPCHOP output format:
    - Column: Target sequence (or sgRNA) - the gRNA sequence  
    - Column: Efficiency - efficiency score
    """
    filepath = None
    gene_lower = gene.lower()
    
    # Search candidates: CHOPCHOP/ subfolder first, then data_dir
    search_dirs = [CHOPCHOP_DIR, data_dir]
    
    for search_dir in search_dirs:
        if not os.path.isdir(search_dir):
            continue
        
        # Exact match
        for candidate in [
            os.path.join(search_dir, f'chopchop_{gene_lower}.tsv'),
            os.path.join(search_dir, f'chopchop_{gene}.tsv'),
        ]:
            if os.path.exists(candidate):
                filepath = candidate
                break
        
        if filepath:
            break
        
        # Handle filenames with stray spaces (e.g., "chopchop_vegfa.tsv")
        pattern = os.path.join(search_dir, f'chopchop_{gene_lower}*.*sv')
        matches = glob.glob(pattern)
        if matches:
            filepath = matches[0]
            break
    
    if filepath is None:
        return None
    
    sequences = []
    scores = {}
    coords = []
    
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq = None
            for col in ['Target sequence', 'sgRNA', 'Guide Sequence', 'target_sequence']:
                if col in row:
                    # Strip PAM (last 3 chars) if the sequence is 23nt (20nt guide + 3nt PAM)
                    raw_seq = row[col].upper().replace(' ', '')
                    if len(raw_seq) == 23:
                        seq = raw_seq[:20]
                    else:
                        seq = raw_seq[:20]
                    break
            
            if seq and len(seq) >= 20:
                sequences.append(seq)
                
                for col in ['Efficiency', 'efficiency', 'Score', 'score']:
                    if col in row:
                        try:
                            scores[seq] = float(row[col])
                        except (ValueError, TypeError):
                            pass
                        break
                
                # Parse genomic location (e.g., 'chr17:7675056')
                loc = row.get('Genomic location', '')
                if ':' in loc:
                    parts = loc.replace('chr', '').split(':')
                    try:
                        chrom = parts[0]
                        coord = int(parts[-1].replace(',', ''))
                        coords.append((chrom, coord, coord + 23))
                    except (ValueError, IndexError):
                        coords.append((None, None, None))
                else:
                    coords.append((None, None, None))
    
    if not sequences:
        return None
    
    return {
        'sequences': sequences,
        'scores': scores,
        'coords': coords,
    }


def extract_chopchop_region(chopchop_data: dict, gene: str) -> tuple:
    """
    Extract the genomic region covered by CHOPCHOP results for a gene.
    
    CHOPCHOP's 'Genomic location' column gives coordinates like 'chr17:7675056'.
    Returns (chrom, min_coord, max_coord) or (None, None, None) if no coordinates found.
    """
    if gene not in chopchop_data:
        return None, None, None
    
    cho_seqs = chopchop_data[gene].get('sequences', [])
    cho_coords = chopchop_data[gene].get('coords', [])
    
    if not cho_coords:
        return None, None, None
    
    valid_coords = [(c, s, e) for c, s, e in cho_coords if s is not None]
    if not valid_coords:
        return None, None, None
    
    chrom = valid_coords[0][0]
    starts = [s for _, s, _ in valid_coords]
    ends = [e for _, _, e in valid_coords]
    
    return chrom, min(starts), max(ends)


def filter_caspins_to_region(caspins_gene: dict, chrom: str, 
                              region_start: int, region_end: int,
                              margin: int = 30) -> dict:
    """
    Filter CasPINS results to a specific genomic region (e.g., CHOPCHOP's target).
    
    Args:
        caspins_gene: CasPINS data dict for one gene
        chrom: Chromosome
        region_start: Start coordinate
        region_end: End coordinate
        margin: Extra bp margin around region
    
    Returns:
        New dict with same structure but only gRNAs within the region
    """
    filtered_seqs = []
    filtered_scores = {}
    filtered_doench = {}
    
    sequences = caspins_gene['sequences']
    coords = caspins_gene['coords']
    scores = caspins_gene['scores']
    doench_scores = caspins_gene['doench_scores']
    
    for i, seq in enumerate(sequences):
        if i < len(coords):
            c, s, e = coords[i]
            if c is not None and s is not None:
                # Normalize chromosome for comparison
                c_norm = c.replace('chr', '')
                chrom_norm = chrom.replace('chr', '') if chrom else ''
                if c_norm == chrom_norm and (region_start - margin) <= s <= (region_end + margin):
                    filtered_seqs.append(seq)
                    if seq in scores:
                        filtered_scores[seq] = scores[seq]
                    if seq in doench_scores:
                        filtered_doench[seq] = doench_scores[seq]
    
    return {
        'sequences': filtered_seqs,
        'scores': filtered_scores,
        'doench_scores': filtered_doench,
    }


def calculate_overlap(list1: list, list2: list, top_n: int = 10) -> dict:
    """Calculate overlap between two ranked gRNA lists."""
    set1 = set(list1[:top_n])
    set2 = set(list2[:top_n])
    
    overlap = set1 & set2
    
    return {
        'top_n': top_n,
        'overlap_count': len(overlap),
        'overlap_fraction': len(overlap) / top_n if top_n > 0 else 0,
        'overlapping_sequences': list(overlap),
        'unique_to_first': list(set1 - set2),
        'unique_to_second': list(set2 - set1)
    }


def calculate_rank_correlation(list1: list, list2: list) -> float:
    """
    Calculate Spearman rank correlation for shared gRNAs.
    
    Uses the Pearson correlation of ranks (general formula) rather than the
    shorthand d-squared formula, because the ranks come from lists of
    different sizes and don't necessarily span 1..n.
    """
    shared = set(list1) & set(list2)
    if len(shared) < 3:
        return None  # Not enough shared gRNAs for meaningful correlation
    
    shared = list(shared)
    ranks1 = {seq: i+1 for i, seq in enumerate(list1)}
    ranks2 = {seq: i+1 for i, seq in enumerate(list2)}
    
    r1 = np.array([ranks1[s] for s in shared], dtype=float)
    r2 = np.array([ranks2[s] for s in shared], dtype=float)
    
    # Pearson correlation of the rank vectors (general Spearman's rho)
    mean1, mean2 = np.mean(r1), np.mean(r2)
    std1, std2 = np.std(r1, ddof=0), np.std(r2, ddof=0)
    
    if std1 == 0 or std2 == 0:
        return None
    
    cov = np.mean((r1 - mean1) * (r2 - mean2))
    rho = cov / (std1 * std2)
    
    return float(rho)


def calculate_full_set_overlap(list1: list, list2: list) -> dict:
    """Calculate overlap between two full gRNA lists (not top-N limited)."""
    set1 = set(list1)
    set2 = set(list2)
    overlap = set1 & set2
    smaller_set = min(len(set1), len(set2))
    
    return {
        'set1_size': len(set1),
        'set2_size': len(set2),
        'overlap_count': len(overlap),
        'overlap_pct_of_smaller': (len(overlap) / smaller_set * 100) if smaller_set > 0 else 0,
        'overlap_pct_of_set1': (len(overlap) / len(set1) * 100) if len(set1) > 0 else 0,
        'overlap_pct_of_set2': (len(overlap) / len(set2) * 100) if len(set2) > 0 else 0,
    }


def validate_crispor_data(data_dir: str, genes: list) -> dict:
    """
    Validate CRISPOR data by checking if all files contain distinct gene data.
    Returns dict with 'valid' bool and 'message' string.
    """
    try:
        import xlrd
    except ImportError:
        return {'valid': False, 'message': 'xlrd not installed'}
    
    # Extract the # Position metadata from each CRISPOR file
    positions = {}
    for gene in genes:
        for prefix in ['crispror_', 'crispor_']:
            filepath = os.path.join(CRISPOR_DIR, f'{prefix}{gene.lower()}.xls')
            if os.path.exists(filepath):
                try:
                    wb = xlrd.open_workbook(filepath)
                    sheet = wb.sheet_by_index(0)
                    for r in range(min(10, sheet.nrows)):
                        key = str(sheet.cell_value(r, 0)).strip()
                        val = str(sheet.cell_value(r, 1)).strip()
                        if key == '# Position':
                            positions[gene] = val
                            break
                except Exception:
                    pass
                break
    
    if not positions:
        return {'valid': False, 'message': 'No CRISPOR files found'}
    
    unique_positions = set(positions.values())
    if len(unique_positions) == 1 and len(positions) > 1:
        pos = list(unique_positions)[0]
        return {
            'valid': False,
            'message': (
                f'ALL {len(positions)} CRISPOR files target the SAME genomic region: {pos}\n'
                f'  This means the same input sequence was submitted for every gene.\n'
                f'  The CRISPOR data must be re-downloaded with the correct gene for each file.\n'
                f'  See MANUAL_BENCHMARKING_GUIDE.md for CRISPOR re-query instructions.'
            )
        }
    
    return {'valid': True, 'message': f'{len(positions)} genes with distinct genomic positions'}


def generate_comparison_report(caspins: dict, crispor: dict, chopchop: dict, 
                                output_dir: str, crispor_valid: bool = True):
    """Generate comprehensive comparison report with region-aware analysis."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("gRNA Design Benchmark: CasPINS vs CRISPOR vs CHOPCHOP")
    report_lines.append("=" * 70)
    report_lines.append("")
    
    # Summary table for manuscript
    summary_data = []
    
    all_genes = set(caspins.keys())
    if chopchop:
        all_genes |= set(chopchop.keys())
    if crispor and crispor_valid:
        all_genes |= set(crispor.keys())
    
    # Aggregate stats for headline summary
    all_region_overlaps = []      # Region-aware: % of CHOPCHOP gRNAs found by CasPINS in same region
    all_full_overlaps = []        # Full-set: % of CHOPCHOP gRNAs found anywhere in CasPINS
    all_region_rhos = []
    
    for gene in sorted(all_genes):
        report_lines.append(f"\n{'='*60}")
        report_lines.append(f"Gene: {gene}")
        report_lines.append(f"{'='*60}")
        
        cas_data = caspins.get(gene, {})
        cas_seqs = cas_data.get('sequences', [])
        cas_total = cas_data.get('n_total_found', len(cas_seqs))
        cho_seqs = chopchop.get(gene, {}).get('sequences', []) if chopchop and gene in chopchop else []
        cri_seqs = crispor.get(gene, {}).get('sequences', []) if crispor and crispor_valid and gene in crispor else []
        
        report_lines.append(f"  CasPINS gRNAs returned: {len(cas_seqs)} (out of {cas_total} found in full gene)")
        if cri_seqs:
            report_lines.append(f"  CRISPOR gRNAs: {len(cri_seqs)}")
        report_lines.append(f"  CHOPCHOP gRNAs: {len(cho_seqs)}")
        
        gene_summary = {'gene': gene, 'caspins_total': len(cas_seqs), 
                        'caspins_found': cas_total}
        
        # === CasPINS vs CHOPCHOP (primary comparison) ===
        if cas_seqs and cho_seqs and gene in chopchop:
            # Extract CHOPCHOP target region
            cho_region = extract_chopchop_region(chopchop, gene)
            chrom, region_start, region_end = cho_region
            
            if chrom is not None:
                region_span = region_end - region_start
                report_lines.append(f"\n  CHOPCHOP target region: chr{chrom}:{region_start}-{region_end} "
                                  f"({region_span:,} bp)")
                
                # Filter CasPINS to CHOPCHOP region
                cas_in_region = filter_caspins_to_region(cas_data, chrom, region_start, region_end)
                cas_region_seqs = cas_in_region['sequences']
                
                report_lines.append(f"  CasPINS gRNAs in CHOPCHOP region: {len(cas_region_seqs)}")
                gene_summary['caspins_in_region'] = len(cas_region_seqs)
                gene_summary['chopchop_count'] = len(cho_seqs)
                gene_summary['chopchop_region_bp'] = region_span
                
                # --- Region-aware overlap (primary metric) ---
                region_overlap = calculate_full_set_overlap(cas_region_seqs, cho_seqs)
                pct = region_overlap['overlap_pct_of_set2']
                report_lines.append(f"\n  ** REGION-AWARE OVERLAP (primary metric) **")
                report_lines.append(f"    CasPINS gRNAs in region: {len(cas_region_seqs)}")
                report_lines.append(f"    CHOPCHOP gRNAs: {len(cho_seqs)}")
                report_lines.append(f"    Shared gRNAs: {region_overlap['overlap_count']}")
                report_lines.append(f"    CasPINS recovers {pct:.1f}% of CHOPCHOP gRNAs "
                                  f"({region_overlap['overlap_count']}/{len(cho_seqs)})")
                if len(cas_region_seqs) > 0:
                    pct_cas = region_overlap['overlap_pct_of_set1']
                    report_lines.append(f"    CHOPCHOP covers {pct_cas:.1f}% of CasPINS-in-region gRNAs "
                                      f"({region_overlap['overlap_count']}/{len(cas_region_seqs)})")
                
                all_region_overlaps.append(pct)
                gene_summary['region_overlap_count'] = region_overlap['overlap_count']
                gene_summary['region_overlap_pct_of_chopchop'] = round(pct, 1)
                
                # Region-aware rank correlation
                if cas_region_seqs and cho_seqs:
                    rho = calculate_rank_correlation(cas_region_seqs, cho_seqs)
                    if rho is not None:
                        n_shared = len(set(cas_region_seqs) & set(cho_seqs))
                        report_lines.append(f"    Spearman rank correlation: {rho:.3f} "
                                          f"(on {n_shared} shared gRNAs)")
                        gene_summary['region_spearman_rho'] = round(rho, 3)
                        all_region_rhos.append(rho)
                
                # --- Top-N region-aware overlap ---
                for top_n in [10, 20, 50]:
                    if len(cas_region_seqs) >= top_n and len(cho_seqs) >= top_n:
                        overlap = calculate_overlap(cas_region_seqs, cho_seqs, top_n)
                        report_lines.append(f"    Top-{top_n} overlap (region): "
                                          f"{overlap['overlap_count']}/{top_n} "
                                          f"({overlap['overlap_fraction']*100:.0f}%)")
            
            # --- Full-set overlap (for reference) ---
            full_overlap = calculate_full_set_overlap(cas_seqs, cho_seqs)
            report_lines.append(f"\n  Full-set overlap (CasPINS full gene vs CHOPCHOP):")
            report_lines.append(f"    Shared gRNAs: {full_overlap['overlap_count']} "
                              f"(of {full_overlap['set2_size']} CHOPCHOP = "
                              f"{full_overlap['overlap_pct_of_set2']:.1f}%)")
            all_full_overlaps.append(full_overlap['overlap_pct_of_set2'])
            gene_summary['full_overlap_count'] = full_overlap['overlap_count']
            gene_summary['full_overlap_pct_of_chopchop'] = round(full_overlap['overlap_pct_of_set2'], 1)
        
        # === CasPINS vs CRISPOR (if valid) ===
        if cas_seqs and cri_seqs:
            full_overlap = calculate_full_set_overlap(cas_seqs, cri_seqs)
            report_lines.append(f"\n  CasPINS vs CRISPOR - Full Set Overlap:")
            report_lines.append(f"    Shared gRNAs: {full_overlap['overlap_count']} "
                              f"(of {full_overlap['set2_size']} CRISPOR = "
                              f"{full_overlap['overlap_pct_of_set2']:.1f}%)")
            gene_summary['crispor_full_overlap_count'] = full_overlap['overlap_count']
            gene_summary['crispor_full_overlap_pct'] = round(full_overlap['overlap_pct_of_set2'], 1)
        
        summary_data.append(gene_summary)
    
    # === Headline Summary ===
    report_lines.insert(3, "")
    headline_idx = 4
    headlines = []
    headlines.append("HEADLINE SUMMARY")
    headlines.append("=" * 50)
    headlines.append("")
    headlines.append("Region-Aware Comparison (CasPINS vs CHOPCHOP, same target region):")
    headlines.append("-" * 50)
    if all_region_overlaps:
        avg_region = np.mean(all_region_overlaps)
        headlines.append(f"  Mean overlap: {avg_region:.1f}% of CHOPCHOP gRNAs recovered by CasPINS")
        headlines.append(f"  Per-gene: {', '.join(f'{v:.0f}%' for v in all_region_overlaps)}")
    if all_region_rhos:
        avg_rho = np.mean(all_region_rhos)
        headlines.append(f"  Mean Spearman rho: {avg_rho:.3f}")
    headlines.append(f"  Genes compared: {len(all_region_overlaps)}")
    headlines.append("")
    headlines.append("Full-Gene Comparison (CasPINS full gene vs CHOPCHOP target region):")
    headlines.append("-" * 50)
    if all_full_overlaps:
        avg_full = np.mean(all_full_overlaps)
        headlines.append(f"  Mean overlap: {avg_full:.1f}% of CHOPCHOP gRNAs found in CasPINS full output")
    headlines.append("")
    
    for i, line in enumerate(headlines):
        report_lines.insert(headline_idx + i, line)
    
    # Write report
    report_path = os.path.join(output_dir, 'grna_comparison_report.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    print(f"\nSaved comparison report to: {report_path}")
    
    # Write summary CSV for manuscript table
    summary_path = os.path.join(output_dir, 'grna_comparison_summary.csv')
    fieldnames = [
        'gene', 'caspins_total', 'caspins_found',
        'chopchop_count', 'chopchop_region_bp',
        'caspins_in_region', 'region_overlap_count', 'region_overlap_pct_of_chopchop',
        'region_spearman_rho',
        'full_overlap_count', 'full_overlap_pct_of_chopchop',
    ]
    if crispor_valid:
        fieldnames += ['crispor_full_overlap_count', 'crispor_full_overlap_pct']
    
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in summary_data:
            writer.writerow(row)
    print(f"Saved summary CSV to: {summary_path}")
    
    return summary_data


def main():
    print("Loading CasPINS results...")
    caspins = load_caspins_results(RESULTS_DIR)
    
    # Validate CRISPOR data before loading
    crispor_validation = validate_crispor_data(DATA_DIR, list(caspins.keys()))
    crispor_valid = crispor_validation['valid']
    
    crispor = {}
    if crispor_valid:
        print("Loading CRISPOR results...")
        for gene in caspins.keys():
            result = load_crispor_results(DATA_DIR, gene)
            if result:
                crispor[gene] = result
                print(f"  Loaded CRISPOR results for {gene}: {len(result['sequences'])} gRNAs")
            else:
                print(f"  No CRISPOR results found for {gene}")
    else:
        print(f"\nWARNING: CRISPOR data is INVALID - skipping CRISPOR comparison.")
        print(f"  Reason: {crispor_validation['message']}")
        print(f"  The comparison will proceed with CHOPCHOP only.\n")
    
    print("Loading CHOPCHOP results...")
    chopchop = {}
    for gene in caspins.keys():
        result = load_chopchop_results(DATA_DIR, gene)
        if result:
            chopchop[gene] = result
            print(f"  Loaded CHOPCHOP results for {gene}: {len(result['sequences'])} gRNAs")
        else:
            print(f"  No CHOPCHOP results found for {gene}")
    
    if not crispor and not chopchop:
        print("\nWARNING: No CRISPOR or CHOPCHOP data found!")
        print("Please follow instructions in:")
        print(f"  {os.path.join(RESULTS_DIR, 'CRISPOR_CHOPCHOP_QUERY_INSTRUCTIONS.md')}")
        print("\nGenerating CasPINS-only report...")
    
    output_dir = os.path.join(RESULTS_DIR, 'comparison')
    summary = generate_comparison_report(caspins, crispor, chopchop, output_dir, 
                                          crispor_valid=crispor_valid)
    
    print("\n" + "=" * 60)
    print("COMPARISON COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    main()
