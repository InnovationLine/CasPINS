"""
AB1 Analyzer Module
Performs detailed analysis of AB1 files and generates reports
"""

from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from.ab1_parser import get_quality_scores
from.sequence_analysis import find_grna_in_sequence


def analyze_ab1_details(control_file, edited_file, grna_sequences, mrna_seq, gene_name):
    """
    Perform detailed AB1 analysis and generate report.
    
    Args:
        control_file: Path to control AB1
        edited_file: Path to edited AB1
        grna_sequences: List of gRNA sequences
        mrna_seq: mRNA reference sequence
        gene_name: Gene name
        
    Returns:
        str: Formatted analysis report
    """
    analysis_lines = []
    analysis_lines.append(f"DETAILED AB1 ANALYSIS FOR {gene_name.upper()}")
    analysis_lines.append("=" * 60)
    analysis_lines.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    analysis_lines.append("")
    
    # Analyze control AB1
    analysis_lines.append("CONTROL AB1 FILE:")
    analysis_lines.append("-" * 40)
    control_seq = _analyze_single_ab1(control_file, analysis_lines)
    
    analysis_lines.append("")
    
    # Analyze edited AB1
    analysis_lines.append("EDITED AB1 FILE:")
    analysis_lines.append("-" * 40)
    edited_seq = _analyze_single_ab1(edited_file, analysis_lines)
    
    analysis_lines.append("")
    
    # Compare sequences if both loaded successfully
    if control_seq and edited_seq:
        _compare_sequences(control_seq, edited_seq, analysis_lines)
    
    analysis_lines.append("")
    
    # Check gRNAs in AB1 sequences
    _analyze_grnas_in_sequences(grna_sequences, control_seq, edited_seq, analysis_lines)
    
    # mRNA reference analysis
    _analyze_mrna_reference(grna_sequences, mrna_seq, control_seq, analysis_lines)
    
    return '\n'.join(analysis_lines)


def _analyze_single_ab1(file_path, analysis_lines):
    """Analyze a single AB1 file and append results to analysis_lines."""
    try:
        record = SeqIO.read(file_path, "abi")
        sequence = str(record.seq)
        
        analysis_lines.append(f"  Sequence length: {len(sequence)} bp")
        analysis_lines.append(f"  First 50 bp: {sequence[:50]}")
        analysis_lines.append(f"  Last 50 bp: {sequence[-50:]}")
        
        # Check sequence quality
        qualities = get_quality_scores(file_path)
        if qualities:
            avg_quality = sum(qualities) / len(qualities)
            analysis_lines.append(f"  Average quality score: {avg_quality:.1f}")
            analysis_lines.append(f"  Quality of first 20bp: {qualities[:20]}")
            analysis_lines.append(f"  Quality of last 20bp: {qualities[-20:]}")
            
            # Quality distribution
            low_quality = sum(1 for q in qualities if q < 20)
            med_quality = sum(1 for q in qualities if 20 <= q < 30)
            high_quality = sum(1 for q in qualities if q >= 30)
            
            analysis_lines.append(f"  Quality distribution:")
            analysis_lines.append(f"    Low (<20): {low_quality} bases ({low_quality/len(qualities)*100:.1f}%)")
            analysis_lines.append(f"    Medium (20-30): {med_quality} bases ({med_quality/len(qualities)*100:.1f}%)")
            analysis_lines.append(f"    High (>30): {high_quality} bases ({high_quality/len(qualities)*100:.1f}%)")
        
        return sequence
        
    except Exception as e:
        analysis_lines.append(f"  ERROR reading file: {e}")
        return None


def _compare_sequences(control_seq, edited_seq, analysis_lines):
    """Compare control and edited sequences."""
    analysis_lines.append("SEQUENCE COMPARISON:")
    analysis_lines.append("-" * 40)
    analysis_lines.append(f"  Length difference: {len(edited_seq) - len(control_seq)} bp")
    
    # Check similarity
    min_len = min(len(control_seq), len(edited_seq))
    if min_len > 0:
        matches = sum(1 for i in range(min_len) if control_seq[i] == edited_seq[i])
        similarity = matches / min_len * 100
        analysis_lines.append(f"  Overall similarity: {similarity:.1f}%")
        
        # Find first difference
        for i in range(min_len):
            if control_seq[i] != edited_seq[i]:
                analysis_lines.append(f"  First difference at position: {i}")
                context_start = max(0, i-10)
                context_end = min(min_len, i+20)
                analysis_lines.append(f"    Control:...{control_seq[context_start:context_end]}...")
                analysis_lines.append(f"    Edited: ...{edited_seq[context_start:context_end]}...")
                break
        else:
            if len(control_seq) == len(edited_seq):
                analysis_lines.append("  Sequences are IDENTICAL!")
            else:
                analysis_lines.append(f"  Sequences are identical up to position {min_len}")


def _analyze_grnas_in_sequences(grna_sequences, control_seq, edited_seq, analysis_lines):
    """Check if gRNAs are present in AB1 sequences."""
    analysis_lines.append("gRNA ANALYSIS IN AB1 SEQUENCES:")
    analysis_lines.append("-" * 40)
    
    for i, grna in enumerate(grna_sequences, 1):
        original_grna = grna
        analysis_lines.append(f"\nGuide RNA {i}: {original_grna}")
        
        # Handle non-standard lengths
        if len(grna) > 20:
            grna_20 = grna[:20]
            analysis_lines.append(f"  Trimmed to 20bp: {grna_20}")
        else:
            grna_20 = grna
        
        # Check in control sequence
        if control_seq:
            pos, strand, cut_site = find_grna_in_sequence(control_seq, grna_20, allow_mismatch=True)
            
            if pos is not None:
                analysis_lines.append(f"  ✓ Found in CONTROL on {strand} strand at position {pos}")
                context = control_seq[max(0,pos-10):min(len(control_seq),pos+len(grna_20)+10)]
                analysis_lines.append(f"    Context:...{context}...")
            else:
                analysis_lines.append(f"  ✗ Not found in CONTROL sequence")
                
                # Check for partial matches
                best_match = 0
                best_pos = -1
                for j in range(len(control_seq) - len(grna_20) + 1):
                    matches = sum(1 for k in range(len(grna_20)) if control_seq[j+k] == grna_20[k])
                    if matches > best_match and matches >= 15:  # At least 75% match
                        best_match = matches
                        best_pos = j
                
                if best_pos >= 0:
                    analysis_lines.append(f"    ~ Partial match at position {best_pos} ({best_match}/{len(grna_20)} matches)")
                    analysis_lines.append(f"      AB1: {control_seq[best_pos:best_pos+len(grna_20)]}")
                    analysis_lines.append(f"      gRNA: {grna_20}")


def _analyze_mrna_reference(grna_sequences, mrna_seq, control_seq, analysis_lines):
    """Analyze gRNA positions in mRNA reference."""
    analysis_lines.append("")
    analysis_lines.append("mRNA REFERENCE ANALYSIS:")
    analysis_lines.append("-" * 40)
    analysis_lines.append(f"  mRNA length: {len(mrna_seq)} bp")
    
    # Check where gRNAs are expected
    for i, grna in enumerate(grna_sequences, 1):
        if len(grna) > 20:
            grna_20 = grna[:20]
        else:
            grna_20 = grna
        
        pos, strand, cut_site = find_grna_in_sequence(mrna_seq, grna_20)
        
        if pos is not None:
            analysis_lines.append(f"  gRNA {i} expected at position {pos} (cut site: {cut_site})")
            
            # Check if this position exists in AB1
            if control_seq and cut_site > len(control_seq):
                analysis_lines.append(f"    WARNING: Cut site is beyond AB1 sequence length ({len(control_seq)} bp)!")
                analysis_lines.append(f"    Need to sequence at least {cut_site + 50} bp to cover this gRNA") 