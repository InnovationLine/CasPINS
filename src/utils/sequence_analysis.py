"""
Sequence Analysis Module
Handles sequence alignment, similarity calculation, and divergence detection
"""

from Bio import pairwise2
from Bio.Seq import Seq
import warnings
warnings.filterwarnings('ignore')


def align_sequences(seq1, seq2):
    """
    Perform global alignment of two sequences to find exact indel locations.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        tuple: (position, aligned_seq1, aligned_seq2)
    """
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


def calculate_similarity(seq1, seq2):
    """
    Calculate percentage similarity between two sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        float: Similarity percentage (0-100)
    """
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0
    
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    return (matches / min_len) * 100


def find_divergence_point(seq1, seq2, window_size=10):
    """
    Find the point where two sequences start to diverge significantly.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        window_size: Size of window to check for mismatches
        
    Returns:
        int: Position of divergence or -1 if no clear divergence
    """
    min_len = min(len(seq1), len(seq2))
    
    # First check if sequences are identical
    if seq1[:min_len] == seq2[:min_len]:
        return -1  # No divergence
    
    # Check overall similarity
    similarity = calculate_similarity(seq1, seq2)
    
    # If sequences are too different, return -1
    if similarity < 50:
        return -1
    
    # Look for divergence point
    for i in range(min_len - window_size):
        window1 = seq1[i:i+window_size]
        window2 = seq2[i:i+window_size]
        mismatches = sum(1 for a, b in zip(window1, window2) if a != b)
        
        if mismatches >= 3:
            # Check if this is real divergence or just noise
            if i + window_size + 10 < min_len:
                next_window1 = seq1[i+5:i+15]
                next_window2 = seq2[i+5:i+15]
                next_mismatches = sum(1 for a, b in zip(next_window1, next_window2) if a != b)
                
                if next_mismatches >= 3:
                    return i
            else:
                return i
    
    return -1


def find_grna_in_sequence(sequence, grna, allow_mismatch=False):
    """
    Find gRNA in a sequence, checking both strands.
    
    Args:
        sequence: Target sequence to search in
        grna: Guide RNA sequence
        allow_mismatch: Whether to allow partial matches
        
    Returns:
        tuple: (position, strand, cut_site) or (None, None, None) if not found
    """
    # Handle non-standard gRNA lengths
    if len(grna) > 20:
        # Try removing potential PAM sequences
        if grna.endswith('GG') or grna.endswith('CC'):
            grna_20bp = grna[:20]
        else:
            grna_20bp = grna[:20]
    else:
        grna_20bp = grna
    
    # Check forward strand
    if grna_20bp in sequence:
        pos = sequence.find(grna_20bp)
        cut_site = pos + 17  # 3bp before PAM
        return pos, 'forward', cut_site
    
    # Check reverse strand
    grna_rc = str(Seq(grna_20bp).reverse_complement())
    if grna_rc in sequence:
        pos = sequence.find(grna_rc)
        cut_site = pos + 3
        return pos, 'reverse', cut_site
    
    # Try partial match if allowed
    if allow_mismatch and len(grna_20bp) >= 19:
        partial_grna = grna_20bp[1:]
        if partial_grna in sequence:
            pos = sequence.find(partial_grna)
            return pos, 'partial', pos + 16
    
    return None, None, None 