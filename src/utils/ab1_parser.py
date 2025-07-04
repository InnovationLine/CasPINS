"""
AB1 File Parser Module
Handles parsing of .ab1 chromatogram files
"""

from Bio import SeqIO


def parse_ab1(file_path):
    """
    Parse an .ab1 file to extract sequence and chromatogram traces.
    
    Args:
        file_path: Path to the .ab1 file
        
    Returns:
        tuple: (sequence, traces_dict)
    """
    record = SeqIO.read(file_path, "abi")
    sequence = str(record.seq)
    traces = {
        'A': record.annotations['abif_raw']['DATA9'],
        'C': record.annotations['abif_raw']['DATA10'],
        'G': record.annotations['abif_raw']['DATA11'],
        'T': record.annotations['abif_raw']['DATA12'],
    }
    return sequence, traces


def get_quality_scores(file_path):
    """
    Extract quality scores from AB1 file.
    
    Args:
        file_path: Path to the .ab1 file
        
    Returns:
        list: Phred quality scores or None if not available
    """
    try:
        record = SeqIO.read(file_path, "abi")
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            return record.letter_annotations['phred_quality']
    except Exception:
        pass
    return None 