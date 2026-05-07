"""
AB1 File Parser Module
Handles parsing of.ab1 chromatogram files
"""

from Bio import SeqIO


def parse_ab1(file_path):
    """
    Parse an.ab1 file to extract sequence and chromatogram traces.
    
    Args:
        file_path: Path to the.ab1 file
        
    Returns:
        tuple: (sequence, traces_dict)
    """
    record = SeqIO.read(file_path, "abi")
    sequence = str(record.seq)
    
    # Check if sequence is all N's and try alternative sources
    if set(sequence.upper()) == {'N'} and 'abif_raw' in record.annotations:
        raw_data = record.annotations['abif_raw']
        
        # Try PBAS1 (called bases) which often contains the actual sequence
        if 'PBAS1' in raw_data:
            pbas1_seq = raw_data['PBAS1'].decode('utf-8') if isinstance(raw_data['PBAS1'], bytes) else raw_data['PBAS1']
            if pbas1_seq and set(pbas1_seq.upper()) != {'N'}:
                sequence = pbas1_seq
        
        # If still no good sequence, try PBAS2
        elif 'PBAS2' in raw_data and set(sequence.upper()) == {'N'}:
            pbas2_seq = raw_data['PBAS2'].decode('utf-8') if isinstance(raw_data['PBAS2'], bytes) else raw_data['PBAS2']
            if pbas2_seq and set(pbas2_seq.upper()) != {'N'}:
                sequence = pbas2_seq
    
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
        file_path: Path to the.ab1 file
        
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