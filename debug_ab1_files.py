#!/usr/bin/env python
"""
Debug script to analyze AB1 files and understand why graphs aren't showing correctly
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def analyze_ab1_file(file_path, label):
    """Analyze an AB1 file and print key information."""
    print(f"\n{label}:")
    print("-" * 40)
    
    try:
        record = SeqIO.read(file_path, "abi")
        sequence = str(record.seq)
        
        print(f"  Sequence length: {len(sequence)} bp")
        print(f"  First 50 bp: {sequence[:50]}")
        print(f"  Last 50 bp: {sequence[-50:]}")
        
        # Check sequence quality
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            qualities = record.letter_annotations['phred_quality']
            avg_quality = sum(qualities) / len(qualities)
            print(f"  Average quality score: {avg_quality:.1f}")
            print(f"  Quality of first 20bp: {qualities[:20]}")
        
        return sequence
        
    except Exception as e:
        print(f"  ERROR reading file: {e}")
        return None

def check_grna_in_ab1(ab1_seq, grna, grna_name):
    """Check if gRNA is present in AB1 sequence."""
    print(f"\n  Checking {grna_name}: {grna}")
    
    grna_rc = str(Seq(grna).reverse_complement())
    
    # Try different lengths if gRNA is >20bp
    if len(grna) > 20:
        grna_20 = grna[:20]
        print(f"    Trimmed to 20bp: {grna_20}")
    else:
        grna_20 = grna
    
    # Search for exact matches
    if grna_20 in ab1_seq:
        pos = ab1_seq.find(grna_20)
        print(f"    ✓ Found on forward strand at position {pos}")
        print(f"    Context: ...{ab1_seq[max(0,pos-10):pos+len(grna_20)+10]}...")
        return True
    elif str(Seq(grna_20).reverse_complement()) in ab1_seq:
        pos = ab1_seq.find(str(Seq(grna_20).reverse_complement()))
        print(f"    ✓ Found on reverse strand at position {pos}")
        print(f"    Context: ...{ab1_seq[max(0,pos-10):pos+len(grna_20)+10]}...")
        return True
    else:
        print(f"    ✗ Not found in AB1 sequence")
        
        # Try partial matches
        for i in range(len(ab1_seq) - 15):
            substr = ab1_seq[i:i+20]
            matches = sum(1 for a, b in zip(substr, grna_20) if a == b)
            if matches >= 17:  # 85% match
                print(f"    ~ Partial match at position {i} ({matches}/20 matches)")
                print(f"    AB1: {substr}")
                print(f"    gRNA: {grna_20}")
                break
        
        return False

def main():
    """Main debug function."""
    print("="*60)
    print("AB1 FILE DEBUG ANALYSIS")
    print("="*60)
    
    genes = ['vmat1', 'vmat2', 'ddc']
    
    for gene in genes:
        print(f"\n{'='*60}")
        print(f"ANALYZING {gene.upper()}")
        print(f"{'='*60}")
        
        # File paths
        control_file = f"data/{gene}/control.ab1"
        edited_file = f"data/{gene}/edited.ab1"
        grna_file = f"data/{gene}/grna.txt"
        mrna_file = f"data/{gene}/mrna.txt"
        
        # Check if files exist
        if not all(os.path.exists(f) for f in [control_file, edited_file, grna_file]):
            print(f"Missing files for {gene}")
            continue
        
        # Analyze AB1 files
        control_seq = analyze_ab1_file(control_file, "Control AB1")
        edited_seq = analyze_ab1_file(edited_file, "Edited AB1")
        
        if control_seq and edited_seq:
            # Compare sequences
            print(f"\nSequence Comparison:")
            print(f"  Length difference: {len(edited_seq) - len(control_seq)} bp")
            
            # Check similarity
            min_len = min(len(control_seq), len(edited_seq))
            matches = sum(1 for i in range(min_len) if control_seq[i] == edited_seq[i])
            similarity = matches / min_len * 100
            print(f"  Overall similarity: {similarity:.1f}%")
            
            # Find first difference
            for i in range(min_len):
                if control_seq[i] != edited_seq[i]:
                    print(f"  First difference at position: {i}")
                    print(f"    Control: ...{control_seq[max(0,i-5):i+15]}...")
                    print(f"    Edited:  ...{edited_seq[max(0,i-5):i+15]}...")
                    break
            else:
                if len(control_seq) == len(edited_seq):
                    print("  Sequences are IDENTICAL!")
                else:
                    print(f"  Sequences are identical up to position {min_len}")
        
        # Read gRNAs
        with open(grna_file, 'r') as f:
            grnas = [line.strip().upper() for line in f if line.strip()]
        
        print(f"\nChecking gRNAs in AB1 sequences:")
        for i, grna in enumerate(grnas, 1):
            if control_seq:
                check_grna_in_ab1(control_seq, grna, f"gRNA {i}")
        
        # Read expected positions from mRNA
        if os.path.exists(mrna_file):
            with open(mrna_file, 'r') as f:
                mrna_seq = ''.join(line.strip() for line in f).upper().replace(' ', '')
            
            print(f"\nmRNA reference length: {len(mrna_seq)} bp")
            
            # Find where gRNAs are in mRNA
            for i, grna in enumerate(grnas, 1):
                grna_20 = grna[:20] if len(grna) > 20 else grna
                if grna_20 in mrna_seq:
                    pos = mrna_seq.find(grna_20)
                    print(f"  gRNA {i} in mRNA at position {pos} (cut site: {pos + 17})")
                    
                    # Check if this position exists in AB1
                    if control_seq and pos + 20 > len(control_seq):
                        print(f"    WARNING: This position is beyond AB1 sequence length!")

if __name__ == '__main__':
    main() 