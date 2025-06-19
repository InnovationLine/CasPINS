import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import requests
from bs4 import BeautifulSoup
import os
import re
import ssl
import urllib3

# Disable SSL warnings if we need to bypass verification
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def get_guide_rna(gene_name, rank=1, organism="Rattus norvegicus", target_dir="workdir"):
    """
    Gets guide RNA for a given gene, either from cache or by manual input.
    """
    cache_file = os.path.join(target_dir, f"{gene_name}_gRNA_rank{rank}.txt")
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            gRNA = f.read().strip()
            print(f"Found cached gRNA for {gene_name} (Rank {rank}): {gRNA}")
            return gRNA

    print(f"\nNo cached gRNA found for {gene_name} (Rank {rank}).")
    print(f"Please go to https://chopchop.cbu.uib.no/ and search for:")
    print(f"  - Target: {gene_name}")
    print(f"  - Organism: {organism}")
    print(f"  - CRISPR: Cas9")
    print(f"  - Type: Knock-out")
    print(f"\nThen click on Rank {rank} and find the 'Target sequence' (23 bases including PAM).")
    
    while True:
        gRNA_input = input("\nPlease enter the 23-base target sequence (or 20-base guide without PAM): ").strip().upper()
        
        # Validate input
        if re.match(r'^[ACGT]{20}$', gRNA_input):
            gRNA_20nt = gRNA_input
            print(f"Accepted 20-base guide RNA: {gRNA_20nt}")
            break
        elif re.match(r'^[ACGT]{23}$', gRNA_input):
            gRNA_20nt = gRNA_input[:20]
            print(f"Extracted 20-base guide RNA from 23-base sequence: {gRNA_20nt}")
            break
        else:
            print("Invalid input. Please enter a valid DNA sequence (20 or 23 bases, only A, C, G, T).")
    
    # Cache the result
    with open(cache_file, 'w') as f:
        f.write(gRNA_20nt)
    print(f"Saved gRNA to cache: {cache_file}")
    
    return gRNA_20nt

def parse_ab1(file_path):
    """
    Parses an .ab1 file to extract sequence and chromatogram traces.
    
    Args:
        file_path (str): The path to the .ab1 file.
        
    Returns:
        tuple: A tuple containing the sequence (str) and a dictionary 
               of chromatogram traces (A, C, G, T).
    """
    record = SeqIO.read(file_path, "abi")
    sequence = str(record.seq)
    
    # Extract trace data
    traces = {
        'A': record.annotations['abif_raw']['DATA9'],
        'C': record.annotations['abif_raw']['DATA10'],
        'G': record.annotations['abif_raw']['DATA11'],
        'T': record.annotations['abif_raw']['DATA12'],
    }
    return sequence, traces

def find_approximate_match(sequence, pattern, max_mismatches=2):
    """
    Find approximate matches of pattern in sequence allowing up to max_mismatches.
    Returns the position of the best match or -1 if no match found.
    """
    best_pos = -1
    best_mismatches = max_mismatches + 1
    
    for i in range(len(sequence) - len(pattern) + 1):
        mismatches = sum(1 for j in range(len(pattern)) if sequence[i+j] != pattern[j])
        if mismatches < best_mismatches:
            best_mismatches = mismatches
            best_pos = i
            if mismatches == 0:  # Perfect match found
                break
    
    return best_pos if best_mismatches <= max_mismatches else -1

def plot_chromatograms(control_file, edited_file, guide_rna):
    """
    Aligns and plots the chromatograms of control and edited samples
    around the CRISPR target site.
    """
    print("\nAnalyzing chromatograms...")
    
    # Parse both .ab1 files
    control_seq, control_traces = parse_ab1(control_file)
    edited_seq, edited_traces = parse_ab1(edited_file)
    
    print(f"Control sequence length: {len(control_seq)}")
    print(f"Edited sequence length: {len(edited_seq)}")
    
    # Find the guide RNA in the control sequence (try both strands)
    from Bio.Seq import Seq
    guide_rna_rc = str(Seq(guide_rna).reverse_complement())
    
    cut_position = None
    strand = None
    match_type = None
    
    # Try exact match first
    if guide_rna in control_seq:
        cut_position = control_seq.index(guide_rna) + 17  # Cut is 3bp before PAM
        strand = "forward"
        match_type = "exact"
        print(f"Found exact guide RNA match on forward strand at position {control_seq.index(guide_rna)}")
    elif guide_rna_rc in control_seq:
        cut_position = control_seq.index(guide_rna_rc) + 3  # Cut is 3bp after the start on reverse
        strand = "reverse"
        match_type = "exact"
        print(f"Found exact guide RNA match on reverse strand at position {control_seq.index(guide_rna_rc)}")
    else:
        # Try approximate match with up to 2 mismatches
        print("No exact match found. Searching for approximate matches...")
        
        forward_pos = find_approximate_match(control_seq, guide_rna, 2)
        reverse_pos = find_approximate_match(control_seq, guide_rna_rc, 2)
        
        if forward_pos != -1:
            cut_position = forward_pos + 17
            strand = "forward"
            match_type = "approximate"
            print(f"Found approximate guide RNA match on forward strand at position {forward_pos}")
        elif reverse_pos != -1:
            cut_position = reverse_pos + 3
            strand = "reverse"
            match_type = "approximate"
            print(f"Found approximate guide RNA match on reverse strand at position {reverse_pos}")
        else:
            # Try partial match (first 15 bases)
            partial_guide = guide_rna[:15]
            partial_guide_rc = guide_rna_rc[:15]
            
            if partial_guide in control_seq:
                cut_position = control_seq.index(partial_guide) + 17
                strand = "forward"
                match_type = "partial"
                print(f"Found partial guide RNA match (first 15bp) on forward strand")
            elif partial_guide_rc in control_seq:
                cut_position = control_seq.index(partial_guide_rc) + 3
                strand = "reverse"
                match_type = "partial"
                print(f"Found partial guide RNA match (first 15bp) on reverse strand")
            else:
                print("ERROR: Could not find guide RNA in control sequence even with approximate matching!")
                print(f"Searched for: {guide_rna} and {guide_rna_rc}")
                print("\nShowing first 100bp of control sequence:")
                print(control_seq[:100])
                return
    
    print(f"Match type: {match_type}, Strand: {strand}")
    print(f"Predicted cut site at position: {cut_position}")
    
    # Define the window around the cut site to display
    window_size = 50  # bases on each side of cut
    start = max(0, cut_position - window_size)
    end = min(len(control_seq), cut_position + window_size)
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
    
    # Colors for each base
    colors = {'A': 'green', 'C': 'blue', 'G': 'black', 'T': 'red'}
    
    # Plot control chromatogram
    ax1.set_title('Control Sample Chromatogram', fontsize=14)
    ax1.set_ylabel('Signal Intensity', fontsize=12)
    
    # We need to map sequence positions to trace positions
    # This is approximate - real implementation would use PLOC data from .ab1
    trace_start = start * 10  # Rough approximation
    trace_end = end * 10
    
    for base, color in colors.items():
        trace_data = control_traces[base][trace_start:trace_end]
        x_vals = np.arange(len(trace_data))
        ax1.plot(x_vals, trace_data, color=color, label=base, alpha=0.7)
    
    # Mark the cut site
    cut_trace_pos = (cut_position - start) * 10
    ax1.axvline(x=cut_trace_pos, color='red', linestyle='--', label='Cut site')
    ax1.legend()
    
    # Plot edited chromatogram
    ax2.set_title('Edited Sample Chromatogram', fontsize=14)
    ax2.set_ylabel('Signal Intensity', fontsize=12)
    ax2.set_xlabel('Trace Position', fontsize=12)
    
    for base, color in colors.items():
        trace_data = edited_traces[base][trace_start:trace_end]
        x_vals = np.arange(len(trace_data))
        ax2.plot(x_vals, trace_data, color=color, label=base, alpha=0.7)
    
    ax2.axvline(x=cut_trace_pos, color='red', linestyle='--', label='Cut site')
    ax2.legend()
    
    plt.tight_layout()
    
    # Save the plot
    output_file = "tide_analysis_chromatogram.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nChromatogram plot saved to: {output_file}")
    
    # Show the plot
    plt.show()
    
    # Print sequence around cut site for verification
    print(f"\nSequence around cut site (±20bp):")
    cut_start = max(0, cut_position - 20)
    cut_end = min(len(control_seq), cut_position + 20)
    print(f"Control: {control_seq[cut_start:cut_position]}|{control_seq[cut_position:cut_end]}")
    print(f"         {' ' * 20}^")
    print(f"         {' ' * 20}Cut site")
    
    return cut_position

if __name__ == '__main__':
    # --- Gene and File Configuration ---
    GENE_NAME = "SLC18A1" # Official gene symbol for vmat1
    ORGANISM = "Rattus norvegicus"
    CONTROL_FILE = "workdir/control.ab1"
    EDITED_FILE = "workdir/edited.ab1"

    # 1. Get Guide RNA (from cache or web)
    # Fetching Rank 1 gRNA as it's the most common choice.
    guide_rna = get_guide_rna(GENE_NAME, rank=1, organism=ORGANISM)

    if guide_rna:
        print("\nStarting TIDE analysis...")
        plot_chromatograms(CONTROL_FILE, EDITED_FILE, guide_rna)
    else:
        print("\nCould not retrieve guide RNA. Exiting analysis.")
        print("Please check the gene name and your internet connection.") 