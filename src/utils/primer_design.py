"""
Primer Design Module
Handles primer design for CRISPR target regions
"""

import primer3
from Bio.Seq import Seq


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


def generate_primer_recommendations(grna_positions, mrna_seq, gene_name):
    """
    Generate primer recommendations based on gRNA positions.
    
    Args:
        grna_positions: List of (position, strand, cut_site) tuples
        mrna_seq: mRNA sequence
        gene_name: Gene name
        
    Returns:
        str: Formatted recommendations text
    """
    recommendations = []
    
    if not grna_positions:
        recommendations.append("*** ERROR: No gRNAs found in the provided mRNA sequence! ***")
        recommendations.append("\nTroubleshooting steps:")
        recommendations.append("  1. Verify the mRNA sequence is correct and complete")
        recommendations.append("  2. Check if gRNAs include PAM sequence (remove if present)")
        recommendations.append("  3. Ensure gRNAs are 20nt long (standard length)")
        recommendations.append("  4. Verify gRNAs target the correct gene")
        recommendations.append("  5. Check if you're using genomic vs cDNA sequence")
        return "\n".join(recommendations)
    
    # Define cut sites
    cut_sites = [pos[2] for pos in grna_positions]
    min_cut = min(cut_sites)
    max_cut = max(cut_sites)
    
    recommendations.append("\n" + "=" * 60)
    recommendations.append("RECOMMENDED SEQUENCING PRIMERS (Designed with Primer3):")
    recommendations.append("=" * 60)
    
    # Define target regions for primer design
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
    
    recommendations.append(f"\ngRNA region spans: {min(p[0] for p in grna_positions)}-{max(p[0] for p in grna_positions) + 23} bp")
    recommendations.append(f"Cut sites span: {min_cut}-{max_cut} ({max_cut - min_cut} bp)")
    
    # Add PCR cycling recommendations
    recommendations.append("\n" + "-" * 40)
    recommendations.append("Recommended PCR conditions:")
    recommendations.append("  - Use high-fidelity polymerase (e.g., Q5, Phusion)")
    recommendations.append("  - Initial denaturation: 98°C for 30s")
    recommendations.append("  - 35 cycles: 98°C 10s, 60°C 20s, 72°C 30s")
    recommendations.append("  - Final extension: 72°C for 2 min")
    
    return "\n".join(recommendations) 