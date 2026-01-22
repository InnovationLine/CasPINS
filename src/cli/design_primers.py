#!/usr/bin/env python
"""
Advanced Primer Design Tool for CRISPR Validation
=================================================

A comprehensive command-line tool to design PCR primers with:
- Direct integration with NCBI, Ensembl, and UniProt databases
- Advanced primer analysis (hairpins, dimers, specificity)
- Detailed quality metrics and recommendations
- Two-tier primer design strategy:
  * PCR I: For genomic DNA amplification (larger products)
  * PCR II: For gene-specific sequencing (optimal for indel analysis)

Usage:
    python -m cli.design_primers GENE_NAME [options]
    
Examples:
    python -m cli.design_primers TP53
    python -m cli.design_primers DDC --species human --check-databases
    python -m cli.design_primers BRCA1 --species mouse --output-format detailed
"""

import sys
import os
import argparse
import json
from datetime import datetime

# Handle imports whether run as module or directly
if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.primer_design import (
    generate_primer_recommendations, 
    generate_enhanced_primer_report,
    PrimerDesigner
)
from utils.sequence_analysis import find_grna_in_sequence

# Import config for data directory
try:
    from config.settings import get_data_directory
except ImportError:
    def get_data_directory():
        """Fallback: check environment or return None."""
        return os.environ.get('CRISPR_DATA_DIR', None)


def get_configured_data_dir():
    """Get data directory from config, environment, or prompt user."""
    data_dir = get_data_directory()
    if data_dir and os.path.isdir(data_dir):
        return data_dir
    
    # Check environment variable
    env_dir = os.environ.get('CRISPR_DATA_DIR', '')
    if env_dir and os.path.isdir(env_dir):
        return env_dir
    
    return None


def validate_inputs(gene_folder, gene_name, data_dir):
    """Validate that required input files exist."""
    grna_file = os.path.join(gene_folder, 'grna.txt')
    mrna_file = os.path.join(gene_folder, 'mrna.txt')
    
    if not os.path.exists(gene_folder):
        print(f"Error: Gene folder '{gene_folder}' not found")
        print(f"\nTo use this tool:")
        print(f"1. Set data directory via:")
        print(f"   - GUI: Settings panel in sidebar")
        print(f"   - CLI: --data-dir /path/to/data")
        print(f"   - Environment: export CRISPR_DATA_DIR=/path/to/data")
        print(f"2. Create folder: <data_dir>/{gene_name.lower()}/")
        print(f"3. Add gRNA sequences to: <data_dir>/{gene_name.lower()}/grna.txt")
        print(f"4. Add mRNA sequence to: <data_dir>/{gene_name.lower()}/mrna.txt")
        
        if data_dir and os.path.isdir(data_dir):
            available_genes = [d.upper() for d in os.listdir(data_dir) 
                             if os.path.isdir(os.path.join(data_dir, d)) and not d.startswith('.')]
            if available_genes:
                print(f"\nAvailable genes in {data_dir}: {', '.join(available_genes)}")
        return None, None
    
    if not os.path.exists(grna_file):
        print(f"Error: gRNA file not found at {grna_file}")
        print("Add your gRNA sequences (one per line) to this file")
        return None, None
    
    if not os.path.exists(mrna_file):
        print(f"Error: mRNA file not found at {mrna_file}")
        print("Add the mRNA sequence to this file")
        return None, None
    
    return grna_file, mrna_file


def read_sequences(grna_file, mrna_file):
    """Read gRNA and mRNA sequences from files."""
    # Read gRNA sequences
    with open(grna_file, 'r') as f:
        grna_sequences = [line.strip().upper() for line in f 
                         if line.strip() and not line.startswith('#')]
    
    # Read mRNA sequence
    with open(mrna_file, 'r') as f:
        mrna_seq = ''.join(line.strip() for line in f 
                          if not line.startswith('>')).upper()
    
    return grna_sequences, mrna_seq


def main():
    parser = argparse.ArgumentParser(
        description='Design PCR primers for CRISPR validation experiments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m cli.design_primers TP53                    # Design primers for TP53
  python -m cli.design_primers DDC --check-databases   # With database checking
  python -m cli.design_primers BRCA1 --output-format detailed

Set data directory via:
  - CLI: --data-dir /path/to/data
  - Environment: export CRISPR_DATA_DIR=/path/to/data
  - GUI: Settings panel in sidebar
        """
    )
    
    parser.add_argument('gene', help='Gene name (e.g., TP53, DDC, BRCA1)')
    parser.add_argument('--data-dir', 
                       help='Path to data directory containing gene folders')
    parser.add_argument('--species', default='human', 
                       choices=['human', 'mouse', 'rat'],
                       help='Species (default: human)')
    parser.add_argument('--check-databases', action='store_true',
                       help='Check online databases for gene information')
    parser.add_argument('--output-format', 
                       choices=['standard', 'detailed', 'both'],
                       default='standard',
                       help='Output format (default: standard)')
    parser.add_argument('--primer-set', 
                       choices=['both', 'pcr1', 'pcr2'], 
                       default='both',
                       help='Which primer set to design (default: both)')
    parser.add_argument('--no-variants', action='store_true',
                       help='Skip variant checking')
    parser.add_argument('--tm-range', nargs=2, type=float,
                       metavar=('MIN', 'MAX'),
                       help='Tm range for primers (default: 57-63)')
    
    args = parser.parse_args()
    
    # Get data directory
    data_dir = args.data_dir or get_configured_data_dir()
    
    if not data_dir:
        print("Error: Data directory not configured!")
        print("\nPlease set data directory using one of:")
        print("  1. --data-dir /path/to/data")
        print("  2. Environment variable CRISPR_DATA_DIR")
        print("  3. GUI Settings panel")
        return 1
    
    # Construct gene folder path
    gene_folder = os.path.join(data_dir, args.gene.lower())
    
    # Validate inputs
    grna_file, mrna_file = validate_inputs(gene_folder, args.gene, data_dir)
    if not grna_file or not mrna_file:
        return 1
    
    # Read sequences
    grna_sequences, mrna_seq = read_sequences(grna_file, mrna_file)
    
    print(f"\n{'='*60}")
    print(f"Designing primers for {args.gene} ({args.species})")
    print(f"{'='*60}")
    print(f"\nFound {len(grna_sequences)} gRNA sequence(s)")
    print(f"mRNA sequence length: {len(mrna_seq)} bp")
    
    # Find gRNA positions
    grna_positions = []
    for i, grna in enumerate(grna_sequences):
        pos, strand, cut_site = find_grna_in_sequence(mrna_seq, grna)
        if pos is not None:
            grna_positions.append({
                'sequence': grna,
                'position': pos,
                'strand': strand,
                'cut_site': cut_site
            })
            print(f"  gRNA {i+1}: Found at position {pos} ({strand} strand), cut site: {cut_site}")
        else:
            print(f"  gRNA {i+1}: Not found in mRNA sequence")
    
    if not grna_positions:
        print("\nError: No gRNAs found in mRNA sequence")
        print("Please verify that:")
        print("  1. gRNA sequences are correct (20bp, no PAM)")
        print("  2. mRNA sequence is correct")
        return 1
    
    # Generate primer recommendations
    print("\n" + "-"*60)
    print("Designing primers...")
    
    recommendations = generate_enhanced_primer_report(
        grna_positions, mrna_seq, args.gene, args.species
    )
    
    # Save results
    output_file = os.path.join(gene_folder, 'primer_recommendations.txt')
    with open(output_file, 'w') as f:
        f.write(recommendations)
    
    print(f"\n✓ Primer recommendations saved to: {output_file}")
    
    # Also save detailed JSON if requested
    if args.output_format in ['detailed', 'both']:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        json_file = os.path.join(gene_folder, f'primer_design_detailed_{timestamp}.json')
        
        detailed_data = {
            'gene': args.gene,
            'species': args.species,
            'timestamp': timestamp,
            'grna_sequences': grna_sequences,
            'grna_positions': grna_positions,
            'mrna_length': len(mrna_seq)
        }
        
        with open(json_file, 'w') as f:
            json.dump(detailed_data, f, indent=2)
        
        print(f"✓ Detailed JSON saved to: {json_file}")
    
    # Print summary
    print("\n" + "="*60)
    print("Next steps:")
    print("  1. Review primer recommendations in the output file")
    print("  2. Validate primers using NCBI Primer-BLAST")
    print("  3. Order primers from your preferred vendor")
    print("  4. Run PCR and sequencing")
    print("  5. Use run.py analysis to analyze AB1 files")
    print("="*60)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
