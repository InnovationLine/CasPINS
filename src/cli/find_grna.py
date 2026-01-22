#!/usr/bin/env python
"""
Find gRNA Tool
==============

A simple command-line tool to find optimal gRNAs for your target genes.
This tool helps you find gRNAs that you can then use with run_analysis.py

Usage:
    python -m cli.find_grna GENE_NAME [options]
    
Examples:
    python -m cli.find_grna TP53
    python -m cli.find_grna DDC --top 5 --save-to data/ddc/grna.txt
    python -m cli.find_grna BRCA1 --cas-type SpCas9 --gc-min 45 --gc-max 55
    python -m cli.find_grna TP53 --species rat --cas-type Cas9-VQR
"""

import argparse
import sys
import os

# Handle imports whether run as module or directly
if __name__ == '__main__':
    # Add parent src directory to path when run directly
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from grna_design import GRNADesigner


def main():
    parser = argparse.ArgumentParser(
        description='Find optimal gRNAs for your target gene',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m cli.find_grna TP53                    # Find gRNAs for human TP53
  python -m cli.find_grna TP53 --species rat      # Find gRNAs for rat TP53
  python -m cli.find_grna DDC --top 5             # Get top 5 gRNAs for DDC
  python -m cli.find_grna BRCA1 --save-to data/brca1/grna.txt  # Save directly
  
Supported species: human, mouse, rat, zebrafish
Supported Cas types: SpCas9, SaCas9, Cas12a, SpCas9-NG, Cas9-VQR, Cas9-EQR, Cas9-VRER, xCas9

After finding gRNAs, copy them to data/<gene>/grna.txt and run:
  python -m cli.run_analysis
        """
    )
    
    # Required arguments
    parser.add_argument('gene', help='Gene name (e.g., TP53, DDC, BRCA1) or Ensembl ID')
    
    # Species and genome parameters
    parser.add_argument('--species', default='human',
                       choices=['human', 'mouse', 'rat', 'zebrafish', 'drosophila', 'c_elegans'],
                       help='Target species (default: human)')
    parser.add_argument('--assembly', 
                       help='Genome assembly version (e.g., GRCh38, GRCm39, Rnor_6.0)')
    parser.add_argument('--transcript', 
                       help='Specific transcript ID to target')
    
    # CRISPR parameters
    parser.add_argument('--cas-type', default='SpCas9',
                       choices=['SpCas9', 'SaCas9', 'Cas12a', 'SpCas9-NG', 'Cas9-VQR', 
                               'Cas9-EQR', 'Cas9-VRER', 'xCas9'],
                       help='CRISPR-Cas system type (default: SpCas9)')
    parser.add_argument('--grna-length', type=int, default=20,
                       help='Length of guide RNA (default: 20)')
    
    # Target region parameters
    parser.add_argument('--target-region', default='all',
                       choices=['all', 'exons', 'cds', '5utr', '3utr', 'introns'],
                       help='Region to target (default: all)')
    parser.add_argument('--target-exons', 
                       help='Specific exons to target (e.g., "1,2,3" or "1-3")')
    parser.add_argument('--avoid-exons', 
                       help='Exons to avoid (e.g., "4,5")')
    
    # Filtering parameters
    parser.add_argument('--top', type=int, default=10, 
                       help='Number of top gRNAs to return (default: 10)')
    parser.add_argument('--gc-min', type=float, default=40,
                       help='Minimum GC content %% (default: 40)')
    parser.add_argument('--gc-max', type=float, default=60,
                       help='Maximum GC content %% (default: 60)')
    parser.add_argument('--homopolymer-max', type=int, default=4,
                       help='Maximum homopolymer length (default: 4)')
    parser.add_argument('--no-poly-t-seed', action='store_true',
                       help='Remove gRNAs with poly-T in seed region')
    
    # Off-target parameters
    parser.add_argument('--check-off-targets', action='store_true',
                       help='Check for off-target sites (requires genome download)')
    parser.add_argument('--max-off-targets', type=int, default=10,
                       help='Maximum allowed off-targets with <= 3 mismatches (default: 10)')
    
    # Output parameters
    parser.add_argument('--save-to', help='Save gRNAs directly to file (e.g., data/gene/grna.txt)')
    parser.add_argument('--export-csv', help='Export detailed results to CSV file')
    parser.add_argument('--show-all-scores', action='store_true',
                       help='Show all scoring metrics')
    parser.add_argument('--output-format', default='sequence', 
                       choices=['sequence', 'fasta', 'detailed'],
                       help='Output format for saved gRNAs (default: sequence)')
    
    # Advanced options
    parser.add_argument('--use-local-only', action='store_true',
                       help='Only use local mRNA files, do not fetch from databases')
    parser.add_argument('--expand-search', type=int, default=0,
                       help='Expand search region upstream/downstream by N bases')
    
    args = parser.parse_args()
    
    # Header
    print("\n" + "="*60)
    print(f"Finding gRNAs for {args.gene} ({args.species})")
    print("="*60)
    
    # Initialize designer with species and assembly
    print(f"\nInitializing {args.cas_type} gRNA designer for {args.species}...")
    designer = GRNADesigner(
        species=args.species,
        assembly=args.assembly,
        cas_type=args.cas_type
    )
    
    # Prepare filters
    filters = {
        'gc_min': args.gc_min,
        'gc_max': args.gc_max,
        'homopolymer_max': args.homopolymer_max,
        'remove_poly_t_seed': args.no_poly_t_seed
    }
    
    # Design gRNAs
    print(f"\nSearching for {args.cas_type} gRNAs in {args.target_region} region...")
    
    # Prepare design parameters
    design_params = {
        'n_results': args.top,
        'filters': filters,
        'include_off_targets': args.check_off_targets
    }
    
    # If using local only, set target type
    if args.use_local_only:
        design_params['target_type'] = 'gene'  # Will check local files first
    
    results = designer.design_grnas(
        target=args.gene,
        **design_params
    )
    
    # Check for errors
    if 'error' in results:
        print(f"\nError: {results['error']}")
        print("\nTroubleshooting:")
        print("1. Make sure the gene name is correct for your species")
        print(f"2. For {args.species}, check if mRNA sequence exists in data/<gene>/mrna.txt")
        print("3. Try using an Ensembl ID (e.g., ENSG00000141510 for TP53)")
        print("4. Use --use-local-only flag if you have local sequence files")
        return 1
    
    # Display results
    grnas = results.get('grnas', [])
    if not grnas:
        print("\nNo suitable gRNAs found with current filters.")
        print("Try adjusting:")
        print("  - GC content range (--gc-min, --gc-max)")
        print("  - Cas type (--cas-type)")
        print("  - Target region (--target-region)")
        return 1
    
    # Display metadata
    metadata = results.get('metadata', {})
    if metadata:
        print(f"\nGene Information:")
        print(f"  Name: {metadata.get('gene_name', 'N/A')}")
        if metadata.get('gene_symbol'):
            print(f"  Symbol: {metadata.get('gene_symbol')}")
        if metadata.get('ensembl_id'):
            print(f"  Ensembl ID: {metadata.get('ensembl_id')}")
        print(f"  Chromosome: {metadata.get('chromosome', 'N/A')}")
        
        # Display strand properly
        strand = metadata.get('strand')
        if strand == 1:
            strand_display = '+ (sense)'
        elif strand == -1:
            strand_display = '- (antisense)'
        else:
            strand_display = 'N/A'
        print(f"  Strand: {strand_display}")
        
        if metadata.get('start') and metadata.get('end'):
            print(f"  Location: {metadata.get('chromosome', '?')}:{metadata.get('start')}-{metadata.get('end')}")
        if metadata.get('assembly'):
            print(f"  Assembly: {metadata.get('assembly')}")
        if metadata.get('source') == 'local_file':
            print(f"  Source: Local file with genomic metadata")
    
    print(f"\nFound {len(grnas)} high-quality gRNAs:")
    print("-" * 80)
    
    # Show gRNAs
    grna_sequences = []
    for i, grna in enumerate(grnas):
        print(f"\n{i+1}. gRNA Sequence: {grna['sequence']}")
        print(f"   PAM: {grna.get('pam', 'N/A')}")
        print(f"   Position: {grna.get('position', 'N/A')} ({grna.get('strand', 'N/A')} strand)")
        print(f"   GC Content: {grna.get('gc_content', 'N/A')}%")
        print(f"   Score: {grna['scores']['composite']:.3f}")
        
        if args.show_all_scores:
            print(f"   Detailed Scores:")
            print(f"     - Doench 2016: {grna['scores']['doench_2016']:.3f}")
            print(f"     - Moreno-Mateos: {grna['scores']['moreno_mateos']:.3f}")
            print(f"     - Xu: {grna['scores']['xu']:.3f}")
        
        if args.check_off_targets and 'off_targets' in grna:
            print(f"   Off-targets: {len(grna['off_targets'])} sites with <= 3 mismatches")
        
        grna_sequences.append(grna['sequence'])
    
    # Save to file if requested
    if args.save_to:
        try:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(args.save_to), exist_ok=True)
            
            with open(args.save_to, 'w') as f:
                if args.output_format == 'sequence':
                    # Simple format: one sequence per line
                    for seq in grna_sequences:
                        f.write(seq + '\n')
                elif args.output_format == 'fasta':
                    # FASTA format
                    for i, grna in enumerate(grnas):
                        f.write(f">{args.gene}_gRNA_{i+1} score={grna['scores']['composite']:.3f}\n")
                        f.write(grna['sequence'] + '\n')
                elif args.output_format == 'detailed':
                    # Detailed format with metadata
                    f.write(f"# gRNAs for {args.gene} ({args.species})\n")
                    f.write(f"# Cas type: {args.cas_type}\n")
                    f.write(f"# Generated by CasPINS - Cas-Primer-Indel Suite\n\n")
                    for i, grna in enumerate(grnas):
                        f.write(f"# gRNA {i+1}\n")
                        f.write(f"# Position: {grna.get('position', 'N/A')}, Strand: {grna.get('strand', 'N/A')}\n")
                        f.write(f"# PAM: {grna.get('pam', 'N/A')}, GC: {grna.get('gc_content', 'N/A')}%\n")
                        f.write(f"# Score: {grna['scores']['composite']:.3f}\n")
                        f.write(grna['sequence'] + '\n\n')
            
            print(f"\n✓ Saved {len(grna_sequences)} gRNAs to {args.save_to}")
            print(f"\nYou can now run: python run.py analysis")
            
        except Exception as e:
            print(f"\nError saving file: {e}")
    
    # Export detailed CSV if requested
    if args.export_csv:
        try:
            export_file = designer.export_results(results, format='csv', output_file=args.export_csv)
            print(f"\n✓ Detailed results exported to {export_file}")
        except Exception as e:
            print(f"\nError exporting CSV: {e}")
    
    # Instructions for next steps
    if not args.save_to:
        print("\n" + "-"*80)
        print("Next steps:")
        print(f"1. Copy the gRNA sequences above")
        print(f"2. Set your data directory (GUI Settings or CRISPR_DATA_DIR env var)")
        print(f"3. Save to: <data_dir>/{args.gene.lower()}/grna.txt")
        print(f"4. Run: python run.py analysis --data-dir <your_data_dir>")
        print("\nOr use --save-to option to save directly:")
        print(f"   python run.py grna {args.gene} --species {args.species} --save-to <data_dir>/{args.gene.lower()}/grna.txt")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
