#!/usr/bin/env python
"""
Benchmarking Script: CasPINS gRNA Design vs CRISPOR/CHOPCHOP

This script generates CasPINS gRNA designs for benchmark genes and outputs
structured data for comparison with CRISPOR and CHOPCHOP results.

Addresses NAR criticism (iii): "A valid experimental comparison with existing
technologies would be necessary."

WORKFLOW:
=========
1. Run this script to generate CasPINS gRNA designs for benchmark genes
2. Manually query CRISPOR (http://crispor.tefor.net/) and CHOPCHOP (https://chopchop.cbu.uib.no/)
   with the SAME input sequences
3. Save CRISPOR/CHOPCHOP results to benchmarking/data/
4. Run benchmark_compare_grna.py to generate comparison plots and statistics

Usage:
    python benchmark_grna_design.py [--genes TP53 ATE1 VEGFA DBH EMX1] [--species human]
"""

import sys
import os
import json
import csv
import argparse
from datetime import datetime

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from src.grna_design.grna_designer import GRNADesigner


# Benchmark gene set: genes with available CRISPOR and CHOPCHOP comparison data
# ATE1 and DBH replace BRCA1 and PCSK9 (CHOPCHOP lacks data for BRCA1/PCSK9)
# ATE1 and DBH also appear in the CasPINS manuscript figures
BENCHMARK_GENES = {
    'TP53': {
        'description': 'Tumor protein p53 - most frequently mutated gene in cancer',
        'ensembl_id': 'ENSG00000141510',
        'refseq_id': 'NM_000546',
        'reason': 'Standard CRISPR benchmark, extensively validated gRNAs available'
    },
    'ATE1': {
        'description': 'Arginyltransferase 1',
        'ensembl_id': 'ENSG00000107669',
        'refseq_id': '',
        'reason': 'Matches manuscript figures, CRISPOR and CHOPCHOP data available'
    },
    'VEGFA': {
        'description': 'Vascular endothelial growth factor A',
        'ensembl_id': 'ENSG00000112715',
        'refseq_id': 'NM_001025366',
        'reason': 'Commonly used in CRISPR specificity studies'
    },
    'DBH': {
        'description': 'Dopamine beta-hydroxylase',
        'ensembl_id': 'ENSG00000123454',
        'refseq_id': '',
        'reason': 'Matches manuscript figures, CRISPOR and CHOPCHOP data available'
    },
    'EMX1': {
        'description': 'Empty spiracles homeobox 1',
        'ensembl_id': 'ENSG00000170370',
        'refseq_id': 'NM_004097',
        'reason': 'Most commonly used CRISPR positive control gene'
    }
}


def run_caspins_grna_design(gene_symbol: str, species: str = 'human', 
                             cas_type: str = 'SpCas9', n_results: int = 20) -> dict:
    """
    Run CasPINS gRNA design for a single gene and return structured results.
    
    Args:
        gene_symbol: Gene symbol (e.g., 'TP53')
        species: Species name (default: 'human')
        cas_type: Cas variant (default: 'SpCas9')
        n_results: Number of top gRNAs to return
    
    Returns:
        Dictionary with gRNA results and metadata
    """
    print(f"\n{'='*60}")
    print(f"Designing gRNAs for {gene_symbol} ({species}, {cas_type})")
    print(f"{'='*60}")
    
    designer = GRNADesigner(
        species=species,
        cas_type=cas_type
    )
    
    try:
        results = designer.design_grnas(
            target=gene_symbol,
            n_results=n_results,
            target_type='auto',
            filters={
                'gc_min': 20,  # Use wide filters to get more candidates
                'gc_max': 80,
                'homopolymer_max': 5
            },
            include_off_targets=False
        )
        
        # Structure output for comparison
        output = {
            'gene': gene_symbol,
            'species': species,
            'cas_type': cas_type,
            'timestamp': datetime.now().isoformat(),
            'tool': 'CasPINS',
            'tool_version': '1.0.0',
            'metadata': results.get('metadata', {}),
            'n_total_found': len(results.get('grnas', [])),
            'grnas': []
        }
        
        for i, grna in enumerate(results.get('grnas', []), 1):
            grna_entry = {
                'rank': i,
                'sequence': grna.get('sequence', ''),
                'pam': grna.get('pam', ''),
                'strand': grna.get('strand', ''),
                'position': grna.get('position', 0),
                'gc_content': grna.get('gc_content', 0),
                'genomic_location': grna.get('location_simple', ''),
                'scores': grna.get('scores', {}),
                'composite_score': grna.get('scores', {}).get('composite', 0),
                'doench_2016_score': grna.get('scores', {}).get('doench_2016', 0),
                'moreno_mateos_score': grna.get('scores', {}).get('moreno_mateos', 0),
                'xu_score': grna.get('scores', {}).get('xu', 0),
            }
            output['grnas'].append(grna_entry)
        
        print(f"  Found {output['n_total_found']} gRNAs for {gene_symbol}")
        for g in output['grnas'][:5]:
            print(f"    #{g['rank']}: {g['sequence']} | PAM: {g['pam']} | "
                  f"Composite: {g['composite_score']:.3f} | "
                  f"Doench: {g['doench_2016_score']:.3f}")
        
        return output
        
    except Exception as e:
        print(f"  ERROR designing gRNAs for {gene_symbol}: {e}")
        return {
            'gene': gene_symbol,
            'species': species,
            'cas_type': cas_type,
            'error': str(e),
            'timestamp': datetime.now().isoformat(),
            'tool': 'CasPINS',
            'grnas': []
        }


def export_for_comparison(results: list, output_dir: str):
    """Export results in formats suitable for comparison with CRISPOR/CHOPCHOP."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Full JSON results
    json_path = os.path.join(output_dir, 'caspins_grna_benchmark_results.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved JSON results to: {json_path}")
    
    # 2. CSV for easy spreadsheet comparison
    csv_path = os.path.join(output_dir, 'caspins_grna_benchmark_results.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Gene', 'Rank', 'gRNA_Sequence', 'PAM', 'Strand', 
            'GC_Content', 'Genomic_Location',
            'Composite_Score', 'Doench_2016', 'Moreno_Mateos', 'Xu_Score'
        ])
        for gene_result in results:
            for grna in gene_result.get('grnas', []):
                writer.writerow([
                    gene_result['gene'],
                    grna['rank'],
                    grna['sequence'],
                    grna['pam'],
                    grna['strand'],
                    f"{grna['gc_content']:.1f}" if grna['gc_content'] else '',
                    grna['genomic_location'],
                    f"{grna['composite_score']:.4f}",
                    f"{grna['doench_2016_score']:.4f}",
                    f"{grna['moreno_mateos_score']:.4f}",
                    f"{grna['xu_score']:.4f}"
                ])
    print(f"Saved CSV results to: {csv_path}")
    
    # 3. Generate input sequences for CRISPOR/CHOPCHOP manual queries
    instructions_path = os.path.join(output_dir, 'CRISPOR_CHOPCHOP_QUERY_INSTRUCTIONS.md')
    with open(instructions_path, 'w') as f:
        f.write("# Instructions for Querying CRISPOR and CHOPCHOP\n\n")
        f.write("## Purpose\n")
        f.write("Run the same genes through CRISPOR and CHOPCHOP to compare gRNA rankings.\n\n")
        f.write("## Steps\n\n")
        f.write("### CRISPOR (http://crispor.tefor.net/)\n")
        f.write("1. Go to http://crispor.tefor.net/\n")
        f.write("2. For each gene below, paste the gene symbol or Ensembl ID\n")
        f.write("3. Select: Genome = Human (hg38), PAM = NGG (SpCas9)\n")
        f.write("4. Download results as TSV/CSV\n")
        f.write("5. Save to: benchmarking/data/crispor_{gene}.tsv\n\n")
        f.write("### CHOPCHOP (https://chopchop.cbu.uib.no/)\n")
        f.write("1. Go to https://chopchop.cbu.uib.no/\n")
        f.write("2. For each gene, enter the gene symbol\n")
        f.write("3. Select: Genome = hg38, Method = Cas9\n")
        f.write("4. Download results\n")
        f.write("5. Save to: benchmarking/data/chopchop_{gene}.tsv\n\n")
        f.write("## Genes to Query\n\n")
        f.write("| Gene | Ensembl ID | RefSeq ID | Reason |\n")
        f.write("|------|-----------|-----------|--------|\n")
        for gene, info in BENCHMARK_GENES.items():
            f.write(f"| {gene} | {info['ensembl_id']} | {info['refseq_id']} | {info['reason']} |\n")
        
        f.write("\n## Comparison Metrics\n\n")
        f.write("After collecting all results, compare:\n")
        f.write("1. **Sequence overlap**: How many of CasPINS top-10 gRNAs appear in CRISPOR/CHOPCHOP top-20?\n")
        f.write("2. **Rank correlation**: Spearman correlation of shared gRNA rankings\n")
        f.write("3. **Score concordance**: Correlation of Doench 2016 scores between tools\n")
        f.write("4. **Unique finds**: gRNAs found by one tool but not others\n\n")
        f.write("## Expected Results\n\n")
        f.write("Since CasPINS implements the same published scoring algorithms (Doench 2016, \n")
        f.write("Moreno-Mateos, Xu), we expect:\n")
        f.write("- High overlap (>70%) in top-10 gRNAs between tools\n")
        f.write("- Strong rank correlation (Spearman rho > 0.7) for shared gRNAs\n")
        f.write("- Minor differences due to implementation details and filter differences\n")
    
    print(f"Saved query instructions to: {instructions_path}")
    
    # 4. Summary statistics
    summary_path = os.path.join(output_dir, 'caspins_benchmark_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("CasPINS gRNA Design Benchmark Summary\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n")
        f.write(f"{'='*50}\n\n")
        for gene_result in results:
            gene = gene_result['gene']
            grnas = gene_result.get('grnas', [])
            if grnas:
                scores = [g['composite_score'] for g in grnas]
                f.write(f"Gene: {gene}\n")
                f.write(f"  gRNAs found: {gene_result['n_total_found']}\n")
                f.write(f"  Top score: {max(scores):.4f}\n")
                f.write(f"  Mean score (top {len(scores)}): {sum(scores)/len(scores):.4f}\n")
                f.write(f"  Score range: {min(scores):.4f} - {max(scores):.4f}\n")
                f.write(f"  Top gRNA: {grnas[0]['sequence']} (PAM: {grnas[0]['pam']})\n\n")
            else:
                f.write(f"Gene: {gene} - ERROR: {gene_result.get('error', 'Unknown')}\n\n")
    
    print(f"Saved summary to: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Benchmark CasPINS gRNA design against CRISPOR/CHOPCHOP'
    )
    parser.add_argument('--genes', nargs='+', default=list(BENCHMARK_GENES.keys()),
                       help='Genes to benchmark (default: TP53 ATE1 VEGFA DBH EMX1)')
    parser.add_argument('--species', default='human', help='Species (default: human)')
    parser.add_argument('--cas-type', default='SpCas9', help='Cas variant (default: SpCas9)')
    parser.add_argument('--n-results', type=int, default=50000, 
                       help='Number of top gRNAs per gene (default: 50000, effectively returns all)')
    parser.add_argument('--output-dir', default=None,
                       help='Output directory (default: benchmarking/results/grna_design/)')
    
    args = parser.parse_args()
    
    if args.output_dir is None:
        args.output_dir = os.path.join(
            os.path.dirname(__file__), '..', 'results', 'grna_design'
        )
    
    print("=" * 60)
    print("CasPINS gRNA Design Benchmarking")
    print("=" * 60)
    print(f"Genes: {', '.join(args.genes)}")
    print(f"Species: {args.species}")
    print(f"Cas type: {args.cas_type}")
    print(f"Top N gRNAs: {args.n_results}")
    
    all_results = []
    for gene in args.genes:
        result = run_caspins_grna_design(
            gene_symbol=gene,
            species=args.species,
            cas_type=args.cas_type,
            n_results=args.n_results
        )
        all_results.append(result)
    
    export_for_comparison(all_results, args.output_dir)
    
    print(f"\n{'='*60}")
    print("BENCHMARKING COMPLETE")
    print(f"{'='*60}")
    print(f"\nNext steps:")
    print(f"1. Review results in: {args.output_dir}")
    print(f"2. Follow instructions in CRISPOR_CHOPCHOP_QUERY_INSTRUCTIONS.md")
    print(f"3. Run benchmark_compare_grna.py after collecting CRISPOR/CHOPCHOP data")


if __name__ == '__main__':
    main()
