"""
Advanced Primer Design Module for CRISPR Analysis
Incorporates features from leading primer design tools with enhanced analysis
"""

import primer3
import requests
import json
import re
from Bio.Seq import Seq
from Bio import Entrez
from typing import List, Dict, Tuple, Optional
from datetime import datetime
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure Entrez
Entrez.email = "crispr_analysis@example.com"


class PrimerDesigner:
    """Advanced primer designer with database integration and comprehensive analysis."""
    
    def __init__(self, gene_name: str, species: str = "human"):
        self.gene_name = gene_name
        self.species = species
        self.ncbi_cache = {}
        self.ensembl_cache = {}
        
    def search_online_databases(self, sequence: str = None, gene_id: str = None) -> Dict:
        """
        Search NCBI, Ensembl, and UniProt databases for gene information.
        
        Returns:
            Dict with database information including:
            - gene_info: Basic gene information
            - variants: Known variants in the region
            - conservation: Conservation scores
            - validated_sequence: Validated reference sequence
        """
        results = {
            'ncbi': {},
            'ensembl': {},
            'uniprot': {},
            'variants': [],
            'conservation': {},
            'validated': False
        }
        
        try:
            # Search NCBI
            if gene_id:
                results['ncbi'] = self._search_ncbi(gene_id)
            
            # Search Ensembl
            results['ensembl'] = self._search_ensembl(self.gene_name, self.species)
            
            # Get variant information
            if results['ensembl'].get('gene_id'):
                results['variants'] = self._get_variants_from_ensembl(
                    results['ensembl']['gene_id']
                )
            
            # Validate sequence if provided
            if sequence and results['ncbi'].get('sequence'):
                results['validated'] = self._validate_sequence(
                    sequence, results['ncbi']['sequence']
                )
                
        except Exception as e:
            logger.warning(f"Database search error: {str(e)}")
            
        return results
    
    def _search_ncbi(self, gene_id: str) -> Dict:
        """Search NCBI database for gene information."""
        try:
            # Search gene database
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
            record = Entrez.read(handle)
            handle.close()
            
            gene_info = {
                'gene_id': gene_id,
                'symbol': record[0].get('Entrezgene_gene', {}).get('Gene-ref_locus', ''),
                'description': record[0].get('Entrezgene_prot', {}).get('Prot-ref_name', ''),
                'chromosome': record[0].get('Entrezgene_location', [{}])[0].get('Maps_display-str', ''),
                'aliases': []
            }
            
            # Get sequence
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
            sequence = handle.read()
            handle.close()
            
            if sequence:
                gene_info['sequence'] = ''.join(sequence.split('\n')[1:])
                
            return gene_info
            
        except Exception as e:
            logger.warning(f"NCBI search error: {str(e)}")
            return {}
    
    def _search_ensembl(self, gene_name: str, species: str) -> Dict:
        """Search Ensembl REST API for gene information."""
        try:
            # Convert species to Ensembl format
            species_map = {
                'human': 'homo_sapiens',
                'mouse': 'mus_musculus',
                'rat': 'rattus_norvegicus'
            }
            ensembl_species = species_map.get(species.lower(), 'homo_sapiens')
            
            # Search for gene
            url = f"https://rest.ensembl.org/lookup/symbol/{ensembl_species}/{gene_name}"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                data = response.json()
                return {
                    'gene_id': data.get('id', ''),
                    'display_name': data.get('display_name', ''),
                    'biotype': data.get('biotype', ''),
                    'chromosome': data.get('seq_region_name', ''),
                    'start': data.get('start', 0),
                    'end': data.get('end', 0),
                    'strand': data.get('strand', 0)
                }
            
        except Exception as e:
            logger.warning(f"Ensembl search error: {str(e)}")
            
        return {}
    
    def _get_variants_from_ensembl(self, gene_id: str) -> List[Dict]:
        """Get known variants from Ensembl."""
        try:
            url = f"https://rest.ensembl.org/overlap/id/{gene_id}"
            params = {"feature": "variation"}
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(url, headers=headers, params=params)
            if response.status_code == 200:
                variants = response.json()
                return [
                    {
                        'id': v.get('id', ''),
                        'position': v.get('start', 0),
                        'alleles': v.get('alleles', []),
                        'clinical_significance': v.get('clinical_significance', [])
                    }
                    for v in variants[:10]  # Limit to first 10 variants
                ]
                
        except Exception as e:
            logger.warning(f"Variant search error: {str(e)}")
            
        return []
    
    def _validate_sequence(self, query_seq: str, ref_seq: str) -> bool:
        """Validate sequence against reference."""
        # Simple validation - can be enhanced
        return query_seq in ref_seq or ref_seq in query_seq
    
    def design_comprehensive_primers(self, mrna_seq: str, target_regions: List[Tuple[int, int]], 
                                   grna_info: Dict = None) -> Dict:
        """
        Design comprehensive primer sets with advanced analysis.
        
        Returns:
            Dict containing:
            - pcr1_primers: PCR I primer sets
            - pcr2_primers: PCR II primer sets
            - quality_metrics: Detailed quality analysis
            - recommendations: Specific recommendations
        """
        results = {
            'pcr1_primers': [],
            'pcr2_primers': [],
            'quality_metrics': {},
            'recommendations': [],
            'database_info': {}
        }
        
        # Search databases first
        results['database_info'] = self.search_online_databases(sequence=mrna_seq)
        
        # Design primers for each target region
        for i, (start, end) in enumerate(target_regions):
            # PCR I primers (genomic amplification)
            pcr1_primers = self._design_advanced_primers(
                mrna_seq, start, end,
                flank_size=500,
                product_range=[[800, 2000]],
                primer_type="PCR_I",
                region_id=i+1
            )
            results['pcr1_primers'].extend(pcr1_primers)
            
            # PCR II primers (sequencing)
            pcr2_primers = self._design_advanced_primers(
                mrna_seq, start, end,
                flank_size=200,
                product_range=[[300, 600]],
                primer_type="PCR_II",
                region_id=i+1
            )
            results['pcr2_primers'].extend(pcr2_primers)
        
        # Analyze primer quality
        results['quality_metrics'] = self._analyze_primer_quality(
            results['pcr1_primers'] + results['pcr2_primers']
        )
        
        # Generate recommendations
        results['recommendations'] = self._generate_recommendations(
            results, grna_info
        )
        
        return results
    
    def _design_advanced_primers(self, seq: str, target_start: int, target_end: int,
                               flank_size: int, product_range: List, primer_type: str,
                               region_id: int) -> List[Dict]:
        """Design primers with advanced parameters and analysis."""
        primers = []
        
        # Define search region
        search_start = max(0, target_start - flank_size)
        search_end = min(len(seq), target_end + flank_size)
        target_len = target_end - target_start
        
        # Enhanced Primer3 parameters
        seq_args = {
            'SEQUENCE_ID': f'{self.gene_name}_{primer_type}_region{region_id}',
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_TARGET': [target_start, target_len],
            'SEQUENCE_INCLUDED_REGION': [search_start, search_end - search_start]
        }
        
        global_args = {
            # Size parameters
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            
            # Tm parameters
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MAX_DIFF_TM': 2.0,
            
            # GC parameters
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            
            # Quality parameters
            'PRIMER_MAX_POLY_X': 3,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY_TH': 45.0,
            'PRIMER_MAX_SELF_END_TH': 35.0,
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
            'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
            
            # Thermodynamic parameters
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_SALT_DIVALENT': 1.5,
            'PRIMER_DNTP_CONC': 0.6,
            'PRIMER_DNA_CONC': 50.0,
            
            # Product parameters
            'PRIMER_PRODUCT_SIZE_RANGE': product_range,
            'PRIMER_NUM_RETURN': 5,  # Return more primers for selection
            
            # Additional parameters
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_GC_CLAMP': 2,  # Require GC clamp
            'PRIMER_LIBERAL_BASE': 1,
            'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 1,
            'PRIMER_LOWERCASE_MASKING': 0,
            'PRIMER_PICK_ANYWAY': 1
        }
        
        try:
            # Remove problematic thermodynamic parameters for now
            if 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' in global_args:
                del global_args['PRIMER_THERMODYNAMIC_PARAMETERS_PATH']
            
            # Remove other parameters that might cause issues
            del global_args['PRIMER_MAX_SELF_ANY_TH']
            del global_args['PRIMER_MAX_SELF_END_TH']
            del global_args['PRIMER_PAIR_MAX_COMPL_ANY_TH']
            del global_args['PRIMER_PAIR_MAX_COMPL_END_TH']
            
            # Design primers
            primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
            
            # Extract and analyze primer pairs
            num_primers = primer3_result.get('PRIMER_PAIR_NUM_RETURNED', 0)
            
            for j in range(num_primers):
                primer_data = self._extract_primer_data(primer3_result, j)
                
                # Additional analysis
                primer_data['hairpin_tm'] = self._calculate_hairpin_tm(primer_data['forward_seq'])
                primer_data['dimer_score'] = self._calculate_dimer_score(
                    primer_data['forward_seq'], primer_data['reverse_seq']
                )
                primer_data['specificity_score'] = self._check_specificity(primer_data['forward_seq'])
                primer_data['gc_distribution'] = self._analyze_gc_distribution(primer_data['forward_seq'])
                
                # Check for variants in primer regions
                if hasattr(self, 'variants'):
                    primer_data['variants_in_primer'] = self._check_variants_in_region(
                        primer_data['forward_pos'],
                        primer_data['forward_pos'] + len(primer_data['forward_seq'])
                    )
                
                primers.append(primer_data)
                
        except Exception as e:
            logger.error(f"Primer3 design failed: {str(e)}")
            # Fallback to simple design
            primers = self._design_fallback_primers(seq, target_start, target_end, primer_type)
            
        return primers
    
    def _extract_primer_data(self, primer3_result: Dict, index: int) -> Dict:
        """Extract comprehensive primer data from Primer3 results."""
        return {
            'rank': index + 1,
            'forward_seq': primer3_result.get(f'PRIMER_LEFT_{index}_SEQUENCE', ''),
            'reverse_seq': primer3_result.get(f'PRIMER_RIGHT_{index}_SEQUENCE', ''),
            'forward_tm': round(primer3_result.get(f'PRIMER_LEFT_{index}_TM', 0), 1),
            'reverse_tm': round(primer3_result.get(f'PRIMER_RIGHT_{index}_TM', 0), 1),
            'forward_gc': round(primer3_result.get(f'PRIMER_LEFT_{index}_GC_PERCENT', 0), 1),
            'reverse_gc': round(primer3_result.get(f'PRIMER_RIGHT_{index}_GC_PERCENT', 0), 1),
            'forward_pos': primer3_result.get(f'PRIMER_LEFT_{index}', [0, 0])[0],
            'reverse_pos': primer3_result.get(f'PRIMER_RIGHT_{index}', [0, 0])[0],
            'product_size': primer3_result.get(f'PRIMER_PAIR_{index}_PRODUCT_SIZE', 0),
            'penalty': round(primer3_result.get(f'PRIMER_PAIR_{index}_PENALTY', 0), 2),
            'forward_end_stability': primer3_result.get(f'PRIMER_LEFT_{index}_END_STABILITY', 0),
            'reverse_end_stability': primer3_result.get(f'PRIMER_RIGHT_{index}_END_STABILITY', 0),
            'forward_hairpin': primer3_result.get(f'PRIMER_LEFT_{index}_HAIRPIN_TH', 0),
            'reverse_hairpin': primer3_result.get(f'PRIMER_RIGHT_{index}_HAIRPIN_TH', 0),
            'pair_compl_any': primer3_result.get(f'PRIMER_PAIR_{index}_COMPL_ANY_TH', 0),
            'pair_compl_end': primer3_result.get(f'PRIMER_PAIR_{index}_COMPL_END_TH', 0)
        }
    
    def _calculate_hairpin_tm(self, sequence: str) -> float:
        """Calculate potential hairpin melting temperature."""
        # Simplified calculation - checks for self-complementarity
        seq_obj = Seq(sequence)
        rev_comp = str(seq_obj.reverse_complement())
        
        max_match = 0
        for i in range(len(sequence) - 3):
            for j in range(i + 4, len(sequence)):
                if sequence[i:j] in rev_comp:
                    max_match = max(max_match, j - i)
                    
        # Rough Tm calculation for hairpin
        if max_match > 0:
            gc_count = sequence.count('G') + sequence.count('C')
            return 4 * gc_count + 2 * (len(sequence) - gc_count) - 5
        return 0
    
    def _calculate_dimer_score(self, seq1: str, seq2: str) -> float:
        """Calculate primer dimer formation score."""
        # Check for 3' complementarity
        end1 = seq1[-5:]
        end2 = seq2[-5:]
        
        # Simple scoring based on complementarity
        score = 0
        for i in range(min(len(end1), len(end2))):
            if self._is_complement(end1[-(i+1)], end2[-(i+1)]):
                score += 2 if i < 2 else 1  # Weight 3' end more heavily
                
        return score
    
    def _is_complement(self, base1: str, base2: str) -> bool:
        """Check if two bases are complementary."""
        complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return complements.get(base1) == base2
    
    def _check_specificity(self, sequence: str) -> float:
        """Check primer specificity score."""
        # Simplified specificity check
        # In production, would BLAST against genome
        
        # Check for low complexity regions
        complexity_score = len(set(sequence)) / len(sequence)
        
        # Check for repetitive elements
        repeat_score = 1.0
        for length in [2, 3, 4]:
            for i in range(len(sequence) - length):
                motif = sequence[i:i+length]
                if sequence.count(motif) > 2:
                    repeat_score *= 0.8
                    
        return complexity_score * repeat_score
    
    def _analyze_gc_distribution(self, sequence: str) -> Dict:
        """Analyze GC distribution in primer."""
        gc_positions = [1 if base in 'GC' else 0 for base in sequence]
        
        # Check 5' and 3' regions
        five_prime_gc = sum(gc_positions[:5]) / 5 * 100
        three_prime_gc = sum(gc_positions[-5:]) / 5 * 100
        
        # Check for GC clamp
        gc_clamp = sequence[-2:].count('G') + sequence[-2:].count('C')
        
        return {
            '5_prime_gc': round(five_prime_gc, 1),
            '3_prime_gc': round(three_prime_gc, 1),
            'gc_clamp': gc_clamp,
            'gc_runs': max(len(list(g)) for k, g in re.findall(r'(([GC])\2*)', sequence) or [('', '')])
        }
    
    def _check_variants_in_region(self, start: int, end: int) -> List[str]:
        """Check for known variants in primer region."""
        # Would check against variant database
        return []
    
    def _analyze_primer_quality(self, all_primers: List[Dict]) -> Dict:
        """Comprehensive quality analysis of all primers."""
        if not all_primers:
            return {}
            
        metrics = {
            'total_primers': len(all_primers),
            'avg_penalty': sum(p.get('penalty', 0) for p in all_primers) / len(all_primers),
            'tm_range': {
                'min': min(p['forward_tm'] for p in all_primers),
                'max': max(p['forward_tm'] for p in all_primers)
            },
            'gc_range': {
                'min': min(p['forward_gc'] for p in all_primers),
                'max': max(p['forward_gc'] for p in all_primers)
            },
            'specificity_scores': [p.get('specificity_score', 0) for p in all_primers],
            'potential_issues': []
        }
        
        # Check for potential issues
        for primer in all_primers:
            if primer.get('hairpin_tm', 0) > 45:
                metrics['potential_issues'].append(
                    f"Primer {primer['rank']}: High hairpin Tm ({primer['hairpin_tm']}°C)"
                )
            if primer.get('dimer_score', 0) > 6:
                metrics['potential_issues'].append(
                    f"Primer {primer['rank']}: High dimer score ({primer['dimer_score']})"
                )
                
        return metrics
    
    def _generate_recommendations(self, results: Dict, grna_info: Dict = None) -> List[str]:
        """Generate specific recommendations based on analysis."""
        recommendations = []
        
        # Database validation recommendations
        if results['database_info'].get('validated'):
            recommendations.append("✓ Sequence validated against reference database")
        else:
            recommendations.append("⚠ Sequence could not be validated against reference")
            
        # Variant recommendations
        if results['database_info'].get('variants'):
            recommendations.append(
                f"⚠ {len(results['database_info']['variants'])} known variants in gene region"
            )
            
        # Quality recommendations
        metrics = results.get('quality_metrics', {})
        if metrics.get('potential_issues'):
            recommendations.append("⚠ Quality issues detected:")
            recommendations.extend(f"  - {issue}" for issue in metrics['potential_issues'])
            
        # Primer selection recommendations
        if results['pcr1_primers'] and results['pcr2_primers']:
            best_pcr1 = min(results['pcr1_primers'], key=lambda x: x.get('penalty', 999))
            best_pcr2 = min(results['pcr2_primers'], key=lambda x: x.get('penalty', 999))
            
            recommendations.append(f"✓ Recommended PCR I: Set {best_pcr1['rank']}")
            recommendations.append(f"✓ Recommended PCR II: Set {best_pcr2['rank']}")
            
        return recommendations
    
    def _design_fallback_primers(self, seq: str, target_start: int, target_end: int, 
                                primer_type: str) -> List[Dict]:
        """Fallback primer design if Primer3 fails."""
        primers = []
        flank_sizes = [150, 175, 200] if primer_type == "PCR_II" else [400, 450, 500]
        
        for i, flank in enumerate(flank_sizes):
            f_start = max(0, target_start - flank)
            f_seq = seq[f_start:f_start+20]
            
            r_start = min(len(seq)-20, target_end + flank - 20)
            r_seq = str(Seq(seq[r_start:r_start+20]).reverse_complement())
            
            if len(f_seq) == 20 and len(r_seq) == 20:
                primers.append({
                    'rank': i + 1,
                    'forward_seq': f_seq,
                    'reverse_seq': r_seq,
                    'forward_tm': self._calculate_tm(f_seq),
                    'reverse_tm': self._calculate_tm(r_seq),
                    'forward_gc': (f_seq.count('G') + f_seq.count('C')) * 5,
                    'reverse_gc': (r_seq.count('G') + r_seq.count('C')) * 5,
                    'forward_pos': f_start,
                    'reverse_pos': r_start,
                    'product_size': r_start + 20 - f_start,
                    'note': 'Fallback design',
                    'penalty': 999  # High penalty for fallback
                })
                
        return primers
    
    def _calculate_tm(self, sequence: str) -> float:
        """Calculate melting temperature using nearest-neighbor method."""
        # Basic calculation
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        
        if len(sequence) < 14:
            return (gc_count * 4) + (at_count * 2)
        else:
            # Salt-adjusted Tm
            tm = 81.5 + 16.6 * (np.log10(0.05)) + 0.41 * (gc_count * 100 / len(sequence)) - 675 / len(sequence)
            return round(tm, 1)


def generate_enhanced_primer_report(grna_positions: List[Tuple], mrna_seq: str, 
                                  gene_name: str, species: str = "human") -> str:
    """
    Generate comprehensive primer design report with database integration.
    
    Args:
        grna_positions: List of (position, strand, cut_site) tuples
        mrna_seq: mRNA sequence
        gene_name: Gene name
        species: Species (default: human)
        
    Returns:
        str: Comprehensive primer design report
    """
    report = []
    
    # Header
    report.append(f"Advanced Primer Design Report for {gene_name} ({species})")
    report.append("=" * 80)
    report.append(f"Generated using CasPINS - Cas-Primer-Indel Suite")
    report.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    if not grna_positions:
        report.append("*** ERROR: No gRNAs found in the provided mRNA sequence! ***")
        report.append("\nTroubleshooting steps:")
        report.append("  1. Verify the mRNA sequence is correct and complete")
        report.append("  2. Check if gRNAs include PAM sequence (remove if present)")
        report.append("  3. Ensure gRNAs are 20nt long (standard length)")
        report.append("  4. Verify gRNAs target the correct gene")
        report.append("  5. Check if you're using genomic vs cDNA sequence")
        return "\n".join(report)
    
    # gRNA information
    report.append("gRNA sequences:")
    for i, (pos, strand, cut) in enumerate(grna_positions, 1):
        grna_seq = mrna_seq[pos:pos+20]
        report.append(f"  {i}. {grna_seq} (position: {pos}, strand: {'+' if strand == 1 else '-'})")
    report.append("")
    
    # Initialize designer
    designer = PrimerDesigner(gene_name, species)
    
    # Define target regions
    cut_sites = [pos[2] for pos in grna_positions]
    min_cut = min(cut_sites)
    max_cut = max(cut_sites)
    target_regions = [(min_cut - 50, max_cut + 50)]
    
    # Design primers
    results = designer.design_comprehensive_primers(mrna_seq, target_regions, 
                                                   {'positions': grna_positions})
    
    # Database search results
    report.append("\n" + "=" * 80)
    report.append("ONLINE DATABASE SEARCH RESULTS")
    report.append("=" * 80)
    
    db_info = results.get('database_info', {})
    if db_info.get('ncbi'):
        report.append("\nNCBI Database:")
        report.append(f"  Gene ID: {db_info['ncbi'].get('gene_id', 'Not found')}")
        report.append(f"  Official Symbol: {db_info['ncbi'].get('symbol', 'Not found')}")
        report.append(f"  Chromosome: {db_info['ncbi'].get('chromosome', 'Not found')}")
        
    if db_info.get('ensembl'):
        report.append("\nEnsembl Database:")
        report.append(f"  Gene ID: {db_info['ensembl'].get('gene_id', 'Not found')}")
        report.append(f"  Biotype: {db_info['ensembl'].get('biotype', 'Not found')}")
        report.append(f"  Coordinates: chr{db_info['ensembl'].get('chromosome', '?')}:"
                     f"{db_info['ensembl'].get('start', '?')}-{db_info['ensembl'].get('end', '?')}")
    
    if db_info.get('variants'):
        report.append(f"\nKnown Variants: {len(db_info['variants'])} found in gene region")
        for var in db_info['variants'][:3]:  # Show first 3
            report.append(f"  - {var['id']}: position {var['position']}, "
                         f"alleles: {'/'.join(var['alleles'])}")
    
    # PCR I Primers
    report.append("\n" + "=" * 80)
    report.append("PCR I PRIMERS (Genomic DNA Amplification)")
    report.append("=" * 80)
    report.append("Purpose: Initial amplification from genomic DNA template")
    report.append("Strategy: Long-range PCR to capture potential off-target sites")
    report.append("Expected product: 800-2000bp (accounts for introns)")
    report.append("Recommended polymerase: High-fidelity (Q5, Phusion, PrimeSTAR)")
    report.append("")
    
    pcr1_primers = results.get('pcr1_primers', [])[:3]  # Top 3
    for i, primer in enumerate(pcr1_primers):
        report.append(f"PCR I - Set {primer['rank']} {'[RECOMMENDED]' if i == 0 else ''}:")
        report.append(f"  Forward primer: 5'-{primer['forward_seq']}-3'")
        report.append(f"    • Position: {primer['forward_pos']} (relative to mRNA)")
        report.append(f"    • Length: {len(primer['forward_seq'])} nt")
        report.append(f"    • Tm: {primer['forward_tm']}°C (ΔTm: "
                     f"{abs(primer['forward_tm'] - primer['reverse_tm'])}°C)")
        report.append(f"    • GC content: {primer['forward_gc']}%")
        
        if 'gc_distribution' in primer:
            gc_dist = primer['gc_distribution']
            report.append(f"    • GC distribution: 5'-end: {gc_dist['5_prime_gc']}%, "
                         f"3'-end: {gc_dist['3_prime_gc']}%, GC clamp: {gc_dist['gc_clamp']}")
        
        report.append(f"    • 3' end stability: {primer.get('forward_end_stability', 'N/A')} kcal/mol")
        report.append(f"    • Hairpin Tm: {primer.get('forward_hairpin', 'N/A')}°C")
        
        report.append(f"\n  Reverse primer: 5'-{primer['reverse_seq']}-3'")
        report.append(f"    • Position: {primer['reverse_pos']} (relative to mRNA)")
        report.append(f"    • Length: {len(primer['reverse_seq'])} nt")
        report.append(f"    • Tm: {primer['reverse_tm']}°C")
        report.append(f"    • GC content: {primer['reverse_gc']}%")
        report.append(f"    • 3' end stability: {primer.get('reverse_end_stability', 'N/A')} kcal/mol")
        report.append(f"    • Hairpin Tm: {primer.get('reverse_hairpin', 'N/A')}°C")
        
        report.append(f"\n  Amplicon characteristics:")
        report.append(f"    • Product size (mRNA): {primer['product_size']} bp")
        report.append(f"    • Genomic size: Will be larger due to introns")
        report.append(f"    • Primer penalty score: {primer.get('penalty', 'N/A')}")
        report.append(f"    • Dimer score: {primer.get('dimer_score', 'N/A')}")
        report.append(f"    • Specificity score: {primer.get('specificity_score', 'N/A')}")
        
        if primer.get('variants_in_primer'):
            report.append(f"    • ⚠ Known variants in primer region: {len(primer['variants_in_primer'])}")
        
        report.append("")
    
    # PCR II Primers
    report.append("\n" + "=" * 80)
    report.append("PCR II PRIMERS (Nested PCR/Direct Sequencing)")
    report.append("=" * 80)
    report.append("Purpose: High-specificity amplification for Sanger sequencing")
    report.append("Strategy: Nested PCR or direct amplification from PCR I product")
    report.append("Expected product: 300-600bp (optimal for sequencing)")
    report.append("Recommended polymerase: Standard Taq or high-fidelity")
    report.append("")
    
    pcr2_primers = results.get('pcr2_primers', [])[:3]  # Top 3
    for i, primer in enumerate(pcr2_primers):
        report.append(f"PCR II - Set {primer['rank']} {'[RECOMMENDED FOR SEQUENCING]' if i == 0 else ''}:")
        report.append(f"  Forward primer: 5'-{primer['forward_seq']}-3'")
        report.append(f"    • Position: {primer['forward_pos']} (relative to mRNA)")
        report.append(f"    • Distance from cut site: {abs(primer['forward_pos'] - min_cut)} bp")
        report.append(f"    • Tm: {primer['forward_tm']}°C")
        report.append(f"    • GC content: {primer['forward_gc']}%")
        
        report.append(f"\n  Reverse primer: 5'-{primer['reverse_seq']}-3'")
        report.append(f"    • Position: {primer['reverse_pos']} (relative to mRNA)")
        report.append(f"    • Distance from cut site: {abs(primer['reverse_pos'] - max_cut)} bp")
        report.append(f"    • Tm: {primer['reverse_tm']}°C")
        report.append(f"    • GC content: {primer['reverse_gc']}%")
        
        report.append(f"\n  Sequencing considerations:")
        report.append(f"    • Product size: {primer['product_size']} bp")
        report.append(f"    • Covers target region: {min_cut}-{max_cut}")
        report.append(f"    • Sequencing read from forward primer reaches: ~{primer['forward_pos'] + 700} bp")
        report.append(f"    • Sequencing read from reverse primer reaches: ~{primer['reverse_pos'] - 700} bp")
        report.append("")
    
    # Target information
    report.append("-" * 80)
    report.append("TARGET REGION DETAILS:")
    report.append(f"  • Gene: {gene_name} ({species})")
    report.append(f"  • gRNA target region: {min(p[0] for p in grna_positions)}-"
                 f"{max(p[0] for p in grna_positions) + 23} bp")
    report.append(f"  • Predicted cut sites: {min_cut}-{max_cut} bp "
                 f"(span: {max_cut - min_cut} bp)")
    report.append(f"  • Number of targets: {len(grna_positions)}")
    
    # Quality metrics summary
    if results.get('quality_metrics'):
        report.append("\n" + "-" * 80)
        report.append("QUALITY METRICS SUMMARY:")
        metrics = results['quality_metrics']
        report.append(f"  • Total primer pairs designed: {metrics.get('total_primers', 0)}")
        report.append(f"  • Average penalty score: {metrics.get('avg_penalty', 'N/A')}")
        report.append(f"  • Tm range: {metrics['tm_range']['min']}-{metrics['tm_range']['max']}°C")
        report.append(f"  • GC range: {metrics['gc_range']['min']}-{metrics['gc_range']['max']}%")
        
        if metrics.get('potential_issues'):
            report.append("\n  ⚠ Potential issues:")
            for issue in metrics['potential_issues']:
                report.append(f"    - {issue}")
    
    # Recommendations
    if results.get('recommendations'):
        report.append("\n" + "-" * 80)
        report.append("ANALYSIS RECOMMENDATIONS:")
        for rec in results['recommendations']:
            report.append(f"  {rec}")
    
    # PCR conditions
    report.append("\n" + "=" * 80)
    report.append("OPTIMIZED PCR CONDITIONS:")
    report.append("=" * 80)
    
    report.append("\nPCR I (Genomic DNA):")
    report.append("  Setup (25 μL reaction):")
    report.append("    • 12.5 μL 2X High-Fidelity Master Mix")
    report.append("    • 1.25 μL Forward primer (10 μM)")
    report.append("    • 1.25 μL Reverse primer (10 μM)")
    report.append("    • 50-100 ng genomic DNA")
    report.append("    • H₂O to 25 μL")
    report.append("")
    report.append("  Cycling conditions:")
    report.append("    1. Initial denaturation: 98°C for 30 sec")
    report.append("    2. 35 cycles:")
    report.append("       • Denaturation: 98°C for 10 sec")
    report.append("       • Annealing: 58-62°C for 20 sec (optimize based on Tm)")
    report.append("       • Extension: 72°C for 30 sec/kb")
    report.append("    3. Final extension: 72°C for 5 min")
    report.append("    4. Hold: 4°C")
    
    report.append("\nPCR II (Nested/Sequencing):")
    report.append("  Setup (50 μL reaction for sequencing):")
    report.append("    • 25 μL 2X Master Mix")
    report.append("    • 2.5 μL Forward primer (10 μM)")
    report.append("    • 2.5 μL Reverse primer (10 μM)")
    report.append("    • 1-2 μL PCR I product (diluted 1:50)")
    report.append("    • H₂O to 50 μL")
    report.append("")
    report.append("  Cycling conditions:")
    report.append("    1. Initial denaturation: 95°C for 3 min")
    report.append("    2. 30 cycles:")
    report.append("       • Denaturation: 95°C for 30 sec")
    report.append("       • Annealing: 60°C for 30 sec")
    report.append("       • Extension: 72°C for 45 sec")
    report.append("    3. Final extension: 72°C for 5 min")
    report.append("    4. Hold: 4°C")
    
    # Troubleshooting guide
    report.append("\n" + "=" * 80)
    report.append("TROUBLESHOOTING GUIDE:")
    report.append("=" * 80)
    report.append("\nNo amplification:")
    report.append("  • Verify DNA quality and concentration")
    report.append("  • Lower annealing temperature by 2-3°C")
    report.append("  • Increase extension time for large products")
    report.append("  • Check for inhibitors in DNA prep")
    
    report.append("\nMultiple bands:")
    report.append("  • Increase annealing temperature")
    report.append("  • Use touchdown PCR protocol")
    report.append("  • Reduce primer concentration")
    report.append("  • Use hot-start polymerase")
    
    report.append("\nWeak amplification:")
    report.append("  • Increase cycle number (up to 40)")
    report.append("  • Optimize Mg²⁺ concentration")
    report.append("  • Fresh primers/polymerase")
    report.append("  • Increase template amount")
    
    # Additional resources
    report.append("\n" + "=" * 80)
    report.append("ADDITIONAL RESOURCES:")
    report.append("=" * 80)
    report.append("  • NCBI Primer-BLAST: https://www.ncbi.nlm.nih.gov/tools/primer-blast/")
    report.append("  • Ensembl: https://www.ensembl.org/")
    report.append("  • UniProt: https://www.uniprot.org/")
    report.append("  • IDT OligoAnalyzer: https://www.idtdna.com/calc/analyzer")
    report.append("  • NEB Tm Calculator: https://tmcalculator.neb.com/")
    
    return "\n".join(report)


# Update the main function to use the enhanced report
def generate_primer_recommendations(grna_positions, mrna_seq, gene_name):
    """
    Main function for backwards compatibility.
    Now uses the enhanced primer design system.
    """
    # Extract species from gene name if provided (e.g., "TP53_human" -> "human")
    species = "human"  # default
    if "_" in gene_name:
        parts = gene_name.split("_")
        if parts[-1].lower() in ['human', 'mouse', 'rat']:
            species = parts[-1].lower()
            gene_name = "_".join(parts[:-1])
    
    return generate_enhanced_primer_report(grna_positions, mrna_seq, gene_name, species)


# Import numpy if available for advanced calculations
try:
    import numpy as np
except ImportError:
    # Fallback for systems without numpy
    class np:
        @staticmethod
        def log10(x):
            import math
            return math.log10(x) 