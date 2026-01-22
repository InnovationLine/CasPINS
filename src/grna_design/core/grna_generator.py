"""
gRNA Generator Module
Handles identification and filtering of potential guide RNAs
"""

import re
from typing import List, Dict, Tuple, Optional
from Bio.Seq import Seq
import numpy as np


class GRNAGenerator:
    """Generate and filter potential guide RNAs from sequences."""
    
    # PAM sequences for different Cas proteins
    PAM_SEQUENCES = {
        'SpCas9': ['NGG', 'NAG'],
        'SaCas9': ['NNGRRT', 'NNGRR'],
        'Cas12a': ['TTTV'],
        'SpCas9-NG': ['NG'],
        'xCas9': ['NG', 'GAA', 'GAT'],
        'Cas9-VQR': ['NGA'],
        'Cas9-EQR': ['NGAG'],
        'Cas9-VRER': ['NGCG']
    }
    
    def __init__(self, cas_type: str = 'SpCas9', grna_length: int = 20):
        """
        Initialize gRNA generator.
        
        Args:
            cas_type: Type of Cas protein (default: SpCas9)
            grna_length: Length of guide RNA (default: 20)
        """
        self.cas_type = cas_type
        self.grna_length = grna_length
        self.pam_sequences = self.PAM_SEQUENCES.get(cas_type, ['NGG'])
        
    def find_all_grnas(self, sequence: str, strand: str = 'both') -> List[Dict]:
        """
        Find all potential gRNAs in a sequence.
        
        Args:
            sequence: DNA sequence to search
            strand: Which strand to search ('forward', 'reverse', or 'both')
            
        Returns:
            List of gRNA dictionaries with positions and sequences
        """
        sequence = sequence.upper()
        grnas = []
        
        # Search forward strand
        if strand in ['forward', 'both']:
            grnas.extend(self._find_grnas_on_strand(sequence, 'forward'))
            
        # Search reverse strand
        if strand in ['reverse', 'both']:
            rev_comp = str(Seq(sequence).reverse_complement())
            rev_grnas = self._find_grnas_on_strand(rev_comp, 'reverse')
            # Adjust positions for reverse strand
            for grna in rev_grnas:
                grna['position'] = len(sequence) - grna['position'] - len(grna['sequence']) - len(grna['pam'])
            grnas.extend(rev_grnas)
            
        return grnas
    
    def _find_grnas_on_strand(self, sequence: str, strand: str) -> List[Dict]:
        """Find gRNAs on a specific strand."""
        grnas = []
        
        for pam in self.pam_sequences:
            # Convert PAM pattern to regex
            pam_regex = self._pam_to_regex(pam)
            
            # Find all PAM occurrences
            for match in re.finditer(pam_regex, sequence):
                pam_pos = match.start()
                
                # Extract gRNA sequence (upstream of PAM for Cas9)
                if self.cas_type.startswith('Cas12'):
                    # Cas12a has PAM upstream
                    grna_start = pam_pos + len(pam)
                    grna_end = grna_start + self.grna_length
                else:
                    # Cas9 variants have PAM downstream
                    grna_start = pam_pos - self.grna_length
                    grna_end = pam_pos
                
                # Check if gRNA is within sequence bounds
                if grna_start >= 0 and grna_end <= len(sequence):
                    grna_seq = sequence[grna_start:grna_end]
                    
                    # Calculate cut site (3-4 bp upstream of PAM for Cas9)
                    if self.cas_type.startswith('Cas12'):
                        cut_site = grna_start + 18  # Cas12a cuts ~18-23 bp from PAM
                    else:
                        cut_site = pam_pos - 3  # Cas9 cuts 3 bp upstream of PAM
                    
                    grna_info = {
                        'sequence': grna_seq,
                        'pam': match.group(),
                        'position': grna_start,
                        'strand': strand,
                        'cut_site': cut_site,
                        'full_sequence': grna_seq + match.group() if not self.cas_type.startswith('Cas12') else match.group() + grna_seq
                    }
                    
                    # Add only if not duplicate
                    if not any(g['sequence'] == grna_seq and g['position'] == grna_start for g in grnas):
                        grnas.append(grna_info)
                        
        return grnas
    
    def _pam_to_regex(self, pam: str) -> str:
        """Convert PAM sequence to regex pattern."""
        # IUPAC nucleotide codes
        iupac = {
            'N': '[ATCG]',
            'R': '[AG]',
            'Y': '[CT]',
            'S': '[GC]',
            'W': '[AT]',
            'K': '[GT]',
            'M': '[AC]',
            'B': '[CGT]',
            'D': '[AGT]',
            'H': '[ACT]',
            'V': '[ACG]'
        }
        
        regex = ''
        for char in pam:
            regex += iupac.get(char, char)
        return regex
    
    def filter_by_gc_content(self, grnas: List[Dict], min_gc: float = 40, max_gc: float = 60) -> List[Dict]:
        """
        Filter gRNAs by GC content.
        
        Args:
            grnas: List of gRNA dictionaries
            min_gc: Minimum GC content (%)
            max_gc: Maximum GC content (%)
            
        Returns:
            Filtered list of gRNAs
        """
        filtered = []
        for grna in grnas:
            gc_count = grna['sequence'].count('G') + grna['sequence'].count('C')
            gc_percent = (gc_count / len(grna['sequence'])) * 100
            
            if min_gc <= gc_percent <= max_gc:
                grna['gc_content'] = round(gc_percent, 1)
                filtered.append(grna)
                
        return filtered
    
    def check_homopolymers(self, grna: str, max_length: int = 4) -> bool:
        """
        Check if gRNA contains homopolymer runs.
        
        Args:
            grna: gRNA sequence
            max_length: Maximum allowed homopolymer length
            
        Returns:
            True if homopolymer exceeds max_length
        """
        for base in 'ATCG':
            if base * (max_length + 1) in grna:
                return True
        return False
    
    def filter_homopolymers(self, grnas: List[Dict], max_length: int = 4) -> List[Dict]:
        """Filter out gRNAs with long homopolymer runs."""
        return [g for g in grnas if not self.check_homopolymers(g['sequence'], max_length)]
    
    def predict_secondary_structure(self, grna: str, threshold: float = -3.0) -> Tuple[float, bool]:
        """
        Predict secondary structure using simple energy calculation.
        
        Args:
            grna: gRNA sequence
            threshold: Energy threshold for filtering (kcal/mol)
            
        Returns:
            Tuple of (folding energy, passes threshold)
        """
        # Simplified energy calculation
        # In production, use RNAfold or similar
        gc_count = grna.count('G') + grna.count('C')
        at_count = grna.count('A') + grna.count('T')
        
        # Simple approximation
        energy = -2.0 * gc_count - 1.0 * at_count
        energy += self._count_gc_runs(grna) * 0.5  # Penalty for GC runs
        
        return energy, energy > threshold
    
    def _count_gc_runs(self, sequence: str) -> int:
        """Count runs of GC bases."""
        runs = 0
        in_run = False
        
        for base in sequence:
            if base in 'GC':
                if not in_run:
                    runs += 1
                    in_run = True
            else:
                in_run = False
                
        return runs
    
    def filter_by_secondary_structure(self, grnas: List[Dict], threshold: float = -3.0) -> List[Dict]:
        """Filter gRNAs by predicted secondary structure."""
        filtered = []
        
        for grna in grnas:
            energy, passes = self.predict_secondary_structure(grna['sequence'], threshold)
            if passes:
                grna['folding_energy'] = round(energy, 2)
                filtered.append(grna)
                
        return filtered
    
    def check_seed_region(self, grna: str) -> Dict[str, any]:
        """
        Analyze the seed region (positions 1-10 from PAM).
        
        Args:
            grna: gRNA sequence
            
        Returns:
            Dictionary with seed region properties
        """
        # Seed region is the PAM-proximal 8-12 nt
        seed_region = grna[-12:]  # Last 12 nt (closest to PAM)
        
        seed_gc = (seed_region.count('G') + seed_region.count('C')) / len(seed_region) * 100
        
        # Check for problematic patterns in seed
        has_poly_t = 'TTTT' in seed_region
        has_poly_g = 'GGGG' in seed_region
        
        return {
            'seed_gc_content': round(seed_gc, 1),
            'has_poly_t': has_poly_t,
            'has_poly_g': has_poly_g,
            'seed_sequence': seed_region
        }
    
    def comprehensive_filter(self, grnas: List[Dict], filters: Dict = None) -> List[Dict]:
        """
        Apply comprehensive filtering to gRNAs.
        
        Args:
            grnas: List of gRNA dictionaries
            filters: Dictionary of filter parameters
            
        Returns:
            Filtered list of gRNAs
        """
        if filters is None:
            filters = {
                'gc_min': 40,
                'gc_max': 60,
                'homopolymer_max': 4,
                'secondary_structure_threshold': -3.0,
                'remove_poly_t_seed': True
            }
        
        # Apply GC content filter
        if 'gc_min' in filters and 'gc_max' in filters:
            grnas = self.filter_by_gc_content(grnas, filters['gc_min'], filters['gc_max'])
        
        # Apply homopolymer filter
        if 'homopolymer_max' in filters:
            grnas = self.filter_homopolymers(grnas, filters['homopolymer_max'])
        
        # Apply secondary structure filter
        if 'secondary_structure_threshold' in filters:
            grnas = self.filter_by_secondary_structure(grnas, filters['secondary_structure_threshold'])
        
        # Apply seed region filter
        if filters.get('remove_poly_t_seed', False):
            filtered = []
            for grna in grnas:
                seed_info = self.check_seed_region(grna['sequence'])
                if not seed_info['has_poly_t']:
                    grna['seed_info'] = seed_info
                    filtered.append(grna)
            grnas = filtered
        
        return grnas 