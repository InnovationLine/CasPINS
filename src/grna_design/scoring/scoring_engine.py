"""
Scoring Engine Module
Implements multiple scoring algorithms for gRNA effectiveness prediction
"""

import numpy as np
from typing import List, Dict, Tuple
import math


class ScoringEngine:
    """Calculate various scores for gRNA effectiveness."""
    
    def __init__(self):
        """Initialize scoring engine with algorithm weights."""
        self.default_weights = {
            'doench_2016': 0.4,
            'moreno_mateos': 0.3,
            'xu': 0.2,
            'cfd': 0.1
        }
        
    def calculate_all_scores(self, grna: str, context: str = None, pam: str = 'NGG') -> Dict[str, float]:
        """
        Calculate all available scores for a gRNA.
        
        Args:
            grna: 20nt guide RNA sequence
            context: Extended sequence context (30nt upstream + grna + pam + 30nt downstream)
            pam: PAM sequence
            
        Returns:
            Dictionary of scores
        """
        scores = {}
        
        # Calculate individual scores
        scores['doench_2016'] = self.calculate_doench_2016_score(grna, context)
        scores['moreno_mateos'] = self.calculate_moreno_mateos_score(grna)
        scores['xu'] = self.calculate_xu_score(grna, pam)
        scores['gc_content'] = self.calculate_gc_score(grna)
        
        # Calculate composite score
        scores['composite'] = self.calculate_composite_score(scores)
        
        return scores
    
    def calculate_doench_2016_score(self, grna: str, context: str = None) -> float:
        """
        Calculate Doench et al. 2016 on-target score.
        
        This is a simplified version. Full implementation would include:
        - Position-specific nucleotide features
        - Di-nucleotide features
        - GC content features
        - Melting temperature
        
        Args:
            grna: 20nt guide RNA sequence
            context: Extended context (4nt + 20nt grna + 3nt PAM + 3nt)
            
        Returns:
            Score between 0 and 1
        """
        if len(grna) != 20:
            return 0.0
            
        score = 0.5  # Base score
        
        # Position-specific scoring (simplified)
        # Preferred nucleotides at specific positions
        position_scores = {
            # Position: (preferred_base, weight)
            0: ('A', 0.05),
            2: ('C', 0.03),
            3: ('C', 0.04),
            7: ('G', 0.03),
            19: ('G', 0.05)
        }
        
        for pos, (base, weight) in position_scores.items():
            if grna[pos] == base:
                score += weight
                
        # GC content scoring
        gc_content = (grna.count('G') + grna.count('C')) / len(grna)
        if 0.4 <= gc_content <= 0.6:
            score += 0.1
        elif 0.3 <= gc_content <= 0.7:
            score += 0.05
            
        # Penalize poly-T
        if 'TTTT' in grna:
            score -= 0.2
            
        # Penalize poly-G in seed region
        seed = grna[-12:]
        if 'GGGG' in seed:
            score -= 0.15
            
        # Di-nucleotide features (simplified)
        good_dinucleotides = ['CC', 'CT', 'GG']
        bad_dinucleotides = ['GT', 'TC']
        
        for i in range(len(grna) - 1):
            dinuc = grna[i:i+2]
            if dinuc in good_dinucleotides:
                score += 0.02
            elif dinuc in bad_dinucleotides:
                score -= 0.02
                
        return max(0, min(1, score))
    
    def calculate_moreno_mateos_score(self, grna: str) -> float:
        """
        Calculate Moreno-Mateos et al. score (CRISPRscan).
        
        Simplified version focusing on:
        - Nucleotide preferences
        - GC content
        - Self-folding energy
        
        Args:
            grna: 20nt guide RNA sequence
            
        Returns:
            Score between 0 and 1
        """
        if len(grna) != 20:
            return 0.0
            
        score = 0.5
        
        # Nucleotide composition scores
        # Based on Moreno-Mateos supplementary data
        base_scores = {
            'G': {'weight': 0.8, 'positions': [4, 6, 10, 14, 16, 17, 18]},
            'C': {'weight': 0.7, 'positions': [3, 5, 7, 8]},
            'A': {'weight': 0.6, 'positions': [9, 12, 15]},
            'T': {'weight': -0.3, 'positions': [1, 2, 11, 13]}
        }
        
        for base, info in base_scores.items():
            for pos in info['positions']:
                if pos < len(grna) and grna[pos] == base:
                    score += info['weight'] * 0.01
                    
        # GC content optimization
        gc_content = (grna.count('G') + grna.count('C')) / len(grna)
        if 0.5 <= gc_content <= 0.65:
            score += 0.15
        elif 0.4 <= gc_content <= 0.7:
            score += 0.08
            
        # Penalize extreme GC
        if gc_content < 0.3 or gc_content > 0.8:
            score -= 0.2
            
        # Self-complementarity penalty (simplified)
        if self._has_self_complementarity(grna):
            score -= 0.1
            
        return max(0, min(1, score))
    
    def calculate_xu_score(self, grna: str, pam: str = 'NGG') -> float:
        """
        Calculate Xu et al. score.
        
        Simplified version based on:
        - PAM preference
        - Position-specific nucleotides
        - GC content
        
        Args:
            grna: 20nt guide RNA sequence
            pam: PAM sequence
            
        Returns:
            Score between 0 and 1
        """
        if len(grna) != 20:
            return 0.0
            
        score = 0.5
        
        # PAM scoring
        pam_scores = {
            'NGG': 1.0,
            'NAG': 0.7,
            'NCG': 0.5,
            'NGA': 0.5
        }
        
        for pam_pattern, pam_weight in pam_scores.items():
            if self._matches_pam(pam, pam_pattern):
                score *= pam_weight
                break
                
        # Position-specific features from Xu et al.
        # Positions are 0-indexed
        preferred_bases = {
            3: {'C': 0.05},
            5: {'G': 0.04},
            7: {'G': 0.03},
            9: {'C': 0.03},
            16: {'G': 0.05},
            17: {'G': 0.05}
        }
        
        for pos, base_weights in preferred_bases.items():
            for base, weight in base_weights.items():
                if grna[pos] == base:
                    score += weight
                    
        # GC content preference
        gc_content = (grna.count('G') + grna.count('C')) / len(grna)
        if 0.45 <= gc_content <= 0.55:
            score += 0.1
            
        return max(0, min(1, score))
    
    def calculate_gc_score(self, grna: str) -> float:
        """
        Calculate GC content score.
        
        Args:
            grna: Guide RNA sequence
            
        Returns:
            Score between 0 and 1 based on GC content
        """
        gc_content = (grna.count('G') + grna.count('C')) / len(grna)
        
        # Optimal GC content is 40-60%
        if 0.4 <= gc_content <= 0.6:
            return 1.0
        elif 0.3 <= gc_content <= 0.7:
            return 0.8
        elif 0.2 <= gc_content <= 0.8:
            return 0.5
        else:
            return 0.2
    
    def calculate_cfd_score(self, grna: str, off_target: str) -> float:
        """
        Calculate Cutting Frequency Determination (CFD) score.
        
        Args:
            grna: On-target sequence
            off_target: Off-target sequence
            
        Returns:
            CFD score between 0 and 1
        """
        if len(grna) != len(off_target):
            return 0.0
            
        # Position-specific mismatch penalties (simplified)
        # Real CFD uses experimentally determined values
        position_penalties = {
            # Distance from PAM: penalty
            1: 0.0,
            2: 0.0,
            3: 0.014,
            4: 0.0,
            5: 0.0,
            6: 0.395,
            7: 0.317,
            8: 0.0,
            9: 0.389,
            10: 0.079,
            11: 0.445,
            12: 0.508,
            13: 0.613,
            14: 0.851,
            15: 0.732,
            16: 0.828,
            17: 0.615,
            18: 0.804,
            19: 0.685,
            20: 0.583
        }
        
        score = 1.0
        
        for i in range(len(grna)):
            if grna[i] != off_target[i]:
                # Position from PAM (1-indexed)
                pos_from_pam = len(grna) - i
                penalty = position_penalties.get(pos_from_pam, 0.5)
                score *= (1 - penalty)
                
        return score
    
    def calculate_composite_score(self, scores: Dict[str, float], weights: Dict[str, float] = None) -> float:
        """
        Calculate weighted composite score.
        
        Args:
            scores: Dictionary of individual scores
            weights: Custom weights for each score
            
        Returns:
            Composite score between 0 and 1
        """
        if weights is None:
            weights = self.default_weights
            
        composite = 0.0
        total_weight = 0.0
        
        for score_type, weight in weights.items():
            if score_type in scores:
                composite += scores[score_type] * weight
                total_weight += weight
                
        if total_weight > 0:
            composite /= total_weight
            
        return min(1.0, max(0.0, composite))
    
    def _has_self_complementarity(self, sequence: str, min_length: int = 4) -> bool:
        """Check if sequence has self-complementarity."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rev_comp = ''.join(complement.get(base, 'N') for base in sequence[::-1])
        
        # Check for complementary regions
        for i in range(len(sequence) - min_length + 1):
            for j in range(len(rev_comp) - min_length + 1):
                if sequence[i:i+min_length] == rev_comp[j:j+min_length]:
                    return True
                    
        return False
    
    def _matches_pam(self, pam: str, pattern: str) -> bool:
        """Check if PAM matches pattern with IUPAC codes."""
        if len(pam) != len(pattern):
            return False
            
        iupac = {
            'N': 'ATCG',
            'R': 'AG',
            'Y': 'CT',
            'S': 'GC',
            'W': 'AT',
            'K': 'GT',
            'M': 'AC'
        }
        
        for p, t in zip(pam, pattern):
            if t in iupac:
                if p not in iupac[t]:
                    return False
            elif p != t:
                return False
                
        return True
    
    def rank_grnas(self, grnas: List[Dict], score_type: str = 'composite') -> List[Dict]:
        """
        Rank gRNAs by score.
        
        Args:
            grnas: List of gRNA dictionaries
            score_type: Which score to use for ranking
            
        Returns:
            Sorted list of gRNAs
        """
        # Calculate scores for all gRNAs
        for grna in grnas:
            if 'scores' not in grna:
                grna['scores'] = self.calculate_all_scores(
                    grna['sequence'],
                    grna.get('context'),
                    grna.get('pam', 'NGG')
                )
        
        # Sort by specified score
        return sorted(grnas, key=lambda x: x['scores'].get(score_type, 0), reverse=True) 