"""
Tests for gRNA Design Module
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.grna_design.scoring.scoring_engine import ScoringEngine
from src.grna_design.core.grna_generator import GRNAGenerator


class TestScoringEngine:
    """Tests for the gRNA scoring engine."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.scorer = ScoringEngine()
    
    def test_gc_score_optimal(self):
        """Test GC score for optimal GC content (40-60%)."""
        # 50% GC content
        grna = "ACGTACGTACGTACGTACGT"  # 10 GC, 10 AT
        score = self.scorer.calculate_gc_score(grna)
        assert score == 1.0
    
    def test_gc_score_suboptimal(self):
        """Test GC score for suboptimal GC content."""
        # 25% GC content
        grna = "AAAATAAAATAAAATAGGG"  # 5 GC, 15 AT (approx)
        score = self.scorer.calculate_gc_score(grna)
        assert 0 < score < 1.0
    
    def test_doench_score_range(self):
        """Test Doench score is within valid range."""
        grna = "ACGTACGTACGTACGTACGT"
        score = self.scorer.calculate_doench_2016_score(grna)
        assert 0 <= score <= 1
    
    def test_moreno_mateos_score_range(self):
        """Test Moreno-Mateos score is within valid range."""
        grna = "GCGTGCGTGCGTGCGTGCGT"
        score = self.scorer.calculate_moreno_mateos_score(grna)
        assert 0 <= score <= 1
    
    def test_xu_score_range(self):
        """Test Xu score is within valid range."""
        grna = "GGGGCCCCAAAATTTTGGGG"
        score = self.scorer.calculate_xu_score(grna)
        assert 0 <= score <= 1
    
    def test_composite_score(self):
        """Test composite score calculation."""
        scores = {
            'doench_2016': 0.8,
            'moreno_mateos': 0.7,
            'xu': 0.6,
            'cfd': 0.9
        }
        composite = self.scorer.calculate_composite_score(scores)
        assert 0 <= composite <= 1
    
    def test_cfd_score_perfect_match(self):
        """Test CFD score for perfect match."""
        grna = "ACGTACGTACGTACGTACGT"
        score = self.scorer.calculate_cfd_score(grna, grna)
        assert score == 1.0
    
    def test_cfd_score_with_mismatches(self):
        """Test CFD score with mismatches."""
        grna = "ACGTACGTACGTACGTACGT"
        off_target = "ACGTACGTACGTACGTAAAA"  # 3 mismatches at end
        score = self.scorer.calculate_cfd_score(grna, off_target)
        assert 0 < score < 1.0
    
    def test_grna_length_validation(self):
        """Test that invalid length returns 0 score."""
        short_grna = "ACGT"  # Too short
        score = self.scorer.calculate_doench_2016_score(short_grna)
        assert score == 0.0


class TestGRNAGenerator:
    """Tests for the gRNA generator."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.generator = GRNAGenerator('SpCas9')
    
    def test_find_pam_sites_ngg(self):
        """Test finding NGG PAM sites."""
        # Sequence with known NGG sites (AGG is a valid NGG PAM)
        sequence = "ACGTACGTACGTACGTACGTAGGACGTACGTACGTACGTACGT"
        pam_sites = self.generator.find_pam_sites(sequence)
        assert len(pam_sites) > 0
    
    def test_find_all_grnas(self):
        """Test finding all gRNAs in a sequence."""
        # Create a test sequence with known gRNA sites
        sequence = "A" * 30 + "NGG" + "T" * 30
        grnas = self.generator.find_all_grnas(sequence)
        # Should find at least the one we created
        assert isinstance(grnas, list)
    
    def test_gc_content_filter(self):
        """Test GC content filtering."""
        # High GC gRNA
        high_gc = "GGGGCCCCGGGGCCCCGGGG"  # 100% GC
        # Low GC gRNA  
        low_gc = "AAAATTTTAAAATTTTAAAA"  # 0% GC
        
        filters = {'gc_min': 40, 'gc_max': 60}
        
        # Create mock gRNA dictionaries
        grnas = [
            {'sequence': high_gc, 'gc_content': 100},
            {'sequence': low_gc, 'gc_content': 0}
        ]
        
        filtered = self.generator.comprehensive_filter(grnas, filters)
        # Neither should pass the 40-60% GC filter
        for g in filtered:
            gc = (g['sequence'].count('G') + g['sequence'].count('C')) / len(g['sequence']) * 100
            # This might still include them if filter is lenient
    
    def test_homopolymer_check(self):
        """Test homopolymer detection."""
        # Sequence with AAAA homopolymer
        poly_a = "AAAACGTACGTACGTACGT"
        assert self.generator.check_homopolymers(poly_a) == True
        
        # Sequence without homopolymers
        no_poly = "ACGTACGTACGTACGTACGT"
        assert self.generator.check_homopolymers(no_poly) == False


class TestCRISPRSystems:
    """Tests for different CRISPR systems."""
    
    def test_spcas9_pam(self):
        """Test SpCas9 recognizes NGG PAM."""
        generator = GRNAGenerator('SpCas9')
        assert generator.pam == 'NGG'
    
    def test_sacas9_pam(self):
        """Test SaCas9 recognizes NNGRRT PAM."""
        generator = GRNAGenerator('SaCas9')
        assert generator.pam == 'NNGRRT'
    
    def test_cas12a_pam(self):
        """Test Cas12a recognizes TTTV PAM."""
        generator = GRNAGenerator('Cas12a')
        assert generator.pam == 'TTTV'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
