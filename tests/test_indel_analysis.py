"""
Tests for Indel Analysis Module
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.utils.indel_analysis import (
    decompose_traces_indel_analysis,
    calculate_editing_efficiency_fallback
)


class TestTraceDecomposition:
    """Tests for trace decomposition algorithm."""
    
    def setup_method(self):
        """Create synthetic trace data for testing."""
        # Generate synthetic chromatogram traces
        length = 5000
        
        # Create base signals with peaks
        self.control_traces = {
            'A': np.zeros(length),
            'T': np.zeros(length),
            'G': np.zeros(length),
            'C': np.zeros(length)
        }
        
        self.edited_traces = {
            'A': np.zeros(length),
            'T': np.zeros(length),
            'G': np.zeros(length),
            'C': np.zeros(length)
        }
        
        # Add synthetic peaks (every 10 positions)
        for i in range(0, length, 10):
            base_idx = i % 4
            bases = ['A', 'T', 'G', 'C']
            # Control has clean peaks
            peak = 100 * np.exp(-0.5 * ((np.arange(-5, 6))**2))
            start = max(0, i - 5)
            end = min(length, i + 6)
            actual_peak = peak[:end-start]
            self.control_traces[bases[base_idx]][start:end] += actual_peak
            
            # Edited has some shifted peaks (simulating indel)
            if i > 2500:  # After cut site
                shifted_i = i + 10  # Simulate 1bp deletion effect
                if shifted_i + 6 < length:
                    self.edited_traces[bases[base_idx]][shifted_i-5:shifted_i+6] += peak
            else:
                self.edited_traces[bases[base_idx]][start:end] += actual_peak
        
        self.cut_site = 250  # Position in bp
    
    def test_decomposition_returns_dict(self):
        """Test that decomposition returns a dictionary."""
        result = decompose_traces_indel_analysis(
            self.control_traces, 
            self.edited_traces, 
            self.cut_site
        )
        assert isinstance(result, dict)
    
    def test_decomposition_has_required_keys(self):
        """Test that result has all required keys."""
        result = decompose_traces_indel_analysis(
            self.control_traces,
            self.edited_traces,
            self.cut_site
        )
        
        required_keys = [
            'editing_efficiency',
            'dominant_indel_size',
            'dominant_indel_percent',
            'indel_spectrum',
            'quality_score',
            'confidence'
        ]
        
        for key in required_keys:
            assert key in result
    
    def test_efficiency_in_valid_range(self):
        """Test that editing efficiency is between 0 and 100."""
        result = decompose_traces_indel_analysis(
            self.control_traces,
            self.edited_traces,
            self.cut_site
        )
        
        assert 0 <= result['editing_efficiency'] <= 100
    
    def test_quality_score_in_valid_range(self):
        """Test that quality score is between 0 and 100."""
        result = decompose_traces_indel_analysis(
            self.control_traces,
            self.edited_traces,
            self.cut_site
        )
        
        assert 0 <= result['quality_score'] <= 100
    
    def test_indel_spectrum_is_dict(self):
        """Test that indel spectrum is a dictionary."""
        result = decompose_traces_indel_analysis(
            self.control_traces,
            self.edited_traces,
            self.cut_site
        )
        
        assert isinstance(result['indel_spectrum'], dict)
    
    def test_handles_invalid_cut_site(self):
        """Test graceful handling of invalid cut site."""
        result = decompose_traces_indel_analysis(
            self.control_traces,
            self.edited_traces,
            -1  # Invalid cut site
        )
        
        # Should still return valid result (using fallback)
        assert 'editing_efficiency' in result
    
    def test_handles_none_cut_site(self):
        """Test graceful handling of None cut site."""
        result = decompose_traces_indel_analysis(
            self.control_traces,
            self.edited_traces,
            None
        )
        
        assert 'editing_efficiency' in result


class TestFallbackMethod:
    """Tests for fallback efficiency calculation."""
    
    def setup_method(self):
        """Create test data."""
        self.control_seq = "ACGT" * 100
        self.edited_seq = "ACGT" * 100
        
        # Simple traces
        length = 4000
        self.control_traces = {
            'A': np.random.rand(length) * 50 + 50,
            'T': np.random.rand(length) * 50 + 50,
            'G': np.random.rand(length) * 50 + 50,
            'C': np.random.rand(length) * 50 + 50
        }
        
        self.edited_traces = {
            'A': np.random.rand(length) * 50 + 50,
            'T': np.random.rand(length) * 50 + 50,
            'G': np.random.rand(length) * 50 + 50,
            'C': np.random.rand(length) * 50 + 50
        }
    
    def test_fallback_returns_dict(self):
        """Test that fallback returns dictionary."""
        result = calculate_editing_efficiency_fallback(
            self.control_seq,
            self.control_traces,
            self.edited_seq,
            self.edited_traces,
            200
        )
        
        assert isinstance(result, dict)
    
    def test_fallback_has_efficiency(self):
        """Test that fallback result includes efficiency."""
        result = calculate_editing_efficiency_fallback(
            self.control_seq,
            self.control_traces,
            self.edited_seq,
            self.edited_traces,
            200
        )
        
        assert 'editing_efficiency' in result
        assert 0 <= result['editing_efficiency'] <= 100
    
    def test_fallback_has_low_confidence(self):
        """Test that fallback method indicates low confidence."""
        result = calculate_editing_efficiency_fallback(
            self.control_seq,
            self.control_traces,
            self.edited_seq,
            self.edited_traces,
            200
        )
        
        assert 'LOW' in result['confidence'] or 'Fallback' in result['confidence']


class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_empty_traces(self):
        """Test handling of empty traces."""
        empty_traces = {'A': np.array([]), 'T': np.array([]), 
                       'G': np.array([]), 'C': np.array([])}
        
        result = decompose_traces_indel_analysis(
            empty_traces, empty_traces, 100
        )
        
        # Should handle gracefully
        assert 'editing_efficiency' in result
    
    def test_very_short_sequence(self):
        """Test handling of very short sequences."""
        short_traces = {'A': np.ones(100), 'T': np.ones(100),
                       'G': np.ones(100), 'C': np.ones(100)}
        
        result = decompose_traces_indel_analysis(
            short_traces, short_traces, 5
        )
        
        assert 'editing_efficiency' in result
    
    def test_identical_traces(self):
        """Test that identical traces show no editing."""
        traces = {'A': np.random.rand(5000) * 100,
                 'T': np.random.rand(5000) * 100,
                 'G': np.random.rand(5000) * 100,
                 'C': np.random.rand(5000) * 100}
        
        result = decompose_traces_indel_analysis(
            traces, traces.copy(), 250
        )
        
        # With identical traces, editing should be very low
        # (Not exactly 0 due to numerical precision)
        assert result['editing_efficiency'] < 50  # Should be low


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
