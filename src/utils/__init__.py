"""
Utils package for CasPINS - Cas-Primer-Indel Suite

Includes:
- AB1 file parsing
- Indel analysis algorithms
- Primer design utilities
- Visualization (standard and TIDE-style)
- File management
"""

from.ab1_parser import parse_ab1
from.indel_analysis import decompose_traces_indel_analysis, calculate_editing_efficiency_fallback
from.primer_design import PrimerDesigner, generate_primer_recommendations, generate_enhanced_primer_report
from.visualization import *
from.visualization_multi import plot_indel_analysis_multi, create_summary_report
from.tide_visualization import (
    create_tide_style_figure,
    display_results_streamlit,
    display_indel_spectrum_streamlit,
    figure_to_base64,
    TIDE_COLORS
)
from.file_management import ensure_output_dir, archive_output_files

__all__ = [
    # AB1 parsing
    'parse_ab1',
    
    # Indel analysis
    'decompose_traces_indel_analysis',
    'calculate_editing_efficiency_fallback',
    
    # Primer design
    'PrimerDesigner',
    'generate_primer_recommendations',
    'generate_enhanced_primer_report',
    
    # Visualization
    'plot_indel_analysis_multi',
    'create_summary_report',
    
    # TIDE-style visualization
    'create_tide_style_figure',
    'display_results_streamlit',
    'display_indel_spectrum_streamlit',
    'figure_to_base64',
    'TIDE_COLORS',
    
    # File management
    'ensure_output_dir',
    'archive_output_files',
] 