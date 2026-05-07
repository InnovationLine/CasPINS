"""
gRNA Design System
A comprehensive tool for CRISPR guide RNA design with advanced scoring and off-target prediction

Supports:
- Multiple input types: Gene symbols, Ensembl IDs, RefSeq IDs, Genomic coordinates
- Multiple editing modes: Knockout, Knock-in (HDR), CRISPRa, CRISPRi, Base editing
- TALEN design
- 90+ species
"""

__version__ = "2.0.0"
__author__ = "CRISPR Analysis Team"

from.core.grna_generator import GRNAGenerator
from.core.sequence_analyzer import SequenceAnalyzer
from.scoring.scoring_engine import ScoringEngine
from.database.genome_manager import GenomeManager
from.grna_designer import GRNADesigner, EditingMode, TargetType
from.hdr_designer import HDRDesigner, HDRTemplate, EditType, design_knockin_template
from.talen_designer import TALENDesigner, TALENPair, find_talen_pairs

__all__ = [
    # Core modules
    'GRNAGenerator',
    'SequenceAnalyzer',
    'ScoringEngine',
    'GenomeManager',
    'GRNADesigner',
    
    # Editing modes and target types
    'EditingMode',
    'TargetType',
    
    # HDR Designer
    'HDRDesigner',
    'HDRTemplate',
    'EditType',
    'design_knockin_template',
    
    # TALEN Designer
    'TALENDesigner',
    'TALENPair',
    'find_talen_pairs',
] 