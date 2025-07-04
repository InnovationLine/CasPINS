#!/usr/bin/env python
"""
Run CRISPR TIDE Analysis Pipeline for Multiple Samples

This script runs the batch analysis on all gene folders with input directories.
Each input folder should contain:
- control.ab1: Control sample
- edited*.ab1: Multiple edited samples (clonally expanded cell lines)
"""

import sys
import os
import gc

# Configure matplotlib for non-interactive mode (prevents hanging)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Import and run the main pipeline
from tide_batch_analysis_multi import main

if __name__ == '__main__':
    print("Starting CRISPR TIDE Analysis Pipeline (Multi-Sample Mode)...")
    print("Processing clonally expanded CRISPR edited cell lines")
    print("-" * 60)
    
    try:
        main()
    finally:
        # Clean up resources
        plt.close('all')
        gc.collect() 