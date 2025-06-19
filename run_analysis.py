#!/usr/bin/env python
"""
Run CRISPR TIDE Analysis Pipeline

This script runs the batch analysis on all samples in data/samples/
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Import and run the main pipeline
from tide_batch_analysis import main

if __name__ == '__main__':
    print("Starting CRISPR TIDE Analysis Pipeline...")
    print("-" * 60)
    main() 