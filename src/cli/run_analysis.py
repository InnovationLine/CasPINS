#!/usr/bin/env python
"""
Run CRISPR Indel Analysis Pipeline for Multiple Samples

This is the CLI entry point for analyzing CRISPR-edited samples.
It processes all genes with input folders automatically.

Usage:
    python -m cli.run_analysis [options]
    
Examples:
    python -m cli.run_analysis                          # Analyze all genes
    python -m cli.run_analysis --data-dir /path/to/data # Specify data directory
    python -m cli.run_analysis --skip-genes gene1 gene2 # Skip specific genes
"""

import sys
import os

# Handle imports whether run as module or directly
if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from indel_analysis_multi import main as run_analysis_main


def main():
    """Entry point for indel analysis."""
    print("\n" + "="*60)
    print("Starting CRISPR Indel Analysis Pipeline (Multi-Sample Mode)...")
    print("="*60)
    
    try:
        run_analysis_main()
    except Exception as e:
        print(f"\nError running analysis: {str(e)}")
        print("\nTroubleshooting:")
        print("1. Make sure you have activated the virtual environment")
        print("2. Check that all required files are in place")
        print("3. Ensure input AB1 files are in data/<gene>/input/")
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
