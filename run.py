#!/usr/bin/env python
"""
CasPINS - Cas-Primer-Indel Suite - Main Entry Point
=========================================

A unified command-line interface for all CRISPR analysis tools.

Usage:
    python run.py gui                          # Launch GUI
    python run.py grna GENE_NAME [options]     # Find gRNAs
    python run.py primers GENE_NAME [options]  # Design primers
    python run.py analysis [options]           # Run indel analysis

Examples:
    python run.py gui
    python run.py grna TP53 --species human --top 10
    python run.py primers TP53
    python run.py analysis --data-dir./data
"""

import sys
import os
import subprocess

# Add src directory to path
_root_dir = os.path.dirname(os.path.abspath(__file__))
_src_dir = os.path.join(_root_dir, 'src')
sys.path.insert(0, _src_dir)


def print_help():
    """Print usage information."""
    help_text = """
CasPINS - Cas-Primer-Indel Suite - Unified Command Interface
==================================================

Usage: python run.py <command> [arguments]

Commands:
  gui                Launch the graphical user interface (recommended for new users)
  grna GENE          Find optimal gRNAs for a target gene
  primers GENE       Design PCR and sequencing primers
  analysis           Run indel analysis on AB1 files

Examples:
  python run.py gui                                    # Start GUI
  python run.py grna TP53 --species human --top 10    # Find gRNAs
  python run.py primers TP53                          # Design primers  
  python run.py analysis --data-dir /path/to/data     # Run analysis

For detailed help on each command:
  python run.py grna --help
  python run.py primers --help
  python run.py analysis --help

Documentation:
  See README.md and docs/ folder for detailed guides.
"""
    print(help_text)


def run_gui():
    """Launch the Streamlit GUI."""
    print("Starting CasPINS - Cas-Primer-Indel Suite GUI...")
    print("=" * 60)
    print("The GUI will open in your default web browser.")
    print("Press Ctrl+C to stop the server.")
    print("=" * 60)
    
    gui_path = os.path.join(_src_dir, 'gui', 'app.py')
    subprocess.run([sys.executable, "-m", "streamlit", "run", gui_path])


def run_grna(args):
    """Run gRNA finder."""
    from cli.find_grna import main
    sys.argv = ['find_grna'] + args
    return main()


def run_primers(args):
    """Run primer designer."""
    from cli.design_primers import main
    sys.argv = ['design_primers'] + args
    return main()


def run_analysis(args):
    """Run indel analysis."""
    from cli.run_analysis import main
    sys.argv = ['run_analysis'] + args
    return main()


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print_help()
        return 0
    
    command = sys.argv[1].lower()
    args = sys.argv[2:]
    
    if command in ['help', '-h', '--help']:
        print_help()
        return 0
    
    elif command == 'gui':
        run_gui()
        return 0
    
    elif command == 'grna':
        return run_grna(args)
    
    elif command == 'primers':
        return run_primers(args)
    
    elif command == 'analysis':
        return run_analysis(args)
    
    else:
        print(f"Unknown command: {command}")
        print("Run 'python run.py --help' for usage information.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
