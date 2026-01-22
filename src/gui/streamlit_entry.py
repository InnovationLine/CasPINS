#!/usr/bin/env python
"""
Streamlit Cloud Entry Point for CRISPR Analysis Suite

This file is the entry point for Streamlit Cloud deployment.
It imports and runs the main GUI application.

Deployment Instructions:
1. Push this repository to GitHub
2. Go to https://share.streamlit.io
3. Connect your GitHub account
4. Select this repository
5. Set the main file path to: src/gui/streamlit_entry.py
6. Deploy!

The app will be available at:
https://[your-app-name].streamlit.app
"""

import sys
import os

# Add parent directories to path
_current_dir = os.path.dirname(os.path.abspath(__file__))
_src_dir = os.path.dirname(_current_dir)
_root_dir = os.path.dirname(_src_dir)

# Add paths
sys.path.insert(0, _src_dir)
sys.path.insert(0, _root_dir)

# Import and run the main GUI
from gui.app import main

if __name__ == "__main__":
    main()
