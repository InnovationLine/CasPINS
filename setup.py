"""
Setup script for CasPINS (Cas-Primer-Indel Suite)

This enables installation via pip and creates command-line entry points.
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="caspins",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="CasPINS: An integrated platform for CRISPR/TALEN gRNA design, primer design, and indel analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/raju1stnov/CasPINS",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements + ["streamlit>=1.28.0"],
    entry_points={
        "console_scripts": [
            # GUI entry point
            "caspins-gui=caspins_gui:cli_entry",
            # CLI entry points
            "caspins-grna=find_grna:main",
            "caspins-primers=design_primers:main",
            "caspins-analyze=run_analysis:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.txt", "*.md", "*.json", "*.bat"],
        "src": ["**/*.py"],
        "docs": ["*.md"],
    },
    zip_safe=False,
) 