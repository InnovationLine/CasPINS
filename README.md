# CasPINS - Cas-Primer-Indel Suite

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18370068.svg)](https://doi.org/10.5281/zenodo.18370068)
[![CI](https://github.com/InnovationLine/CasPINS/actions/workflows/ci.yml/badge.svg)](https://github.com/InnovationLine/CasPINS/actions/workflows/ci.yml)

**CasPINS** is a comprehensive, integrated platform for CRISPR/TALEN gRNA design, primer design, and indel analysis with both **GUI** and **command-line** interfaces.

## 🎯 Overview

CasPINS provides an integrated workflow for genome editing experiments:
1. **gRNA Design**: Find optimal guide RNAs for your target genes
2. **Primer Design**: Generate primers for validation experiments  
3. **Indel Analysis**: Quantify editing efficiency from AB1 sequencing files

The tool is **completely generic** - it works with any gene from any species. Simply provide your data in the correct folder structure and the tool will process it automatically.

## 🚀 Quick Start

### Option 1: Graphical User Interface (GUI) - Recommended

```bash
# Install dependencies
pip install -r requirements.txt

# Launch the GUI
python run.py gui

# Or on Windows: double-click run_gui.bat
```

The GUI opens in your web browser and provides an intuitive interface for all features.

### Option 2: Command Line Interface (CLI)

```bash
# Setup
pip install -r requirements.txt

# View all available commands
python run.py --help

# Find gRNAs for your gene
python run.py grna TP53 --species human --top 10

# Design validation primers
python run.py primers TP53

# Run indel analysis
python run.py analysis --data-dir ./data
```

## 🖥️ GUI Features

The graphical interface provides:

- **🔍 gRNA Design Tab**
  - Search for optimal guide RNAs
  - Support for 90+ species and multiple Cas types
  - Interactive results table with filtering
  - One-click export to gene folders

- **🧪 Primer Design Tab**
  - Automatic primer generation for validation
  - PCR I and PCR II primer sets
  - Database integration (NCBI, Ensembl)
  - Direct export functionality

- **📊 Indel Analysis Tab**
  - Single sample or batch analysis modes
  - Drag-and-drop AB1 file upload
  - TIDE-style visualization
  - Real-time progress tracking

- **📚 Documentation Tab**
  - Built-in help and tutorials
  - Quick reference guides
  - Troubleshooting tips

## 📋 Features

### Module 1: gRNA Design
- **Multi-species support**: 90+ species including human, mouse, rat, zebrafish
- **Multiple Cas variants**: SpCas9, SaCas9, Cas12a, SpCas9-NG, and more
- **Advanced scoring**: Doench 2016, Moreno-Mateos, and Xu algorithms
- **Comprehensive filtering**: GC content, homopolymers, off-targets

### Module 2: Primer Design
- **Dual primer sets**: PCR I for genomic DNA, PCR II for sequencing
- **Automated design**: Using Primer3 with optimized parameters
- **CRISPR-aware**: Primers flank cut sites appropriately
- **Database integration**: NCBI, Ensembl sequence validation

### Module 3: Indel Analysis
- **Trace Decomposition**: Quantify editing efficiency with NNLS algorithm
- **Multi-sample support**: Analyze multiple clones simultaneously
- **Visual reports**: Publication-ready plots and summaries
- **Quality metrics**: Signal quality and confidence scoring

## 📁 Project Structure

```
CasPINS/
├── run.py                   # Main unified entry point
├── setup.py                 # Package installation
├── requirements.txt         # Python dependencies
├── README.md
├── LICENSE
├── docs/                    # Documentation
│   ├── workflow_guide.md
│   ├── grna_finding_guide.md
│   └── ...
├── tests/                   # Test suite
└── src/                     # Source code
    ├── cli/                 # Command-line tools
    │   ├── find_grna.py
    │   ├── design_primers.py
    │   └── run_analysis.py
    ├── gui/                 # GUI application
    │   ├── app.py
    │   └── streamlit_entry.py
    ├── grna_design/         # gRNA design system
    ├── utils/               # Analysis utilities
    └── config/              # Configuration
```

## 📂 Data Directory Structure

The tool works with **any gene** - just organize your data like this:

```
your_data_directory/         # Set via GUI Settings or CRISPR_DATA_DIR env var
├── gene_name/               # Any gene (e.g., tp53, brca1, myod1)
│   ├── grna.txt             # gRNA sequences (one per line)
│   ├── mrna.txt             # mRNA reference sequence
│   └── input/               # AB1 sequencing files
│       ├── control.ab1      # Control/wild-type sample
│       └── edited*.ab1      # Edited samples (any number)
└── another_gene/
    └── ...
```

## 🔧 Usage Examples

### Find gRNAs
```bash
# Find gRNAs for human TP53
python run.py grna TP53 --species human --top 10

# Save directly to data folder
python run.py grna TP53 --save-to data/tp53/grna.txt
```

### Design Primers
```bash
# Design primers for TP53
python run.py primers TP53 --data-dir ./data
```

### Run Indel Analysis
```bash
# Analyze all genes with input data
python run.py analysis --data-dir ./data

# The pipeline automatically:
# 1. Finds all gene folders with control.ab1 and edited*.ab1 files
# 2. Runs trace decomposition analysis
# 3. Generates plots and JSON results
```

## 📊 Output Files

Each analysis generates:
- `indel_analysis_gene_sample_timestamp.png` - Individual sample analysis plot
- `indel_analysis_summary_gene_timestamp.png` - Summary of all samples
- `indel_analysis_gene_timestamp.json` - Detailed results data
- `primer_recommendations.txt` - PCR and sequencing primers

## 🧬 Adding New Genes

```bash
# 1. Create gene folder in your data directory
mkdir data/your_gene

# 2. Add gRNA sequence
echo "ACGTACGTACGTACGTACGT" > data/your_gene/grna.txt

# 3. Add mRNA sequence (from NCBI/Ensembl)
# Save to: data/your_gene/mrna.txt

# 4. Add your AB1 sequencing files
# data/your_gene/input/control.ab1
# data/your_gene/input/edited*.ab1

# 5. Run analysis
python run.py analysis --data-dir ./data
```

## 📖 Documentation

- [User Guide](docs/USER_GUIDE.md) - Complete usage instructions and workflows
- [Developer Guide](docs/DEVELOPER_GUIDE.md) - Technical documentation and architecture
- [Deployment Guide](DEPLOYMENT_GUIDE.md) - Installation and deployment options

## ⚖️ Algorithm Details

Our indel analysis uses **Non-Negative Least Squares (NNLS)** decomposition:
- Open-source, public domain algorithm
- Decomposes edited traces into shifted control traces
- Suitable for clonal cell lines (50-100% efficiency)
- TIDE-style visualization

## 📦 Installation Options

### Local Installation (Recommended)
```bash
# Clone the repository
git clone https://github.com/InnovationLine/CasPINS.git
cd CasPINS

# Create virtual environment (optional but recommended)
python -m venv caspins_env
source caspins_env/bin/activate  # Linux/Mac
# or: caspins_env\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Run GUI
python run.py gui
```

### PyPI Installation (Planned)

> **Note:** PyPI package publication is planned for a future release.

```bash
# Coming soon:
pip install caspins
caspins-gui
```

## 📚 Citation

If you use this software in your research, please cite:

```bibtex
@software{caspins,
  author = {Dasgupta, Rinki and Das, Kaushik},
  title = {CasPINS: An Integrated Platform for CRISPR/TALEN gRNA Design, Primer Generation, and Indel Analysis},
  year = {2026},
  publisher = {Zenodo},
  url = {https://github.com/InnovationLine/CasPINS},
  doi = {10.5281/zenodo.18370068}
}
```

## 🤝 Contributing

Contributions are welcome! Please read our [Contributing Guidelines](CONTRIBUTING.md) before submitting issues or pull requests.

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## 🙏 Acknowledgments

- BioPython for sequence handling
- Primer3 for primer design
- NumPy/SciPy for signal processing
- Streamlit for the GUI framework
- The CRISPR community for continued innovation

## 📧 Contact

For questions or support, please open an issue on GitHub.
