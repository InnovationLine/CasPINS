# CRISPR TIDE Analysis Pipeline

A comprehensive Python pipeline for analyzing CRISPR editing efficiency using TIDE (Tracking of Indels by Decomposition) analysis on Sanger sequencing data.

## 🚀 Quick Start

```bash
# Clone repository
git clone <repository-url>
cd crispr_analysis

# Setup environment (Windows)
python -m venv crispr_env
.\crispr_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run analysis
python run_analysis.py
```

## 📋 Features

- **Batch Processing**: Analyze multiple CRISPR samples simultaneously
- **TIDE Analysis**: Quantify editing efficiency with confidence scores
- **Smart Validation**: Automatic sequence quality checks with configurable thresholds
- **gRNA Detection**: Automated identification in mRNA sequences with mismatch tolerance
- **Primer Design**: Professional primer recommendations using Primer3
- **Modular Architecture**: Clean, maintainable code structure
- **Flexible Control**: Command-line arguments for customized analysis
- **Comprehensive Output**: Detailed reports, visualizations, and quality metrics
- **Archive System**: Automatic organization of previous analysis runs

## 📁 Project Structure

```
crispr_analysis/
├── data/                    # Input data directory
│   ├── vmat1/              # Gene-specific folders
│   ├── vmat2/              
│   └── ddc/                
├── src/                    # Source code
│   ├── tide_batch_analysis.py  # Main orchestrator
│   └── utils/              # Modular components
├── docs/                   # Documentation
│   ├── modular_architecture.md
│   └── project_structure.md
├── run_analysis.py         # Entry point
├── requirements.txt        # Dependencies
├── USAGE.md               # Comprehensive usage guide
└── README.md              # This file
```

## 📖 Documentation

- **[USAGE.md](USAGE.md)** - Comprehensive usage guide with examples
- **[docs/modular_architecture.md](docs/modular_architecture.md)** - Technical architecture details
- **[docs/project_structure.md](docs/project_structure.md)** - Detailed project organization

## 🔧 Requirements

- Python 3.7+
- Windows/Linux/Mac OS
- Dependencies listed in `requirements.txt`:
  - biopython
  - numpy
  - matplotlib
  - scipy
  - pandas
  - primer3-py

## 💡 Basic Usage

### Standard Analysis
```bash
python run_analysis.py
```

### Force Plotting Despite Issues
```bash
python run_analysis.py --force-plot
```

### Custom Similarity Threshold
```bash
python run_analysis.py --similarity-threshold 50
```

### Process Specific Genes
```bash
python run_analysis.py --genes vmat1 ddc
```

See [USAGE.md](USAGE.md) for comprehensive examples and options.

## 📊 Input File Requirements

Each gene folder must contain:
- `control.ab1` - Control sample chromatogram
- `edited.ab1` - Edited sample chromatogram
- `grna.txt` - Guide RNA sequences (one per line)
- `mrna.txt` - mRNA reference sequence

## 📈 Output Files

- **ab1_analysis_gene.txt** - Detailed sequence quality analysis
- **recommendations_gene_timestamp.txt** - TIDE results and primer recommendations
- **tide_analysis_gene_timestamp.png** - 4-panel visualization plot

## 🛠️ Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--force-plot` | Generate plots even with validation failures | False |
| `--similarity-threshold` | Minimum sequence similarity (%) | 70 |
| `--skip-validation` | Skip all validation checks | False |
| `--interactive` | Ask before processing each gene | False |
| `--no-archive` | Don't archive previous outputs | False |
| `--genes` | Specific genes to analyze | all |

## 🐛 Troubleshooting

Common issues:
1. **Low sequence similarity** - Use `--force-plot` to investigate
2. **gRNA not found** - Check for PAM sequences or length issues
3. **Cut site out of bounds** - Need longer sequencing reads

See [USAGE.md](USAGE.md#troubleshooting) for detailed solutions.

## 📝 License

This project is for research use only.

## 👥 Contributors

[Add contributors here]

## 📧 Contact

For questions or issues, please contact the project maintainer.

---

For detailed usage instructions, see [USAGE.md](USAGE.md) 