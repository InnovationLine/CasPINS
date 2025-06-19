# Project Structure Guide

## Overview

This project follows a standard GitHub repository structure for better organization and maintainability.

## Directory Structure

```
crisper_analysis/
│
├── src/                          # Source code
│   └── tide_batch_analysis.py    # Main batch processing pipeline
│
├── scripts/                      # Utility scripts
│   ├── tide_analysis.py          # Original individual analysis
│   ├── tide_analysis_demo.py     # Demo/test script
│   └── tide_analysis_pipeline.py # Interactive pipeline version
│
├── data/                         # Data directory with gene folders
│   ├── vmat1/
│   │   ├── control.ab1           # Wild-type sequencing
│   │   ├── edited.ab1            # CRISPR-edited sequencing
│   │   ├── grna.txt              # Guide RNA sequences
│   │   ├── mrna.txt              # Full mRNA sequence (REQUIRED)
│   │   └── output/               # Analysis outputs folder
│   │       ├── tide_analysis_vmat1_20250618_201234.png
│   │       └── recommendations_vmat1_20250618_201234.txt
│   ├── vmat2/                    # (to be added)
│   │   └── output/               # (created automatically)
│   └── ddc/                      # (to be added)
│       └── output/               # (created automatically)
│
├── results/                      # General analysis outputs
│   └── tide_analysis_demo.png    # Demo analysis result
│
├── docs/                         # Documentation
│   └── project_structure.md      # This file
│
├── crispr_env/                   # Python virtual environment
├── gRNA/                         # Original gRNA data (to be removed)
├── Results_1267065/              # Original sequencing data (to be removed)
├── 2025-06-17_14C_1267196/       # Original sequencing data (to be removed)
│
├── run_analysis.py               # Main entry point
├── requirements.txt              # Python dependencies
├── README.md                     # Project documentation
└── .gitignore                    # Git ignore rules
```

## File Naming Standards

### Required files in each gene folder (e.g., `data/vmat1/`):
- `control.ab1` - Control/wild-type sample chromatogram
- `edited.ab1` - CRISPR-edited sample chromatogram
- `grna.txt` - Guide RNA sequences (one per line)
- `mrna.txt` - Full mRNA sequence (REQUIRED for primer design)

### Generated output files in `output/` subfolder:
- `tide_analysis_[gene]_[YYYYMMDD_HHMMSS].png` - Analysis plot
- `recommendations_[gene]_[YYYYMMDD_HHMMSS].txt` - Primer recommendations

The timestamp format (YYYYMMDD_HHMMSS) ensures that multiple analyses can be run without overwriting previous results.

## Usage

From the project root directory:
```bash
python run_analysis.py
```

This will:
1. Process all gene folders in `data/`
2. Create an `output/` subfolder in each gene directory (if it doesn't exist)
3. Generate timestamped results in the output folder

## Adding New Genes

1. Create a new folder: `data/[gene_name]/`
2. Add ALL required files:
   - `control.ab1`
   - `edited.ab1`
   - `grna.txt`
   - `mrna.txt` (REQUIRED)
3. Run `python run_analysis.py`
4. Find results in `data/[gene_name]/output/`

## Best Practices

1. Keep all sample data directly in `data/[gene_name]/` folders
2. Use standardized filenames (control.ab1, edited.ab1, grna.txt, mrna.txt)
3. Store guide RNAs in `grna.txt` (one sequence per line)
4. Always include `mrna.txt` with the full mRNA sequence
5. Generated outputs are automatically organized in the `output/` subfolder
6. Use `results/` for cross-gene analyses or summary reports

## Why Output Subfolders?

The `output/` subfolder organization provides:
- **Clean separation**: Input files stay separate from generated outputs
- **Easy cleanup**: You can delete the entire output folder to remove all results
- **Better organization**: All analysis results are in one place
- **Version control friendly**: You can easily exclude output folders from git

## Why Timestamps?

Timestamped output files allow you to:
- Run multiple analyses with different parameters
- Keep a history of your analysis results
- Compare results from different runs
- Never accidentally overwrite important results 