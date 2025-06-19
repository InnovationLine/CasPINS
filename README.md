# CRISPR TIDE Analysis Pipeline

A Python-based pipeline for analyzing CRISPR editing efficiency using TIDE (Tracking of Indels by Decomposition) methodology.

## Overview

This pipeline processes Sanger sequencing chromatogram files (.ab1) to visualize and analyze CRISPR editing outcomes. It compares control and edited samples to identify divergence points and generates primer recommendations for proper TIDE analysis based on the target gene's mRNA sequence.

## Project Structure

```
crisper_analysis/
├── src/
│   └── tide_batch_analysis.py    # Main analysis pipeline
├── scripts/
│   ├── tide_analysis.py          # Individual analysis script
│   ├── tide_analysis_demo.py     # Demo analysis script
│   └── tide_analysis_pipeline.py # Interactive pipeline
├── data/
│   ├── vmat1/
│   │   ├── control.ab1           # Control sample chromatogram
│   │   ├── edited.ab1            # Edited sample chromatogram
│   │   ├── grna.txt              # Guide RNA sequences
│   │   ├── mrna.txt              # mRNA sequence (REQUIRED)
│   │   └── output/               # All generated files go here
│   │       ├── tide_analysis_vmat1_20250618_212231.png
│   │       └── recommendations_vmat1_20250618_212231.txt
│   ├── vmat2/
│   └── ddc/
│       └── output/                   # Created automatically when processed
├── results/                      # Analysis outputs
├── docs/                         # Documentation
├── requirements.txt              # Python dependencies
└── README.md                     # This file
```

## Requirements

- Python 3.7+
- BioPython
- NumPy
- Matplotlib
- BeautifulSoup4
- Requests

Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### 1. Prepare Your Data

For each gene you want to analyze, create a folder in `data/` with ALL required files:
- `control.ab1` - Sanger sequencing file from wild-type sample
- `edited.ab1` - Sanger sequencing file from CRISPR-edited sample
- `grna.txt` - Text file with guide RNA sequences (one per line)
- `mrna.txt` - Full mRNA sequence for primer design (REQUIRED)

Example `grna.txt`:
```
ACGATAACACCCCCAGTAGA
TGTCTATGCCATCGCCGATG
```

### 2. Run the Analysis

```bash
python run_analysis.py
```

The pipeline will:
1. Process all gene folders in `data/`
2. Create an `output/` subfolder in each gene directory
3. Generate timestamped TIDE analysis plots in the output folder
4. Create timestamped primer recommendations in the output folder

### 3. Outputs

For each gene, the pipeline generates timestamped files in `data/[gene]/output/`:
- `tide_analysis_[gene]_[YYYYMMDD_HHMMSS].png` - Chromatogram comparison plot
- `recommendations_[gene]_[YYYYMMDD_HHMMSS].txt` - Sequencing primer recommendations

The timestamp format ensures you can run multiple analyses without overwriting previous results.

## Features

- **Batch Processing**: Analyze multiple genes in one run
- **Automatic gRNA Detection**: Finds guide RNA positions in mRNA sequences
- **Primer Design**: Generates sequencing primers optimized for TIDE analysis
- **Visual Analysis**: Three-panel plots showing control, edited, and overlay views
- **Timestamped Outputs**: All results are timestamped to preserve analysis history
- **Organized Outputs**: Results are saved in dedicated output folders

## Adding New Genes

1. Create a new folder in `data/` (e.g., `data/mygene/`)
2. Add ALL required files:
   - `control.ab1`
   - `edited.ab1`
   - `grna.txt`
   - `mrna.txt`
3. Run `python run_analysis.py`
4. Find results in `data/mygene/output/`

## Important Notes

- **mrna.txt is REQUIRED**: The pipeline needs the full mRNA sequence to design proper sequencing primers
- **Timestamped outputs**: Each run creates new files with timestamps, preserving your analysis history
- **Output organization**: All generated files are saved in the `output/` subfolder of each gene directory
- **File naming**: Use exactly the names specified (control.ab1, edited.ab1, grna.txt, mrna.txt)

## Troubleshooting

- **"Missing required files"**: Ensure ALL four files (control.ab1, edited.ab1, grna.txt, mrna.txt) are present
- **"Could not find guide RNA in control sequence"**: Your sequencing doesn't cover the gRNA target region. Use the recommended primers to re-sequence
- **Noisy chromatograms**: Check DNA quality and sequencing conditions. Poly-G/T regions often cause issues

## License

This project is for research use only.

## Contact

For questions or issues, please contact the project maintainer. 