# CRISPR TIDE Analysis Pipeline - Usage Guide

## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#installation)
- [Input File Requirements](#input-file-requirements)
- [Basic Usage](#basic-usage)
- [Command Line Options](#command-line-options)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Best Practices](#best-practices)
- [Multi-Sample Mode](#multi-sample-mode)

## Quick Start

```bash
# Activate virtual environment (Windows)
.\crispr_env\Scripts\activate

# Find gRNAs for rat TP53 with standard SpCas9
python run.py grna TP53 --species rat --top 10 --save-to data/tp53/grna.txt

# Find gRNAs with SpCas9-VQR (for NGA PAMs like CHOPCHOP)
python run.py grna TP53 --species rat --cas-type Cas9-VQR --gc-min 65 --gc-max 75

python run.py grna TP53 --save-to data/tp53/grna_new.txt --species rat --cas-type Cas9-VQR --gc-min 65 --gc-max 75 --top 5

# Run analysis on all genes with default settings
python run.py analysis

# Run with forced plotting despite validation issues
python run.py analysis --force-plot

# Run on specific gene with custom threshold
python run.py analysis --genes tp53 --similarity-threshold 50
```

## Installation

1. **Clone the repository**

   ```bash
   git clone <repository-url>
   cd crispr_analysis
   ```
2. **Create virtual environment**

   ```bash
   # Windows
   python -m venv crispr_env
   .\crispr_env\Scripts\activate

   # Linux/Mac
   python -m venv crispr_env
   source crispr_env/bin/activate
   ```
3. **Install dependencies**

   ```bash
   pip install -r requirements.txt
   ```

## Input File Requirements

Each gene folder must contain these four files:

```
data/
└── gene_name/
    ├── control.ab1      # Control sample chromatogram
    ├── edited.ab1       # Edited sample chromatogram
    ├── grna.txt         # Guide RNA sequences (one per line)
    └── mrna.txt         # mRNA reference sequence
```

### File Formats

**grna.txt** - Plain text, one gRNA per line:

```
ACGTCCAGTGTACCCTGACG
TGTCTATGCCATCGCCGATG
```

**mrna.txt** - Plain text, mRNA sequence (can be multi-line):

```
ATGGCGGCGGTGGCGGCGGCGGCGGCCGCCTCGGCCCCGGCGGCCGCC...
```

## Basic Usage

### Run all genes with default settings

```bash
python run.py analysis
```

- Uses 70% similarity threshold
- Creates error plots for validation failures
- Archives previous outputs
- Processes all genes in the data directory

### Run from src directory

```bash
python src/tide_batch_analysis.py
```

## Command Line Options

### `--force-plot`

Generate plots even when validation fails. Plots will include warning boxes.

```bash
python run.py analysis --force-plot
```

### `--similarity-threshold <value>`

Set minimum sequence similarity threshold (default: 70%)

```bash
python run.py analysis --similarity-threshold 50
```

### `--skip-validation`

Skip all validation checks (NOT RECOMMENDED - use only for debugging)

```bash
python run.py analysis --skip-validation
```

### `--interactive`

Ask for confirmation before processing each gene

```bash
python run.py analysis --interactive
```

### `--no-archive`

Do not archive previous output files

```bash
python run.py analysis --no-archive
```

**Note about archiving**: The pipeline archives all files from the output folder at the START of each run. Files are moved to `archive/run_YYYYMMDD_HHMMSS/` folders. This ensures the output folder only contains results from the current analysis.

### `--genes <gene1> [<gene2> ...]`

Process specific genes only

```bash
python run.py analysis --genes tp53 ddc
```

### `--help`

Show help message with all options

```bash
python run.py analysis --help
```

## Usage Examples

### Example 1: Standard Analysis

Run analysis on all genes with default settings:

```bash
python run.py analysis
```

Output:

```
Starting CRISPR TIDE Analysis Pipeline...
------------------------------------------------------------

============================================================
TIDE ANALYSIS BATCH PIPELINE
============================================================
Settings:
  Similarity threshold: 70.0%
  Force plot: False
  Skip validation: False
  Interactive mode: False
  Archive previous files: True

Processing tp53...
  Found 2 gRNA sequences
  [SUCCESS] AB1 analysis saved: data\tp53\output\ab1_analysis_tp53.txt
  ...
```

### Example 2: Force Plotting with Low Similarity

When your sequences have low similarity but you still want to see the analysis:

```bash
python run.py analysis --force-plot --similarity-threshold 25
```

This will:

- Generate plots even with 25% similarity
- Include warning messages in the plots
- Show which validation rules failed

### Example 3: Interactive Mode for Specific Genes

Process only selected genes with confirmation:

```bash
python run.py analysis --interactive --genes tp53 brca1
```

Output:

```
Processing tp53...
  Process tp53? (y/n): y
  Found 2 gRNA sequences
  ...
  
Processing brca1...
  Process brca1? (y/n): n
  Skipping brca1
```

### Example 4: Quick Re-run Without Archiving

For rapid iteration during development:

```bash
python run.py analysis --no-archive --genes ddc
```

### Example 5: Debug Mode with All Information

Skip validation to see what happens (use carefully):

```bash
python run.py analysis --skip-validation --force-plot
```

### Example 6: Custom Pipeline for Specific Requirements

Low similarity threshold, specific gene, force plotting:

```bash
python run.py analysis --genes tp53 --similarity-threshold 30 --force-plot
```

## Output Files

### Directory Structure

```
data/gene_name/
├── archive/                           # Previous runs
│   └── run_20240620_143822/         # Timestamped folders
│       ├── ab1_analysis_gene.txt
│       ├── recommendations_gene_timestamp.txt
│       └── tide_analysis_gene_timestamp.png
└── output/                           # Current run
    ├── ab1_analysis_gene.txt
    ├── recommendations_gene_timestamp.txt
    └── tide_analysis_gene_timestamp.png
```

### Output File Types

#### 1. `ab1_analysis_gene.txt`

Detailed AB1 file analysis including:

- Sequence lengths
- Quality scores (Phred)
- First/last 50bp of sequences
- Quality distribution
- Sequence comparison
- gRNA presence analysis

Example content:

```
DETAILED AB1 ANALYSIS FOR TP53
============================================================
Analysis Date: 2024-06-20 14:38:22

CONTROL AB1 FILE:
----------------------------------------
  Sequence length: 594 bp
  First 50 bp: GGGGGTTTTTAATCAAAACATCAACGGATGTTGGATCCCACTTTGTTCCG
  Average quality score: 13.7
  Quality distribution:
    Low (<20): 522 bases (87.9%)
    Medium (20-30): 52 bases (8.8%)
    High (>30): 20 bases (3.4%)
```

#### 2. `recommendations_gene_timestamp.txt`

Comprehensive analysis report with:

- gRNA summary and positions
- TIDE analysis results
- Editing efficiency
- Primer recommendations
- PCR conditions

Example content:

```
============================================================
CRISPR ANALYSIS REPORT - TP53
============================================================
Generated: 2024-06-20 14:38:22

gRNA SUMMARY:
----------------------------------------
gRNA 1: ACGATAACACCCCCAGTAGA
  Position in mRNA: 1326
  Strand: reverse
  Cut site: 1329

TIDE ANALYSIS RESULTS:
----------------------------------------
Editing Efficiency: 75.3%
Dominant Indel: -3bp (45.2%)
Quality Score: 82%
Confidence: HIGH - Clear editing pattern
```

#### 3. `tide_analysis_gene_timestamp.png`

Four-panel visualization showing:

- **Top Left**: Chromatogram overlay (control vs edited)
- **Top Right**: Sequence alignment visualization
- **Bottom Left**: Indel spectrum or signal analysis
- **Bottom Right**: Editing efficiency pie chart

File naming indicates status:

- `tide_analysis_gene_timestamp.png` - Normal analysis
- `tide_analysis_gene_timestamp_WARNING.png` - Analysis with warnings
- `tide_analysis_gene_timestamp_ERROR.png` - Failed validation

## Troubleshooting

### Common Issues and Solutions

#### 1. "Missing required files"

```
ERROR: Missing required files in data/tp53:
  - mrna.txt (mRNA reference sequence)
```

**Solution**: Ensure all four required files exist in the gene folder

#### 2. "gRNA not found in mRNA sequence"

```
ERROR: gRNA 1 (ACGTCCAGTGTACCCTGACG) not found in mRNA sequence!
```

**Solutions**:

- Check if gRNA includes PAM sequence (remove if present)
- Verify gRNA is exactly 20bp
- Ensure correct mRNA sequence
- Check if using genomic vs cDNA sequence

#### 3. "Low sequence similarity"

```
Sequence similarity: 27.6%
ERROR: Control and edited sequences have very low similarity
```

**Solutions**:

- Verify control and edited samples are from same region
- Check primer design
- Use `--force-plot` to see analysis anyway
- Lower similarity threshold: `--similarity-threshold 25`

#### 4. "Cut site beyond sequence length"

```
Expected cut site (1329) is beyond sequence length (594)
```

**Solution**: Need longer sequencing reads that cover the gRNA target region

#### 5. Import errors

```
ImportError: No module named 'Bio'
```

**Solution**: Activate virtual environment and install requirements:

```bash
.\crispr_env\Scripts\activate  # Windows
pip install -r requirements.txt
```

## Best Practices

### 1. Data Validation

- Always check `ab1_analysis_*.txt` files first
- Verify sequence quality scores
- Ensure gRNAs are found in reference

### 2. Analysis Strategy

- Start with default settings
- Use `--force-plot` to investigate issues
- Adjust similarity threshold based on your data
- Archive files for analysis history

### 3. Quality Control

- Aim for >70% sequence similarity
- Check quality scores in AB1 files
- Verify cut sites are within sequenced region
- Use multiple gRNAs for validation

### 4. Workflow Recommendations

```bash
# 1. First run - check data quality
python run.py analysis

# 2. If validation fails, investigate
python run.py analysis --force-plot

# 3. Adjust thresholds if needed
python run.py analysis --similarity-threshold 50 --force-plot

# 4. Process specific genes after fixes
python run.py analysis --genes tp53 --no-archive
```

### 5. Interpreting Results

- **High Confidence**: Clear editing pattern, good quality scores
- **Medium Confidence**: Some validation warnings, but analysis possible
- **Low Confidence**: Fallback methods used, interpret with caution

## Advanced Usage

### Running from different directories

```bash
# From project root
python run.py analysis

# From src directory
cd src
python tide_batch_analysis.py

# With full path
C:\path\to\python.exe C:\path\to\run_analysis.py
```

### Batch processing with logging

```bash
# Windows
python run.py analysis > analysis_log.txt 2>&1

# Linux/Mac
python run.py analysis 2>&1 | tee analysis_log.txt
```

### Automated pipeline

```bash
# Process all genes without interaction
python run.py analysis --skip-validation --force-plot --no-archive
```

## Support

For issues or questions:

1. Check error messages in console output
2. Review `ab1_analysis_*.txt` files
3. Examine plot warning messages
4. Verify input file formats
5. Check GitHub issues or create new one

---

*Last updated: December 2024*

# CRISPR TIDE Analysis Pipeline - Multi-Sample Mode

## Overview

This pipeline analyzes multiple clonally expanded CRISPR-edited cell lines against a control sample using the TIDE (Tracking of Indels by DEcomposition) algorithm. It generates comprehensive visualizations similar to the TIDE web tool (https://apps.datacurators.nl/tide/).

## Directory Structure

The pipeline expects the following directory structure:

```
data/
├── gene1/
│   ├── input/
│   │   ├── control.ab1       # Control sample (unedited)
│   │   ├── editedA2.ab1      # Edited clone A2
│   │   ├── editedA3.ab1      # Edited clone A3
│   │   └── edited*.ab1       # Any number of edited clones
│   ├── grna.txt              # gRNA sequence (20bp)
│   └── mrna.txt              # mRNA reference sequence
├── gene2/
│   └── ...
└── gene3/
    └── ...
```

### Required Files

1. **input/** folder containing:

   - `control.ab1`: Control sample chromatogram
   - `edited*.ab1`: One or more edited sample chromatograms
2. **grna.txt**: Contains the 20bp guide RNA sequence (without PAM)

   ```
   GCCGTCCCGAGTACAGCCAG
   ```
3. **mrna.txt**: Contains the mRNA reference sequence

## Running the Analysis

### Basic Usage

```bash
# Activate the virtual environment
.\crispr_env\Scripts\activate  # Windows
source crispr_env/bin/activate  # Linux/Mac

# Run the analysis
python run.py analysis
```

### Command Line Options

```bash
# Skip specific genes
python src/tide_batch_analysis_multi.py --skip-genes gene1

# Force analysis even without input folder
python src/tide_batch_analysis_multi.py --force

# Don't archive previous results
python src/tide_batch_analysis_multi.py --no-archive
```

## Output Files

For each gene with an input folder, the pipeline generates:

### Individual Sample Analysis

- `tide_[gene]_[sample]_[timestamp].png`: TIDE analysis plot for each edited sample
  - Input visualization with alignment and decomposition windows
  - Alignment window detail (before cut site)
  - Decomposition window detail (after cut site)
  - Indel spectrum bar chart
  - Editing efficiency pie chart

### Summary Reports

- `tide_summary_[gene]_[timestamp].png`: Overview of all samples with efficiency pie charts
- `tide_results_[gene]_[timestamp].json`: Detailed analysis results in JSON format

### Example Output Structure

```
data/
└── gene1/
    └── output/
        ├── tide_gene1_editedA2_20250703_194117.png
        ├── tide_gene1_editedA3_20250703_194117.png
        ├── tide_summary_gene1_20250703_194117.png
        └── tide_results_gene1_20250703_194117.json
```

## Understanding the Results

### TIDE Analysis Plot Components

1. **Input Visualization (Top Panel)**

   - Shows the full chromatogram trace
   - Green shaded area: Alignment window (used for aligning control and edited sequences)
   - Orange shaded area: Decomposition window (used for TIDE analysis)
   - Red dashed line: Expected cut site
   - Purple bar: gRNA location
2. **Alignment Window (Middle Left)**

   - Detailed view of the sequence before the cut site
   - Blue trace: Control sample
   - Red trace: Edited sample
   - Green bars: Sequence match indicator
3. **Decomposition Window (Middle Right)**

   - Detailed view of the sequence after the cut site
   - Shows divergence between control and edited samples
   - Yellow area: Signal difference
4. **Indel Spectrum (Bottom Left)**

   - Bar chart showing frequency of different indel sizes
   - Red bars: Insertions
   - Blue bars: Deletions
   - Gold border: Dominant indel
5. **Efficiency Pie Chart (Bottom Right)**

   - Visual representation of editing efficiency
   - Red: Edited cells
   - Blue: Unedited cells

### JSON Results Format

```json
{
  "gene": "tp53",
  "timestamp": "20250703_194117",
  "grna_sequences": ["GCCGTCCCGAGTACAGCCAG"],
  "num_samples": 12,
  "samples": [
    {
      "editing_efficiency": 75.5,
      "dominant_indel_size": -3,
      "dominant_indel_percent": 45.2,
      "indel_spectrum": {
        "-3": 45.2,
        "-1": 15.3,
        "+1": 10.0
      },
      "quality_score": 85.0,
      "method": "TIDE",
      "sequence_similarity": 95.5,
      "sample_name": "editedA2",
      "expected_cut_site": 832
    }
  ]
}
```

## Troubleshooting

### No input folder found

- Ensure the `input/` folder exists under the gene directory
- Check that it contains `control.ab1` and at least one `edited*.ab1` file

### Low quality scores

- Verify sequence quality in the AB1 files
- Ensure the gRNA sequence is correct
- Check that control and edited samples are from the same amplicon

### gRNA not found

- Verify the gRNA sequence in `grna.txt` is correct (20bp, no PAM)
- Check that the mRNA sequence contains the target site

## Technical Details

### TIDE Algorithm

The pipeline implements the TIDE decomposition algorithm to:

1. Align control and edited chromatogram traces
2. Decompose the edited trace into a linear combination of shifted control traces
3. Calculate the frequency of each indel size (-10 to +10 bp)
4. Determine overall editing efficiency

### Quality Metrics

- **Editing Efficiency**: Percentage of cells with indels
- **Dominant Indel**: Most frequent indel size and its percentage
- **Quality Score**: Confidence in the decomposition (higher is better)
- **Sequence Similarity**: Alignment quality between control and edited samples

## References

Brinkman et al., "Easy quantitative assessment of genome editing by sequence trace decomposition." Nucleic Acids Research, 2014.

TIDE Web Tool: https://apps.datacurators.nl/tide/
