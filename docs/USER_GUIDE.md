# CRISPR Analysis Suite - User Guide

## Overview

The CRISPR Analysis Suite provides a complete workflow for CRISPR experiments:

1. **gRNA Design** - Find optimal guide RNAs for your target genes
2. **Primer Design** - Generate validation primers  
3. **Indel Analysis** - Quantify editing efficiency from sequencing data

This tool is **completely generic** - it works with any gene from any species.

---

## Quick Start

### Option 1: GUI (Recommended)

```bash
python run.py gui
```

### Option 2: Command Line

```bash
# Find gRNAs
python run.py grna TP53 --species human --top 10

# Design primers
python run.py primers TP53

# Run indel analysis
python run.py analysis --data-dir ./data
```

---

## Module 1: gRNA Design

### Basic Usage

```bash
# Find gRNAs for a gene
python run.py grna TP53

# Save directly to data folder
python run.py grna TP53 --save-to data/tp53/grna.txt

# Customize parameters
python run.py grna BRCA1 --top 5 --gc-min 45 --gc-max 55
```

### Command Options

| Option | Description | Default |
|--------|-------------|---------|
| `--top` | Number of gRNAs to return | 10 |
| `--species` | Target species | human |
| `--cas-type` | Cas protein type | SpCas9 |
| `--gc-min` | Minimum GC content (%) | 40 |
| `--gc-max` | Maximum GC content (%) | 60 |
| `--save-to` | Save directly to file | None |

### Supported Cas Types

| Cas Type | PAM | Notes |
|----------|-----|-------|
| SpCas9 | NGG | Most common |
| SaCas9 | NNGRRT | Smaller, AAV-compatible |
| Cas12a | TTTV | T-rich PAM |
| SpCas9-NG | NG | Relaxed PAM |

### TALEN Design

Select "TALEN" under Nuclease System in the GUI for TALEN pair design with:
- RVD sequence generation
- Off-target analysis
- Restriction site identification

---

## Module 2: Primer Design

### Usage

```bash
python run.py primers GENE_NAME --data-dir ./data
```

Or use the GUI Primer Design tab for interactive design with:
- PCR I primers (genomic DNA amplification)
- PCR II primers (sequencing)
- Database validation (NCBI, Ensembl)

---

## Module 3: Indel Analysis

### Prerequisites

1. Control and edited `.ab1` files in `data/<gene>/input/`
2. gRNA sequences in `data/<gene>/grna.txt`
3. mRNA reference in `data/<gene>/mrna.txt`

### Running Analysis

```bash
# Analyze all genes
python run.py analysis --data-dir ./data

# The pipeline automatically:
# 1. Finds all gene folders with AB1 files
# 2. Runs trace decomposition
# 3. Generates plots and results
```

### Output Files

- `indel_analysis_*.png` - Analysis plots
- `indel_analysis_summary_*.png` - Summary visualization

---

## Adding New Genes

### Step-by-Step

```bash
# 1. Create gene folder
mkdir data/your_gene

# 2. Add gRNA sequence (one per line)
echo "ACGTACGTACGTACGTACGT" > data/your_gene/grna.txt

# 3. Add mRNA sequence from NCBI/Ensembl
# Save to: data/your_gene/mrna.txt

# 4. Add your AB1 sequencing files
# data/your_gene/input/control.ab1
# data/your_gene/input/edited_sample1.ab1

# 5. Run analysis
python run.py analysis --data-dir ./data
```

### Data Directory Structure

```
your_data_directory/
├── gene_name/
│   ├── grna.txt         # gRNA sequences (one per line)
│   ├── mrna.txt         # mRNA reference sequence
│   └── input/           # AB1 sequencing files
│       ├── control.ab1  # Control/wild-type sample
│       └── edited*.ab1  # Edited samples
└── another_gene/
    └── ...
```

---

## Troubleshooting

### gRNA Design

**"Could not retrieve sequence"**
- Check gene name spelling
- Try using Ensembl ID (e.g., ENSG00000165646)
- Add mRNA sequence manually to `data/<gene>/mrna.txt`

**"No suitable gRNAs found"**
- Adjust GC range: `--gc-min 35 --gc-max 65`
- Try different Cas type
- Check sequence for low complexity regions

### Indel Analysis

**Analysis fails**
- Verify AB1 file quality
- Check gRNA is present in mRNA sequence
- Ensure control sample is clean

**Low confidence results**
- Use higher quality sequencing
- Sequence 400-700bp around cut site
- Include proper controls

---

## Best Practices

1. **gRNA Design**: Start with defaults, design 3-5 gRNAs per target
2. **Sequencing**: Use high-quality sequencing (>Q20)
3. **Controls**: Always include proper wild-type controls
4. **Validation**: Scores are predictive - validate experimentally
