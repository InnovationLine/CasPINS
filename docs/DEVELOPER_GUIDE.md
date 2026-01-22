# CasPINS - Developer Guide

## Architecture Overview

CasPINS is organized into modular components:

```
src/
├── cli/                  # Command-line interfaces
│   ├── find_grna.py      # gRNA design CLI
│   ├── design_primers.py # Primer design CLI
│   └── run_analysis.py   # Indel analysis CLI
├── gui/                  # Streamlit GUI
│   ├── app.py            # Main GUI application
│   └── streamlit_entry.py # Cloud deployment entry
├── grna_design/          # gRNA design system
│   ├── grna_designer.py  # Main designer class
│   ├── talen_designer.py # TALEN design
│   ├── hdr_designer.py   # HDR template design
│   ├── core/             # Core functionality
│   ├── scoring/          # Scoring algorithms
│   └── database/         # Genome management
├── utils/                # Analysis utilities
│   ├── indel_analysis.py # NNLS decomposition
│   ├── primer_design.py  # Primer3 integration
│   └── visualization.py  # Plotting functions
└── config/               # Configuration
    ├── settings.py       # Global settings
    └── species.json      # Species database
```

---

## Algorithm Details

### Indel Analysis: NNLS Decomposition

Our implementation uses **Non-Negative Least Squares (NNLS)**:

1. **Build Reference Matrix**: Create shifted versions of control trace (-10bp to +10bp)
2. **Decompose Signal**: Use NNLS to find optimal linear combination
3. **Extract Frequencies**: Coefficients represent indel frequencies
4. **Calculate Efficiency**: 100% - wild-type fraction

```python
# Simplified approach
A = matrix_of_shifted_control_traces
coefficients = nnls(A, edited_trace)
editing_efficiency = (1 - coefficients[wild_type]) * 100
```

### gRNA Scoring System

Composite score combines multiple algorithms:

| Algorithm | Weight | Description |
|-----------|--------|-------------|
| Doench 2016 | 40% | Position-specific preferences |
| Moreno-Mateos | 30% | CRISPRscan algorithm |
| Xu Score | 20% | Machine learning-based |
| GC Content | 10% | Binding stability |

### TALEN Efficiency Scoring

TALEN pairs are scored based on:
- Arm lengths (optimal: 17-18bp)
- Spacer length (optimal: 14-16bp)
- GC content (optimal: 40-60%)
- RVD composition
- Homopolymer avoidance

---

## Key Classes

### GRNADesigner

```python
from grna_design.grna_designer import GRNADesigner

designer = GRNADesigner(species='human', cas_type='SpCas9')
results = designer.design_grnas(
    target='TP53',
    n_results=10,
    include_off_targets=True
)
```

### TALENDesigner

```python
from grna_design.talen_designer import TALENDesigner

designer = TALENDesigner(species='human')
results = designer.design_talens(
    target='ENSG00000123454',
    n_results=10
)
```

### PrimerDesigner

```python
from utils.primer_design import PrimerDesigner

designer = PrimerDesigner(gene_name='TP53', species='human')
primers = designer.design_primers(target_sequence, grna_positions)
```

---

## Fallback Methods

### Primary: NNLS Trace Decomposition
- High confidence, quantitative
- Requires high-quality control trace
- Provides detailed indel spectrum

### Fallback: Signal Decay Analysis
- Used when control trace quality is low
- Based on signal degradation patterns
- Lower confidence but more robust

---

## Adding New Features

### New Cas Type

1. Add PAM pattern to `src/grna_design/grna_designer.py`:
```python
CAS_TYPES = {
    'NewCas': {'pam': 'NNNN', 'pam_position': '3prime'}
}
```

2. Update scoring if needed in `src/grna_design/scoring/scoring_engine.py`

### New Species

1. Add to `src/config/species.json`:
```json
{
    "new_species": {
        "ensembl_name": "genus_species",
        "assembly": "assembly_name"
    }
}
```

---

## Testing

```bash
# Run all tests
python -m pytest tests/

# Run specific test
python -m pytest tests/test_grna_design.py -v
```

---

## Deployment

### Local
```bash
python run.py gui
```

### Streamlit Cloud
Set main file to: `src/gui/streamlit_entry.py`

### Docker
```bash
docker-compose up
```

---

## Dependencies

Core requirements:
- Python 3.8+
- BioPython (sequence handling)
- NumPy/SciPy (signal processing)
- Streamlit (GUI)
- Primer3-py (primer design)

See `requirements.txt` for full list.
