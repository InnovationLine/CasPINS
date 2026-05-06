# Feature Comparison: CasPINS vs Existing Genome Editing Tools

This table directly addresses NAR criticism (iii): "there are no detailed examples and there is no thorough benchmarking."

## Table 1. Feature Comparison of Genome Editing Software Tools

| Feature | CasPINS | CHOPCHOP v3 | CRISPOR | Benchling | Primer3/Primer-BLAST | TIDE | ICE |
|---------|---------|-------------|---------|-----------|---------------------|------|-----|
| **gRNA Design** | | | | | | | |
| gRNA identification | Yes | Yes | Yes | Yes | No | No | No |
| On-target scoring (Doench 2016) | Yes | Yes | Yes | Yes | No | No | No |
| On-target scoring (Moreno-Mateos) | Yes | Yes | Yes | No | No | No | No |
| On-target scoring (Xu) | Yes | No | Yes | No | No | No | No |
| Off-target analysis | Yes | Yes | Yes | Yes | No | No | No |
| Multiple Cas variants | 14 variants | 4 (Cas9, Cpf1, Cas13, Nickase) | 10+ (incl. Cas12Max) | 3 variants | N/A | N/A | N/A |
| Species support | 90+ | 162+ | ~150+ | Limited | N/A | N/A | N/A |
| Cas12a/Cpf1 support | Yes | Yes | Yes | Yes | No | No | No |
| High-fidelity Cas9 variants | Yes (eSpCas9, HF1, Hypa) | No | Yes | No | No | No | No |
| Multiple editing modes (KO, KI, CRISPRa/i, BE, PE) | Yes | Limited | No | Yes | No | No | No |
| **TALEN Design** | | | | | | | |
| TALEN pair identification | Yes | No | No | No | No | No | No |
| RVD sequence generation | Yes | No | No | No | No | No | No |
| TALEN off-target analysis | Yes | No | No | No | No | No | No |
| TALEN restriction site mapping | Yes | No | No | No | No | No | No |
| **Primer Design** | | | | | | | |
| PCR primer design | Yes | No | No | Yes | Yes | No | No |
| Sequencing primer design | Yes | No | No | No | Yes | No | No |
| CRISPR-aware primer positioning | Yes | No | No | No | No | No | No |
| Nested PCR design (PCR I + PCR II) | Yes | No | No | No | No | No | No |
| NCBI/Ensembl database integration | Yes | Ensembl | UCSC | Internal | NCBI | No | No |
| Automated cut-site-relative positioning | Yes | No | No | No | No | No | No |
| **Indel Analysis** | | | | | | | |
| Sanger trace decomposition | Yes (NNLS) | No | No | No | No | Yes (NNLS) | Yes (ML) |
| Editing efficiency quantification | Yes | No | No | Yes* | No | Yes | Yes |
| Indel spectrum characterization | Yes | No | No | No | No | Yes | Yes |
| AB1 file direct input | Yes | No | No | No | No | Yes | Yes |
| Batch/multi-sample analysis | Unlimited | N/A | N/A | N/A | N/A | No (1 pair) | Yes (up to 700) |
| **Workflow Integration** | | | | | | | |
| Integrated gRNA + Primer design | Yes | No | No | Partial | No | No | No |
| Integrated gRNA + Indel analysis | Yes | No | No | No | No | No | No |
| Integrated Primer + Indel analysis | Yes | No | No | No | No | No | No |
| End-to-end workflow (design to analysis) | **Yes** | No | No | Partial | No | No | No |
| Shared data architecture across modules | Yes | N/A | N/A | Yes | N/A | N/A | N/A |
| **Accessibility & Deployment** | | | | | | | |
| Open source | Yes (MIT) | Yes | Yes | No | Yes | Partial (R code on request) | Yes (GitHub) |
| Local deployment | Yes | Yes | Yes | No (cloud) | Yes | No (web) | No (web) |
| Web GUI | Yes (Streamlit) | Yes | Yes | Yes | No | Yes | Yes |
| Command-line interface | Yes | Yes | Yes | No | Yes | No | No |
| Docker support | Yes | Yes | No | N/A | No | No | No |
| No internet required (local mode) | Yes | Partial | No | No | Yes | No | No |
| Batch processing via CLI | Yes | Yes | Limited | No | Yes | No | No |
| Export formats (CSV/JSON/GenBank) | Yes | Yes | Yes | Yes | Text | PDF | CSV |

*Benchling provides NGS-based analysis, not Sanger trace decomposition.

## Key Differentiators

### 1. Workflow Integration (Unique to CasPINS)
No existing open-source tool provides end-to-end integration from gRNA design through primer engineering to indel quantification. CasPINS is the only platform where:
- gRNA sequences selected during design are **automatically available** for primer positioning
- Sample organization established during project setup **propagates** to indel analysis
- All intermediate and final outputs are **preserved in a unified project structure**

### 2. TALEN + CRISPR in One Platform
CasPINS uniquely combines both CRISPR and TALEN design in a single interface, enabling researchers to evaluate both nuclease systems for a given target.

### 3. Unlimited Batch Indel Analysis
While TIDE processes one sample pair at a time and ICE supports up to 96 samples, CasPINS supports unlimited batch analysis with automatic sample detection and summary statistics.

### 4. Open-Source NNLS Decomposition
CasPINS implements fully open, reproducible NNLS-based trace decomposition, eliminating dependency on proprietary algorithms (TIDE is not open source; ICE is free but proprietary).

---

## How to Cite This Table

This feature comparison table is maintained at: https://github.com/InnovationLine/CasPINS/blob/main/benchmarking/feature_comparison_table.md

For the manuscript, this table should be included as **Table 1** or as a **Supplementary Table** depending on journal formatting requirements.
