---
title: 'CasPINS: An Integrated Platform for CRISPR/TALEN gRNA Design, Primer Generation, and Indel Analysis'
tags:
  - Python
  - CRISPR
  - TALEN
  - genome editing
  - gRNA design
  - indel analysis
  - bioinformatics
authors:
  - name: Rinki Dasgupta
    orcid: 0009-0008-5478-840X
    corresponding: true
    affiliation: "1, 2"
  - name: Kaushik Das
    orcid: 0009-0009-8828-2631
    affiliation: 3
affiliations:
  - name: Department of Psychiatry and Neuroscience, Dell Medical Center at University of Texas at Austin, Austin, TX, USA
    index: 1
  - name: Department of Biology, Texas Woman's University, Denton, TX, USA
    index: 2
  - name: Department of Computer Applications, National Institute of Technology Jharkhand, Jamshedpur, India
    index: 3
date: 25 January 2026
bibliography: paper.bib
---

# Summary

CasPINS (Cas-Primer-Indel Suite) is an open-source Python platform that integrates the complete genome editing workflow—gRNA/TALEN design, primer generation, and indel analysis—into a single interactive web application. The platform supports 90+ species, 14 CRISPR-Cas variants, and provides both graphical and command-line interfaces, enabling researchers without programming expertise to efficiently design and validate genome editing experiments.

# Statement of Need

Genome editing with CRISPR-Cas [@doudna2014new; @jinek2012programmable] and TALENs [@miller2011tale] has transformed biological research, yet researchers must navigate multiple disconnected tools: gRNA design platforms (CHOPCHOP [@labun2019chopchop], CRISPOR [@concordet2018crispor]), primer design software (Primer3 [@untergasser2012primer3]), and indel analysis tools (TIDE [@brinkman2014easy]). This fragmented workflow creates inefficiencies, increases error potential during data transfer, and presents barriers for laboratories without bioinformatics support.

CasPINS addresses this gap by consolidating all computational steps into a unified interface. A researcher can progress from target gene identification through publication-ready editing efficiency data without exporting files, reformatting inputs, or learning multiple software interfaces.

# Software Design

CasPINS employs a modular architecture with four core modules (gRNA Design, TALEN Design, Primer Design, Indel Analysis) sharing common data structures for workflow continuity. We chose Streamlit for the web interface to enable interactive visualization without requiring client-side installation—a critical accessibility decision over desktop-only alternatives.

**Build vs. Contribute Justification:** We evaluated contributing to existing tools but found fundamental architectural barriers: CHOPCHOP and CRISPOR are web-only services without extensible codebases for primer integration; TIDE is closed-source; Primer3 lacks genome editing awareness. Rather than attempting incompatible integrations, CasPINS wraps Primer3 and BioPython [@biopython] within a unified workflow layer, extending their capabilities through CRISPR-aware positioning logic. This approach respects existing tool boundaries while providing the integration researchers need.

The NNLS-based indel decomposition was chosen over proprietary algorithms to ensure full reproducibility and avoid licensing restrictions that would limit academic use.

# Research Impact Statement

CasPINS enables researchers to complete genome editing workflows that previously required 3-4 separate tools and manual data transfer. The platform's accessibility—no installation beyond Python, no programming required—directly addresses the reproducibility crisis in genome editing by providing standardized, documented workflows.

**Community Readiness:** CasPINS provides comprehensive documentation (User Guide, Developer Guide, Deployment Guide), automated tests with CI/CD via GitHub Actions, Docker containerization for reproducible deployment, and MIT licensing for unrestricted academic and commercial use. The software is archived on Zenodo (DOI: 10.5281/zenodo.18370068) for long-term preservation.

**Credible Near-term Significance:** CasPINS fills a documented gap—no existing open-source tool integrates gRNA design, CRISPR-aware primer design, and indel quantification. The platform's 90+ species support and 14 Cas variants exceed alternatives, enabling research across model organisms, agriculture, and emerging species.

# Key Features

- **gRNA Design:** 14 Cas variants (SpCas9, SaCas9, Cas12a, etc.) with Doench [@doench2016optimized], Moreno-Mateos [@moreno2015crisprscan], and Xu [@xu2015sequence] scoring algorithms
- **TALEN Design:** RVD sequence generation, off-target analysis, restriction site mapping
- **Primer Design:** NCBI/Ensembl integration, nested PCR primers positioned relative to cut sites
- **Indel Analysis:** NNLS trace decomposition, single and batch modes, publication-ready visualizations

# AI Usage Disclosure

Generative AI tools (Claude, Anthropic) were used for manuscript drafting assistance and code documentation. All AI-generated content was reviewed, verified, and edited by the authors. AI was not used for algorithm development, data analysis, or scientific interpretation. The authors take full responsibility for all content.

# Acknowledgements

We thank the developers of Primer3, BioPython, and Streamlit. We acknowledge NCBI and Ensembl for public genome databases.

# References
